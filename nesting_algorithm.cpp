// nesting_clip2.cpp — NFP‑greedy nesting (актуальный Clipper2 v2.x)
// --------------------------------------------------------------------
//   • DXF‑парсер  (минимум LWPOLYLINE) → Paths64 (+ дырки)
//   • 36 ориентаций (0…350° шаг 10°)
//   • NFP через Clipper2::MinkowskiSum (clipper.minkowski.h)
//   • Greedy corner‑search, размещение в отверстия
//   • CLI:  nest -s 3000x1500 -r 10 *.dxf -o layout.csv
//   • Header‑only зависимость: Clipper2 (https://github.com/AngusJohnson/Clipper2)
// --------------------------------------------------------------------
//  Build (g++ / MSVC):
//      g++ -std=c++17 -O3 nesting_clip2.cpp -I./Clipper3/CPP/Clipper2Lib/include -o nest
// --------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <cstdint>

#define _USE_MATH_DEFINES
#include <math.h>

#include "clipper3/clipper.h"
#include "clipper3/clipper.minkowski.h"
using namespace Clipper2Lib;

// ───────── helpers ─────────
constexpr double SCALE = 1000.0;              // мм → int64
inline int64_t I64(double v){ return llround(v * SCALE); }
inline double  Dbl(int64_t v){ return double(v) / SCALE; }

static Rect64 getBBox(const Paths64& P){
    Rect64 r{INT64_MAX,INT64_MAX,INT64_MIN,INT64_MIN};
    for(const auto&path:P) for(auto pt:path){
        r.left = std::min(r.left, pt.x);
        r.bottom = std::min(r.bottom, pt.y);
        r.right = std::max(r.right, pt.x);
        r.top = std::max(r.top, pt.y);
    }
    return r;
}

static Paths64 unionAll(const Paths64& in){
    return Union(in, FillRule::NonZero, 1.0);   
}

static bool overlap(const Paths64&a,const Paths64&b){
    return !Intersect(a,b,FillRule::NonZero,1.0).empty();
}

// ───────── DXF mini loader ─────────
struct RawPart{ Paths64 rings; double area=0; int id=0; };

static std::vector<RawPart> loadDXF(const std::vector<std::string>& files){
    std::vector<RawPart> parts;
    for(size_t idx=0; idx<files.size(); ++idx){
        std::ifstream fin(files[idx]); if(!fin) throw std::runtime_error("open "+files[idx]);
        std::string c,v; bool inEnt=false; Path64 cur; Paths64 rings;
        auto nxt=[&](){ return bool(std::getline(fin,c)&&std::getline(fin,v)); };
        while(nxt()){
            if(c=="0"&&v=="SECTION"){ nxt(); inEnt=(c=="2"&&v=="ENTITIES"); continue; }
            if(!inEnt||c!="0") continue;
            if(v=="LWPOLYLINE"){
                cur.clear(); bool closed=false;
                while(nxt()){
                    if(c=="0") break;
                    if(c=="70") closed=(std::stoi(v)&1);
                    else if(c=="10"){ double x=std::stod(v); nxt(); double y=std::stod(v); cur.push_back({I64(x),I64(y)}); }
                }
                if(closed) rings.push_back(cur);
            }
        }
        if(rings.empty()) throw std::runtime_error("no rings in "+files[idx]);
        std::sort(rings.begin(),rings.end(),[](auto&a,auto&b){return std::abs(Area(a))>std::abs(Area(b));});
        if(Area(rings[0])<0) ReversePath(rings[0]);
        for(size_t i=1;i<rings.size();++i) if(Area(rings[i])>0) ReversePath(rings[i]);
        parts.push_back({rings,std::abs(Area(rings[0]))/(SCALE*SCALE),int(idx)});
    }
    return parts;
}

// ───────── orientations ─────────
struct Orient{ Paths64 poly; Rect64 bb; int id; int ang; };

static Paths64 rot(const Paths64&src,double s,double c){
    Paths64 out; out.reserve(src.size());
    for(const auto&ring:src){ Path64 r; r.reserve(ring.size());
        for(auto pt:ring){ double x=Dbl(pt.x), y=Dbl(pt.y); r.push_back({I64(x*c - y*s), I64(x*s + y*c)}); }
        out.push_back(std::move(r));
    }
    return out;
}

static std::vector<std::vector<Orient>> makeOrient(const std::vector<RawPart>&parts,int step){
    std::vector<std::vector<Orient>> all;
    for(const auto &p:parts){ std::vector<Orient> ovec;
        for(int ang=0;ang<360;ang+=step){ double rad=ang*M_PI/180.0; auto r=rot(p.rings,sin(rad),cos(rad));
            ovec.push_back({ unionAll(r), getBBox(r), p.id, ang }); }
        all.push_back(std::move(ovec)); }
    return all;
}

// ───────── NFP cache ─────────
using Key=uint64_t; static std::unordered_map<Key,Paths64> NFP;
static Key kfn(int a,int ra,int b,int rb){return ((uint64_t)a<<48)^((uint64_t)ra<<32)^((uint64_t)b<<16)^rb;}
static const Paths64& nfp(const Orient&A,const Orient&B){
    Key k=kfn(A.id,A.ang,B.id,B.ang); auto it=NFP.find(k); if(it!=NFP.end()) return it->second;
    return NFP[k]=MinkowskiSum(A.poly[0],B.poly);  

// ───────── greedy placer ─────────
struct Place{int id; double x,y; int ang;};

static std::vector<Place> greedy(const std::vector<RawPart>&parts,double W,double H,int step){
    auto ord=parts; std::sort(ord.begin(),ord.end(),[](auto&a,auto&b){return a.area>b.area;});
    auto orient=makeOrient(ord,step);
    Path64 sheet={{0,0},{I64(W),0},{I64(W),I64(H)},{0,I64(H)}};
    std::vector<Paths64> placed; std::vector<Place> out;
    for(size_t i=0;i<ord.size();++i){ bool ok=false; Place best{};
        for(auto &op:orient[i]){
            std::vector<Point64> cand={{0,0}}; for(auto&pp:placed) for(auto&pt:pp[0]) cand.push_back(pt);
            for(auto c:cand){ Rect64 bb=op.bb; bb.left+=c.x; bb.right+=c.x; bb.bottom+=c.y; bb.top+=c.y;
                if(bb.right>sheet[2].x||bb.top>sheet[2].y||bb.left<0||bb.bottom<0) continue;
                Paths64 sh; sh.reserve(op.poly.size()); for(auto&ring:op.poly){ Path64 r; r.reserve(ring.size()); for(auto pt:ring) r.push_back({pt.x+c.x, pt.y+c.y}); sh.push_back(std::move(r)); }
                bool clash=false; for(auto&pl:placed){ if(overlap(pl,sh)){ clash=true; break; }} if(clash) continue;
                best={op.id,Dbl(c.x),Dbl(c.y),op.ang}; placed.push_back(std::move(sh)); ok=true; break; }
            if(ok) break; }
        if(ok) out.push_back(best); else std::cerr<<"skip "<<ord[i].id<<"\n";
    }
    return out;
}

// ───────── CLI ─────────
struct CLI{double W=0,H=0;int rot=10;std::string out="layout.csv";std::vector<std::string> files;};
static CLI parse(int ac,char**av){CLI c;for(int i=1;i<ac;++i){std::string a=av[i];if(a=="-s"||a=="--sheet"){auto s=std::string(av[++i]);auto x=s.find('x');c.W=std::stod(s.substr(0,x));c.H=std::stod(s.substr(x+1));}else if(a=="-r"||a=="--rot") c.rot=std::stoi(av[++i]); else if(a=="-o"||a=="--out") c.out=av[++i]; else c.files.push_back(a);} if(c.W==0||c.H==0||c.files.empty()) throw std::runtime_error("usage: nest -s WxH [-r 10] *.dxf -o out"); return c;}

// ───────── main ─────────
int main(int argc,char*argv[]){ try{
        auto cli=parse(argc,argv);
        auto parts=loadDXF(cli.files);
        auto placed=greedy(parts,cli.W,cli.H,cli.rot);
        std::ofstream f(cli.out); f<<"part,x_mm,y_mm,angle\n"; for(auto&p:placed) f<<p.id<<","<<std::fixed<<std::setprecision(3)<<p.x<<","<<p.y<<","<<p.ang<<"\n";
        std::cout<<"placed "<<placed.size()<<"/"<<parts.size()<<" → "<<cli.out<<"\n";
        return 0;
    }catch(const std::exception&e){ std::cerr<<"ERR: "<<e.what()<<"\n"; return 1; }}
