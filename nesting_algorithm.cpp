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
constexpr int    DEFAULT_ROT_STEP = 10;
inline int64_t I64(double v){ return llround(v * SCALE); }
inline double  Dbl(int64_t v){ return double(v) / SCALE; }

// Reverse polygon orientation in place

// Reverse polygon orientation
template <typename T>
void ReversePath(std::vector<T>& p){
    std::reverse(p.begin(), p.end());
}

// Compute bounding box of multiple polygons

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

// Union all rings to simplify part geometry
// Merge all polygons in a set
static Paths64 unionAll(const Paths64& in){
    return Union(in, FillRule::NonZero);
}

// Test for polygon overlap

static bool overlap(const Paths64&a,const Paths64&b){
    return !Intersect(a,b,FillRule::NonZero).empty();
}

// ───────── DXF mini loader ─────────
struct RawPart{ Paths64 rings; double area=0; int id=0; };

// Geometry helpers ------------------------------------------------------

// Approximate circular arc in degrees
static Path64 approxArc(PointD c,double r,double a1,double a2,double step=5.0){
    if(a2<a1) a2+=360.0; int seg=std::max(1,int(std::ceil((a2-a1)/step)));
    Path64 p; p.reserve(seg+1);
    for(int i=0;i<=seg;++i){ double ang=(a1+(a2-a1)*i/seg)*M_PI/180.0;
        p.push_back({I64(c.x+r*std::cos(ang)),I64(c.y+r*std::sin(ang))}); }
    return p;
}

// Approximate ellipse arc using major axis and ratio
static Path64 approxEllipse(PointD c,PointD maj,double ratio,double a1,double a2,double step=5.0){
    double ml=std::hypot(maj.x,maj.y); if(ml==0) return {};
    PointD u{maj.x/ml,maj.y/ml}; double b=ml*ratio;
    if(a2<a1) a2+=360.0; int seg=std::max(1,int(std::ceil((a2-a1)/step)));
    Path64 p; p.reserve(seg+1);
    for(int i=0;i<=seg;++i){ double ang=(a1+(a2-a1)*i/seg)*M_PI/180.0; double ca=std::cos(ang),sa=std::sin(ang);
        double x=c.x+maj.x*ca - b*u.y*sa; double y=c.y+maj.y*ca + b*u.x*sa;
        p.push_back({I64(x),I64(y)}); }
    return p;
}

// Simple Catmull-Rom spline approximation for SPLINE fit points
static Path64 approxSpline(const std::vector<PointD>& pts,int segPer=8){
    if(pts.size()<2) return {};
    Path64 out; auto CR=[&](const PointD&p0,const PointD&p1,const PointD&p2,const PointD&p3,double t){
        double t2=t*t,t3=t2*t; return PointD{
            0.5*((2*p1.x)+(-p0.x+p2.x)*t+(2*p0.x-5*p1.x+4*p2.x-p3.x)*t2+(-p0.x+3*p1.x-3*p2.x+p3.x)*t3),
            0.5*((2*p1.y)+(-p0.y+p2.y)*t+(2*p0.y-5*p1.y+4*p2.y-p3.y)*t2+(-p0.y+3*p1.y-3*p2.y+p3.y)*t3)}; };
    for(size_t i=0;i<pts.size()-1;++i){
        PointD p0=i?pts[i-1]:pts[i], p1=pts[i], p2=pts[i+1], p3=(i+2<pts.size()?pts[i+2]:pts[i+1]);
        for(int s=0;s<segPer;++s){ double t=double(s)/segPer; auto p=CR(p0,p1,p2,p3,t); out.push_back({I64(p.x),I64(p.y)}); }
    }
    out.push_back({I64(pts.back().x),I64(pts.back().y)}); return out;
}

static constexpr int64_t JOIN_TOL=10000; // coordinate tolerance
static bool eqPt(const Point64&a,const Point64&b){ return llabs(a.x-b.x)<=JOIN_TOL && llabs(a.y-b.y)<=JOIN_TOL; }

// Small helper to read DXF code/value pairs with single element lookahead
struct DXFReader{
    std::ifstream& in;
    std::string code,value;
    bool has=false;
    DXFReader(std::ifstream&f):in(f){}
    bool next(){
        if(has){ has=false; return true; }
        return bool(std::getline(in,code) && std::getline(in,value));
    }
    void push(){ has=true; }
};

// Read pair helper
static bool expect(DXFReader&rd,const char* need,int seg,const char* what){
    if(!rd.next()){ std::cerr<<"seg "<<seg<<" unexpected EOF while reading "<<what<<"\n"; return false; }
    if(rd.code!=need) std::cerr<<"seg "<<seg<<" expected "<<need<<" after "<<what<<" got "<<rd.code<<"\n"; 
    return true;
}

// Parse LINE entity
static bool parseLine(DXFReader&rd,Path64&out,int seg){
    double x1=0,y1=0,x2=0,y2=0; bool p1=false,p2=false; out.clear();
    while(rd.next()){
        if(rd.code=="0"){ rd.push(); break; }
        if(rd.code=="10"){ x1=std::stod(rd.value); if(expect(rd,"20",seg,"10")) y1=std::stod(rd.value); p1=true; }
        else if(rd.code=="11"){ x2=std::stod(rd.value); if(expect(rd,"21",seg,"11")) y2=std::stod(rd.value); p2=true; }
    }
    if(!p1||!p2){ std::cerr<<"seg "<<seg<<" LINE missing point\n"; return false; }
    out={{I64(x1),I64(y1)},{I64(x2),I64(y2)}}; return true;
}

// Parse ARC entity
static bool parseArc(DXFReader&rd,Path64&out,int seg){
    PointD cen{0,0}; double radius=0,a1=0,a2=0; bool haveC=false,haveR=false; out.clear();
    while(rd.next()){
        if(rd.code=="0"){ rd.push(); break; }
        if(rd.code=="10"){ cen.x=std::stod(rd.value); if(expect(rd,"20",seg,"10")) cen.y=std::stod(rd.value); haveC=true; }
        else if(rd.code=="40"){ radius=std::stod(rd.value); haveR=true; }
        else if(rd.code=="50") a1=std::stod(rd.value);
        else if(rd.code=="51") a2=std::stod(rd.value);
    }
    if(!haveC||!haveR){ std::cerr<<"seg "<<seg<<" ARC missing data\n"; return false; }
    out=approxArc(cen,radius,a1,a2); return !out.empty();
}

// Parse CIRCLE entity
static bool parseCircle(DXFReader&rd,Path64&out,int seg){
    PointD cen{0,0}; double radius=0; bool haveC=false,haveR=false; out.clear();
    while(rd.next()){
        if(rd.code=="0"){ rd.push(); break; }
        if(rd.code=="10"){ cen.x=std::stod(rd.value); if(expect(rd,"20",seg,"10")) cen.y=std::stod(rd.value); haveC=true; }
        else if(rd.code=="40"){ radius=std::stod(rd.value); haveR=true; }
    }
    if(!haveC||!haveR){ std::cerr<<"seg "<<seg<<" CIRCLE missing data\n"; return false; }
    out=approxArc(cen,radius,0,360); return true;
}

// Parse ELLIPSE entity
static bool parseEllipse(DXFReader&rd,Path64&out,int seg){
    PointD cen{0,0},maj{0,0}; double ratio=1,a1=0,a2=0; bool haveC=false,haveMaj=false; out.clear();
    while(rd.next()){
        if(rd.code=="0"){ rd.push(); break; }
        if(rd.code=="10"){ cen.x=std::stod(rd.value); if(expect(rd,"20",seg,"10")) cen.y=std::stod(rd.value); haveC=true; }
        else if(rd.code=="11"){ maj.x=std::stod(rd.value); if(expect(rd,"21",seg,"11")) maj.y=std::stod(rd.value); haveMaj=true; }
        else if(rd.code=="40") ratio=std::stod(rd.value);
        else if(rd.code=="41") a1=std::stod(rd.value)*180.0/M_PI;
        else if(rd.code=="42") a2=std::stod(rd.value)*180.0/M_PI;
    }
    if(!haveC||!haveMaj){ std::cerr<<"seg "<<seg<<" ELLIPSE missing data\n"; return false; }
    out=approxEllipse(cen,maj,ratio,a1,a2); return !out.empty();
}

// Parse SPLINE entity using fit points only
static bool parseSpline(DXFReader&rd,Path64&out,int seg){
    std::vector<PointD> pts; out.clear();
    while(rd.next()){
        if(rd.code=="0"){ rd.push(); break; }
        if(rd.code=="10"){ double x=std::stod(rd.value); if(expect(rd,"20",seg,"10")) { double y=std::stod(rd.value); pts.push_back({x,y}); } }
    }
    if(pts.empty()){ std::cerr<<"seg "<<seg<<" SPLINE no fit points\n"; return false; }
    out=approxSpline(pts); return !out.empty();
}

// Parse LWPOLYLINE entity
static bool parseLWPolyline(DXFReader&rd,Path64&out,bool&closed,int seg){
    out.clear(); closed=false; bool ok=true;
    while(rd.next()){
        if(rd.code=="0"){ rd.push(); break; }
        if(rd.code=="70") closed=(std::stoi(rd.value)&1);
        else if(rd.code=="10"){ double x=std::stod(rd.value); if(expect(rd,"20",seg,"10")) { double y=std::stod(rd.value); out.push_back({I64(x),I64(y)}); } }
    }
    if(out.empty()){ std::cerr<<"seg "<<seg<<" LWPOLYLINE empty\n"; ok=false; }
    return ok && !out.empty();
}

// Parse classic POLYLINE entity (sequence of VERTEX)
static bool parsePolyline(DXFReader&rd,Path64&out,int seg){
    out.clear();
    while(rd.next()){
        if(rd.code=="0" && rd.value=="VERTEX"){
            double x=0,y=0; bool have=false;
            while(rd.next()){
                if(rd.code=="0"){ rd.push(); break; }
                if(rd.code=="10"){ x=std::stod(rd.value); if(expect(rd,"20",seg,"10")) y=std::stod(rd.value); have=true; }
            }
            if(have) out.push_back({I64(x),I64(y)});
        }else if(rd.code=="0" && rd.value=="SEQEND"){
            break;
        }else if(rd.code=="0"){
            rd.push();
            break;
        }
    }
    if(out.empty()){ std::cerr<<"seg "<<seg<<" POLYLINE empty\n"; return false; }
    return true;
}


// Connect line/arc segments into closed rings
static Paths64 connectSegments(std::vector<Path64> segs){
    Paths64 out;
    std::vector<char> used(segs.size(), 0);

    for(size_t i=0;i<segs.size();++i){
        if(used[i]) continue;
        Path64 cur = segs[i];
        used[i] = 1;
        bool extended = true;

        while(extended){
            extended = false;

            for(size_t j=0;j<segs.size();++j){
                if(used[j]) continue;
                auto &s = segs[j];

                // хвост cur  → голова s
                if(eqPt(cur.back(), s.front())){
                    cur.insert(cur.end(), s.begin()+1, s.end());
                    used[j]=1; extended=true; break;
                }
                // хвост cur  → хвост s (надо развернуть s)
                if(eqPt(cur.back(), s.back())){
                    cur.insert(cur.end(), s.rbegin()+1, s.rend());
                    used[j]=1; extended=true; break;
                }
                // голова cur → хвост s  (дописываем в начало)
                if(eqPt(cur.front(), s.back())){
                    cur.insert(cur.begin(), s.rbegin()+1, s.rend());
                    used[j]=1; extended=true; break;
                }
                // голова cur → голова s  (дописываем в начало, разворачивая s)
                if(eqPt(cur.front(), s.front())){
                    cur.insert(cur.begin(), s.begin()+1, s.end());
                    used[j]=1; extended=true; break;
                }
            }
        }

        // замкнулся?
        if(eqPt(cur.front(), cur.back()) && cur.size() > 3)
            out.push_back(std::move(cur));
    }
    return out;
}


// Parse DXF entities into polygons
static std::vector<RawPart> loadDXF(const std::vector<std::string>& files){
    std::vector<RawPart> parts;
    for(size_t idx=0; idx<files.size(); ++idx){
        std::ifstream fin(files[idx]);
        if(!fin) throw std::runtime_error("open "+files[idx]);

        DXFReader rd(fin);
        bool inEnt=false;
        std::vector<Path64> rings,segs;
        int segNo=0;

        while(rd.next()){
            if(rd.code=="0" && rd.value=="SECTION"){
                if(rd.next() && rd.code=="2") inEnt=(rd.value=="ENTITIES");
                continue;
            }
            if(!inEnt || rd.code!="0") continue;

            Path64 tmp; bool closed=false;
            if(rd.value=="LWPOLYLINE"){
                if(parseLWPolyline(rd,tmp,closed,segNo++))
                    (closed?rings:segs).push_back(std::move(tmp));
            }else if(rd.value=="LINE"){
                if(parseLine(rd,tmp,segNo++)) segs.push_back(std::move(tmp));
            }else if(rd.value=="ARC"){
                if(parseArc(rd,tmp,segNo++)) segs.push_back(std::move(tmp));
            }else if(rd.value=="CIRCLE"){
                if(parseCircle(rd,tmp,segNo++)) segs.push_back(std::move(tmp));
            }else if(rd.value=="ELLIPSE"){
                if(parseEllipse(rd,tmp,segNo++)) segs.push_back(std::move(tmp));
            }else if(rd.value=="SPLINE"){
                if(parseSpline(rd,tmp,segNo++)) segs.push_back(std::move(tmp));
            }else if(rd.value=="POLYLINE"){
                if(parsePolyline(rd,tmp,segNo++)) segs.push_back(std::move(tmp));
            }
        }

        auto joined=connectSegments(segs);
        std::cerr<<"DXF debug  "<<files[idx]<<"  segs="<<segs.size() <<"  joined="<<joined.size()<<"  rings="<<rings.size()<<"\n";
        rings.insert(rings.end(),joined.begin(),joined.end());
        if(rings.empty()) throw std::runtime_error("no rings in "+files[idx]);
        std::sort(rings.begin(),rings.end(),[](const Path64&a,const Path64&b){return std::abs(Area(a))>std::abs(Area(b));});
        if(Area(rings[0])<0) ReversePath(rings[0]);
        for(size_t i=1;i<rings.size();++i) if(Area(rings[i])>0) ReversePath(rings[i]);
        parts.push_back({rings,std::abs(Area(rings[0]))/(SCALE*SCALE),int(idx)});
    }
    return parts;
}

// ───────── orientations ─────────
struct Orient{ Paths64 poly; Rect64 bb; int id; int ang; };

// Rotate polygon set using sine and cosine
static Paths64 rot(const Paths64&src,double s,double c){
    Paths64 out; out.reserve(src.size());
    for(const auto&ring:src){ Path64 r; r.reserve(ring.size());
        for(auto pt:ring){ double x=Dbl(pt.x), y=Dbl(pt.y); r.push_back({I64(x*c - y*s), I64(x*s + y*c)}); }
        out.push_back(std::move(r));
    }
    return out;
}

// Precompute rotated variants for each part
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
static Key kfn(int a, int ra, int b, int rb){
    return ((uint64_t)a<<48)^((uint64_t)ra<<32)^((uint64_t)b<<16)^rb;
}

// Calculate and cache no-fit polygon (NFP) for oriented parts
static const Paths64& nfp(const Orient& A, const Orient& B){
    Key k = kfn(A.id, A.ang, B.id, B.ang);
    auto it = NFP.find(k);
    if (it != NFP.end())
        return it->second;
    return NFP[k] = MinkowskiSum(A.poly[0], B.poly[0], true);
}

// ───────── greedy placer ─────────
// Place a set of parts onto the sheet using a simple greedy algorithm
struct Place{int id; double x,y; int ang;};

static std::vector<Place> greedy(const std::vector<RawPart>&parts,double W,double H,int step){
    auto ord=parts; std::sort(ord.begin(),ord.end(),[](auto&a,auto&b){return a.area>b.area;});
    auto orient=makeOrient(ord,step);
    Path64 sheet;
    sheet.reserve(4);
    sheet.emplace_back(int64_t(0), int64_t(0));
    sheet.emplace_back(I64(W), int64_t(0));
    sheet.emplace_back(I64(W), I64(H));
    sheet.emplace_back(int64_t(0), I64(H));
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
struct CLI{
    double W = 0, H = 0;
    int rot = DEFAULT_ROT_STEP;
    std::string out = "layout.csv";
    std::vector<std::string> files;
};

// Display usage information
static void printHelp(const char* exe){
    std::cout << "Usage: " << exe << " -s WxH [options] files...\n"
              << "Options:\n"
              << "  -s, --sheet WxH    Sheet size in mm\n"
              << "  -r, --rot N       Rotation step in degrees (default 10)\n"
              << "  -o, --out FILE    Output CSV file (default layout.csv)\n"
              << "  -h, --help        Show this help message\n";
}

// Parse command line arguments
static CLI parse(int ac, char** av){
    CLI c;
    for(int i=1;i<ac;++i){
        std::string a = av[i];
        if(a=="-h"||a=="--help"){ printHelp(av[0]); std::exit(0); }
        else if(a=="-s"||a=="--sheet"){
            if(i+1>=ac) throw std::runtime_error("missing value for --sheet");
            std::string s = av[++i]; auto x = s.find('x');
            if(x==std::string::npos) throw std::runtime_error("sheet format WxH");
            c.W = std::stod(s.substr(0,x));
            c.H = std::stod(s.substr(x+1));
        }else if(a=="-r"||a=="--rot"){
            if(i+1>=ac) throw std::runtime_error("missing value for --rot");
            c.rot = std::stoi(av[++i]);
        }else if(a=="-o"||a=="--out"){
            if(i+1>=ac) throw std::runtime_error("missing value for --out");
            c.out = av[++i];
        }else if(a.size() && a[0]=='-'){
            throw std::runtime_error("unknown option " + a);
        }else{
            c.files.push_back(a);
        }
    }
    if(c.W<=0 || c.H<=0 || c.files.empty())
        throw std::runtime_error("use --help for usage");
    return c;
}

// ───────── main ─────────
// Entry point
int main(int argc,char*argv[]){
    try{
        auto cli=parse(argc,argv);
        auto parts=loadDXF(cli.files);
        auto placed=greedy(parts,cli.W,cli.H,cli.rot);
        std::ofstream f(cli.out); f<<"part,x_mm,y_mm,angle\n"; for(auto&p:placed) f<<p.id<<","<<std::fixed<<std::setprecision(3)<<p.x<<","<<p.y<<","<<p.ang<<"\n";
        std::cout<<"placed "<<placed.size()<<"/"<<parts.size()<<" → "<<cli.out<<"\n";
        return 0;
    }catch(const std::exception&e){ std::cerr<<"ERR: "<<e.what()<<"\n"; return 1; }}