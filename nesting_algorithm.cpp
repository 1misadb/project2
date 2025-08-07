// nesting_clip2.cpp — NFP‑greedy nesting (актуальный Clipper2 v2.x)
// --------------------------------------------------------------------
//   • DXF‑парсер  (минимум LWPOLYLINE) → Paths64 (+ дырки)
//   • 36 ориентаций (0…350° шаг 10°)
//   • NFP через Clipper2::MinkowskiSum (clipper.minkowski.h)
//   • Greedy corner‑search, размещение в отверстия
//   • CLI:  nest -s 2000x2000 -r 10 *.dxf -o layout.csv
//   • Header‑only зависимость: Clipper2 (https://github.com/AngusJohnson/Clipper2)
// --------------------------------------------------------------------
//  Build (g++ / MSVC):
//      cl /EHsc /std:c++17 /O2 /openmp clipper3\src\clipper.engine.cpp clipper3\src\clipper.offset.cpp clipper3\src\clipper.rectclip.cpp nesting_algorithm.cpp /Iclipper3\include /I. /Fe:nest.exe
// --------------------------------------------------------------------
// для нестинга полная для раскладки через питон парсер полная схема лежит в редми используюет джейсон пакет предпологает что координаты будут лежать в массиве 
// в виде poliny так как перебирать багованые координаты через C++ это дело не состоятельное и лучшен использовать готовую провереную итоговую питон библиотеку  
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cctype>
#include <map>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cstdint>
#include <array>
#include <deque>
#include <tuple>
#include <random>
#include <mutex>
#include "cxxopts.hpp"
#include <chrono>
#include <atomic>
#include <thread>
#include <set>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif
#define _USE_MATH_DEFINES
#include <math.h>
#include "clipper3/clipper.h"
#include "clipper3/clipper.minkowski.h"
#include <json.hpp>
#include <spdlog/spdlog.h>
#if defined(__clang__)
# define GUARDED_BY(x) __attribute__((guarded_by(x)))
#else
# define GUARDED_BY(x)
#endif
#include "lru_cache.h"
double gUnit = 1.0;
bool gVerbose = false;
extern double gKerf;
extern double gGap;
using namespace Clipper2Lib;

static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

// ────────── types ──────────
// ───────── helpers ─────────
// Rect64: 64-bit integer rectangle
// Масштаб преобразования мм → int64 (Clipper scale)
double getDxfUnitFactor(int code) {
    switch (code) {
        case 0: return 1.0;       // Не указано, по умолчанию 1
        case 1: return 25.4;      // Дюймы → мм
        case 2: return 25.4 * 12; // Футы → мм
        case 3: return 25.4 * 12 * 5280; // Мили → мм
        case 4: return 1.0;       // Миллиметры
        case 5: return 10.0;      // Сантиметры → мм
        case 6: return 1000.0;    // Метры → мм
        case 7: return 1000000.0; // Километры → мм
        case 8: return 0.001;     // Микроны → мм
        case 9: return 0.0254;    // Милли (mil) → мм
        default: return 1.0;
    }
}

int detectDxfUnits(const std::string& dxf_file) {
    std::ifstream fin(dxf_file);
    if (!fin) return 0;
    std::string code, value;
    while (std::getline(fin, code) && std::getline(fin, value)) {
        if (trim(code) == "9" && trim(value) == "$INSUNITS") {
            if (std::getline(fin, code) && std::getline(fin, value) && trim(code) == "70") {
                std::string val = trim(value);
                try { return std::stoi(val); }
                catch (...) { return 0; }
            }
        }
        if (trim(code) == "0" && trim(value) == "ENDSEC") break;
    }
    return 0;
}
using Key = uint64_t;
static Key kfn(int a, int ra, int b, int rb){
    return ((uint64_t)a<<48) ^ ((uint64_t)ra<<32) ^ ((uint64_t)b<<16) ^ (uint64_t)rb;
}
// --- Overlap cache key structure ---
struct OverlapKey {
    int part_a, ang_a;
    int64_t x_a, y_a;
    int part_b, ang_b;
    int64_t x_b, y_b;
    bool operator==(const OverlapKey& o) const {
        return part_a==o.part_a && ang_a==o.ang_a && x_a==o.x_a && y_a==o.y_a &&
               part_b==o.part_b && ang_b==o.ang_b && x_b==o.x_b && y_b==o.y_b;
    }
};
static inline OverlapKey makeKey(int pa,int aa,int64_t xa,int64_t ya,
                                 int pb,int ab,int64_t xb,int64_t yb){
    OverlapKey k{pa,aa,xa,ya,pb,ab,xb,yb};
    if(std::tie(k.part_a,k.ang_a,k.x_a,k.y_a) >
       std::tie(k.part_b,k.ang_b,k.x_b,k.y_b)){
        std::swap(k.part_a,k.part_b);
        std::swap(k.ang_a,k.ang_b);
        std::swap(k.x_a,k.x_b);
        std::swap(k.y_a,k.y_b);
    }
    return k;
}
namespace std {
    template<>
    struct hash<OverlapKey> {
        size_t operator()(const OverlapKey& k) const {
            size_t h = hash<int>()(k.part_a) ^ (hash<int>()(k.ang_a)<<1);
            h ^= (hash<int64_t>()(k.x_a)<<2) ^ (hash<int64_t>()(k.y_a)<<3);
            h ^= (hash<int>()(k.part_b)<<4) ^ (hash<int>()(k.ang_b)<<5);
            h ^= (hash<int64_t>()(k.x_b)<<6) ^ (hash<int64_t>()(k.y_b)<<7);
            return h;
        }
    };
}

// --- Global overlap cache ---
#include <tbb/concurrent_unordered_map.h>
static tbb::concurrent_unordered_map<OverlapKey, bool> overlapCache GUARDED_BY(overlapMutex);
static std::deque<OverlapKey> overlapOrder GUARDED_BY(overlapMutex);
static std::mutex overlapMutex;
constexpr size_t OVERLAP_CACHE_MAX = 20000;
static LRUCache<Key, Paths64> sharedNFP(20000);
std::mutex output_mutex;
constexpr double SCALE = 1e4;
constexpr int DEFAULT_ROT_STEP = 10;
constexpr double TOL_MM = 1.0;       
// ───────── CLI ─────────
struct CLI{
    double W = 0, H = 0;
    int rot = DEFAULT_ROT_STEP;
    double grid = 5;
    std::string out = "layout.csv";
    std::string dxf = "layout.dxf"; 
    std::vector<std::string> files;
    std::vector<int> nums;
    int num = 1;
    bool joinSegments = false;
    int iter = 1;
    int runs = 1;
    int generations = 50;
    int pop_size = 50;
    int polish = 0;
    std::string strategy = "area";
    int cand_limit = 2000;
    int nfp_limit = 20000;
    int overlap_limit = 20000;
    double time_limit = 3000.0;
    bool verbose = false;
    int run = 0;            // run index for multi strategy runs
    bool fill_gaps = false; // enable gap nesting
};
// base tolerance in mm
inline int64_t tol_mm() { return llround(TOL_MM * SCALE); }
// Преобразование double → int64 с учетом SCALE
inline int64_t I64(double v) {          // из мм → int64
    return llround(v * gUnit * SCALE);  // ← добавили gUnit
}

inline double Dbl(int64_t v) {          // из int64 → мм
    return double(v) / SCALE; // ← делим тоже на gUnit
}
// Reverse polygon orientation in place
template <typename T>
void ReversePath(std::vector<T>& p) {
    std::reverse(p.begin(), p.end());
}
// Compute bounding box of multiple polygons

static Rect64 getBBox(const Paths64& P){
    if (P.empty() || P[0].empty())
        return Rect64{0,0,0,0};
    auto pt0 = P[0][0];
    Rect64 r{pt0.x, pt0.y, pt0.x, pt0.y};
    for(const auto&path:P)
        for(auto pt:path){
            r.left = std::min(r.left, pt.x);
            r.bottom = std::min(r.bottom, pt.y);
            r.right = std::max(r.right, pt.x);
            r.top = std::max(r.top, pt.y);
        }
    return r;
}

// Compute bounding box of a single polygon
static Rect64 getBBox(const Path64& p)
{
    if(p.empty()) return Rect64{0,0,0,0};
    Rect64 r{p[0].x,p[0].y,p[0].x,p[0].y};
    for(const auto&pt : p){
        r.left   = std::min(r.left, pt.x);
        r.bottom = std::min(r.bottom, pt.y);
        r.right  = std::max(r.right, pt.x);
        r.top    = std::max(r.top, pt.y);
    }
    return r;
}

// Path length in mm
static double pathLength(const Path64& p)
{
    if(p.size() < 2) return 0.0;
    double len = 0.0;
    for(size_t i=1;i<p.size();++i)
        len += std::hypot(double(p[i].x - p[i-1].x),
                          double(p[i].y - p[i-1].y));
    return len / SCALE;
}

// Union all rings to simplify part geometry
// Merge all polygons in a set
static Paths64 unionAll(const Paths64& in){
    return Union(in, FillRule::NonZero);
}
static std::string rectToStr(const Rect64& r)
{
    std::ostringstream ss;
    ss << "(" << Dbl(r.left) << "," << Dbl(r.bottom) << ")-("
       << Dbl(r.right) << "," << Dbl(r.top) << ")";
    return ss.str();
}

static std::string debugPaths64(const Paths64& p)
{
    std::ostringstream ss;
    ss << "{";
    for(size_t i=0;i<p.size();++i){
        if(i) ss << ",";
        ss << "[";
        for(size_t j=0;j<p[i].size();++j){
            if(j) ss << ",";
            ss << "(" << p[i][j].x << "," << p[i][j].y << ")";
        }
        ss << "]";
    }
    ss << "}";
    return ss.str();
}

static inline double areaMM2(const Paths64& p)
{
    return std::abs(Area(p)) / (SCALE * SCALE);
}

#include "geometry.h"

// ───────── DXF mini loader ─────────
struct RawPart{
    Paths64 rings;
    double area = 0;
    double perimeter = 0;
    int id = 0;
    Paths64 holes;
};

// New raw entity representation from DXF
struct RawEntity{
    std::string type;                                             // DXF entity type
    std::vector<std::pair<double,double>> pts;                    // polyline points in mm
    std::vector<std::vector<std::pair<double,double>>> extra;     // reserved for blocks etc
};

// Fallback spline → polyline conversion (simply connect control points)
static std::vector<std::pair<double,double>>
splineToPolylineFallback(const std::vector<std::pair<double,double>>& ctrlPts)
{
    return ctrlPts;
}

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

// Point comparison within tolerance
static bool eqPt(const Point64& a, const Point64& b)
{
    int64_t dx = a.x - b.x;
    int64_t dy = a.y - b.y;
    return (dx * dx + dy * dy) <= tol_mm() * tol_mm();
}

// Check if a path is closed without modifying it 
static bool isClosed(const Path64& path) 
{ 
    return path.size() > 1 && eqPt(path.front(), path.back()); 
}

// Euclidean distance between two points
inline int64_t distance(const Point64& a, const Point64& b)
{
    double dx = double(a.x) - double(b.x);
    double dy = double(a.y) - double(b.y);
    return static_cast<int64_t>(std::sqrt(dx*dx + dy*dy));
}

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
static bool parseLine(DXFReader& rd, Path64& out, int seg) {
    double x1=0, y1=0, x2=0, y2=0;
    bool have_x1 = false, have_y1 = false, have_x2 = false, have_y2 = false;
    out.clear();

    while (rd.next()) {
        if (rd.code == "0") { rd.push(); break; }

        try {
            if (rd.code == "10") { x1 = std::stod(rd.value); have_x1 = true; }
            else if (rd.code == "20") { y1 = std::stod(rd.value); have_y1 = true; }
            else if (rd.code == "11") { x2 = std::stod(rd.value); have_x2 = true; }
            else if (rd.code == "21") { y2 = std::stod(rd.value); have_y2 = true; }
        } catch(const std::exception& e) {
            std::cerr << "seg " << seg << " LINE parse error: " << e.what() << "\n";
            return false;
        }
        // остальные коды просто игнорируем!
    }

    if (!have_x1 || !have_y1 || !have_x2 || !have_y2) {
        std::cerr << "seg " << seg << " LINE missing coordinate(s): "
                  << (!have_x1 ? "x1 " : "") << (!have_y1 ? "y1 " : "")
                  << (!have_x2 ? "x2 " : "") << (!have_y2 ? "y2 " : "") << "\n";
        return false;
    }

    out.push_back({I64(x1), I64(y1)});
    out.push_back({I64(x2), I64(y2)});
    return true;
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
// Функция разбора SPLINE из DXF
static bool parseSpline(DXFReader& rd, Path64& out, int segNo)
{
    struct Vec2 { double x, y; };

    // Де Бура для B-spline
    auto deBoor = [](int k, int degree, double t,
                     const std::vector<double>& knots,
                     const std::vector<Vec2>& ctrl)
    {
        std::vector<Vec2> d;
        d.reserve(degree + 1);
        for (int j = 0; j <= degree; ++j)
            d.push_back(ctrl[k - degree + j]);

        for (int r = 1; r <= degree; ++r)
            for (int j = degree; j >= r; --j)
            {
                int i = k - degree + j;
                double den = knots[i + degree - r + 1] - knots[i];
                double alpha;
                if (den == 0) alpha = 0;
                else          alpha = (t - knots[i]) / den;
                d[j].x = (1 - alpha) * d[j - 1].x + alpha * d[j].x;
                d[j].y = (1 - alpha) * d[j - 1].y + alpha * d[j].y;
            }
        return d[degree];
    };

    // Аппроксимация по B-spline (samples ~ 100)
    auto approximateSpline = [&](const std::vector<Vec2>& ctrl,
                                 const std::vector<double>& knots,
                                 int degree, int samples)
    {
        std::vector<Point64> poly;
        int n = (int)ctrl.size() - 1;
        double t0 = knots[degree];
        double t1 = knots[n + 1];
        for (int s = 0; s < samples; ++s)
        {
            double t = t0 + (t1 - t0) * s / (samples - 1);
            int k = degree;
            while (k < n + 1 && !(t >= knots[k] && t < knots[k + 1])) ++k;
            if (k == n + 1) k = n;

            Vec2 p = deBoor(k, degree, t, knots, ctrl);
            if (std::isnan(p.x) || std::isnan(p.y)) continue;
            poly.push_back(Point64(I64(p.x), I64(p.y)));
        }
        return poly;
    };

    // --- Чтение DXF полей ---
    std::vector<std::array<double,3>> fitPts, ctrlPts;
    std::vector<double> knotVec;
    out.clear();

    bool closed = false;
    int degree = 3;

    double fx=0, fy=0, fz=0; bool fxOk=false, fyOk=false, fzOk=false;
    double cx=0, cy=0, cz=0; bool cxOk=false, cyOk=false, czOk=false;

    while (rd.next())
    {
        if (rd.code[0]=='0') { rd.push(); break; }
        int code = std::stoi(rd.code);
        switch (code)
        {
        // fit-пойнты
        case 11: fx = std::stod(rd.value); fxOk = true; if(fxOk&&fyOk&&fzOk) { fitPts.push_back({fx,fy,fz}); fxOk=fyOk=fzOk=false; } break;
        case 21: fy = std::stod(rd.value); fyOk = true; if(fxOk&&fyOk&&fzOk) { fitPts.push_back({fx,fy,fz}); fxOk=fyOk=fzOk=false; } break;
        case 31: fz = std::stod(rd.value); fzOk = true; if(fxOk&&fyOk&&fzOk) { fitPts.push_back({fx,fy,fz}); fxOk=fyOk=fzOk=false; } break;
        // control-points
        case 10: cx = std::stod(rd.value); cxOk = true; break;
        case 20: cy = std::stod(rd.value); cyOk = true; break;
        case 30:
            cz = std::stod(rd.value); czOk = true;
            if(cxOk&&cyOk&&czOk) { ctrlPts.push_back({cx,cy,cz}); cxOk=cyOk=czOk=false; }
            break;
        // knots
        case 40: knotVec.push_back(std::stod(rd.value)); break;
        // flags etc
        // DXF SPLINE flag bit 0 == 1 --> closed spline
        // earlier code erroneously checked bit 3 (planar)
        case 70: closed = (std::stoi(rd.value) & 1); break;
        case 71: degree = std::stoi(rd.value); break;
        default: break;
        }
    }

    // --- Формируем полилинию ---
    if (!fitPts.empty()) {
        for (auto &p : fitPts)
            if (!std::isnan(p[0]) && !std::isnan(p[1]))
                out.push_back(Point64(I64(p[0]), I64(p[1])));
    }
    else if (!ctrlPts.empty() && knotVec.size() >= ctrlPts.size() + degree + 1) {
        std::vector<Vec2> ctrlVec;
        for (auto& p : ctrlPts)
            ctrlVec.push_back({p[0], p[1]});
        auto poly = approximateSpline(ctrlVec, knotVec, degree, 1000);
        out.insert(out.end(), poly.begin(), poly.end());
    }
    if(gVerbose){
        std::cerr << "ctrlPts.size()=" << ctrlPts.size() << ", knotVec.size()=" << knotVec.size() << ", degree=" << degree << '\n';
        std::cerr << "SPLINE closed flag=" << closed << "   isClosed(path)=" << isClosed(out) << '\n';
    }
    // --- Замыкание ---
    if (!out.empty()) {
        if (closed) {
            int64_t dist = distance(out.front(), out.back());
            if (dist > tol_mm()) {
                if(gVerbose)
                    std::cerr << "[FORCE] spline " << segNo
                              << " forcibly closed, big gap=" << Dbl(dist) << " mm\n";
                out.push_back(out.front());
            } else if (dist <= tol_mm() && out.front() != out.back()) {
                out.push_back(out.front());
                if(gVerbose)
                    std::cerr << "[AutoClose] spline " << segNo
                              << " closed (distance=" << Dbl(dist) << " mm)\n";
            } else if (dist > tol_mm()) {
                if(gVerbose)
                    std::cerr << "[Warn] spline " << segNo
                              << " not closed, gap=" << Dbl(dist) << " mm\n";
            }
        }
    }

    // --- Убираем мусорные (0,0) точки на концах ---
    while (!out.empty() && out.back().x == 0 && out.back().y == 0)
        out.pop_back();
    if (!out.empty() && out.front().x == 0 && out.front().y == 0)
        out.erase(out.begin());

    // --- Debug ---
    if(gVerbose)
        std::cerr << " flags: closed="<<closed<<"  degree="<<degree
                  << "  fit="<<fitPts.size()<<"  ctrl="<<ctrlPts.size()
                  << "  knots="<<knotVec.size() << '\n';

    if(out.empty() && gVerbose)
        std::cerr << "  !! spline #" << segNo << " → no points parsed\n";
    return !out.empty();
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
static Paths64 connectSegments(const std::vector<Path64>& segs)
{
    using Point = Point64;
    const int64_t TOL = tol_mm();

    auto key = [&](const Point& p)
    {
        return std::make_pair((p.x + TOL/2)/TOL, (p.y + TOL/2)/TOL);
    };

    std::multimap<std::pair<int64_t,int64_t>, size_t> ends;
    for(size_t i=0;i<segs.size();++i)
    {
        if(segs[i].empty()) continue;
        ends.emplace(key(segs[i].front()), i);
        ends.emplace(key(segs[i].back()),  i);
    }

    std::vector<char> used(segs.size(),0);
    Paths64 out;
    size_t closedCnt = 0;

    for(size_t i=0;i<segs.size();++i)
    {
        if(used[i] || segs[i].empty()) continue;
        Path64 path = segs[i];
        used[i] = 1;
        bool extended = true;

        while(extended)
        {
            extended = false;
            auto r = ends.equal_range(key(path.back()));
            for(auto it=r.first; it!=r.second; ++it)
            {
                size_t j = it->second;
                if(used[j] || j==i) continue;
                if(eqPt(path.back(), segs[j].front()))
                {
                    path.insert(path.end(), segs[j].begin()+1, segs[j].end());
                    used[j]=1; extended=true; break;
                }
                if(eqPt(path.back(), segs[j].back()))
                {
                    path.insert(path.end(), segs[j].rbegin()+1, segs[j].rend());
                    used[j]=1; extended=true; break;
                }
            }

            r = ends.equal_range(key(path.front()));
            for(auto it=r.first; it!=r.second; ++it)
            {
                size_t j = it->second;
                if(used[j] || j==i) continue;
                if(eqPt(path.front(), segs[j].back()))
                {
                    path.insert(path.begin(), segs[j].rbegin()+1, segs[j].rend());
                    used[j]=1; extended=true; break;
                }
                if(eqPt(path.front(), segs[j].front()))
                {
                    path.insert(path.begin(), segs[j].begin()+1, segs[j].end());
                    used[j]=1; extended=true; break;
                }
            }
        }

        if(path.size()>2 && !eqPt(path.front(), path.back()))
        {
            int64_t dist = distance(path.front(), path.back());
            if(dist <= tol_mm())
            {
                path.push_back(path.front());
                if(gVerbose)
                    std::cerr << "[AutoClose] path " << out.size()
                              << " closed (distance=" << Dbl(dist) << " mm)\n";
            }
            else
            {
                if(gVerbose)
                    std::cerr << "[Warn] open path " << out.size()
                              << " gap=" << Dbl(dist) << " mm\n";
            }
        }

        bool closed = !path.empty() && eqPt(path.front(), path.back());
        if(closed) closedCnt++;

        if(path.size()>1)
        {
            out.push_back(std::move(path));
            Rect64 bb = getBBox(out.back());
            double len = pathLength(out.back());
            if(gVerbose)
                std::cerr << " [path " << out.size()-1 << "] len=" << len << " mm bbox="
                          << Dbl(bb.right-bb.left) << "x" << Dbl(bb.top-bb.bottom)
                          << " closed=" << (closed?"yes":"no") << "\n";
        }
    }

    std::map<std::pair<int64_t,int64_t>, int> ptcnt;
    for(size_t i=0;i<segs.size();++i)
    {
        if(segs[i].empty()) continue;
        ptcnt[key(segs[i].front())]++;
        ptcnt[key(segs[i].back())]++;
    }
    if(gVerbose){
        std::cerr << "Dangling ends: ";
        for (auto& p : ptcnt)
            if (p.second % 2 != 0)
                std::cerr << "(" << Dbl(p.first.first*TOL) << "," << Dbl(p.first.second*TOL) << ") ";
        std::cerr << "\n";
        std::cerr << "Segments=" << segs.size() << " Paths=" << out.size()
                  << " Closed=" << closedCnt << "\n";
    }
    return out;
}

// ---------------------------------------------------------------------------
// Minimal DXF entity parser returning RawEntity list
static std::vector<RawEntity> parseEntities(std::ifstream& fin){
    std::vector<RawEntity> out;
    DXFReader rd(fin);
    bool inEnt = false;
    while(rd.next()){
        if(trim(rd.code)=="0" && rd.value=="SECTION"){
            if(rd.next() && trim(rd.code)=="2")
                inEnt = (rd.value=="ENTITIES");
            continue;
        }
        if(trim(rd.code)=="0" && rd.value=="ENDSEC"){ inEnt=false; continue; }
        if(!inEnt || trim(rd.code)!="0") continue;

        std::string etype = rd.value;
        RawEntity ent; ent.type = etype; Path64 tmp; bool closed=false;
        try{
            if(etype=="LWPOLYLINE"){
                if(parseLWPolyline(rd,tmp,closed,0))
                    for(auto&p:tmp) ent.pts.push_back({Dbl(p.x),Dbl(p.y)});
            }else if(etype=="LINE"){
                if(parseLine(rd,tmp,0))
                    for(auto&p:tmp) ent.pts.push_back({Dbl(p.x),Dbl(p.y)});
            }else if(etype=="ARC"){
                if(parseArc(rd,tmp,0))
                    for(auto&p:tmp) ent.pts.push_back({Dbl(p.x),Dbl(p.y)});
            }else if(etype=="CIRCLE"){
                if(parseCircle(rd,tmp,0))
                    for(auto&p:tmp) ent.pts.push_back({Dbl(p.x),Dbl(p.y)});
            }else if(etype=="ELLIPSE"){
                if(parseEllipse(rd,tmp,0))
                    for(auto&p:tmp) ent.pts.push_back({Dbl(p.x),Dbl(p.y)});
            }else if(etype=="SPLINE"){
                if(parseSpline(rd,tmp,0))
                    for(auto&p:tmp) ent.pts.push_back({Dbl(p.x),Dbl(p.y)});
            }else if(etype=="POLYLINE"){
                if(parsePolyline(rd,tmp,0))
                    for(auto&p:tmp) ent.pts.push_back({Dbl(p.x),Dbl(p.y)});
            }else if(etype=="INSERT"){
                std::cerr<<"[WARN] INSERT entity encountered, not expanded.\n";
            }
        }catch(...){
            std::cerr << "[WARN] Exception parsing "<<etype<<"\n";
        }
        if(!ent.pts.empty()) out.push_back(std::move(ent));
    }
    return out;
}

// Convert RawEntity → RawPart (single closed polyline)
// ───────── helpers ──────────────────────────────────────────────────
//  int64 ← мм  (SCALE уже учтён, gUnit - нет, т.к. RawEntity.pts = мм)
static inline std::int64_t I64mm(double mm)
{
    return static_cast<std::int64_t>( std::llround(mm * SCALE) );
}

// ───────── RawEntity → RawPart ──────────────────────────────────────
static RawPart entityToPart(const RawEntity& e)
{
    RawPart part;

    Path64 outer; outer.reserve(e.pts.size());
    for (auto [x,y] : e.pts)             // mm → int64 (без gUnit!)
        outer.push_back({ llround(x * SCALE), llround(y * SCALE) });

    // bbox only for outer ring
    Rect64 bb = getBBox(outer);
    const int64_t dx = -bb.left, dy = -bb.bottom;

    auto apply_shift = [&](Path64& p){
        if (dx || dy)
            for (auto& pt : p){ pt.x += dx; pt.y += dy; }
    };

    apply_shift(outer);
    if (outer.size() > 1 && outer.front() != outer.back())
        outer.push_back(outer.front());
    part.rings.push_back(std::move(outer));

    // convert holes if any (not shifted individually)
    for (const auto& h : e.extra){
        Path64 hole; hole.reserve(h.size());
        for (auto [x,y] : h)
            hole.push_back({ llround(x * SCALE), llround(y * SCALE) });
        apply_shift(hole);
        if (hole.size() > 1 && hole.front() != hole.back())
            hole.push_back(hole.front());
        part.holes.push_back(std::move(hole));
    }

    part.area = std::abs(Area(part.rings[0])) / (SCALE * SCALE);
    part.id = 0;
    return part;
}

static std::vector<RawPart> entitiesToParts(const std::vector<RawEntity>& ents)
{
    std::vector<RawPart> parts;
    parts.reserve(ents.size());

    for (const auto& e : ents)
        parts.emplace_back( entityToPart(e) );

    return parts;
}

// ---------------------------------------------------------------
// Load parts from JSON file produced by extract_polylines.py
static std::vector<RawPart> loadPartsFromJson(const std::string& filename)
{
    spdlog::info("[JSON] parsing {}", filename);
    std::ifstream fin(filename);
    if(!fin){
        spdlog::warn("[JSON] cannot open {}", filename);
        return {};
    }

    std::stringstream buf; buf << fin.rdbuf();
    std::string data = buf.str();
    auto pos = data.find_first_not_of(" \t\r\n");
    if(pos == std::string::npos){
        spdlog::warn("[JSON] empty file {}", filename);
        return {};
    }
    char first = data[pos];
    if(first != '[' && first != '{'){
        spdlog::warn("[JSON] skipping non JSON {}", filename);
        return {};
    }

    nlohmann::json j;
    try{
        j = nlohmann::json::parse(data);
    }catch(const std::exception& e){
        spdlog::error("[JSON] parse {}: {}", filename, e.what());
        return {};
    }
    if(!j.is_array()){
        spdlog::error("[JSON] invalid JSON format in {}", filename);
        return {};
    }

    auto parse_one = [&](const nlohmann::json& arr)->RawPart {
        RawPart p;
        bool first = true;
        if(!arr.is_array()) return p;
        for(const auto& path : arr){
            if(!path.is_array()){
                spdlog::warn("[JSON] path is not array in {}", filename);
                continue;
            }
            Path64 ring;
            for(const auto& pt : path){
                if(!pt.is_array() || pt.size() < 2 || !pt[0].is_number() || !pt[1].is_number()){
                    spdlog::warn("[JSON] invalid point in {}", filename);
                    continue;
                }
                double x = pt[0].get<double>();
                double y = pt[1].get<double>();
                if(!std::isfinite(x) || !std::isfinite(y)){
                    spdlog::warn("[JSON] non finite coordinate in {}", filename);
                    continue;
                }
                ring.emplace_back(I64mm(x), I64mm(y));
            }
            if(ring.size() < 3) {
                spdlog::warn("[JSON] skipping degenerate ring in {}", filename);
                continue;
            }
            if(ring.front() != ring.back())
                ring.push_back(ring.front());
            if(first){
                p.rings.push_back(ring);
                first = false;
            }else{
                p.holes.push_back(ring);
            }
        }
        if(p.rings.empty())
            return p;

        Rect64 bb = getBBox(p.rings[0]);
        int64_t dx = -bb.left, dy = -bb.bottom;
        if(dx || dy){
            for(auto& pt : p.rings[0]){ pt.x += dx; pt.y += dy; }
            for(auto& hole : p.holes)
                for(auto& pt : hole){ pt.x += dx; pt.y += dy; }
        }

        p.area = std::abs(Area(p.rings[0])) / (SCALE * SCALE);
        p.perimeter = pathLength(p.rings[0]);
        return p;
    };

    bool multi = false;
    if(j.size() > 0 && j[0].is_array()){
        if(j[0].size() > 0 && j[0][0].is_array()){
            if(j[0][0].size() > 0 && j[0][0][0].is_array())
                multi = true;
        }
    }

    std::vector<RawPart> parts;
    if(multi){
        for(const auto& detail : j){
            RawPart part = parse_one(detail);
            if(!part.rings.empty()){ 
                part.id = static_cast<int>(parts.size());
                parts.push_back(std::move(part));
            }
        }
    }else{
        RawPart part = parse_one(j);
        if(!part.rings.empty()){
            part.id = 0;
            parts.push_back(std::move(part));
        }
    }
    spdlog::info("[JSON] loaded {} part(s) from {}", parts.size(), filename);
    return parts;
}



// ───────── orientations ─────────
struct Orient{ Paths64 poly; Rect64 bb; int id; int ang; };

// Rotate polygon set using sine and cosine
static Paths64 rot(const Paths64& src, double s, double c) {
    Paths64 out; out.reserve(src.size());
    for (const auto& ring : src) {
        Path64 r; r.reserve(ring.size());
        for (auto pt : ring) {
            double x = Dbl(pt.x);     // мм
            double y = Dbl(pt.y);     // мм
            // mm → int64   (без gUnit!)
            r.push_back({ I64mm(x * c - y * s),
                          I64mm(x * s + y * c) });
        }
        out.push_back(std::move(r));
    }
    return out;
}
// Precompute rotated variants for each part
static std::vector<std::vector<Orient>> makeOrient(const std::vector<RawPart>& parts, int step){
    std::vector<std::vector<Orient>> all(parts.size());
    for (const auto& p : parts) {
        std::vector<Orient> ovec;
        for (int ang = 0; ang < 360; ang += step) {
            double rad = ang * M_PI / 180.0;
            auto r = rot(p.rings, sin(rad), cos(rad));
            ovec.push_back({ unionAll(r), getBBox(r), p.id, ang });
        }
        all[p.id] = std::move(ovec);
    }
    return all;
}

// Precompute NFPs on GPU in batches for convex pairs
static void precomputeNFPGPU(const std::vector<std::vector<Orient>>& orients,
                             LRUCache<Key, Paths64>& cache) {
#ifdef USE_CUDA
    bool use_gpu = cuda_available();
    auto t0 = std::chrono::steady_clock::now();
    const size_t BATCH = NFP_BATCH_SIZE;
    std::vector<Path64> batchA, batchB;
    std::vector<Key> keys;
    size_t gpu_cnt = 0, cpu_cnt = 0;
    if (use_gpu) {
        for (const auto& oa : orients) {
            for (const auto& a : oa) {
                for (const auto& ob : orients) {
                    for (const auto& b : ob) {
                        Key k = kfn(a.id, a.ang, b.id, b.ang);
                        Paths64 tmp;
                        if (cache.get(k, tmp)) continue;
                        Path64 pa = SimplifyPath(a.poly[0], 0.05 * SCALE);
                        Path64 pb = SimplifyPath(b.poly[0], 0.05 * SCALE);
                        batchA.push_back(pa); batchB.push_back(pb); keys.push_back(k);
                        ++gpu_cnt;
                        if (batchA.size() >= BATCH) {
                            auto sols = minkowskiBatchGPU(batchA, batchB);
                            for (size_t i = 0; i < sols.size(); ++i)
                                cache.put(keys[i], sols[i]);
                            batchA.clear(); batchB.clear(); keys.clear();
                        }
                    }
                }
            }
        }
        if (!batchA.empty()) {
            auto sols = minkowskiBatchGPU(batchA, batchB);
            for (size_t i = 0; i < sols.size(); ++i)
                cache.put(keys[i], sols[i]);
        }
    } else {
        for (const auto& oa : orients) {
            for (const auto& a : oa) {
                for (const auto& ob : orients) {
                    for (const auto& b : ob) {
                        Key k = kfn(a.id, a.ang, b.id, b.ang);
                        Paths64 tmp;
                        if (cache.get(k, tmp)) continue;
                        Path64 pa = SimplifyPath(a.poly[0], 0.05 * SCALE);
                        Path64 pb = SimplifyPath(b.poly[0], 0.05 * SCALE);
                        auto res = MinkowskiSum(pa, pb, true);
                        cache.put(k, res);
                        ++cpu_cnt;
                    }
                }
            }
        }
    }
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                  std::chrono::steady_clock::now() - t0).count();
    std::cerr << "[NFP-GPU] batchGPU: " << gpu_cnt
              << " pairs, CPU fallback: " << cpu_cnt
              << " pairs, time: " << ms << " ms" << std::endl;
#else
    (void)orients; (void)cache;
#endif
}


// ───────── NFP cache ─────────
static const Paths64& nfp(
    const Orient& A, const Orient& B,
    std::unordered_map<Key,Paths64>& localNFP,
    LRUCache<Key,Paths64>& globalCache)
{
    Key k = kfn(A.id, A.ang, B.id, B.ang);
    auto it = localNFP.find(k);
    if(it != localNFP.end()) return it->second;

    // 1) try global cache
    Paths64 cached;
    if(globalCache.get(k, cached)) {
        localNFP[k] = cached;
        return localNFP[k];
    }

    // 2) Считаем NFP без блокировки
    Path64 a = SimplifyPath(A.poly[0], 0.05 * SCALE);
    Path64 b = SimplifyPath(B.poly[0], 0.05 * SCALE);
    if (std::abs(Area(a)) / (SCALE * SCALE) < 1.0) {
        Rect64 bb = getBBox(a);
        a = { {bb.left, bb.bottom}, {bb.right, bb.bottom}, {bb.right, bb.top}, {bb.left, bb.top}, {bb.left, bb.bottom} };
    }
    if (std::abs(Area(b)) / (SCALE * SCALE) < 1.0) {
        Rect64 bb = getBBox(b);
        b = { {bb.left, bb.bottom}, {bb.right, bb.bottom}, {bb.right, bb.top}, {bb.left, bb.top}, {bb.left, bb.bottom} };
    }
    Paths64 res = MinkowskiSum(a, b, true);
    // kerf compensation: expand NFP by 0.1mm
    const double kerf_mm = 0.1;
    if(kerf_mm > 0.0) {
        res = InflatePaths(res, kerf_mm * SCALE, JoinType::Miter, EndType::Polygon);
    }

    // 3) store in global cache
    globalCache.put(k, res);
    localNFP.emplace(k, res);
    return localNFP[k];
}


static std::vector<Point64> candidates(const std::vector<Paths64>& placed,
                                       const Paths64& current,
                                       int grid,
                                       const Point64& sheetBR) {
    constexpr size_t MAX_GRID_PTS = 100000;
    std::vector<Point64> out;

    for (const auto& shape : placed) {
        Rect64 bb = getBBox(shape);
        out.push_back({bb.left, bb.bottom});
        out.push_back({bb.right, bb.bottom});
        out.push_back({bb.left, bb.top});
        out.push_back({bb.right, bb.top});
    }


    int64_t step = static_cast<int64_t>(std::llround(grid * SCALE));
    if (step <= 0) step = SCALE; // safety
    for (int64_t y = 0; y <= sheetBR.y && out.size() < MAX_GRID_PTS; y += step) {
        for (int64_t x = 0; x <= sheetBR.x && out.size() < MAX_GRID_PTS; x += step) {
            out.push_back({x, y});
        }
    }
    if (out.size() >= MAX_GRID_PTS) {
        std::lock_guard<std::mutex> lock(output_mutex);
        if(gVerbose)
            std::cerr << "[WARN] grid candidate limit reached" << std::endl;
    }
    return out;
}
// ───────── greedy placer ─────────
// Place a set of parts onto the sheet using a simple greedy algorithm
struct Place{int id; double x,y; int ang;};

struct LayoutStats{ double bboxArea=0.0; double placedArea=0.0; };

static LayoutStats computeStats(const std::vector<Place>& layout, const std::vector<RawPart>& parts){
    LayoutStats st{}; if(layout.empty()) return st; Paths64 all; for(const auto& pl:layout){
        double rad=pl.ang*M_PI/180.0; auto r=rot(parts[pl.id].rings,sin(rad),cos(rad));
        Paths64 moved=movePaths(r,I64mm(pl.x),I64mm(pl.y)); all.insert(all.end(),moved.begin(),moved.end()); }
    Rect64 bb=getBBox(all); st.bboxArea=Dbl(bb.right-bb.left)*Dbl(bb.top-bb.bottom);
    auto uni=Union(all,FillRule::NonZero); st.placedArea=areaMM2(uni); return st; }

static double computeArea(const std::vector<Place>& layout,const std::vector<RawPart>& parts){
    return computeStats(layout, parts).bboxArea;
}
const int MAX_OVERLAP_CHECKS = 100;          // Максимум overlap-проверок по умолчанию
const double MAX_TIME_PER_PART = 0.5;        // Секунд на одну деталь по умолчанию

auto part_start = std::chrono::steady_clock::now();
int overlap_checks = 0;
// Функция раскладки для одной случайной перестановки (ОДНА итерация!)
static std::vector<Place> greedy(
    const std::vector<RawPart>& parts,
    const std::vector<std::vector<Orient>>& all_orients,
    double W, double H, int grid,
    unsigned int seed,
    std::unordered_map<Key, Paths64>& localNFP,
    LRUCache<Key, Paths64>& sharedNFP,
    const CLI* cli_ptr
) {
    size_t CAND_LIMIT = 1e6;
    size_t NFP_LIMIT = 50000;
    int MAX_OVERLAP_CHECKS = 100;
    double MAX_TIME_PER_PART = 0.5;
    if (cli_ptr) {
        CAND_LIMIT = cli_ptr->cand_limit;
        NFP_LIMIT = cli_ptr->nfp_limit;
        MAX_OVERLAP_CHECKS = cli_ptr->overlap_limit;
        MAX_TIME_PER_PART = cli_ptr->time_limit;
    }
    auto mm2i = [](double mm) { return static_cast<int64_t>(std::llround(mm * SCALE)); };
    Point64 sheetBR{ mm2i(W), mm2i(H) };

    std::vector<RawPart> ord = parts;
    if (seed != 0) {
        std::mt19937 rng(seed);
        std::shuffle(ord.begin(), ord.end(), rng);
    }

    std::vector<Paths64> placedShapes;
    std::vector<Orient> placedOrient;
    std::vector<Place> layout;
    size_t total_checks_gpu = 0;
    size_t total_checks_all = 0;

    for (size_t i = 0; i < ord.size(); ++i) {
        overlap_checks = 0;                      // ① обнулить перед деталью
        auto part_start = std::chrono::steady_clock::now();   // ② время – локально
        spdlog::info("Placing part {} / {} (ID {})", i+1, ord.size(), ord[i].id);
        if (ord[i].rings.empty()) {
            std::cerr << "[ERROR] Деталь " << i << " пуста! Пропускаю.\n";
            continue;
        }
        auto t_start = std::chrono::steady_clock::now();
        bool ok = false; Place best{}; Orient bestO{};
        const auto& oriSet = all_orients[ ord[i].id ];
        size_t invalidBB_total = 0;
        bool firstTooBig = (i == 0);

        // --- Параллелизация по ориентациям ---
        int bestOrientIdx = -1;
        double bestCost = 1e100;
        Place bestPl{};
        Orient bestOrient{};

#pragma omp parallel
        {   
            std::unordered_map<Key, Paths64> localNFP;
            int localBestOrientIdx = -1;
            double localBestCost = 1e100;
            Place localBestPl{};
            Orient localBestOrient{};
#pragma omp for nowait
            for (int opIdx = 0; opIdx < (int)oriSet.size(); ++opIdx) {
                const auto& op = oriSet[opIdx];
                if (op.bb.right <= sheetBR.x && op.bb.top <= sheetBR.y)
                    if(i==0) firstTooBig = false;
                if (op.bb.right > sheetBR.x || op.bb.top > sheetBR.y)
                    continue;

                const size_t CAND_LIMIT = 100000;   // Было 10000, увеличено для большего перебора
                const size_t NFP_LIMIT = 20000;     // Было 2000, увеличено для большего перебора
                std::vector<Point64> cand;
                cand.reserve(CAND_LIMIT);
                cand.push_back({0,0});
                for(const auto& ps : placedShapes){
                    if (cand.size() >= CAND_LIMIT) break;
                    Rect64 bb = getBBox(ps);
                    cand.push_back({bb.left, bb.bottom});
                    if (cand.size() >= CAND_LIMIT) break;
                    cand.push_back({bb.right, bb.bottom});
                    if (cand.size() >= CAND_LIMIT) break;
                    cand.push_back({bb.left, bb.top});
                    if (cand.size() >= CAND_LIMIT) break;
                    cand.push_back({bb.right, bb.top});
                }
                int64_t step = mm2i(grid);
                for(int64_t y=0; y<=sheetBR.y && cand.size()<CAND_LIMIT; y+=step) {
                    for(int64_t x=0; x<=sheetBR.x && cand.size()<CAND_LIMIT; x+=step) {
                        cand.push_back({x,y});
                    }
                }
                size_t nfp_count = 0;
                for(const auto& po : placedOrient) {
                    if (cand.size() >= CAND_LIMIT) break;
                    const Paths64& nfpp = nfp(po, op, localNFP, sharedNFP);
                    for(const auto& pth : nfpp) {
                        for(const auto& pt : pth) {
                            if (cand.size() >= CAND_LIMIT || nfp_count >= NFP_LIMIT) break;
                            cand.push_back(pt);
                            nfp_count++;
                        }
                        if (cand.size() >= CAND_LIMIT || nfp_count >= NFP_LIMIT) break;
                    }
                    if (cand.size() >= CAND_LIMIT) break;
                }
                // --- сортировка кандидатов по (y, x) для "лево-низ" ---
                if (cand.size() > 1) {
                    std::sort(cand.begin(), cand.end(), [](const Point64& a, const Point64& b) {
                        return a.y == b.y ? a.x < b.x : a.y < b.y;
                    });
                    cand.erase(std::unique(cand.begin(), cand.end(), [](const Point64& a, const Point64& b) {
                        return a.x == b.x && a.y == b.y;
                    }), cand.end());
                }
                if (cand.size() > CAND_LIMIT)
                    cand.resize(CAND_LIMIT);

                size_t valid = 0;
                for (size_t idx = 0; idx < cand.size(); ++idx) {
                    if (overlap_checks >= MAX_OVERLAP_CHECKS) break;
                    // исправлено: объявить now здесь
                    auto now = std::chrono::steady_clock::now();
                    double elapsed = std::chrono::duration<double>(now - part_start).count();
                    if (elapsed > MAX_TIME_PER_PART) break;
                    const auto& c = cand[idx];
                    Rect64 bb = op.bb; bb.left += c.x; bb.right += c.x; bb.bottom += c.y; bb.top += c.y;
                    if (bb.left < 0 || bb.bottom < 0 || bb.right > sheetBR.x || bb.top > sheetBR.y){
                        ++invalidBB_total;
                        continue;
                    }
                    Paths64 moved = movePaths(op.poly, c.x, c.y);
                    bool clash = false;
                    std::vector<Paths64> gpuBatch;
                    std::vector<OverlapKey> gpuKeys;
                    for (size_t pi = 0; pi < placedShapes.size(); ++pi) {
                        const auto& pl = placedShapes[pi];
                        Rect64 bbPl = getBBox(pl);
                        Rect64 bbMoved = getBBox(moved);
                        if (bbPl.right < bbMoved.left || bbPl.left > bbMoved.right ||
                            bbPl.top   < bbMoved.bottom || bbPl.bottom > bbMoved.top)
                            continue;
                        overlap_checks++;
                        total_checks_all++;
                        const auto& po = placedOrient[pi];
                        int64_t cx = I64mm(c.x), cy = I64mm(c.y);
                        int64_t xb = bbPl.left;
                        int64_t yb = bbPl.bottom;
                        OverlapKey key = makeKey(
                            op.id, op.ang, cx, cy,
                            po.id, po.ang, xb, yb);
                        bool found = false; bool cached = false;
                        {
                            std::lock_guard<std::mutex> lock(overlapMutex);
                            auto it = overlapCache.find(key);
                            if (it != overlapCache.end()) {
                                cached = it->second;
                                found = true;
                            }
                        }
                        if (found) {
                            if (cached) { clash = true; break; }
                            continue;
                        }
                        gpuBatch.push_back(pl);
                        gpuKeys.push_back(key);
                    }
                    if (!clash && !gpuBatch.empty()) {
                        std::vector<bool> res;
#ifdef USE_CUDA
                        if (cuda_available()) {
                            res = overlapBatchGPU(moved, gpuBatch);
                            total_checks_gpu += gpuBatch.size();
                        } else
#endif
                        {
                            res.resize(gpuBatch.size());
                            for (size_t bi = 0; bi < gpuBatch.size(); ++bi)
                                res[bi] = overlap(moved, gpuBatch[bi]);
                        }
                        for (size_t bi = 0; bi < res.size(); ++bi) {
                            bool ov = res[bi];
                            {
                                std::lock_guard<std::mutex> lock(overlapMutex);
                                overlapCache[gpuKeys[bi]] = ov;
                                overlapOrder.push_back(gpuKeys[bi]);
                                if (overlapOrder.size() > OVERLAP_CACHE_MAX) {
                                    OverlapKey old = overlapOrder.front();
                                    overlapOrder.pop_front();
                                    overlapCache.unsafe_erase(old);
                                }
                            }
                            if (ov) { clash = true; break; }
                        }
                    }
                    if (clash) continue;
                    ++valid;
                    double cost = Dbl(bb.top);
                    // --- Улучшение: если cost одинаковый, выбирать левее ---
                    if (cost < localBestCost ||
                        (std::abs(cost - localBestCost) < 1e-6 && c.x < localBestPl.x)) {
                        localBestCost   = cost;
                        localBestPl     = { op.id, Dbl(c.x), Dbl(c.y), op.ang };
                        localBestOrient = op;
                        localBestOrientIdx = opIdx;
                    }
                }
            }
            // --- Критическая секция для выбора глобального лучшего ---
#pragma omp critical
            {
                if (localBestCost < bestCost ||
                    (std::abs(localBestCost - bestCost) < 1e-6 && localBestPl.x < bestPl.x)) {
                    bestCost = localBestCost;
                    bestPl = localBestPl;
                    bestOrient = localBestOrient;
                    bestOrientIdx = localBestOrientIdx;
                }
            }
        }
        // --- после параллельного цикла ---
        if (bestCost < 1e100 && bestOrientIdx >= 0) {
            auto& op = oriSet[bestOrientIdx];
            auto& c = bestPl;
            Paths64 moved = movePaths(op.poly, I64mm(c.x), I64mm(c.y));
            placedShapes.push_back(std::move(moved));
            placedOrient.push_back(op);
            layout.push_back(c);
            if(gVerbose){
                std::lock_guard<std::mutex> lock(output_mutex);
                std::cerr << "[PLACED] Part " << i << " (ID " << op.id << ") at (x=" << c.x << ", y=" << c.y << ", angle=" << c.ang << ")\n";
            }
            if(i==0){
                Rect64 fb = op.bb;
                fb.left   += I64mm(c.x);
                fb.right  += I64mm(c.x);
                fb.bottom += I64mm(c.y);
                fb.top    += I64mm(c.y);
                std::lock_guard<std::mutex> lock(output_mutex);
                if(gVerbose)
                    std::cerr << "[FIRST] bbox=" << rectToStr(fb)
                              << " pos=(" << c.x << "," << c.y << ")\n";
            }
            ok = true;
        }
        auto t_end = std::chrono::steady_clock::now();
        double t_sec = std::chrono::duration<double>(t_end - t_start).count();
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            if(gVerbose){
                std::cerr << "[INFO] Время размещения детали " << i << ": " << t_sec << " сек\n";
                std::cerr << "[INFO] Невалидных кандидатов по bbox для детали " << i
                          << ": " << invalidBB_total << "\n";
                if(i==0 && firstTooBig)
                    std::cerr << "Первый part слишком велик\n";
            }
        }
        if (!ok && gVerbose) {
            std::cerr << "[WARN] Part " << i << " не удалось разместить – пропускаю\n";
        }
    }
    spdlog::info("[GREEDY] overlap checks {} (GPU {})", total_checks_all, total_checks_gpu);
    return layout;
}


// Display usage information
static void printHelp(const char* exe){
    std::cout << "Usage: " << exe << " -s WxH [options] files...\n"
              << "Options:\n"
              << "  -s, --sheet WxH    Sheet size in mm\n"
             << "  -r, --rot N       Rotation step in degrees (default 45)\n"
             << "  -o, --out FILE    Output CSV file (default layout.csv)\n"
             << "  --dxf FILE        Output DXF layout (default layout.dxf)\n"
             << "  -j, --join-segments  Join broken segments into loops\n"
             << "  -i, --iter N      Number of random iterations (default 1)\n"
             << "  -g, --grid N      Grid step in mm (default 10)\n"
             << "  -h, --help        Show this help message\n"
             << "  --cand N          Candidate positions per part (default 1000)\n"
            << "  --nfp N           NFPs per part (default 200)\n"
            << "  --overlap N       Overlap checks per part (default 500)\n"
            << "  --time N          Time limit per part in seconds (default 2.0)\n"
             << "  -v, --verbose     Enable verbose output\n"
             << "  --run N           Run index for multi runs\n"
             << "  --fill-gaps       Enable gap nesting\n";
}

static CLI parse(int ac, char** av){
    CLI c;
    cxxopts::Options options(av[0], "Nesting algorithm");
    options.add_options()
        ("s,sheet", "sheet WxH", cxxopts::value<std::string>())
        ("r,rot", "rotation step", cxxopts::value<int>())
        ("o,out", "output csv", cxxopts::value<std::string>())
        ("dxf", "output dxf", cxxopts::value<std::string>())
        ("iter", "iterations", cxxopts::value<int>())
        ("pop", "population size", cxxopts::value<int>())
        ("gen", "generations", cxxopts::value<int>())
        ("strategy", "strategy", cxxopts::value<std::string>())
        ("kerf", "kerf width", cxxopts::value<double>())
        ("gap", "gap between parts", cxxopts::value<double>())
        ("verbose", "verbose", cxxopts::value<bool>()->default_value("false"))
        ("fill-gaps", "fill gaps")
        ("help", "print help");

    auto result = options.parse(ac, av);
    if(result.count("help") || ac==1){
        std::cout << options.help() << "\n";
        std::exit(0);
    }

    // sheet WxH
    if(result.count("sheet")) {
        auto sheet = result["sheet"].as<std::string>();
        auto x = sheet.find('x');
        c.W = std::stod(sheet.substr(0,x));
        c.H = std::stod(sheet.substr(x+1));
    }

    // Заполнение с дефолтами 
    c.rot         = result.count("rot")      ? result["rot"].as<int>()           : DEFAULT_ROT_STEP;
    c.out         = result.count("out")      ? result["out"].as<std::string>()   : "layout.csv";
    c.dxf         = result.count("dxf")      ? result["dxf"].as<std::string>()   : "layout.dxf";
    c.iter        = result.count("iter")     ? result["iter"].as<int>()          : 1;
    c.pop_size    = result.count("pop")      ? result["pop"].as<int>()           : 50;
    c.generations = result.count("gen")      ? result["gen"].as<int>()           : 50;
    c.strategy    = result.count("strategy") ? result["strategy"].as<std::string>() : "area";
    gKerf         = result.count("kerf")     ? result["kerf"].as<double>()       : 0.0;
    gGap          = result.count("gap")      ? result["gap"].as<double>()        : 0.0;
    c.verbose     = result["verbose"].as<bool>();
    c.fill_gaps   = result.count("fill-gaps");

    // Файлы (только .json)
    for(int i=1;i<ac;++i){
        std::string a = av[i];
        if(a.size() && a[0] != '-' && a.find(".json") != std::string::npos)
            c.files.push_back(a);
    }
    if(c.W<=0 || c.H<=0 || c.files.empty())
        throw std::runtime_error("use --help for usage");
    return c;
}


std::string PolylineToDXF(const Path64& path) {
    std::ostringstream ss;
    ss << "0\nLWPOLYLINE\n8\n0\n90\n" << path.size() << "\n70\n1\n"; // 70=1 -> closed polyline
    for (const auto& p : path) {
        ss << "10\n" << Dbl(p.x) << "\n20\n" << Dbl(p.y) << "\n";
    }
    return ss.str();
}
void ExportToDXF(const std::string& filename, const std::vector<RawPart>& parts) {
    std::ofstream f(filename);
    f << "0\nSECTION\n2\nENTITIES\n";

    for (const auto& part : parts) {
        for (const auto& ring : part.rings) {  
            f << PolylineToDXF(ring);
        }
        for (const auto& hole : part.holes) {
            f << PolylineToDXF(hole);
        }
    }

    f << "0\nENDSEC\n0\nEOF\n";
}
std::string PlacedPolylineToDXF(const Path64& path, double x, double y, double angle_deg) {
    std::ostringstream ss;
    double angle_rad = angle_deg * M_PI / 180.0;
    double s = sin(angle_rad), c = cos(angle_rad);
    ss << "0\nLWPOLYLINE\n8\n0\n90\n" << path.size() << "\n70\n1\n";
    for (const auto& p : path) {
        // Применяем поворот и сдвиг!
        double px = Dbl(p.x), py = Dbl(p.y);
        double nx = c * px - s * py + x;
        double ny = s * px + c * py + y;
        ss << "10\n" << nx << "\n20\n" << ny << "\n";
    }
    return ss.str();
}

// Экспорт всей раскладки в DXF:
void ExportPlacedToDXF(const std::string& filename,
                       const std::vector<RawPart>& parts,
                       const std::vector<Place>& placed,
                       double W, double H)

{
    std::ofstream f(filename);
    f << "0\nSECTION\n2\nENTITIES\n";
    // Рисуем границу листа
    Path64 sheet;
    sheet.push_back(Point64(int64_t(0), int64_t(0)));
    sheet.push_back(Point64(I64mm(W), int64_t(0)));
    sheet.push_back(Point64(I64mm(W), I64mm(H)));
    sheet.push_back(Point64(int64_t(0), I64mm(H)));
    sheet.push_back(Point64(int64_t(0), int64_t(0)));
    f << PolylineToDXF(sheet);

    for (const auto& pl : placed) {
        const RawPart& part = parts[pl.id];
        for (const auto& ring : part.rings) {
            f << PlacedPolylineToDXF(ring, pl.x, pl.y, pl.ang);
        }
        for (const auto& hole : part.holes) {
            f << PlacedPolylineToDXF(hole, pl.x, pl.y, pl.ang);
        }
    }
    f << "0\nENDSEC\n0\nEOF\n";
}

void ExportPlacedToSVG(const std::string& filename,
                       const std::vector<RawPart>& parts,
                       const std::vector<Place>& placed,
                       double W, double H)
{
    std::ofstream f(filename);
    f << "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 " << W << ' ' << H << "'>\n";
    f << "<rect x='0' y='0' width='" << W << "' height='" << H << "' fill='none' stroke='black'/>\n";
    for (const auto& pl : placed) {
        const RawPart& part = parts[pl.id];
        double rad = pl.ang * M_PI / 180.0;
        double s = sin(rad), c = cos(rad);
        auto emit = [&](const Path64& ring) {
            f << "<polyline fill='none' stroke='red' points='";
            for (const auto& p : ring) {
                double px = Dbl(p.x); double py = Dbl(p.y);
                double nx = c*px - s*py + pl.x;
                double ny = s*px + c*py + pl.y;
                f << nx << ',' << ny << ' ';
            }
            f << "'/>\n";
        };
        for (const auto& ring : part.rings) emit(ring);
        for (const auto& hole : part.holes) emit(hole);
    }
    f << "</svg>\n";
}

// --- Genetic Algorithm Structures ---
struct Genome {
    std::vector<int> order;      // порядок деталей
    std::vector<int> angles;     // ориентация (индекс в all_orients) для каждой детали
    double fitness = 1e100;      // площадь bounding box
    std::vector<Place> layout;   // результат greedy
};

// --- Генерация случайной особи ---
Genome random_genome(size_t n, size_t n_orients, std::mt19937& rng) {
    Genome g;
    g.order.resize(n);
    std::iota(g.order.begin(), g.order.end(), 0);
    std::shuffle(g.order.begin(), g.order.end(), rng);
    g.angles.resize(n);
    for (size_t i = 0; i < n; ++i)
        g.angles[i] = rng() % n_orients;
    return g;
}

// --- Кроссовер (однородный) ---
Genome crossover(const Genome& a, const Genome& b, std::mt19937& rng) {
    size_t n = a.order.size();
    Genome child;
    child.order = a.order;
    // PMX crossover для порядка
    std::uniform_int_distribution<size_t> dist(0, n-1);
    size_t l = dist(rng), r = dist(rng);
    if (l > r) std::swap(l, r);
    std::set<int> taken;
    for (size_t i = l; i <= r; ++i) {
        child.order[i] = b.order[i];
        taken.insert(b.order[i]);
    }
    size_t ai = 0;
    for (size_t i = 0; i < n; ++i) {
        if (i >= l && i <= r) continue;
        while (taken.count(a.order[ai])) ++ai;
        child.order[i] = a.order[ai++];
    }
    // Ориентации — случайно от одного из родителей
    child.angles.resize(n);
    for (size_t i = 0; i < n; ++i)
        child.angles[i] = (rng() & 1) ? a.angles[i] : b.angles[i];
    return child;
}

// --- Мутация ---
void mutate(Genome& g, size_t n_orients, std::mt19937& rng, double rate=0.1) {
    std::uniform_real_distribution<double> prob(0.0, 1.0);
    if (prob(rng) < rate) {
        // swap two parts
        size_t i = rng() % g.order.size(), j = rng() % g.order.size();
        std::swap(g.order[i], g.order[j]);
    }
    for (size_t i = 0; i < g.angles.size(); ++i)
        if (prob(rng) < rate)
            g.angles[i] = rng() % n_orients;
}

// --- Оценка приспособленности (fitness) ---
double evaluate_genome(
    Genome& g,
    const std::vector<RawPart>& parts,
    const std::vector<std::vector<Orient>>& all_orients,
    double W, double H, double grid,
    LRUCache<Key, Paths64>& sharedNFP,
    const CLI* cli_ptr
) {
    // Собираем перестановку деталей в соответствии с порядком из генома
    std::vector<RawPart> perm(parts.size());
    std::vector<std::vector<Orient>> perm_orients(parts.size());
    for (size_t i = 0; i < parts.size(); ++i) {
        perm[i] = parts[g.order[i]];
        perm_orients[i] = all_orients[g.order[i]];
    }

    // Вызываем greedy с полной таблицей ориентаций (perm_orients)
    std::unordered_map<Key, Paths64> localNFP;
    auto layout = greedy(perm, perm_orients, W, H, grid, 0, localNFP, sharedNFP, cli_ptr);

    double area = computeArea(layout, perm);
    g.fitness = area;
    g.layout = layout;
    return area;
}

static std::vector<Place> local_search(const std::vector<Place>& layout,
                                       const std::vector<RawPart>& parts,
                                       const std::vector<std::vector<Orient>>& all_orients,
                                       double W, double H, double grid,
                                       size_t max_passes){
    auto best = layout;
    double bestA = computeArea(best, parts);
    double step = grid*2.0;
    for(size_t pass=0; pass<max_passes; ++pass){
        bool improved=false;
        for(size_t i=0;i<best.size();++i){
            Place orig = best[i];
            Place bestLocal = orig;
            double localBest = bestA;
            for(const auto& o : all_orients[orig.id]){
                for(int dx=-1;dx<=1;++dx) for(int dy=-1;dy<=1;++dy){
                    double nx=orig.x+dx*step, ny=orig.y+dy*step;
                    if(nx<0||ny<0||nx>W||ny>H) continue;
                    Place cand{orig.id,nx,ny,o.ang};
                    auto temp=best; temp[i]=cand; bool clash=false;
                    for(size_t j=0;j<temp.size();++j){ if(i==j) continue;
                        double r1=cand.ang*M_PI/180.0; auto p1=movePaths(rot(parts[cand.id].rings,sin(r1),cos(r1)),I64mm(cand.x),I64mm(cand.y));
                        double r2=temp[j].ang*M_PI/180.0; auto p2=movePaths(rot(parts[temp[j].id].rings,sin(r2),cos(r2)),I64mm(temp[j].x),I64mm(temp[j].y));
                        if(overlap(p1,p2)){clash=true;break;} }
                    if(clash) continue;
                    double a=computeArea(temp, parts);
                    if(a<localBest){localBest=a; bestLocal=cand;}
                }
            }
            if(localBest < bestA){ best[i]=bestLocal; bestA=localBest; improved=true; }
        }
        // pairwise swap search
        for(size_t a=0;a<best.size();++a) for(size_t b=a+1;b<best.size();++b){
            auto temp=best; std::swap(temp[a], temp[b]);
            double aA=computeArea(temp, parts);
            if(aA < bestA){ best=temp; bestA=aA; improved=true; }
        }
        if(!improved){
            if(step<=grid/2) break;
            step/=2.0;
        }
    }
    return best;
}

// --- Gap fill: try to place remaining parts into free areas ---
static std::vector<Place> fillGaps(std::vector<Place> layout,
                                   const std::vector<RawPart>& parts,
                                   const std::vector<std::vector<Orient>>& all_orients,
                                   double W, double H, double grid){
    std::set<int> placedIds;
    for(const auto& pl: layout) placedIds.insert(pl.id);

    Paths64 placedPaths;
    for(const auto& pl: layout){
        double rad = pl.ang*M_PI/180.0;
        Paths64 r = movePaths(rot(parts[pl.id].rings,sin(rad),cos(rad)),I64mm(pl.x),I64mm(pl.y));
        placedPaths.insert(placedPaths.end(), r.begin(), r.end());
    }
    Path64 sheet;
    sheet.push_back(Point64(int64_t(0),int64_t(0)));
    sheet.push_back(Point64(I64mm(W),int64_t(0)));
    sheet.push_back(Point64(I64mm(W),I64mm(H)));
    sheet.push_back(Point64(int64_t(0),I64mm(H)));
    sheet.push_back(Point64(int64_t(0),int64_t(0)));
    Paths64 free = Difference(Paths64{sheet}, Union(placedPaths,FillRule::NonZero), FillRule::NonZero);

    for(size_t pid=0; pid<parts.size(); ++pid){
        if(placedIds.count(pid)) continue;
        const auto& oris = all_orients[pid];
        bool done=false;
        for(const auto& gap : free){
            Rect64 gb = getBBox(Paths64{gap});
            for(const auto& o: oris){
                if(o.bb.right-o.bb.left > gb.right-gb.left ||
                   o.bb.top-o.bb.bottom > gb.top-gb.bottom) continue;
                Place pl{(int)pid, Dbl(gb.left), Dbl(gb.bottom), o.ang};
                Paths64 moved = movePaths(o.poly, gb.left, gb.bottom);
                bool clash=false;
                for(const auto& pp: placedPaths){
                    if(overlap(Paths64{pp}, moved)){ clash=true; break; }
                }
                if(!clash){
                    layout.push_back(pl);
                    placedPaths.insert(placedPaths.end(), moved.begin(), moved.end());
                    free = Difference(Paths64{sheet}, Union(placedPaths,FillRule::NonZero), FillRule::NonZero);
                    done=true; break;
                }
            }
            if(done) break;
        }
    }
    return layout;
}


// --- Генетический алгоритм основной цикл ---
// Добавлен параметр seed_layout для передачи стартовой особи
std::vector<Place> genetic_nesting(
    const std::vector<RawPart>& parts,
    const std::vector<std::vector<Orient>>& all_orients,
    double W, double H, double grid,
    LRUCache<Key, Paths64>& sharedNFP,
    const CLI* cli_ptr,
    int generations = 4,
    int pop_size = 4,
    const std::vector<Place>& seed_layout = {} 
) {
    std::mt19937 rng(std::random_device{}());
    size_t n = parts.size();
    size_t n_orients = all_orients[0].size();
    std::vector<Genome> pop(pop_size);
    // Инициализация популяции
    for (int i = 0; i < pop_size; ++i)
        pop[i] = random_genome(n, n_orients, rng);

    // Если есть seed_layout, используем его как первую особь
    if (!seed_layout.empty()) {
        // Восстановить порядок и ориентации из seed_layout
        Genome g;
        g.order.resize(seed_layout.size());
        g.angles.resize(seed_layout.size());
        for (size_t i = 0; i < seed_layout.size(); ++i) {
            g.order[i] = seed_layout[i].id;
            // Найти индекс ориентации по углу
            int ang_idx = 0;
            for (size_t j = 0; j < all_orients[g.order[i]].size(); ++j) {
                if (all_orients[g.order[i]][j].ang == seed_layout[i].ang) {
                    ang_idx = int(j);
                    break;
                }
            }
            g.angles[i] = ang_idx;
        }
        evaluate_genome(g, parts, all_orients, W, H, grid, sharedNFP, cli_ptr);
        pop[0] = g;
    }

    // ...existing code for GA generations...
    for (int gen = 0; gen < generations; ++gen) {
        // Сортировка по fitness
        std::sort(pop.begin(), pop.end(), [](const Genome& a, const Genome& b) {
            return a.fitness < b.fitness;
        });
        // Элитизм: сохраняем лучших
        int elite = pop_size / 5;
        std::vector<Genome> next(pop.begin(), pop.begin() + elite);
        // Новое поколение
        std::vector<Genome> children(pop_size - elite);

        auto tournament = [&](std::mt19937& trng) {
            int best = trng() % pop_size;
            for(int t=1;t<3;++t){
                int cand = trng() % pop_size;
                if(pop[cand].fitness < pop[best].fitness) best = cand;
            }
            return best;
        };

        #pragma omp parallel for
        for (int k = 0; k < int(children.size()); ++k) {
            std::mt19937 thread_rng(rng());
            int i1 = tournament(thread_rng);
            int i2 = tournament(thread_rng);
            Genome child = crossover(pop[i1], pop[i2], thread_rng);
            mutate(child, n_orients, thread_rng, 0.2);
            evaluate_genome(child, parts, all_orients, W, H, grid, sharedNFP, cli_ptr);
            children[k] = std::move(child);
        }

        // потом добавь их к next (после elits):
        for (auto& child : children)
            next.push_back(std::move(child));
        pop = std::move(next);
        spdlog::info("GA generation {} best area: {} mm2 placed {}/{}",
                     gen+1, pop[0].fitness, pop[0].layout.size(), n);
    }
    // Вернуть лучшую раскладку
    return pop[0].layout;
}


// Сконцентрированная версия функции main с «слиянием» логики greedy‑итераций и GA‑улучшения.
// * greedy‑раскладка теперь считается один раз для каждой «итерации»
// * лучший результат greedy передаётся в genetic_nesting как стартовая особь
// * итоговый bestLayout определяется как min(area) среди greedy & GA
// -----------------------------------------------------------------------------
#ifndef NEST_UNIT_TEST
int main(int argc, char* argv[])
{
#ifdef _OPENMP
    std::cout << "OpenMP threads: " << omp_get_max_threads() << '\n';
#else
    std::cout << "OpenMP DISABLED\n";
#endif
#ifdef USE_CUDA
    cudaDeviceProp prop{};
    bool gpu = cuda_available();
    if (gpu) cudaGetDeviceProperties(&prop, 0);
    std::cout << "[DEBUG] NEST build " << __DATE__ << ' ' << __TIME__
              << ", CUDA=" << (gpu ? 1 : 0)
              << ", device=" << (gpu ? prop.name : "none") << "\n";
    if (gpu)
        std::cout << "[INFO] CUDA GPU detected and will be used for overlap checks\n";
    else
        std::cout << "[INFO] CUDA NOT detected, using CPU only\n";
#else
    std::cout << "[DEBUG] NEST build " << __DATE__ << ' ' << __TIME__
              << ", CUDA=0, device=none\n";
    std::cout << "[INFO] CUDA NOT compiled in this binary, only CPU\n";
#endif
    


    try {
        CLI cli = parse(argc, argv);
        gVerbose = cli.verbose;
        spdlog::set_pattern("[%H:%M:%S] %^%l%$ %v");
        spdlog::set_level(gVerbose ? spdlog::level::info : spdlog::level::warn);
        // ─── Debug‑вывод параметров ───────────────────────────────────────────
        if(gVerbose)
            spdlog::info("Sheet: {}x{} mm, grid={} mm, rot step={}°, iter={}",
                          cli.W, cli.H, cli.grid, cli.rot, cli.iter);

        // ─── Загрузка деталей из JSON ────────────────────────────────────────
        spdlog::info("Loading parts from JSON files...");
        std::vector<RawPart> parts;
        for (size_t idx = 0; idx < cli.files.size(); ++idx) {
            auto batch = loadPartsFromJson(cli.files[idx]);
            int copies = (idx < cli.nums.size() ? cli.nums[idx] : 1);
            for (int c = 0; c < copies; ++c)
                for (const auto& base : batch) {
                    RawPart p = base;
                    p.id = int(parts.size());
                    parts.push_back(std::move(p));
                }
        }
        if (parts.empty())
            throw std::runtime_error("no parts loaded");
        spdlog::info("Loaded {} part instances", parts.size());

        // ─── Предрасчёт поворотов ────────────────────────────────────────────
        spdlog::info("Precomputing orientations...");
        const auto all_orients = makeOrient(parts, cli.rot);

        // Precompute NFP cache using GPU where possible
        spdlog::info("Precomputing NFP cache...");
        precomputeNFPGPU(all_orients, sharedNFP);

        // ─── Greedy‑итерации (параллельно) ───────────────────────────────────
        const int total_iter = std::max(1, cli.iter);
        spdlog::info("Starting greedy iterations: {}", total_iter);
        std::vector<std::vector<Place>> greedyLayouts(total_iter);
        std::vector<double>                greedyArea(total_iter, 1e100);

        std::atomic<int> done{0};
        std::vector<std::thread> workers;
        for (int it = 0; it < total_iter; ++it) {
            workers.emplace_back([&, it]() {
                std::unordered_map<Key, Paths64> localNFP;
                unsigned int seed = (it == 0 ? 0u :                         // первая итерация детерминирована
                                      unsigned(std::chrono::steady_clock::now()
                                               .time_since_epoch().count()) + it);
                auto layout = greedy(parts, all_orients, cli.W, cli.H, cli.grid,
                                      seed, localNFP, sharedNFP, &cli);
                double area  = computeArea(layout, parts);
                greedyLayouts[it] = std::move(layout);
                greedyArea[it]   = area;
                if ((++done % 2 == 0 || done == total_iter) && gVerbose)
                    std::cerr << "[GREEDY] finished " << done << "/" << total_iter << "\n";
            });
        }
        for (auto& t : workers) t.join();

        // ─── Лучший greedy‑результат ─────────────────────────────────────────
        int bestGreedyIdx = int(std::min_element(greedyArea.begin(), greedyArea.end()) - greedyArea.begin());
        const auto& seedLayout = greedyLayouts[bestGreedyIdx];
        double bestGreedyArea  = greedyArea[bestGreedyIdx];
        if(gVerbose)
            spdlog::info("Greedy best area {} mm² (iter {})", bestGreedyArea, bestGreedyIdx);

        // ─── Genetic‑улучшение ───────────────────────────────────────────────
        spdlog::info("Starting genetic optimisation (gen={}, pop={})", cli.generations, cli.pop_size);
        auto gaLayout = genetic_nesting(parts, all_orients,
                                        cli.W, cli.H, cli.grid,
                                        sharedNFP, &cli,
                                        cli.generations,
                                        cli.pop_size,
                                        seedLayout);
        if(cli.polish>0)
            gaLayout = local_search(gaLayout, parts, all_orients, cli.W, cli.H, cli.grid, cli.polish);
        double gaArea = computeArea(gaLayout, parts);
        if(gVerbose)
            spdlog::info("GA final area {} mm²", gaArea);

        // ─── Выбор окончательного расклада ───────────────────────────────────
        auto bestLayout = (gaArea < bestGreedyArea ? gaLayout : seedLayout);
        if(cli.fill_gaps)
            bestLayout = fillGaps(bestLayout, parts, all_orients, cli.W, cli.H, cli.grid);

        // ─── Экспорт результатов ─────────────────────────────────────────────
        spdlog::info("Exporting results to {} and {}", cli.out, cli.dxf);
        std::ofstream csv(cli.out);
        LayoutStats st = computeStats(bestLayout, parts);
        double fill = st.placedArea/(cli.W*cli.H)*100.0;
        double waste = cli.W*cli.H - st.placedArea;
        csv << "run,strategy,part,x_mm,y_mm,angle,fill_pct,waste_mm2\n";
        for (const auto& pl : bestLayout)
            csv << cli.run << ',' << cli.strategy << ',' << pl.id << ','
                << std::fixed << std::setprecision(3)
                << pl.x << ',' << pl.y << ',' << pl.ang << ','
                << fill << ',' << waste << "\n";
        std::cout << "placed " << bestLayout.size() << '/' << parts.size()
                  << "  ->  " << cli.out << "\n";
        ExportPlacedToDXF(cli.dxf, parts, bestLayout, cli.W, cli.H);
        ExportPlacedToSVG("layout.svg", parts, bestLayout, cli.W, cli.H);
        spdlog::info("DXF exported to {}", cli.dxf);
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "ERR: " << e.what() << "\n";
        return 1;
    }
}
#endif // NEST_UNIT_TEST