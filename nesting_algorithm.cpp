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
//      cl /EHsc /std:c++17 /O2 clipper3\src\clipper.engine.cpp clipper3\src\clipper.offset.cpp clipper3\src\clipper.rectclip.cpp nesting_algorithm.cpp /Iclipper3\include /I. /Fe:nest.exe
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
#include <map>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cstdint>
#include <array>
#include <random>
#include <mutex>
#include <chrono>
#include <atomic>
#include <thread>
#define _USE_MATH_DEFINES
#include <math.h>

#include "clipper3/clipper.h"
#include "clipper3/clipper.minkowski.h"
#include <json.hpp>
double gUnit = 1.0;
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
static std::unordered_map<Key, Paths64> gNfpCache;
static std::mutex gNfpMutex;
std::mutex result_mutex;
std::mutex output_mutex;
std::mutex nfp_mutex;
std::unordered_map<Key, Paths64> sharedNFP;
constexpr double SCALE = 1e4;
constexpr int DEFAULT_ROT_STEP = 45;
constexpr double TOL_MM = 1.0;                 // base tolerance in mm
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

// Test for polygon overlap

bool overlap(const Paths64& a, const Paths64& b) {
    Paths64 ua = Union(a, FillRule::NonZero);
    Paths64 ub = Union(b, FillRule::NonZero);
    Paths64 isect = Intersect(ua, ub, FillRule::NonZero);
    return std::abs(Area(isect)) > 0.0;
}
// Translate polygon set by (dx,dy)
static Paths64 movePaths(const Paths64& src, int64_t dx, int64_t dy){
    Paths64 out; out.reserve(src.size());
    for(const auto& path: src){
        Path64 p; p.reserve(path.size());
        for(auto pt: path) p.push_back({pt.x + dx, pt.y + dy});
        out.push_back(std::move(p));
    }
    return out;
}

// ───────── DXF mini loader ─────────
struct RawPart{ Paths64 rings; double area=0; int id=0; Paths64 holes; };

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

        if (rd.code == "10") { x1 = std::stod(rd.value); have_x1 = true; }
        else if (rd.code == "20") { y1 = std::stod(rd.value); have_y1 = true; }
        else if (rd.code == "11") { x2 = std::stod(rd.value); have_x2 = true; }
        else if (rd.code == "21") { y2 = std::stod(rd.value); have_y2 = true; }
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
    std::cerr << "ctrlPts.size()=" << ctrlPts.size() << ", knotVec.size()=" << knotVec.size() << ", degree=" << degree << '\n';
    std::cerr << "SPLINE closed flag=" << closed << "   isClosed(path)=" << isClosed(out) << '\n';
    // --- Замыкание ---
    if (!out.empty()) {
        if (closed) {
            int64_t dist = distance(out.front(), out.back());
            if (dist > tol_mm()) {
                std::cerr << "[FORCE] spline " << segNo
                        << " forcibly closed, big gap=" << Dbl(dist) << " mm\n";
                out.push_back(out.front());
            } else if (dist <= tol_mm() && out.front() != out.back()) {
                out.push_back(out.front());
                std::cerr << "[AutoClose] spline " << segNo
                        << " closed (distance=" << Dbl(dist) << " mm)\n";
            } else if (dist > tol_mm()) {
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
    std::cerr << " flags: closed="<<closed<<"  degree="<<degree
              << "  fit="<<fitPts.size()<<"  ctrl="<<ctrlPts.size()
              << "  knots="<<knotVec.size() << '\n';

    if(out.empty())
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
                std::cerr << "[AutoClose] path " << out.size()
                          << " closed (distance=" << Dbl(dist) << " mm)\n";
            }
            else
            {
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
    std::cerr << "Dangling ends: ";
    for (auto& p : ptcnt)
        if (p.second % 2 != 0)
            std::cerr << "(" << Dbl(p.first.first*TOL) << "," << Dbl(p.first.second*TOL) << ") ";
    std::cerr << "\n";
    std::cerr << "Segments=" << segs.size() << " Paths=" << out.size()
              << " Closed=" << closedCnt << "\n";
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
    std::ifstream fin(filename);
    if(!fin)
        throw std::runtime_error("open "+filename);

    nlohmann::json j;
    fin >> j;
    if(!j.is_array())
        throw std::runtime_error("invalid JSON format in "+filename);

    auto parse_one = [](const nlohmann::json& arr)->RawPart {
        RawPart p;
        bool first = true;
        for(const auto& path : arr){
            if(!path.is_array())
                continue;
            Path64 ring;
            for(const auto& pt : path){
                if(pt.size() < 2) continue;
                double x = pt[0].get<double>();
                double y = pt[1].get<double>();
                ring.emplace_back(I64mm(x), I64mm(y));
            }
            if(ring.size() < 2) continue;
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


// ───────── NFP cache ─────────
static const Paths64& nfp(
    const Orient& A, const Orient& B,
    std::unordered_map<Key,Paths64>& localNFP,
    std::unordered_map<Key,Paths64>& gNfpCache,
    std::mutex& gNfpMutex)
{
    Key k = kfn(A.id, A.ang, B.id, B.ang);
    auto it = localNFP.find(k);
    if(it != localNFP.end()) return it->second;

    {
        std::lock_guard<std::mutex> lock(gNfpMutex);
        auto git = gNfpCache.find(k);
        if(git != gNfpCache.end()) {
            localNFP[k] = git->second;
            return localNFP[k];
        }
    }

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

    {
        std::lock_guard<std::mutex> lock(gNfpMutex);
        gNfpCache[k] = res;
        localNFP[k] = res;
    }
    return localNFP[k];
}


static std::vector<Point64> candidates(const std::vector<Paths64>& placed, const Paths64& current, int grid) {
    std::vector<Point64> out;

    for (const auto& shape : placed) {
        Rect64 bb = getBBox(shape);
        out.push_back({bb.left, bb.bottom});
        out.push_back({bb.right, bb.bottom});
        out.push_back({bb.left, bb.top});
        out.push_back({bb.right, bb.top});
    }

    // добавим сетку по всему листу
    int64_t step = static_cast<int64_t>(std::llround(grid * SCALE));
    for (int64_t y = 0; y <= 1e6; y += step) {
        for (int64_t x = 0; x <= 1e6; x += step) {
            out.push_back({x, y});
        }
    }

    return out;
}
// ───────── greedy placer ─────────
// Place a set of parts onto the sheet using a simple greedy algorithm
struct Place{int id; double x,y; int ang;};

static double computeArea(const std::vector<Place>& layout, const std::vector<RawPart>& parts){
    if(layout.empty()) return 0.0;
    Paths64 all;
    for(const auto& pl : layout){
        double rad = pl.ang * M_PI / 180.0;
        auto r = rot(parts[pl.id].rings, sin(rad), cos(rad));
        Paths64 moved = movePaths(r, I64mm(pl.x), I64mm(pl.y));
        all.insert(all.end(), moved.begin(), moved.end());
    }
    Rect64 bb = getBBox(all);
    return Dbl(bb.right - bb.left) * Dbl(bb.top - bb.bottom);
}
// Функция раскладки для одной случайной перестановки (ОДНА итерация!)
static std::vector<Place> greedy(
    const std::vector<RawPart>& parts,
    const std::vector<std::vector<Orient>>& all_orients,
    double W, double H, int grid,
    unsigned int seed,
    std::unordered_map<Key, Paths64>& localNFP,
    std::unordered_map<Key, Paths64>& sharedNFP,
    std::mutex& nfp_mutex
) {
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

    for (size_t i = 0; i < ord.size(); ++i) {
        std::cerr << "[STEP] Размещаю деталь " << i << " (ID " << ord[i].id << ") ...\n";
        if (ord[i].rings.empty()) {
            std::cerr << "[ERROR] Деталь " << i << " пуста! Пропускаю.\n";
            continue;
        }
        bool ok = false; Place best{}; Orient bestO{};
        const auto& oriSet = all_orients[ ord[i].id ];
        size_t invalidBB_total = 0;
        bool firstTooBig = (i == 0);

        for (const auto& op : oriSet) {
            if (op.bb.right <= sheetBR.x && op.bb.top <= sheetBR.y)
                if(i==0) firstTooBig = false;
            if (op.bb.right > sheetBR.x || op.bb.top > sheetBR.y)
                continue;
        std::vector<Point64> cand;  // просто объявление (БЕЗ = ...)

        auto generated = candidates(placedShapes, op.poly, grid);
        cand.insert(cand.end(), generated.begin(), generated.end());

        for(const auto& po : placedOrient) {
            const Paths64& nfpp = nfp(po, op, localNFP, sharedNFP, nfp_mutex);
            for(const auto& pth : nfpp)
                for(const auto& pt : pth)
                    cand.push_back(pt);
        }
        cand.push_back({0,0});
        for(const auto& ps : placedShapes){
            Rect64 bb = getBBox(ps);
            cand.push_back({bb.left, bb.bottom});
            cand.push_back({bb.right, bb.bottom});
            cand.push_back({bb.left, bb.top});
            cand.push_back({bb.right, bb.top});
        }
        int64_t step_i = mm2i(grid);
        for(int64_t y=0; y<=sheetBR.y; y+=step_i)
            for(int64_t x=0; x<=sheetBR.x; x+=step_i)
                cand.push_back({x,y});

        // --- Удалить дубликаты
        std::sort(cand.begin(), cand.end(), [](const Point64& a, const Point64& b) {
            return a.y == b.y ? a.x < b.x : a.y < b.y;
        });
        cand.erase(std::unique(cand.begin(), cand.end(), [](const Point64& a, const Point64& b) {
            return a.x == b.x && a.y == b.y;
        }), cand.end());

        // --- Ограничить число кандидатов
        const size_t CAND_LIMIT = 2000;
        if (cand.size() > CAND_LIMIT) cand.resize(CAND_LIMIT);


            double bestCost = 1e100;
            Place  bestPl{};    // лучший кандидат (id, x, y, ang)
            Orient bestO{};     // лучший поворот
            size_t bestIdx = -1;
            size_t valid = 0;
            for (size_t idx = 0; idx < cand.size(); ++idx) {
                const auto& c = cand[idx];
                Rect64 bb = op.bb; bb.left += c.x; bb.right += c.x; bb.bottom += c.y; bb.top += c.y;
                if (bb.left < 0 || bb.bottom < 0 || bb.right > sheetBR.x || bb.top > sheetBR.y){
                    ++invalidBB_total;
                    continue;
                }
                Paths64 moved = movePaths(op.poly, c.x, c.y);
                bool clash = false;
                for (size_t pi = 0; pi < placedShapes.size(); ++pi) {
                    const auto& pl = placedShapes[pi];
                    Rect64 bbPl = getBBox(pl);
                    Rect64 bbMoved = getBBox(moved);
                    double areaA = areaMM2(pl);
                    double areaB = areaMM2(moved);
                    {
                        std::lock_guard<std::mutex> lock(output_mutex);
                        std::cerr << "[CHECK] Aarea=" << areaA << " bb=" << rectToStr(bbPl)
                                  << " Barea=" << areaB << " bb=" << rectToStr(bbMoved) << "\n";
                    }
                    bool ov = overlap(pl, moved);
                    if (ov && (bbPl.right < bbMoved.left || bbPl.left > bbMoved.right ||
                               bbPl.top < bbMoved.bottom || bbPl.bottom > bbMoved.top)) {
                        std::lock_guard<std::mutex> lock(output_mutex);
                        std::cerr << "[DEBUG] Overlap true but bbox separate ...";
                    }
                    if (bbPl.right < bb.left || bbPl.left > bb.right ||
                        bbPl.top < bb.bottom || bbPl.bottom > bb.top)
                        continue;
                    if (ov) {
                        {
                            std::lock_guard<std::mutex> lock(output_mutex);
                            std::cerr << "[OVERLAP] new part " << ord[i].id
                                    << " at (" << Dbl(c.x) << "," << Dbl(c.y)
                                    << ") ang=" << op.ang
                                    << " intersects placed part "
                                    << placedOrient[pi].id << "\n";
                        }
                        clash = true;
                        break;
                    }
                }
                if (clash) continue;
                ++valid;
                if (valid <= 10 || valid == 100 || valid == 500 || valid == 1000) {
                    std::lock_guard<std::mutex> lock(output_mutex);
                    std::cerr << "  [CAND] #" << idx << ": x=" << Dbl(c.x)
                            << " y=" << Dbl(c.y)
                            << " (valid=" << valid << ")\n";
                }

                // --- Вот здесь критерий: например, минимальный верхний Y (минимальная высота слоя)
                double cost = Dbl(bb.top);
                if (cost < bestCost) {
                    bestCost = cost;
                    bestPl   = { op.id, Dbl(c.x), Dbl(c.y), op.ang };
                    bestO    = op;
                    bestIdx  = idx;
                }
            }
            std::cerr << "[STEP] Всего валидных кандидатов: " << valid << '\n';
            // --- если нашёлся хоть один валидный кандидат — добавляем только его
            if (bestCost < 1e100) {
                auto& c = cand[bestIdx];
                Paths64 moved = movePaths(op.poly, I64mm(bestPl.x), I64mm(bestPl.y));
                placedShapes.push_back(std::move(moved));
                placedOrient.push_back(bestO);
                layout.push_back(bestPl);
                if(i==0){
                    Rect64 fb = bestO.bb;
                    fb.left   += I64mm(bestPl.x);
                    fb.right  += I64mm(bestPl.x);
                    fb.bottom += I64mm(bestPl.y);
                    fb.top    += I64mm(bestPl.y);
                    std::lock_guard<std::mutex> lock(output_mutex);
                    std::cerr << "[FIRST] bbox=" << rectToStr(fb)
                              << " pos=(" << bestPl.x << "," << bestPl.y << ")\n";
                }
                ok = true;
                break;
            }
            
        }
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cerr << "[INFO] Невалидных кандидатов по bbox для детали " << i
                      << ": " << invalidBB_total << "\n";
            if(i==0 && firstTooBig)
                std::cerr << "Первый part слишком велик\n";
        }
            if (!ok) {
                std::cerr << "[WARN] Part " << i << " не удалось разместить – пропускаю\n";
            }
    }
    return layout;
}

// ───────── CLI ─────────
struct CLI{
    double W = 0, H = 0;
    int rot = DEFAULT_ROT_STEP;
    int grid = 100;
    std::string out = "layout.csv";
    std::vector<std::string> files;
    std::vector<int> nums;
    int num = 1;
    bool joinSegments = false;
    int iter = 1; // number of random iterations
};

// Display usage information
static void printHelp(const char* exe){
    std::cout << "Usage: " << exe << " -s WxH [options] files...\n"
              << "Options:\n"
              << "  -s, --sheet WxH    Sheet size in mm\n"
             << "  -r, --rot N       Rotation step in degrees (default 45)\n"
             << "  -o, --out FILE    Output CSV file (default layout.csv)\n"
             << "  -j, --join-segments  Join broken segments into loops\n"
             << "  -i, --iter N      Number of random iterations (default 1)\n"
             << "  -g, --grid N      Grid step in mm (default 10)\n"
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
        }else if(a=="-j"||a=="--join-segments"){
            c.joinSegments = true;
        }else if(a=="-i"||a=="--iter"){
            if(i+1>=ac) throw std::runtime_error("missing value for --iter");
            c.iter = std::stoi(av[++i]);
        }else if(a=="-g"||a=="--grid"){
            if(i+1>=ac) throw std::runtime_error("missing value for --grid");
            c.grid = std::stoi(av[++i]);
            if(c.grid <= 0) throw std::runtime_error("--grid must be > 0");
        }else if(a=="-n"||a=="--num"){
            while (i+1<ac && std::isdigit(av[i+1][0]))
            {
                c.nums.push_back(std::stoi(av[++i]));
            }
            
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

// ───────── main ─────────
int main(int argc, char* argv[]) {
    try {
        auto cli = parse(argc, argv);
        
        std::vector<RawPart> parts;
        for(size_t idx=0; idx<cli.files.size(); ++idx) {
            auto pvec = loadPartsFromJson(cli.files[idx]);
            int cnt = (idx < cli.nums.size() ? cli.nums[idx] : 1);
            for(int r=0; r<cnt; ++r)
                for(const auto& base : pvec) {
                    RawPart p = base;
                    p.id = int(parts.size());
                    parts.push_back(p);
                }
        }
        
        // ВАЖНО: вычисляем повороты только один раз!
        auto all_orients = makeOrient(parts, cli.rot);
        for (size_t i = 0; i < parts.size(); ++i) {
            auto& p = parts[i];
            if (p.rings.empty()) {
                std::cerr << "[WARN] Part #" << i << " пустой!\n";
            } else {
                auto bb = getBBox(p.rings[0]);
                std::cerr << "[CHECK] Part #" << i
                        << " bbox: " << Dbl(bb.right - bb.left) << "x" << Dbl(bb.top - bb.bottom)
                        << " мм, area: " << p.area << " мм²\n";
            }
        }
        int total_iter = cli.iter > 0 ? cli.iter : 1;
        std::vector<std::vector<Place>> allLayouts(total_iter);
        std::cerr << "[INFO] Загружено " << parts.size() << " деталей\n";
        std::cerr << "[INFO] Начинаю " << total_iter << " итераций...\n";
        std::vector<double> allAreas(total_iter, std::numeric_limits<double>::max());

        std::atomic<int> finished_iters{0};
        std::vector<std::thread> threads;
        
        for (int it = 0; it < total_iter; ++it) {
            threads.emplace_back([&, it]() {
                std::unordered_map<Key, Paths64> localNFP;

                unsigned int my_seed = static_cast<unsigned int>(
                    std::chrono::steady_clock::now().time_since_epoch().count() + it * 1337
                );
                if (it == 0) my_seed = 0;
                std::cerr << "[INFO] Итер " << it+1 << "/" << total_iter << "...\n";
                auto layout = greedy(parts, all_orients, cli.W, cli.H, cli.grid, my_seed,
                                    localNFP, sharedNFP, nfp_mutex);
                double area = computeArea(layout, parts);

                {
                    std::lock_guard<std::mutex> lock(result_mutex);
                    allLayouts[it] = layout;
                    allAreas[it] = area;
                }
                {
                    std::lock_guard<std::mutex> out_lock(output_mutex);
                    std::cerr << "[DONE] Итерация " << it+1 << ": размещено "
                            << layout.size() << "/" << parts.size() << " деталей\n" << std::flush;
                }
                int done = ++finished_iters;
                if (done % 2 == 0 || done == total_iter) {
                    std::lock_guard<std::mutex> out_lock(output_mutex);
                    std::cerr << "[PROGRESS] " << done << "/" << total_iter << " итераций завершено\n";
                }
            });
        }
        
        for (auto& th : threads) th.join();
        {
            std::lock_guard<std::mutex> out_lock(output_mutex);
            std::cerr << "[INFO] Все итерации завершены\n";
        }
        // --- ОСТАЛЬНОЕ КАК РАНЬШЕ ---
        size_t bestIdx = 0;
        for (size_t i = 1; i < allAreas.size(); ++i)
            if (allAreas[i] < allAreas[bestIdx])
                bestIdx = i;
        const std::vector<Place>& finalBest = allLayouts[bestIdx];

        // Экспорт CSV
        std::ofstream f(cli.out);
        f << "part,x_mm,y_mm,angle\n";
        for (const auto& p : finalBest)
            f << p.id << "," << std::fixed << std::setprecision(3)
              << p.x << "," << p.y << "," << p.ang << "\n";

        std::cout << "placed " << finalBest.size() << "/" << parts.size()
                  << " → " << cli.out << "\n";
        ExportToDXF("output.dxf", parts);
        std::cout << "DXF exported to output.dxf\n";
        ExportPlacedToDXF("layout.dxf", parts, finalBest, cli.W, cli.H);
        std::cout << "Placed DXF exported to layout.dxf\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "ERR: " << e.what() << "\n";
        return 1;
    }
}
