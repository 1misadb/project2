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
#include <map>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <cstdint>
#include <array>
#define _USE_MATH_DEFINES
#include <math.h>

#include "clipper3/clipper.h"
#include "clipper3/clipper.minkowski.h"
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

constexpr double SCALE = 1e4;
constexpr int DEFAULT_ROT_STEP = 10;
constexpr double TOL_MM = 10.0;                 // base tolerance in mm
inline int64_t tol_mm() { return llround(TOL_MM * gUnit * SCALE); }
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

// Point comparison within tolerance
static bool eqPt(const Point64& a, const Point64& b)
{
    int64_t dx = a.x - b.x;
    int64_t dy = a.y - b.y;
    return (dx * dx + dy * dy) <= tol_mm() * tol_mm();
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
        case 70: closed = (std::stoi(rd.value) & 0x08); break;
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
    // --- Замыкание ---
    if (!out.empty()) {
        if (closed) {
            int64_t dist = distance(out.front(), out.back());
            if (dist <= tol_mm() && out.front() != out.back()) {
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


std::map<std::string, int> dxfStats;
std::vector<std::string> unsupportedEntities;
// Parse DXF entities into polygons
static std::vector<RawPart> loadDXF(const std::vector<std::string>& files,
                                    bool joinSegments){
    std::vector<RawPart> parts;
    for(size_t idx=0; idx<files.size(); ++idx){
        std::ifstream fin(files[idx]);
        if(!fin) throw std::runtime_error("open "+files[idx]);
        
        DXFReader rd(fin);
        bool inEnt = false;
        std::vector<Path64> rings, segs;
        int segNo = 0;
        int splineNo = 0;
        while(rd.next()) {
            dxfStats[rd.value]++;
            // Начало новой секции
            if(trim(rd.code) == "0" && rd.value == "SECTION") {
                if(rd.next() && trim(rd.code) == "2")
                    inEnt = (rd.value == "ENTITIES");
                continue;
            }
            // Конец секции — просто сбрасываем флаг!
            if(trim(rd.code) == "0" && rd.value == "ENDSEC") {
                inEnt = false;
                continue;
            }
            // Обрабатываем только объекты в секции ENTITIES
            if(!inEnt || trim(rd.code) != "0") continue;

            std::cerr << "ENTITY code=" << rd.code << " value=" << rd.value << "\n";
            Path64 tmp; bool closed = false;
            if(rd.value == "LWPOLYLINE") {
                if(parseLWPolyline(rd, tmp, closed, segNo++))
                    (closed ? rings : segs).push_back(std::move(tmp));
            } else if(rd.value == "LINE") {
                if(parseLine(rd, tmp, segNo++)) segs.push_back(std::move(tmp));
            } else if(rd.value == "ARC") {
                if(parseArc(rd, tmp, segNo++)) segs.push_back(std::move(tmp));
            } else if(rd.value == "CIRCLE") {
                if(parseCircle(rd, tmp, segNo++)) segs.push_back(std::move(tmp));
            } else if(rd.value == "ELLIPSE") {
                if(parseEllipse(rd, tmp, segNo++)) segs.push_back(std::move(tmp));
            } else if(rd.value == "SPLINE") {
                std::cerr << "== SPLINE #" << splineNo << " ==\n";
                if(parseSpline(rd, tmp, segNo++)) segs.push_back(std::move(tmp));
                splineNo++;
            } else if(rd.value == "POLYLINE") {
                if(parsePolyline(rd, tmp, segNo++)) segs.push_back(std::move(tmp));
            }
            std::cerr << "\n======= DXF ENTITY STATS =======\n";
        for (const auto& p : dxfStats) {
            std::cerr << "  " << p.first << ": " << p.second << "\n";
        }

        if (!unsupportedEntities.empty()) {
            std::cerr << "\n!!! UNSUPPORTED ENTITIES DETECTED !!!\n";
            std::map<std::string, int> unsupportedCount;
            for (const auto& e : unsupportedEntities) unsupportedCount[e]++;
            for (const auto& p : unsupportedCount)
                std::cerr << "  " << p.first << ": " << p.second << "\n";
        }
        std::cerr << "================================\n";
        }

        Paths64 joined = joinSegments ? connectSegments(segs)
                                      : Paths64(segs.begin(), segs.end());
        std::cerr << "DXF debug  " << files[idx]
                << "  segs=" << segs.size()
                << "  joined=" << joined.size()
                << "  rings=" << rings.size()
                << "  joinSegments=" << (joinSegments?"on":"off")
                << "\n";
        // Диагностика: первая точка первого path
        if (!joined.empty() && !joined[0].empty()) {
            auto pt = joined[0][0];
            std::cerr << " first pt: " << Dbl(pt.x) << ", " << Dbl(pt.y) << "\n";
        }
        // Добавляем все контуры из joined
        rings.insert(rings.end(), joined.begin(), joined.end());

        // Главное: если rings всё ещё пуст, копируем все контуры из joined!
        if (rings.empty() && !joined.empty())
            rings = joined;
        if (rings.empty()) throw std::runtime_error("no rings in " + files[idx]);
        // Диагностика: площадь каждого контура
        for (size_t i = 0; i < rings.size(); ++i)
            std::cerr << "rings[" << i << "].Area=" << Area(rings[i]) << "\n";
        std::sort(rings.begin(), rings.end(),
            [](const Path64 &a, const Path64 &b) {
                return std::abs(Area(a)) > std::abs(Area(b));
            });
        if (Area(rings[0]) < 0) ReversePath(rings[0]);
        for (size_t i = 1; i < rings.size(); ++i)
            if (Area(rings[i]) > 0) ReversePath(rings[i]);
        // Вот здесь выводим точки для каждой фигуры!
        for (size_t i = 0; i < rings.size(); ++i) {
            std::cerr << "rings[" << i << "]:";
            for (const auto& pt : rings[i])
                std::cerr << " (" << Dbl(pt.x) << "," << Dbl(pt.y) << ")";
            std::cerr << "\n";
        }
        // Далее обычное смещение к (0,0)
        Rect64 bb = getBBox(rings);
        int64_t dx = -bb.left, dy = -bb.bottom;
        if (dx || dy)
            for (auto &path : rings)
                for (auto &pt : path) { pt.x += dx; pt.y += dy; }
        // Диагностика: финальный bbox
        Rect64 bb_mm = getBBox(rings);
        std::cerr << " part " << idx
                << " bbox = " << Dbl(bb_mm.right - bb_mm.left)
                << " × " << Dbl(bb_mm.top - bb_mm.bottom) << " mm\n"
                << " bottom=" << Dbl(bb_mm.bottom)
                << " top=" << Dbl(bb_mm.top) << "\n";
        // Добавляем деталь в parts
        parts.push_back({rings, std::abs(Area(rings[0])) / (SCALE * SCALE), int(idx)});
        

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
    std::vector<int> nums;
    int num = 1;
    bool joinSegments = false;
};

// Display usage information
static void printHelp(const char* exe){
    std::cout << "Usage: " << exe << " -s WxH [options] files...\n"
              << "Options:\n"
              << "  -s, --sheet WxH    Sheet size in mm\n"
              << "  -r, --rot N       Rotation step in degrees (default 10)\n"
             << "  -o, --out FILE    Output CSV file (default layout.csv)\n"
             << "  -j, --join-segments  Join broken segments into loops\n"
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
        for (const auto& ring : part.rings) {  // экспорт всех колец каждой детали
            f << PolylineToDXF(ring);
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
                       const std::vector<Place>& placed)
{
    std::ofstream f(filename);
    f << "0\nSECTION\n2\nENTITIES\n";
    for (const auto& pl : placed) {
        const RawPart& part = parts[pl.id];
        for (const auto& ring : part.rings) {
            f << PlacedPolylineToDXF(ring, pl.x, pl.y, pl.ang);
        }
    }
    f << "0\nENDSEC\n0\nEOF\n";
}
// ───────── main ─────────
int main(int argc, char* argv[]) {
    try {
        auto cli = parse(argc, argv);
        int ucode      = detectDxfUnits(cli.files.front());
        std::cerr << "INSUNITS detected: " << ucode << "  gUnit=" << gUnit << std::endl;
        gUnit          = getDxfUnitFactor(ucode);
        std::cerr << "[DXF] INSUNITS=" << ucode
                  << "  (scale=" << gUnit << " mm per DXF unit)\n";
        auto parts = loadDXF(cli.files, cli.joinSegments);
        size_t nFile = cli.files.size();
        for(size_t i = 0; i < nFile; ++i){
            int cnt = (i < cli.nums.size()?cli.nums[i] : 1);
            for(int j = 0; j < cnt; ++j){
                RawPart p = parts[i];
                p.id = int(parts.size());
                parts.push_back(p);
            }
        }
        auto placed = greedy(parts, cli.W, cli.H, cli.rot);

        // Экспорт CSV
        std::ofstream f(cli.out);
        f << "part,x_mm,y_mm,angle\n";
        for (auto& p : placed)
            f << p.id << "," << std::fixed << std::setprecision(3)
              << p.x << "," << p.y << "," << p.ang << "\n";

        std::cout << "placed " << placed.size() << "/" << parts.size()
                  << " → " << cli.out << "\n";

        // Экспорт оригинальных контуров (можно оставить для сравнения)
        ExportToDXF("output.dxf", parts);
        std::cout << "DXF exported to output.dxf\n";

        // Экспорт размещённой раскладки!
        ExportPlacedToDXF("layout.dxf", parts, placed);
        std::cout << "Placed DXF exported to layout.dxf\n";

        return 0;

    } catch (const std::exception& e) {
        std::cerr << "ERR: " << e.what() << "\n";
        return 1;
    }
}