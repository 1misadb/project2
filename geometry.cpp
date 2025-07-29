#include "geometry.h"
#include <numeric>

double gKerf = 0.0;
double gGap = 0.0;
static constexpr double SCALE = 1e4;

Paths64 movePaths(const Paths64& src, int64_t dx, int64_t dy){
    Paths64 out; out.reserve(src.size());
    for(const auto& path: src){
        Path64 p; p.reserve(path.size());
        for(auto pt: path) p.push_back({pt.x + dx, pt.y + dy});
        out.push_back(std::move(p));
    }
    return out;
}

bool overlap(const Paths64& a, const Paths64& b) {
    double delta = (gKerf * 0.5 + gGap) * SCALE;
    Paths64 ea = InflatePaths(a, delta, JoinType::Miter, EndType::Polygon);
    Paths64 eb = InflatePaths(b, delta, JoinType::Miter, EndType::Polygon);
    Paths64 ua = Union(ea, FillRule::NonZero);
    Paths64 ub = Union(eb, FillRule::NonZero);
    Paths64 isect = Intersect(ua, ub, FillRule::NonZero);
    return std::abs(Area(isect)) > 0.0;
}

