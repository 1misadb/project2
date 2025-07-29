#include "geometry.h"
#include <limits>
#include <algorithm>
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
    auto bvhA = buildBVH(ua);
    auto bvhB = buildBVH(ub);
    return overlapBVH(bvhA, ua, bvhB, ub);
}

// --- BVH utilities ---
static Rect64 bbox(const Path64& p){
    if(p.empty()) return Rect64{0,0,0,0};
    Rect64 r{p[0].x,p[0].y,p[0].x,p[0].y};
    for(auto pt: p){
        r.left = std::min(r.left, pt.x);
        r.right = std::max(r.right, pt.x);
        r.bottom = std::min(r.bottom, pt.y);
        r.top = std::max(r.top, pt.y);
    }
    return r;
}

std::vector<BVHNode> buildBVH(const Paths64& paths){
    std::vector<BVHNode> nodes;
    if(paths.empty()) return nodes;
    std::vector<int> idx(paths.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::function<int(int,int)> build = [&](int l,int r){
        int nodeIdx = (int)nodes.size();
        nodes.push_back({});
        if(l==r){
            nodes[nodeIdx].box = bbox(paths[idx[l]]);
            nodes[nodeIdx].index = idx[l];
            return nodeIdx;
        }
        Rect64 box{INT64_MAX,INT64_MAX,INT64_MIN,INT64_MIN};
        for(int i=l;i<=r;++i){
            Rect64 b = bbox(paths[idx[i]]);
            box.left = std::min(box.left,b.left);
            box.right= std::max(box.right,b.right);
            box.bottom=std::min(box.bottom,b.bottom);
            box.top= std::max(box.top,b.top);
        }
        nodes[nodeIdx].box = box;
        int axis = (box.right-box.left > box.top-box.bottom) ? 0:1;
        std::sort(idx.begin()+l, idx.begin()+r+1, [&](int a,int b){
            Rect64 ba=bbox(paths[a]), bb=bbox(paths[b]);
            int64_t ca = axis? (ba.bottom+ba.top):(ba.left+ba.right);
            int64_t cb = axis? (bb.bottom+bb.top):(bb.left+bb.right);
            return ca < cb;
        });
        int m=(l+r)/2;
        nodes[nodeIdx].left = build(l,m);
        nodes[nodeIdx].right = build(m+1,r);
        return nodeIdx;
    };
    build(0, (int)idx.size()-1);
    return nodes;
}

static bool boxesIntersect(const Rect64&a,const Rect64&b){
    return !(a.right < b.left || a.left > b.right || a.top < b.bottom || a.bottom > b.top);
}

static bool overlapBVHRec(const std::vector<BVHNode>& ta,const Paths64& pa,
                          const std::vector<BVHNode>& tb,const Paths64& pb,
                          int ia,int ib){
    const BVHNode& na=ta[ia];
    const BVHNode& nb=tb[ib];
    if(!boxesIntersect(na.box, nb.box)) return false;
    if(na.index>=0 && nb.index>=0){
        Paths64 ua{pa[na.index]};
        Paths64 ub{pb[nb.index]};
        Paths64 is=Intersect(ua,ub,FillRule::NonZero);
        return !is.empty();
    }
    if(na.index>=0){
        return overlapBVHRec(ta,pa,tb,pb,ia,nb.left) || overlapBVHRec(ta,pa,tb,pb,ia,nb.right);
    }
    if(nb.index>=0){
        return overlapBVHRec(ta,pa,tb,pb,na.left,ib) || overlapBVHRec(ta,pa,tb,pb,na.right,ib);
    }
    return overlapBVHRec(ta,pa,tb,pb,na.left,nb.left) ||
           overlapBVHRec(ta,pa,tb,pb,na.left,nb.right) ||
           overlapBVHRec(ta,pa,tb,pb,na.right,nb.left) ||
           overlapBVHRec(ta,pa,tb,pb,na.right,nb.right);
}

bool overlapBVH(const std::vector<BVHNode>& treeA, const Paths64& pa,
                const std::vector<BVHNode>& treeB, const Paths64& pb){
    if(treeA.empty() || treeB.empty()) return false;
    return overlapBVHRec(treeA,pa,treeB,pb,0,0);
}

