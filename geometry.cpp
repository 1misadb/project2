#include "geometry.h"
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif
#include <limits>
#include <algorithm>
#include <numeric>
#include <cstdint>

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

#ifdef USE_CUDA
struct GPUPath { int start; int size; };
struct GPUShape { int start; int size; };
extern "C" void overlapKernelLauncher(const long long* d_xs, const long long* d_ys, const GPUPath* d_paths,
                                      GPUShape d_cand, const GPUShape* d_shapes, int n, int* d_out);
#endif

bool cuda_available(){
#ifdef USE_CUDA
    int cnt=0; return cudaGetDeviceCount(&cnt)==cudaSuccess && cnt>0;
#else
    return false;
#endif
}

std::vector<bool> overlapBatchGPU(const Paths64& cand,
                                  const std::vector<Paths64>& others){
#ifdef USE_CUDA
    if(!cuda_available()){ // fallback
        std::vector<bool> r(others.size());
        for(size_t i=0;i<others.size();++i) r[i]=overlap(cand, others[i]);
        return r;
    }
    std::vector<long long> xs,ys; xs.reserve(1024); ys.reserve(1024);
    std::vector<GPUPath> paths; paths.reserve(256);
    GPUShape candShape{(int)paths.size(), (int)cand.size()};
    for(const auto& p:cand){
        GPUPath gp{(int)xs.size(), (int)p.size()};
        for(auto pt:p){ xs.push_back(pt.x); ys.push_back(pt.y); }
        paths.push_back(gp);
    }
    std::vector<GPUShape> shapes; shapes.reserve(others.size());
    for(const auto& sh:others){
        GPUShape s{(int)paths.size(), (int)sh.size()};
        for(const auto& p:sh){
            GPUPath gp{(int)xs.size(), (int)p.size()};
            for(auto pt:p){ xs.push_back(pt.x); ys.push_back(pt.y); }
            paths.push_back(gp);
        }
        shapes.push_back(s);
    }
    long long* d_xs; long long* d_ys; GPUPath* d_paths; GPUShape* d_shapes; int* d_out;
    size_t pts_sz=xs.size()*sizeof(long long); cudaMalloc(&d_xs,pts_sz); cudaMalloc(&d_ys,pts_sz);
    cudaMemcpy(d_xs,xs.data(),pts_sz,cudaMemcpyHostToDevice);
    cudaMemcpy(d_ys,ys.data(),pts_sz,cudaMemcpyHostToDevice);
    cudaMalloc(&d_paths,paths.size()*sizeof(GPUPath));
    cudaMemcpy(d_paths,paths.data(),paths.size()*sizeof(GPUPath),cudaMemcpyHostToDevice);
    cudaMalloc(&d_shapes,shapes.size()*sizeof(GPUShape));
    cudaMemcpy(d_shapes,shapes.data(),shapes.size()*sizeof(GPUShape),cudaMemcpyHostToDevice);
    cudaMalloc(&d_out,shapes.size()*sizeof(int));
    GPUShape d_cand=candShape; // copy by value

    // <<<--- вот тут вызывем только launcher-обёртку!
    overlapKernelLauncher(d_xs, d_ys, d_paths, d_cand, d_shapes, (int)shapes.size(), d_out);

    std::vector<int> res_raw(others.size());
    cudaMemcpy(res_raw.data(), d_out, others.size()*sizeof(int), cudaMemcpyDeviceToHost);
    std::vector<bool> res(others.size());
    for (size_t i = 0; i < others.size(); ++i)
        res[i] = (res_raw[i] != 0);

    cudaFree(d_xs); cudaFree(d_ys); cudaFree(d_paths); cudaFree(d_shapes); cudaFree(d_out);
    return res;
#else
    std::vector<bool> r(others.size());
    for(size_t i=0;i<others.size();++i) r[i]=overlap(cand, others[i]);
    return r;
#endif
}

// ---- Convex helpers ----
static long long cross64(const Point64& a, const Point64& b, const Point64& c){
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

bool isConvex(const Path64& p){
    if(p.size() < 4) return true;
    long long prev = 0;
    size_t n = p.size();
    for(size_t i=0;i<n;++i){
        long long cr = cross64(p[i], p[(i+1)%n], p[(i+2)%n]);
        if(cr != 0){
            if(prev != 0 && (cr > 0) != (prev > 0)) return false;
            prev = cr;
        }
    }
    return true;
}

Path64 convexHull(const std::vector<Point64>& pts){
    if(pts.size() < 3) return Path64(pts.begin(), pts.end());
    std::vector<Point64> p = pts;
    std::sort(p.begin(), p.end(), [](const Point64& a, const Point64& b){
        return a.x < b.x || (a.x == b.x && a.y < b.y);
    });
    std::vector<Point64> lower, upper;
    for(const auto& pt : p){
        while(lower.size() >= 2 && cross64(lower[lower.size()-2], lower.back(), pt) <= 0)
            lower.pop_back();
        lower.push_back(pt);
    }
    for(auto it=p.rbegin(); it!=p.rend(); ++it){
        while(upper.size() >= 2 && cross64(upper[upper.size()-2], upper.back(), *it) <= 0)
            upper.pop_back();
        upper.push_back(*it);
    }
    lower.pop_back();
    upper.pop_back();
    lower.insert(lower.end(), upper.begin(), upper.end());
    return Path64(lower.begin(), lower.end());
}

#ifdef USE_CUDA
struct GPUPair { GPUPath a; GPUPath b; int start; };
extern "C" void minkowskiKernelLauncher(const long long* ax,const long long* ay,
                                         const long long* bx,const long long* by,
                                         const GPUPair* pairs,int pairCount,
                                         long long* outx,long long* outy);
#endif

std::vector<Paths64> minkowskiBatchGPU(const std::vector<Path64>& A,
                                       const std::vector<Path64>& B){
#ifdef USE_CUDA
    if(!cuda_available()){
        std::vector<Paths64> out(A.size());
        for(size_t i=0;i<A.size();++i) out[i]=MinkowskiSum(A[i], B[i], true);
        return out;
    }
    size_t pairCount = A.size();
    std::vector<long long> ax, ay, bx, by;
    std::vector<GPUPath> ainfo(pairCount), binfo(pairCount);
    for(size_t i=0;i<pairCount;++i){
        ainfo[i] = { (int)ax.size(), (int)A[i].size() };
        for(auto pt:A[i]){ ax.push_back(pt.x); ay.push_back(pt.y); }
    }
    for(size_t i=0;i<pairCount;++i){
        binfo[i] = { (int)bx.size(), (int)B[i].size() };
        for(auto pt:B[i]){ bx.push_back(pt.x); by.push_back(pt.y); }
    }
    std::vector<GPUPair> pairs(pairCount);
    size_t out_sz = 0;
    for(size_t i=0;i<pairCount;++i){
        pairs[i] = { ainfo[i], binfo[i], (int)out_sz };
        out_sz += (size_t)ainfo[i].size * (size_t)binfo[i].size;
    }
    long long *d_ax,*d_ay,*d_bx,*d_by,*d_outx,*d_outy; GPUPair* d_pairs;
    cudaMalloc(&d_ax, ax.size()*sizeof(long long));
    cudaMalloc(&d_ay, ay.size()*sizeof(long long));
    cudaMalloc(&d_bx, bx.size()*sizeof(long long));
    cudaMalloc(&d_by, by.size()*sizeof(long long));
    cudaMemcpy(d_ax, ax.data(), ax.size()*sizeof(long long), cudaMemcpyHostToDevice);
    cudaMemcpy(d_ay, ay.data(), ay.size()*sizeof(long long), cudaMemcpyHostToDevice);
    cudaMemcpy(d_bx, bx.data(), bx.size()*sizeof(long long), cudaMemcpyHostToDevice);
    cudaMemcpy(d_by, by.data(), by.size()*sizeof(long long), cudaMemcpyHostToDevice);
    cudaMalloc(&d_pairs, pairs.size()*sizeof(GPUPair));
    cudaMemcpy(d_pairs, pairs.data(), pairs.size()*sizeof(GPUPair), cudaMemcpyHostToDevice);
    cudaMalloc(&d_outx, out_sz*sizeof(long long));
    cudaMalloc(&d_outy, out_sz*sizeof(long long));

    minkowskiKernelLauncher(d_ax,d_ay,d_bx,d_by,d_pairs,(int)pairCount,d_outx,d_outy);

    std::vector<long long> outx(out_sz), outy(out_sz);
    cudaMemcpy(outx.data(), d_outx, out_sz*sizeof(long long), cudaMemcpyDeviceToHost);
    cudaMemcpy(outy.data(), d_outy, out_sz*sizeof(long long), cudaMemcpyDeviceToHost);
    cudaFree(d_ax); cudaFree(d_ay); cudaFree(d_bx); cudaFree(d_by);
    cudaFree(d_pairs); cudaFree(d_outx); cudaFree(d_outy);

    std::vector<Paths64> res(pairCount);
    for(size_t i=0;i<pairCount;++i){
        size_t count = (size_t)ainfo[i].size * (size_t)binfo[i].size;
        std::vector<Point64> tmp; tmp.reserve(count);
        for(size_t j=0;j<count;++j){
            tmp.push_back({ outx[pairs[i].start + j], outy[pairs[i].start + j] });
        }
        Path64 hull = convexHull(tmp);
        res[i] = { hull };
    }
    return res;
#else
    std::vector<Paths64> out(A.size());
    for(size_t i=0;i<A.size();++i) out[i]=MinkowskiSum(A[i], B[i], true);
    return out;
#endif
}
