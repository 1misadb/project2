#include "geometry.h"
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif
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

#ifdef USE_CUDA
struct GPUPath { int start; int size; };
struct GPUShape { int start; int size; };
extern void overlapKernel(const long long* xs,const long long* ys,
                          const GPUPath* paths, GPUShape cand,
                          const GPUShape* shapes,int numShapes,bool* res);
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
    long long* d_xs; long long* d_ys; GPUPath* d_paths; GPUShape* d_shapes; bool* d_out;
    size_t pts_sz=xs.size()*sizeof(long long); cudaMalloc(&d_xs,pts_sz); cudaMalloc(&d_ys,pts_sz);
    cudaMemcpy(d_xs,xs.data(),pts_sz,cudaMemcpyHostToDevice);
    cudaMemcpy(d_ys,ys.data(),pts_sz,cudaMemcpyHostToDevice);
    cudaMalloc(&d_paths,paths.size()*sizeof(GPUPath));
    cudaMemcpy(d_paths,paths.data(),paths.size()*sizeof(GPUPath),cudaMemcpyHostToDevice);
    cudaMalloc(&d_shapes,shapes.size()*sizeof(GPUShape));
    cudaMemcpy(d_shapes,shapes.data(),shapes.size()*sizeof(GPUShape),cudaMemcpyHostToDevice);
    cudaMalloc(&d_out,shapes.size()*sizeof(bool));
    GPUShape d_cand=candShape; // copy by value
    int threads=128; int blocks=(shapes.size()+threads-1)/threads;
    overlapKernel<<<blocks,threads>>>(d_xs,d_ys,d_paths,d_cand,d_shapes,(int)shapes.size(),d_out);
    cudaDeviceSynchronize();
    std::vector<bool> res(others.size());
    cudaMemcpy(res.data(),d_out,others.size()*sizeof(bool),cudaMemcpyDeviceToHost);
    cudaFree(d_xs); cudaFree(d_ys); cudaFree(d_paths); cudaFree(d_shapes); cudaFree(d_out);
    return res;
#else
    std::vector<bool> r(others.size());
    for(size_t i=0;i<others.size();++i) r[i]=overlap(cand, others[i]);
    return r;
#endif
}

