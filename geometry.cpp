#include "geometry.h"
#ifdef USE_CUDA
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

#define CUDA_CHECK(x)                                                        \
    do {                                                                    \
        cudaError_t err = (x);                                              \
        if (err != cudaSuccess) {                                           \
            fprintf(stderr,"[ERR] %s failed: %s\n", #x,                    \
                    cudaGetErrorString(err));                               \
            exit(1);                                                        \
        }                                                                   \
    } while(0)

#define CHECK_ALLOC(ptr)                                                    \
    do {                                                                    \
        if ((ptr) == nullptr) {                                             \
            fprintf(stderr,"[FATAL] cudaMalloc вернул NULL для %s!\n",    \
                    #ptr);                                                 \
            exit(1);                                                        \
        }                                                                   \
    } while(0)

#define CHECK_PTR(ptr)                                                      \
    do {                                                                    \
        if ((ptr) == nullptr) {                                             \
            fprintf(stderr,"[FATAL] Нельзя использовать NULL pointer: %s == 0x0\n",#ptr);\
            exit(1);                                                        \
        }                                                                   \
    } while(0)

#define FATAL_KERNEL_NULL(ptr)                                              \
    do {                                                                    \
        if ((ptr) == nullptr) {                                             \
            fprintf(stderr,                                               \
                    "[FATAL] Нельзя вызывать ядро с NULL pointer: %s == 0x0\n",\
                    #ptr);                                                 \
            exit(1);                                                        \
        }                                                                   \
    } while(0)
#endif
#include <spdlog/spdlog.h>
#include <limits>
#include <algorithm>
#include <numeric>
#include <cstdint>

// --- внешние символы для Clipper2 (предполагается, что geometry.h тянет clipper3 headers)
using namespace Clipper2Lib;

double gKerf = 0.0;
double gGap = 0.0;
static constexpr double SCALE = 1e4;

// ---- Утилиты перемещения и пересечений ----
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
    int cnt=0; cudaError_t rc = cudaGetDeviceCount(&cnt);
    if(rc != cudaSuccess){
        fprintf(stderr,"[ERR] cudaGetDeviceCount failed: %s\n", cudaGetErrorString(rc));
        return false;
    }
    return cnt>0;
#else
    return false;
#endif
}

// -------------------- FIXED VERSION --------------------
std::vector<bool> overlapBatchGPU(const Paths64& cand,
                                  const std::vector<Paths64>& others){
#ifdef USE_CUDA
    // Нечего проверять — быстрый выход
    if (others.empty()) {
        spdlog::debug("overlapBatchGPU: no shapes to check");
        return {};
    }
    // CUDA недоступна — CPU путь
    if (!cuda_available()){
        spdlog::debug("overlapBatchGPU: CUDA unavailable, using CPU path");
        std::vector<bool> r(others.size());
        for(size_t i=0;i<others.size();++i) r[i]=overlap(cand, others[i]);
        return r;
    }

    spdlog::debug("overlapBatchGPU: cand paths={} others={} total_shapes={}",
                  cand.size(), others.size(), others.size()+1);

    // Flatten координаты и метаданные
    std::vector<long long> xs,ys; xs.reserve(1024); ys.reserve(1024);
    size_t totalPathCount = cand.size();
    for(const auto& sh:others) totalPathCount += sh.size();
    std::vector<GPUPath> paths; paths.reserve(totalPathCount);

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

    // FIX: если после подготовки нечего обрабатывать — не запускаем ядро
    if (shapes.empty()){
        spdlog::debug("overlapBatchGPU: shapes vector is empty, nothing to process");
        return std::vector<bool>(others.size(), false);
    }

    // --- Thread local device buffers reused across calls ---
    struct DevBuf {
        long long* xs = nullptr; size_t xs_cap = 0;
        long long* ys = nullptr; size_t ys_cap = 0;
        GPUPath* paths = nullptr; size_t paths_cap = 0;
        GPUShape* shapes = nullptr; size_t shapes_cap = 0;
        int* out = nullptr; size_t out_cap = 0;
        ~DevBuf(){
            if(xs) cudaFree(xs);
            if(ys) cudaFree(ys);
            if(paths) cudaFree(paths);
            if(shapes) cudaFree(shapes);
            if(out) cudaFree(out);
        }
    };
    static thread_local DevBuf buf;

    auto ensure = [](auto*& ptr,size_t& cap,size_t bytes){
        if(bytes==0) return;                 // важно: не аллоцировать 0 байт
        if(cap < bytes){
            if(ptr) CUDA_CHECK(cudaFree(ptr));
            CUDA_CHECK(cudaMalloc(&ptr, bytes));
            CHECK_ALLOC(ptr);
            cap = bytes;
        }
    };

    size_t pts_sz = xs.size()*sizeof(long long);
    ensure(buf.xs, buf.xs_cap, pts_sz);
    ensure(buf.ys, buf.ys_cap, pts_sz);
    CUDA_CHECK(cudaMemcpy(buf.xs, xs.data(), pts_sz, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(buf.ys, ys.data(), pts_sz, cudaMemcpyHostToDevice));

    ensure(buf.paths,  buf.paths_cap,  paths.size()*sizeof(GPUPath));
    CUDA_CHECK(cudaMemcpy(buf.paths, paths.data(), paths.size()*sizeof(GPUPath), cudaMemcpyHostToDevice));
    ensure(buf.shapes, buf.shapes_cap, shapes.size()*sizeof(GPUShape));
    CUDA_CHECK(cudaMemcpy(buf.shapes, shapes.data(), shapes.size()*sizeof(GPUShape), cudaMemcpyHostToDevice));
    if (shapes.empty()) {
        spdlog::debug("overlapBatchGPU: shapes empty, nothing to process");
        return std::vector<bool>(others.size(), false);
    }

    ensure(buf.out, buf.out_cap, shapes.size() > 0 ? shapes.size()*sizeof(int) : sizeof(int));
    CHECK_ALLOC(buf.out);
    
    ensure(buf.out,    buf.out_cap,    shapes.size()*sizeof(int));
    CHECK_ALLOC(buf.out); // дополнительная гарантия до запуска ядра

    GPUShape d_cand = candShape; // copy by value

    // Запуск GPU
    printf("[call launcher] shapes=%zu cand.size=%zu\n", shapes.size(), cand.size());
    FATAL_KERNEL_NULL(buf.xs); FATAL_KERNEL_NULL(buf.ys); FATAL_KERNEL_NULL(buf.paths);
    FATAL_KERNEL_NULL(buf.shapes); FATAL_KERNEL_NULL(buf.out);

    overlapKernelLauncher(buf.xs, buf.ys, buf.paths, d_cand, buf.shapes, (int)shapes.size(), buf.out);

    std::vector<int> res_raw(others.size());
    CUDA_CHECK(cudaMemcpy(res_raw.data(), buf.out, others.size()*sizeof(int), cudaMemcpyDeviceToHost));
    for(size_t i=0;i<others.size() && i<10;++i)
        printf("[host result] i=%zu val=%d\n", i, res_raw[i]);

    std::vector<bool> res(others.size());
    for(size_t i=0;i<others.size();++i) res[i] = (res_raw[i]!=0);

    // buffers intentionally kept for reuse
    return res;
#else
    std::vector<bool> r(others.size());
    for(size_t i=0;i<others.size();++i) r[i]=overlap(cand, others[i]);
    return r;
#endif
}
// ------------------ /FIXED VERSION ---------------------

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
        spdlog::debug("minkowskiBatchGPU: CUDA unavailable, using CPU path");
        std::vector<Paths64> out(A.size());
        for(size_t i=0;i<A.size();++i) out[i]=MinkowskiSum(A[i], B[i], true);
        return out;
    }
    size_t pairCount = A.size();
    spdlog::debug("minkowskiBatchGPU: pairCount={}", pairCount);
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
    CUDA_CHECK(cudaMalloc(&d_ax, ax.size()*sizeof(long long))); CHECK_ALLOC(d_ax);
    CUDA_CHECK(cudaMalloc(&d_ay, ay.size()*sizeof(long long))); CHECK_ALLOC(d_ay);
    CUDA_CHECK(cudaMalloc(&d_bx, bx.size()*sizeof(long long))); CHECK_ALLOC(d_bx);
    CUDA_CHECK(cudaMalloc(&d_by, by.size()*sizeof(long long))); CHECK_ALLOC(d_by);
    CHECK_PTR(d_ax); CUDA_CHECK(cudaMemcpy(d_ax, ax.data(), ax.size()*sizeof(long long), cudaMemcpyHostToDevice));
    CHECK_PTR(d_ay); CUDA_CHECK(cudaMemcpy(d_ay, ay.data(), ay.size()*sizeof(long long), cudaMemcpyHostToDevice));
    CHECK_PTR(d_bx); CUDA_CHECK(cudaMemcpy(d_bx, bx.data(), bx.size()*sizeof(long long), cudaMemcpyHostToDevice));
    CHECK_PTR(d_by); CUDA_CHECK(cudaMemcpy(d_by, by.data(), by.size()*sizeof(long long), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(&d_pairs, pairs.size()*sizeof(GPUPair))); CHECK_ALLOC(d_pairs);
    CHECK_PTR(d_pairs); CUDA_CHECK(cudaMemcpy(d_pairs, pairs.data(), pairs.size()*sizeof(GPUPair), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMalloc(&d_outx, out_sz*sizeof(long long))); CHECK_ALLOC(d_outx);
    CUDA_CHECK(cudaMalloc(&d_outy, out_sz*sizeof(long long))); CHECK_ALLOC(d_outy);

    FATAL_KERNEL_NULL(d_ax); FATAL_KERNEL_NULL(d_ay); FATAL_KERNEL_NULL(d_bx);
    FATAL_KERNEL_NULL(d_by); FATAL_KERNEL_NULL(d_pairs); FATAL_KERNEL_NULL(d_outx);
    FATAL_KERNEL_NULL(d_outy);
    minkowskiKernelLauncher(d_ax,d_ay,d_bx,d_by,d_pairs,(int)pairCount,d_outx,d_outy);

    std::vector<long long> outx(out_sz), outy(out_sz);
    CHECK_PTR(d_outx); CUDA_CHECK(cudaMemcpy(outx.data(), d_outx, out_sz*sizeof(long long), cudaMemcpyDeviceToHost));
    CHECK_PTR(d_outy); CUDA_CHECK(cudaMemcpy(outy.data(), d_outy, out_sz*sizeof(long long), cudaMemcpyDeviceToHost));
    if(d_ax){ CUDA_CHECK(cudaFree(d_ax)); d_ax=nullptr; }
    if(d_ay){ CUDA_CHECK(cudaFree(d_ay)); d_ay=nullptr; }
    if(d_bx){ CUDA_CHECK(cudaFree(d_bx)); d_bx=nullptr; }
    if(d_by){ CUDA_CHECK(cudaFree(d_by)); d_by=nullptr; }
    if(d_pairs){ CUDA_CHECK(cudaFree(d_pairs)); d_pairs=nullptr; }
    if(d_outx){ CUDA_CHECK(cudaFree(d_outx)); d_outx=nullptr; }
    if(d_outy){ CUDA_CHECK(cudaFree(d_outy)); d_outy=nullptr; }

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
