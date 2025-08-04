#include <cuda_runtime.h>
#include <cassert>
#include <stdexcept>
#include <stdio.h>
#include <cstdlib>
#include "geometry.h"

#define DEBUG_LIMIT 20

#define ASSERT_MSG(cond, msg) do { \
    if(!(cond)){ \
        printf("ASSERT FAIL: %s\n", msg); \
        assert(cond); \
    } \
} while(0)

// ---------------------------------------------------------------------------------
// Safety macro for CUDA API calls
// ---------------------------------------------------------------------------------
#ifndef CUDA_CHECK
#define CUDA_CHECK(x)                                                        \
    do {                                                                    \
        cudaError_t err = (x);                                              \
        if (err != cudaSuccess) {                                           \
            fprintf(stderr,"[ERR] %s failed: %s\n", #x,                    \
                    cudaGetErrorString(err));                               \
            exit(1);                                                        \
        }                                                                   \
    } while(0)
#endif

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
            fprintf(stderr,"[FATAL] Нельзя использовать NULL pointer: %s == 0x0\n",#ptr);
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

// ---------------------------------------------------------------------------------
// Simple geometry structures used on GPU
// ---------------------------------------------------------------------------------
struct GPUPath { int start; int size; };
struct GPUShape { int start; int size; };
struct Pt { long long x; long long y; };

// ---------------------------------------------------------------------------------
// Device utility functions
// ---------------------------------------------------------------------------------
__device__ long long cross(Pt a, Pt b, Pt c){
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

__device__ bool onSeg(Pt a, Pt b, Pt c){
    return b.x >= min(a.x, c.x) && b.x <= max(a.x, c.x) &&
           b.y >= min(a.y, c.y) && b.y <= max(a.y, c.y);
}

__device__ bool segInt(Pt p1, Pt q1, Pt p2, Pt q2){
    long long o1 = cross(p1, q1, p2);
    long long o2 = cross(p1, q1, q2);
    long long o3 = cross(p2, q2, p1);
    long long o4 = cross(p2, q2, q1);
    if(o1 == 0 && onSeg(p1, p2, q1)) return true;
    if(o2 == 0 && onSeg(p1, q2, q1)) return true;
    if(o3 == 0 && onSeg(p2, p1, q2)) return true;
    if(o4 == 0 && onSeg(p2, q1, q2)) return true;
    return ((o1 > 0) != (o2 > 0)) && ((o3 > 0) != (o4 > 0));
}

__device__ bool pointInPoly(Pt p,const long long* xs,const long long* ys,GPUPath poly){
    bool inside = false; int j = poly.size - 1;
    for(int i=0;i<poly.size;++i){
        long long xi = xs[poly.start + i];
        long long yi = ys[poly.start + i];
        long long xj = xs[poly.start + j];
        long long yj = ys[poly.start + j];
        if(yj == yi){ j = i; continue; } // avoid division by zero
        bool inter = ((yi > p.y) != (yj > p.y)) &&
            (p.x < (double)(xj - xi) * (p.y - yi) / (double)(yj - yi) + xi);
        if(inter) inside = !inside; j = i;
    }
    return inside;
}

__device__ bool polyOverlap(const long long* xs,const long long* ys,GPUPath a,GPUPath b){
    for(int i=0;i<a.size;++i){
        Pt a1{ xs[a.start + i], ys[a.start + i] };
        Pt a2{ xs[a.start + (i+1)%a.size], ys[a.start + (i+1)%a.size] };
        for(int j=0;j<b.size;++j){
            Pt b1{ xs[b.start + j], ys[b.start + j] };
            Pt b2{ xs[b.start + (j+1)%b.size], ys[b.start + (j+1)%b.size] };
            if(segInt(a1,a2,b1,b2)) return true;
        }
    }
    Pt p{ xs[a.start], ys[a.start] };
    if(pointInPoly(p,xs,ys,b)) return true;
    p.x = xs[b.start]; p.y = ys[b.start];
    if(pointInPoly(p,xs,ys,a)) return true;
    return false;
}

// ---------------------------------------------------------------------------------
// Kernel performing flat edge-to-edge intersection checks.
// Each thread processes exactly one pair of edges. Results are written atomically
// into res[shapeIdx]. `offsets` array contains prefix sums of edge pair counts
// for each shape in `shapes`.
// ---------------------------------------------------------------------------------
__global__ void overlapEdgesKernel(const long long* xs,const long long* ys,
                                  const GPUPath* paths,GPUShape cand,
                                  const GPUShape* shapes,int numShapes,
                                  const size_t* offsets,int* res){
    assert(offsets && "offsets pointer is null"); // <FIX offsets>
    size_t gid = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = offsets[numShapes];
    if(gid >= total) return;

    if(gid < DEBUG_LIMIT){
        printf("[kstart] tid=%zu total=%zu offsets=%p shapes=%p\n",
               gid, total, offsets, shapes);
    }

    // determine shape index via binary search
    int l=0, r=numShapes;
    while(l < r){
        int m = (l + r) / 2;
        if(gid >= offsets[m+1]) l = m + 1; else r = m;
    }
    int shapeIdx = l;
    if(!(shapeIdx >= 0 && shapeIdx < numShapes)){
        printf("[kassert] invalid shapeIdx=%d gid=%zu numShapes=%d\n", shapeIdx, gid, numShapes);
        assert(shapeIdx >= 0 && shapeIdx < numShapes);
    }

    if(gid < DEBUG_LIMIT)
        printf("[dbg] tid=%zu shapeIdx=%d off=%zu offNext=%zu\n", gid, shapeIdx,
               offsets[shapeIdx], offsets[shapeIdx+1]);

    if(!(gid >= offsets[shapeIdx] && gid < offsets[shapeIdx+1])){
        printf("[kassert] gid=%zu outside [%zu,%zu) for shapeIdx=%d\n",
               gid, offsets[shapeIdx], offsets[shapeIdx+1], shapeIdx);
        assert(gid >= offsets[shapeIdx] && gid < offsets[shapeIdx+1]);
    }

    size_t local = gid - offsets[shapeIdx];
    GPUShape sh = shapes[shapeIdx];

    if(gid < DEBUG_LIMIT)
        printf("[dbg] tid=%zu local=%zu\n", gid, local);

    for(int ia=0; ia<cand.size; ++ia){
        GPUPath pa = paths[cand.start + ia];
        if(pa.size < 2) continue;
        for(int ib=0; ib<sh.size; ++ib){
            GPUPath pb = paths[sh.start + ib];
            if(pb.size < 2) continue;
            size_t cnt = (size_t)pa.size * (size_t)pb.size;
            if(gid < DEBUG_LIMIT)
                printf("[dbg] tid=%zu ia=%d ib=%d local=%zu cnt=%zu\n", gid, ia, ib, local, cnt);
            if(local < cnt){
                int ea = (int)(local / pb.size);
                int eb = (int)(local % pb.size);

                if(!(ea >= 0 && ea < pa.size)){
                    printf("[kassert] ea out of range ea=%d pa.size=%d\n", ea, pa.size);
                    assert(ea >= 0 && ea < pa.size);
                }
                if(!(eb >= 0 && eb < pb.size)){
                    printf("[kassert] eb out of range eb=%d pb.size=%d\n", eb, pb.size);
                    assert(eb >= 0 && eb < pb.size);
                }

                if(gid < DEBUG_LIMIT){
                    printf("[dbg] tid=%zu -> ea=%d eb=%d pa.start=%d pb.start=%d\n",
                           gid, ea, eb, pa.start, pb.start);
                }
                Pt a1{ xs[pa.start + ea], ys[pa.start + ea] };
                Pt a2{ xs[pa.start + (ea+1)%pa.size], ys[pa.start + (ea+1)%pa.size] };
                Pt b1{ xs[pb.start + eb], ys[pb.start + eb] };
                Pt b2{ xs[pb.start + (eb+1)%pb.size], ys[pb.start + (eb+1)%pb.size] };
                if(gid < DEBUG_LIMIT){
                    printf("[coords] a1=(%lld,%lld) a2=(%lld,%lld) b1=(%lld,%lld) b2=(%lld,%lld)\n",
                           a1.x,a1.y,a2.x,a2.y,b1.x,b1.y,b2.x,b2.y);
                }
                bool hit = segInt(a1,a2,b1,b2);
                if(hit) atomicExch(&res[shapeIdx], 1);
                if(shapeIdx < 10 && gid < 10){
                    printf("edgePair dbg pair=%d tid=%zu ea=%d eb=%d hit=%d\n",
                           shapeIdx,gid,ea,eb,(int)hit);
                }
                if(gid < DEBUG_LIMIT)
                    printf("[dbg] tid=%zu done hit=%d\n", gid, (int)hit);
                return;
            }
            local -= cnt;
        }
    }
    if(gid < DEBUG_LIMIT)
        printf("[dbg] tid=%zu no-pair\n", gid);
}

// ---------------------------------------------------------------------------------
// Secondary kernel: for all pairs still marked as 0 after the edge kernel,
// run full polygon-overlap check (handles containment cases).
// ---------------------------------------------------------------------------------
__global__ void overlapInsideKernel(const long long* xs,const long long* ys,
                                   const GPUPath* paths,GPUShape cand,
                                   const GPUShape* shapes,int numShapes,int* res){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < DEBUG_LIMIT){
        printf("[kinside-start] idx=%d shapes=%p res=%p\n", idx, shapes, res);
    }
    if(idx >= numShapes){
        printf("[kinside] idx %d >= numShapes %d\n", idx, numShapes);
        return;
    }
    if(res[idx]) return; // already found intersection

    bool ov = false;
    GPUShape sh = shapes[idx];
    for(int ia=0; ia<cand.size && !ov; ++ia){
        GPUPath pa = paths[cand.start + ia];
        for(int ib=0; ib<sh.size && !ov; ++ib){
            GPUPath pb = paths[sh.start + ib];
            if(polyOverlap(xs,ys,pa,pb)) ov = true;
        }
    }
    if(ov) res[idx] = 1;
    if(idx < 10 && threadIdx.x == 0 && blockIdx.x == 0){
        printf("inside dbg pair=%d ov=%d\n", idx, (int)ov);
    }
}

// ---------------------------------------------------------------------------------
// Host launcher. Computes prefix sums for edge pairs and launches both kernels.
// ---------------------------------------------------------------------------------
#ifdef __cplusplus
extern "C" {
#endif

void overlapKernelLauncher(const long long* d_xs,const long long* d_ys,
                           const GPUPath* d_paths,GPUShape d_cand,
                           const GPUShape* d_shapes,int n, const size_t* d_offset, int* d_out){
    // bring shapes and paths back to host to compute prefix sums
    printf("[launcher] d_xs=%p d_ys=%p d_paths=%p d_shapes=%p d_out=%p n=%d\n",
           d_xs, d_ys, d_paths, d_shapes, d_out, n);
    printf("[launcher] cand start=%d size=%d\n", d_cand.start, d_cand.size);
    CHECK_PTR(d_xs); CHECK_PTR(d_ys); CHECK_PTR(d_paths); CHECK_PTR(d_shapes); CHECK_PTR(d_out);
    ASSERT_MSG(n > 0, "n must be >0");
    std::vector<GPUShape> h_shapes(n);
    CHECK_PTR(d_shapes);
    CUDA_CHECK(cudaMemcpy(h_shapes.data(), d_shapes, n*sizeof(GPUShape), cudaMemcpyDeviceToHost));
    printf("[launcher] copied shapes to host ptr=%p size=%zu\n", h_shapes.data(), h_shapes.size());

    int maxPath = d_cand.start + d_cand.size;
    for(int i=0;i<n;++i){
        int end = h_shapes[i].start + h_shapes[i].size;
        if(end > maxPath) maxPath = end;
    }
    std::vector<GPUPath> h_paths(maxPath);
    CHECK_PTR(d_paths);
    CUDA_CHECK(cudaMemcpy(h_paths.data(), d_paths, maxPath*sizeof(GPUPath), cudaMemcpyDeviceToHost));
    printf("[launcher] copied %d paths to host ptr=%p\n", maxPath, h_paths.data());

    std::vector<size_t> h_off(n+1,0);
    printf("[prefix] host offsets ptr=%p\n", h_off.data());
    size_t total = 0;
    for(int i=0;i<n;++i){
        size_t cnt = 0;
        for(int ia=0; ia<d_cand.size; ++ia){
            GPUPath pa = h_paths[d_cand.start + ia];
            for(int ib=0; ib<h_shapes[i].size; ++ib){
                GPUPath pb = h_paths[h_shapes[i].start + ib];
                cnt += (size_t)pa.size * (size_t)pb.size;
            }
        }
        h_off[i] = total;
        if (cnt == 0) continue; // <FIX offsets>
        total += cnt;
        if(i < DEBUG_LIMIT)
            printf("[prefix] i=%d addr=%p val=%zu cnt=%zu\n", i, &h_off[i], h_off[i], cnt);
        if(i > 0 && h_off[i] < h_off[i-1]){
            printf("[prefix-error] monotonic fail at %d prev=%zu curr=%zu\n", i, h_off[i-1], h_off[i]);
            assert(h_off[i] >= h_off[i-1]);
        }
    }
    h_off[n] = total;
    ASSERT_MSG(total >= h_off[n-1], "total prefix invalid");
    if(n < DEBUG_LIMIT)
        printf("[prefix] i=%d addr=%p val=%zu (total)\n", n, &h_off[n], total);

    // device memory for prefix sums
    size_t* d_off = nullptr;
    CUDA_CHECK(cudaMalloc(&d_off, (n+1)*sizeof(size_t))); // <FIX offsets>
    CHECK_ALLOC(d_off);
    printf("[malloc] d_off=%p bytes=%zu\n", d_off, (n+1)*sizeof(size_t));
    CHECK_PTR(d_off);
    CUDA_CHECK(cudaMemcpy(d_off, h_off.data(), (n+1)*sizeof(size_t),
                          cudaMemcpyHostToDevice));
    printf("[memcpy] host->dev d_off=%p from=%p bytes=%zu\n",
           d_off, h_off.data(), (n+1)*sizeof(size_t));

    // verify that the copy succeeded by reading the array back
    std::vector<size_t> h_check(n+1, 0);
    CHECK_PTR(d_off);
    CUDA_CHECK(cudaMemcpy(h_check.data(), d_off, (n+1)*sizeof(size_t),
                          cudaMemcpyDeviceToHost));
    for(int i=0;i<=n;++i){
        if(i < DEBUG_LIMIT)
            printf("[check] off[%d]=%zu\n", i, h_check[i]);
        if(h_check[i] != h_off[i]){
            fprintf(stderr,"prefix mismatch at %d: host=%zu dev=%zu\n",
                    i,h_off[i],h_check[i]);
        }
        assert(h_check[i] == h_off[i]);
    }

    // zero result array
    printf("[memset] d_out=%p bytes=%zu\n", d_out, n*sizeof(int));
    CHECK_PTR(d_out);
    CUDA_CHECK(cudaMemset(d_out, 0, n*sizeof(int)));
    std::vector<int> zeroCheck(n, -1);
    CHECK_PTR(d_out);
    CUDA_CHECK(cudaMemcpy(zeroCheck.data(), d_out, n*sizeof(int), cudaMemcpyDeviceToHost));
    for(int i=0;i<n && i<DEBUG_LIMIT; ++i){
        printf("[memset-check] i=%d val=%d\n", i, zeroCheck[i]);
        ASSERT_MSG(zeroCheck[i] == 0, "memset failed");
    }

    int t = 256;
    int b = (total + t - 1) / t;
    ASSERT_MSG(t > 0 && b > 0, "invalid launch config");
    printf("[launch edge] blocks=%d threads=%d total=%zu\n", b, t, total);
    printf("[launch edge] d_off=%p d_out=%p\n", d_off, d_out);
    FATAL_KERNEL_NULL(d_xs);
    FATAL_KERNEL_NULL(d_ys);
    FATAL_KERNEL_NULL(d_paths);
    FATAL_KERNEL_NULL(d_shapes);
    FATAL_KERNEL_NULL(d_off);
    FATAL_KERNEL_NULL(d_out);
    overlapEdgesKernel<<<b,t>>>(d_xs,d_ys,d_paths,d_cand,d_shapes,n,d_off,d_out);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<int> afterEdge(n);
    CHECK_PTR(d_out);
    CUDA_CHECK(cudaMemcpy(afterEdge.data(), d_out, n*sizeof(int), cudaMemcpyDeviceToHost));
    for(int i=0;i<n && i<DEBUG_LIMIT; ++i)
        printf("[after edge] i=%d val=%d\n", i, afterEdge[i]);

    int t2 = 128;
    int b2 = (n + t2 - 1) / t2;
    printf("[launch inside] blocks=%d threads=%d n=%d\n", b2, t2, n);
    FATAL_KERNEL_NULL(d_xs);
    FATAL_KERNEL_NULL(d_ys);
    FATAL_KERNEL_NULL(d_paths);
    FATAL_KERNEL_NULL(d_shapes);
    FATAL_KERNEL_NULL(d_out);
    overlapInsideKernel<<<b2,t2>>>(d_xs,d_ys,d_paths,d_cand,d_shapes,n,d_out);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    std::vector<int> finalRes(n);
    CHECK_PTR(d_out);
    CUDA_CHECK(cudaMemcpy(finalRes.data(), d_out, n*sizeof(int), cudaMemcpyDeviceToHost));
    for(int i=0;i<n && i<DEBUG_LIMIT; ++i)
        printf("[final result] i=%d val=%d\n", i, finalRes[i]);

    printf("[free] d_off=%p\n", d_off);
    if(d_off){
        CUDA_CHECK(cudaFree(d_off));
        d_off = nullptr;
    }
}

#ifdef __cplusplus
}
#endif

// ---------------------------------------------------------------------------------
// Quick reference:
// - Modify `t` above to change threads per block for the edge kernel.
// - If CUDA errors are printed, check parameters and GPU memory availability.
// ---------------------------------------------------------------------------------
