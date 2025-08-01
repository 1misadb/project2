#include <cuda_runtime.h>
#include <stdio.h>
#include "geometry.h"

// ---------------------------------------------------------------------------------
// Safety macro for CUDA API calls
// ---------------------------------------------------------------------------------
#ifndef CUDA_CHECK
#define CUDA_CHECK(call)                                                   \
    do {                                                                  \
        cudaError_t err__ = (call);                                       \
        if (err__ != cudaSuccess) {                                       \
            fprintf(stderr,"CUDA error %s in %s at line %d\n",              \
                    cudaGetErrorString(err__), __FILE__, __LINE__);       \
        }                                                                 \
    } while (0)
#endif

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
    size_t gid = blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = offsets[numShapes];
    if(gid >= total) return;

    // determine shape index via binary search
    int l=0, r=numShapes;
    while(l < r){
        int m = (l + r) / 2;
        if(gid >= offsets[m+1]) l = m + 1; else r = m;
    }
    int shapeIdx = l;
    if(shapeIdx >= numShapes) return; // safety

    size_t local = gid - offsets[shapeIdx];
    GPUShape sh = shapes[shapeIdx];

    for(int ia=0; ia<cand.size; ++ia){
        GPUPath pa = paths[cand.start + ia];
        if(pa.size < 2) continue;
        for(int ib=0; ib<sh.size; ++ib){
            GPUPath pb = paths[sh.start + ib];
            if(pb.size < 2) continue;
            size_t cnt = (size_t)pa.size * (size_t)pb.size;
            if(local < cnt){
                int ea = local / pb.size;
                int eb = local % pb.size;
                Pt a1{ xs[pa.start + ea], ys[pa.start + ea] };
                Pt a2{ xs[pa.start + (ea+1)%pa.size], ys[pa.start + (ea+1)%pa.size] };
                Pt b1{ xs[pb.start + eb], ys[pb.start + eb] };
                Pt b2{ xs[pb.start + (eb+1)%pb.size], ys[pb.start + (eb+1)%pb.size] };
                bool hit = segInt(a1,a2,b1,b2);
                if(hit) atomicExch(&res[shapeIdx], 1);
                if(shapeIdx < 10 && gid < 10){
                    printf("edgePair dbg pair=%d tid=%zu ea=%d eb=%d hit=%d\n",
                           shapeIdx,gid,ea,eb,(int)hit);
                }
                return;
            }
            local -= cnt;
        }
    }
}

// ---------------------------------------------------------------------------------
// Secondary kernel: for all pairs still marked as 0 after the edge kernel,
// run full polygon-overlap check (handles containment cases).
// ---------------------------------------------------------------------------------
__global__ void overlapInsideKernel(const long long* xs,const long long* ys,
                                   const GPUPath* paths,GPUShape cand,
                                   const GPUShape* shapes,int numShapes,int* res){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= numShapes) return;
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
                           const GPUShape* d_shapes,int n,int* d_out){
    // bring shapes and paths back to host to compute prefix sums
    std::vector<GPUShape> h_shapes(n);
    CUDA_CHECK(cudaMemcpy(h_shapes.data(), d_shapes, n*sizeof(GPUShape), cudaMemcpyDeviceToHost));

    int maxPath = d_cand.start + d_cand.size;
    for(int i=0;i<n;++i){
        int end = h_shapes[i].start + h_shapes[i].size;
        if(end > maxPath) maxPath = end;
    }
    std::vector<GPUPath> h_paths(maxPath);
    CUDA_CHECK(cudaMemcpy(h_paths.data(), d_paths, maxPath*sizeof(GPUPath), cudaMemcpyDeviceToHost));

    std::vector<size_t> h_off(n+1,0);
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
        total += cnt;
    }
    h_off[n] = total;

    size_t* d_off = nullptr;
    CUDA_CHECK(cudaMalloc(&d_off, (n+1)*sizeof(size_t)));
    CUDA_CHECK(cudaMemcpy(d_off, h_off.data(), (n+1)*sizeof(size_t), cudaMemcpyHostToDevice));

    // zero result array
    CUDA_CHECK(cudaMemset(d_out, 0, n*sizeof(int)));

    int t = 256;
    int b = (total + t - 1) / t;
    overlapEdgesKernel<<<b,t>>>(d_xs,d_ys,d_paths,d_cand,d_shapes,n,d_off,d_out);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    int t2 = 128;
    int b2 = (n + t2 - 1) / t2;
    overlapInsideKernel<<<b2,t2>>>(d_xs,d_ys,d_paths,d_cand,d_shapes,n,d_out);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaFree(d_off));
}

#ifdef __cplusplus
}
#endif

// ---------------------------------------------------------------------------------
// Quick reference:
// - Modify `t` above to change threads per block for the edge kernel.
// - If CUDA errors are printed, check parameters and GPU memory availability.
// ---------------------------------------------------------------------------------
