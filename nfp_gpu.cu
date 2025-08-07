#include <cuda_runtime.h>
#include "geometry.h"
#include <cstdio>
#include <cstdlib>

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

#define FATAL_KERNEL_NULL(ptr)                                              \
    do {                                                                    \
        if ((ptr) == nullptr) {                                             \
            fprintf(stderr,                                               \
                    "[FATAL] Нельзя вызывать ядро с NULL pointer: %s == 0x0\n",\
                    #ptr);                                                 \
            exit(1);                                                        \
        }                                                                   \
    } while(0)

struct GPUPath { int start; int size; };
struct GPUPair { GPUPath a; GPUPath b; int start; };

__global__ void minkowskiPairsKernel(const long long* ax,const long long* ay,
                                     const long long* bx,const long long* by,
                                     const GPUPair* pairs,int pairCount,
                                     long long* outx,long long* outy){
    int pairIdx = blockIdx.x;
    if(pairIdx >= pairCount) return;
    GPUPair p = pairs[pairIdx];
    for(int i = threadIdx.x; i < p.a.size; i += blockDim.x){
        long long axv = ax[p.a.start + i];
        long long ayv = ay[p.a.start + i];
        for(int j=0;j<p.b.size;++j){
            long long bxv = bx[p.b.start + j];
            long long byv = by[p.b.start + j];
            int idx = p.start + i*p.b.size + j;
            outx[idx] = axv + bxv;
            outy[idx] = ayv + byv;
        }
    }
}

#ifdef __cplusplus
extern "C" {
#endif
void minkowskiKernelLauncher(const long long* ax,const long long* ay,
                             const long long* bx,const long long* by,
                             const GPUPair* pairs,int pairCount,
                             long long* outx,long long* outy){
    if(pairCount <= 0){
        fprintf(stderr, "[WARN] minkowskiKernelLauncher called with pairCount=%d\n", pairCount);
        return;
    }
    FATAL_KERNEL_NULL(ax); FATAL_KERNEL_NULL(ay); FATAL_KERNEL_NULL(bx);
    FATAL_KERNEL_NULL(by); FATAL_KERNEL_NULL(pairs); FATAL_KERNEL_NULL(outx);
    FATAL_KERNEL_NULL(outy);
    int threads = 128;
    minkowskiPairsKernel<<<pairCount, threads>>>(ax,ay,bx,by,pairs,pairCount,outx,outy);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
}
#ifdef __cplusplus
}
#endif
