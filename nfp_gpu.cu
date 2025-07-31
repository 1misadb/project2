#include <cuda_runtime.h>
#include "geometry.h"

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
    int threads = 128;
    minkowskiPairsKernel<<<pairCount, threads>>>(ax,ay,bx,by,pairs,pairCount,outx,outy);
    cudaDeviceSynchronize();
}
#ifdef __cplusplus
}
#endif
