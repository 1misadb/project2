#include <cuda_runtime.h>
#include "geometry.h"
#include <vector>
#include <limits>

struct Int2LL { long long x; long long y; };

__device__ long long cross_ll(Int2LL a, Int2LL b, Int2LL c){
    return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}

__device__ bool segments_intersect(Int2LL a1, Int2LL a2, Int2LL b1, Int2LL b2){
    long long d1 = cross_ll(a1,a2,b1);
    long long d2 = cross_ll(a1,a2,b2);
    long long d3 = cross_ll(b1,b2,a1);
    long long d4 = cross_ll(b1,b2,a2);
    if(((d1>0 && d2<0)||(d1<0 && d2>0)) &&
       ((d3>0 && d4<0)||(d3<0 && d4>0))) return true;
    if(d1==0 && min(a1.x,a2.x)<=b1.x && b1.x<=max(a1.x,a2.x) &&
                min(a1.y,a2.y)<=b1.y && b1.y<=max(a1.y,a2.y)) return true;
    if(d2==0 && min(a1.x,a2.x)<=b2.x && b2.x<=max(a1.x,a2.x) &&
                min(a1.y,a2.y)<=b2.y && b2.y<=max(a1.y,a2.y)) return true;
    if(d3==0 && min(b1.x,b2.x)<=a1.x && a1.x<=max(b1.x,b2.x) &&
                min(b1.y,b2.y)<=a1.y && a1.y<=max(b1.y,b2.y)) return true;
    if(d4==0 && min(b1.x,b2.x)<=a2.x && a2.x<=max(b1.x,b2.x) &&
                min(b1.y,b2.y)<=a2.y && a2.y<=max(b1.y,b2.y)) return true;
    return false;
}

__device__ bool point_in_poly(const Int2LL* poly, int n, Int2LL p){
    bool c=false;
    for(int i=0,j=n-1;i<n;j=i++){
        Int2LL pi=poly[i], pj=poly[j];
        if(((pi.y>p.y)!=(pj.y>p.y)) &&
            (p.x < (double)(pj.x-pi.x)*(p.y-pi.y)/(double)(pj.y-pi.y)+pi.x))
            c=!c;
    }
    return c;
}

__device__ bool poly_overlap(const Int2LL* a, int na, const Int2LL* b, int nb){
    for(int i=0;i<na;i++){
        Int2LL a1=a[i];
        Int2LL a2=a[(i+1)%na];
        for(int j=0;j<nb;j++){
            Int2LL b1=b[j];
            Int2LL b2=b[(j+1)%nb];
            if(segments_intersect(a1,a2,b1,b2)) return true;
        }
    }
    if(point_in_poly(a,na,b[0])) return true;
    if(point_in_poly(b,nb,a[0])) return true;
    return false;
}

extern "C" __global__ void overlapKernel(const Int2LL* ptsA,const int* offsA,const int* lensA,
                                          const Int2LL* ptsB,const int* offsB,const int* lensB,
                                          unsigned char* res,int N){
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if(idx>=N) return;
    const Int2LL* pa = ptsA + offsA[idx];
    const Int2LL* pb = ptsB + offsB[idx];
    int na = lensA[idx];
    int nb = lensB[idx];
    res[idx] = poly_overlap(pa,na,pb,nb);
}

static void flatten(const std::vector<Paths64>& shapes,
                    std::vector<Int2LL>& pts,
                    std::vector<int>& offs,
                    std::vector<int>& lens){
    int offset=0;
    for(const auto& p: shapes){
        const Path64& poly = p.empty()? Path64{}: p[0];
        offs.push_back(offset);
        lens.push_back((int)poly.size());
        for(auto pt: poly){ pts.push_back({pt.x, pt.y}); }
        offset += poly.size();
    }
}

bool overlapBatchCUDA(const std::vector<Paths64>& A,const std::vector<Paths64>& B,std::vector<bool>& out){
    int N = (int)A.size();
    out.resize(N);
    if(N==0) return true;
    std::vector<Int2LL> hPtsA, hPtsB;
    std::vector<int> hOffA, hOffB, hLenA, hLenB;
    flatten(A,hPtsA,hOffA,hLenA);
    flatten(B,hPtsB,hOffB,hLenB);
    Int2LL *dPtsA,*dPtsB; int *dOffA,*dOffB,*dLenA,*dLenB; unsigned char *dRes;
    cudaMalloc(&dPtsA,hPtsA.size()*sizeof(Int2LL));
    cudaMalloc(&dPtsB,hPtsB.size()*sizeof(Int2LL));
    cudaMalloc(&dOffA,hOffA.size()*sizeof(int));
    cudaMalloc(&dOffB,hOffB.size()*sizeof(int));
    cudaMalloc(&dLenA,hLenA.size()*sizeof(int));
    cudaMalloc(&dLenB,hLenB.size()*sizeof(int));
    cudaMalloc(&dRes,N*sizeof(unsigned char));
    cudaMemcpy(dPtsA,hPtsA.data(),hPtsA.size()*sizeof(Int2LL),cudaMemcpyHostToDevice);
    cudaMemcpy(dPtsB,hPtsB.data(),hPtsB.size()*sizeof(Int2LL),cudaMemcpyHostToDevice);
    cudaMemcpy(dOffA,hOffA.data(),hOffA.size()*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(dOffB,hOffB.data(),hOffB.size()*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(dLenA,hLenA.data(),hLenA.size()*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(dLenB,hLenB.data(),hLenB.size()*sizeof(int),cudaMemcpyHostToDevice);
    int block=128; int grid=(N+block-1)/block;
    overlapKernel<<<grid,block>>>(dPtsA,dOffA,dLenA,dPtsB,dOffB,dLenB,dRes,N);
    std::vector<unsigned char> hRes(N);
    cudaMemcpy(hRes.data(),dRes,N*sizeof(unsigned char),cudaMemcpyDeviceToHost);
    out.resize(N);
    for(int i=0;i<N;++i) out[i]=hRes[i];
    cudaFree(dPtsA); cudaFree(dPtsB); cudaFree(dOffA); cudaFree(dOffB);
    cudaFree(dLenA); cudaFree(dLenB); cudaFree(dRes);
    return true;
}

bool cudaAvailable(){ return true; }
