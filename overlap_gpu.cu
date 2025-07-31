#include <cuda_runtime.h>
#include "geometry.h"

// --- Структуры ---
struct GPUPath { int start; int size; };
struct GPUShape { int start; int size; };
struct Pt { long long x; long long y; };

// --- device-функции ---
__device__ long long cross(Pt a, Pt b, Pt c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}
__device__ bool onSeg(Pt a, Pt b, Pt c) {
    return b.x >= min(a.x, c.x) && b.x <= max(a.x, c.x) &&
           b.y >= min(a.y, c.y) && b.y <= max(a.y, c.y);
}
__device__ bool segInt(Pt p1, Pt q1, Pt p2, Pt q2) {
    long long o1 = cross(p1, q1, p2);
    long long o2 = cross(p1, q1, q2);
    long long o3 = cross(p2, q2, p1);
    long long o4 = cross(p2, q2, q1);
    if (o1 == 0 && onSeg(p1, p2, q1)) return true;
    if (o2 == 0 && onSeg(p1, q2, q1)) return true;
    if (o3 == 0 && onSeg(p2, p1, q2)) return true;
    if (o4 == 0 && onSeg(p2, q1, q2)) return true;
    return ((o1 > 0) != (o2 > 0)) && ((o3 > 0) != (o4 > 0));
}
__device__ bool pointInPoly(Pt p, const long long* xs, const long long* ys, GPUPath poly) {
    bool inside = false; int j = poly.size - 1;
    for (int i = 0; i < poly.size; i++) {
        long long xi = xs[poly.start + i], yi = ys[poly.start + i];
        long long xj = xs[poly.start + j], yj = ys[poly.start + j];
        bool inter = ((yi > p.y) != (yj > p.y)) &&
            (p.x < (double)(xj - xi) * (p.y - yi) / (double)(yj - yi) + xi);
        if (inter) inside = !inside; j = i;
    }
    return inside;
}
__device__ bool polyOverlap(const long long* xs, const long long* ys, GPUPath a, GPUPath b) {
    for (int i = 0; i < a.size; i++) {
        Pt a1{ xs[a.start + i], ys[a.start + i] };
        Pt a2{ xs[a.start + (i + 1) % a.size], ys[a.start + (i + 1) % a.size] };
        for (int j = 0; j < b.size; j++) {
            Pt b1{ xs[b.start + j], ys[b.start + j] };
            Pt b2{ xs[b.start + (j + 1) % b.size], ys[b.start + (j + 1) % b.size] };
            if (segInt(a1, a2, b1, b2)) return true;
        }
    }
    Pt p{ xs[a.start], ys[a.start] };
    if (pointInPoly(p, xs, ys, b)) return true;
    p.x = xs[b.start]; p.y = ys[b.start];
    if (pointInPoly(p, xs, ys, a)) return true;
    return false;
}

// --- CUDA kernel (ВНЕ extern "C") ---
__global__ void overlapKernel(const long long* xs, const long long* ys,
    const GPUPath* paths, GPUShape cand,
    const GPUShape* shapes, int numShapes, bool* res) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= numShapes) return;
    bool ov = false; GPUShape sh = shapes[idx];
    for (int ia = 0; ia < cand.size && !ov; ++ia) {
        GPUPath pa = paths[cand.start + ia];
        for (int ib = 0; ib < sh.size && !ov; ++ib) {
            GPUPath pb = paths[sh.start + ib];
            if (polyOverlap(xs, ys, pa, pb)) ov = true;
        }
    }
    res[idx] = ov;
}

// --- Обёртка (Launcher) с extern "C" ---
#ifdef __cplusplus
extern "C" {
#endif

void overlapKernelLauncher(const long long* d_xs, const long long* d_ys, const GPUPath* d_paths,
    GPUShape d_cand, const GPUShape* d_shapes, int n, uint8_t* d_out)
{
    int threads = 128;
    int blocks = (n + threads - 1) / threads;
    overlapKernel<<<blocks, threads>>>(d_xs, d_ys, d_paths, d_cand, d_shapes, n, (bool*)d_out);
    cudaDeviceSynchronize();
}

#ifdef __cplusplus
}
#endif
