#pragma once

struct GPUPath {
    int start;
    int size;
};

struct GPUShape {
    int start;
    int size;
};

static_assert(sizeof(GPUPath) == 8 && alignof(GPUPath) == 4, "GPUPath layout mismatch");
static_assert(sizeof(GPUShape) == 8 && alignof(GPUShape) == 4, "GPUShape layout mismatch");

extern "C" void overlapKernelLauncher(const long long*, const long long*, const GPUPath*, GPUShape, const GPUShape*, int, int*);

