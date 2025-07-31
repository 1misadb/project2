#pragma once
#include "clipper3/clipper.h"
#include <vector>
using namespace Clipper2Lib;

extern double gKerf;
extern double gGap;

Paths64 movePaths(const Paths64& src, int64_t dx, int64_t dy);
bool overlap(const Paths64& a, const Paths64& b);
bool cuda_available();
std::vector<bool> overlapBatchGPU(const Paths64& cand,
                                  const std::vector<Paths64>& others);

struct BVHNode {
    Rect64 box;
    int left = -1;
    int right = -1;
    int index = -1;
};

std::vector<BVHNode> buildBVH(const Paths64& paths);
bool overlapBVH(const std::vector<BVHNode>& treeA, const Paths64& pa,
                const std::vector<BVHNode>& treeB, const Paths64& pb);

// --- convex utilities and GPU Minkowski ---
bool isConvex(const Path64& p);
Path64 convexHull(const std::vector<Point64>& pts);
std::vector<Paths64> minkowskiBatchGPU(const std::vector<Path64>& A,
                                       const std::vector<Path64>& B);

