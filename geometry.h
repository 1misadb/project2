#pragma once
#include "clipper3/clipper.h"
using namespace Clipper2Lib;

extern double gKerf;
extern double gGap;

Paths64 movePaths(const Paths64& src, int64_t dx, int64_t dy);
bool overlap(const Paths64& a, const Paths64& b);

struct BVHNode {
    Rect64 box;
    int left = -1;
    int right = -1;
    int index = -1;
};

std::vector<BVHNode> buildBVH(const Paths64& paths);
bool overlapBVH(const std::vector<BVHNode>& treeA, const Paths64& pa,
                const std::vector<BVHNode>& treeB, const Paths64& pb);

