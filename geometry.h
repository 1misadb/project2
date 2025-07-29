#pragma once
#include "clipper3/clipper.h"
#include "bvh.h"
using namespace Clipper2Lib;

extern double gKerf;
extern double gGap;

Paths64 movePaths(const Paths64& src, int64_t dx, int64_t dy);
bool overlap(const Paths64& a, const Paths64& b);
BVH buildBVH(const Paths64& p);
bool overlapBVH(const BVH& A,const BVH& B,const Paths64& pA,const Paths64& pB);

