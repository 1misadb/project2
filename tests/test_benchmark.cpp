#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include "geometry.h"

TEST_CASE("overlap bvh benchmark") {
    Paths64 a = { { {0,0},{1000,0},{1000,1000},{0,1000} } };
    Paths64 b = { { {500,500},{1500,500},{1500,1500},{500,1500} } };
    auto ta = buildBVH(a);
    auto tb = buildBVH(b);
    BENCHMARK("overlapBVH") {
        return overlapBVH(ta, a, tb, b);
    };
}
