#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "geometry.h"
#include "clipper3/clipper.minkowski.h"

TEST_CASE("minkowski gpu matches cpu") {
    Path64 a = { {0,0},{10,0},{10,10},{0,10} };
    Path64 b = { {0,0},{5,0},{5,5},{0,5} };
    std::vector<Path64> As{a};
    std::vector<Path64> Bs{b};
    auto cpu = MinkowskiSum(a, b, true);
    auto gpu = minkowskiBatchGPU(As, Bs);
    REQUIRE(gpu.size() == 1);
    REQUIRE(gpu[0][0].size() == cpu[0].size());
}
