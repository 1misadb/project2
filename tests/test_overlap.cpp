#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "geometry.h"

TEST_CASE("overlap basic") {
    Paths64 square1 = { { {0,0},{10,0},{10,10},{0,10} } };
    Paths64 square2 = { { {5,5},{15,5},{15,15},{5,15} } };
    REQUIRE(overlap(square1, square2) == true);
}

TEST_CASE("overlap none") {
    Paths64 a = { { {0,0},{10,0},{10,10},{0,10} } };
    Paths64 b = { { {20,0},{30,0},{30,10},{20,10} } };
    REQUIRE_FALSE(overlap(a,b));
}

TEST_CASE("overlap touch no kerf") {
    gKerf = 0.0; gGap = 0.0;
    Paths64 a = { { {0,0},{10,0},{10,10},{0,10} } };
    Paths64 b = { { {10,0},{20,0},{20,10},{10,10} } };
    REQUIRE_FALSE(overlap(a,b));
}

TEST_CASE("overlap touch with kerf") {
    gKerf = 1.0; gGap = 0.0;
    Paths64 a = { { {0,0},{10,0},{10,10},{0,10} } };
    Paths64 b = { { {10,0},{20,0},{20,10},{10,10} } };
    REQUIRE(overlap(a,b));
}

