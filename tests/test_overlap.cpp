#define NEST_UNIT_TEST
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "nesting_algorithm.cpp"

TEST_CASE("overlap basic") {
    Paths64 square1 = { { {0,0},{10,0},{10,10},{0,10} } };
    Paths64 square2 = { { {5,5},{15,5},{15,15},{5,15} } };
    REQUIRE(overlap(square1, square2) == true);
}
