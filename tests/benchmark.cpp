#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "geometry.h"
#include "lru_cache.h"

TEST_CASE("NFP LRU") {
    LRUCache<int,int> cache(2);
    cache.put(1,1);cache.put(2,2);
    int v; REQUIRE(cache.get(1,v));
    cache.put(3,3);
    REQUIRE_FALSE(cache.get(2,v));
}

