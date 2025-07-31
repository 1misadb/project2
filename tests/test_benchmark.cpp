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

TEST_CASE("cuda batch vs cpu") {
    Paths64 base = { { {0,0},{100,0},{100,100},{0,100} } };
    std::vector<Paths64> others;
    for(int i=0;i<10000;i++) others.push_back(movePaths(base,i*50,i*50));
    if(!cuda_available()) {
        SUCCEED("cuda not available");
        return;
    }
    BENCHMARK("cpu batch") {
        std::vector<bool> r(others.size());
        for(size_t i=0;i<others.size();++i) r[i]=overlap(base,others[i]);
        return r[0];
    };
    BENCHMARK("cuda batch") {
        auto r = overlapBatchGPU(base, others);
        return r[0];
    };
}

TEST_CASE("overlap throughput") {
    Paths64 base = { { {0,0},{100,0},{100,100},{0,100} } };
    std::vector<Paths64> others;
    for(int i=0;i<10000;i++) others.push_back(movePaths(base,i*30,i*30));
    auto t0 = std::chrono::steady_clock::now();
    size_t cnt = 0;
    for(const auto& o:others){ overlap(base,o); cnt++; }
    auto ms_cpu = std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::steady_clock::now()-t0).count();
    std::cout << "CPU throughput: " << (cnt/1000.0)/(ms_cpu/1000.0)
              << " k-checks/s\n";
    if(cuda_available()) {
        t0 = std::chrono::steady_clock::now();
        overlapBatchGPU(base, others);
        auto ms_gpu = std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::steady_clock::now()-t0).count();
        std::cout << "GPU throughput: " << (cnt/1000.0)/(ms_gpu/1000.0)
                  << " k-checks/s\n";
    } else {
        std::cout << "GPU not available\n";
    }
}
