#ifdef USE_CUDA
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <cuda_runtime.h>
#include <vector>

TEST_CASE("offset prefix copy"){
    const int n = 2;
    std::vector<size_t> off = {0, 16, 32};
    size_t* d_off = nullptr;
    REQUIRE(cudaMalloc(&d_off, (n+1)*sizeof(size_t)) == cudaSuccess);
    REQUIRE(cudaMemcpy(d_off, off.data(), (n+1)*sizeof(size_t), cudaMemcpyHostToDevice) == cudaSuccess);
    std::vector<size_t> back(n+1, 0);
    REQUIRE(cudaMemcpy(back.data(), d_off, (n+1)*sizeof(size_t), cudaMemcpyDeviceToHost) == cudaSuccess);
    REQUIRE(back == off);
    cudaFree(d_off);
}
#else
int main(){ return 0; }
#endif
