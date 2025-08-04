#ifdef USE_CUDA
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include <cuda_runtime.h>
#include <vector>

TEST_CASE("offset prefix copy"){
    const int n = 2;
    std::vector<size_t> off = {0, 16, 32};
    size_t* d_off = nullptr;
    cudaError_t rc = cudaMalloc(&d_off, (n+1)*sizeof(size_t));
    REQUIRE(rc == cudaSuccess);
    REQUIRE(d_off != nullptr);
    rc = cudaMemcpy(d_off, off.data(), (n+1)*sizeof(size_t), cudaMemcpyHostToDevice);
    REQUIRE(rc == cudaSuccess);
    std::vector<size_t> back(n+1, 0);
    rc = cudaMemcpy(back.data(), d_off, (n+1)*sizeof(size_t), cudaMemcpyDeviceToHost);
    REQUIRE(rc == cudaSuccess);
    REQUIRE(back == off);
    if(d_off){
        rc = cudaFree(d_off);
        REQUIRE(rc == cudaSuccess);
        d_off = nullptr;
    }
}
#else
int main(){ return 0; }
#endif
