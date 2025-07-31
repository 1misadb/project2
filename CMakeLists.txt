cmake_minimum_required(VERSION 3.14)
project(nesting)
set(CMAKE_CXX_STANDARD 17)
# Relax TBB version requirement for broader compatibility
find_package(TBB REQUIRED)
find_package(CUDAToolkit)
add_library(clipper3
    clipper3/src/clipper.engine.cpp
    clipper3/src/clipper.offset.cpp
    clipper3/src/clipper.rectclip.cpp
)
target_include_directories(clipper3 PUBLIC clipper3)
add_library(geometry geometry.cpp)
target_include_directories(geometry PUBLIC . clipper3)
target_link_libraries(geometry PUBLIC clipper3)
if(CUDAToolkit_FOUND)
    enable_language(CUDA)
    add_library(cuda_kernels STATIC overlap_gpu.cu nfp_gpu.cu)
    set_target_properties(cuda_kernels PROPERTIES POSITION_INDEPENDENT_CODE ON)
    target_link_libraries(cuda_kernels PUBLIC CUDA::cudart)
    target_compile_definitions(cuda_kernels PUBLIC USE_CUDA=1)
    target_include_directories(cuda_kernels PUBLIC .)
    target_link_libraries(geometry PUBLIC cuda_kernels)
    add_definitions(-DUSE_CUDA)
endif()

add_executable(nest
    nesting_algorithm.cpp
)
target_include_directories(nest PRIVATE . clipper3 spdlog/include)
target_link_libraries(nest PRIVATE geometry TBB::tbb)

if (MSVC)
    target_compile_options(nest PRIVATE /utf-8)
endif()

add_executable(test_overlap tests/test_overlap.cpp)
target_include_directories(test_overlap PRIVATE . clipper3)
target_link_libraries(test_overlap PRIVATE geometry TBB::tbb)

add_executable(test_benchmark tests/test_benchmark.cpp)
target_include_directories(test_benchmark PRIVATE . clipper3)
target_link_libraries(test_benchmark PRIVATE geometry TBB::tbb)

add_executable(test_minkowski tests/test_minkowski_gpu.cpp)
target_include_directories(test_minkowski PRIVATE . clipper3)
target_link_libraries(test_minkowski PRIVATE geometry TBB::tbb)
