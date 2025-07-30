cmake_minimum_required(VERSION 3.14)
project(nesting)
set(CMAKE_CXX_STANDARD 17)
option(USE_CUDA "Build with CUDA" ON)
if(USE_CUDA)
  enable_language(CUDA)
endif()
find_package(TBB REQUIRED)
add_library(clipper3
    clipper3/src/clipper.engine.cpp
    clipper3/src/clipper.offset.cpp
    clipper3/src/clipper.rectclip.cpp
)
target_include_directories(clipper3 PUBLIC clipper3)
add_library(geometry geometry.cpp)
target_include_directories(geometry PUBLIC . clipper3)
target_link_libraries(geometry PUBLIC clipper3)

if(USE_CUDA)
  add_library(cuda_overlap cuda/overlap_cuda.cu)
  target_include_directories(cuda_overlap PUBLIC .)
  set_target_properties(cuda_overlap PROPERTIES CUDA_STANDARD 17 CUDA_SEPARABLE_COMPILATION ON)
  find_package(CUDAToolkit REQUIRED)
  target_link_libraries(cuda_overlap PUBLIC CUDA::cudart)
  target_compile_features(cuda_overlap PUBLIC cxx_std_17)
  target_link_libraries(geometry PUBLIC cuda_overlap)
  target_compile_definitions(geometry PUBLIC USE_CUDA)
endif()

add_executable(nest
    nesting_algorithm.cpp
)
target_include_directories(nest PRIVATE . clipper3 spdlog/include)
target_link_libraries(nest PRIVATE geometry TBB::tbb)

add_executable(test_overlap tests/test_overlap.cpp)
target_include_directories(test_overlap PRIVATE . clipper3)
target_link_libraries(test_overlap PRIVATE geometry TBB::tbb)

add_executable(test_benchmark tests/test_benchmark.cpp)
target_include_directories(test_benchmark PRIVATE . clipper3)
target_link_libraries(test_benchmark PRIVATE geometry TBB::tbb)
