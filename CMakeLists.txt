cmake_minimum_required(VERSION 3.18)
project(nesting LANGUAGES CXX)

# Build with C++20 and emit debug information even in release builds
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE)
endif()

# Enable aggressive warnings for both host and device code.  AddressSanitizer
# is enabled for host objects.  CUDA receives equivalent flags via
# --compiler-options.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -fsanitize=address")
if(CUDAToolkit_FOUND)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -G --generate-line-info --device-debug --compiler-options=-Wall,-Wextra,-Werror,-fsanitize=address")
endif()
link_libraries(-fsanitize=address)

# --- DEPENDENCIES ---
find_package(TBB CONFIG REQUIRED)

find_package(CUDAToolkit)
find_package(OpenMP REQUIRED)  # <-- добавляем OpenMP

# --- CLIPPER3 LIBRARY ---
add_library(clipper3
    clipper3/src/clipper.engine.cpp
    clipper3/src/clipper.offset.cpp
    clipper3/src/clipper.rectclip.cpp
)
target_include_directories(clipper3 PUBLIC clipper3)

# --- GEOMETRY LIBRARY ---
add_library(geometry geometry.cpp)
target_include_directories(geometry PUBLIC . clipper3 spdlog/include)
target_link_libraries(geometry PUBLIC clipper3)

# --- CUDA (OPTIONAL) ---
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

# --- MAIN EXECUTABLE ---
add_executable(nest nesting_algorithm.cpp)
target_include_directories(nest PRIVATE . clipper3 spdlog/include)
target_link_libraries(nest PRIVATE geometry TBB::tbb)
target_link_libraries(nest PRIVATE TBB::tbb)

# --- OpenMP SUPPORT ---
if(OpenMP_CXX_FOUND)
    target_link_libraries(nest PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP enabled")
endif()

# --- OPTIMIZATIONS & MSVC ---
if(MSVC)
    target_compile_options(nest PRIVATE /utf-8 /openmp)
else()
    target_compile_options(nest PRIVATE -O3 -fopenmp)
endif()

# --- TESTS ---
add_executable(test_overlap tests/test_overlap.cpp)
target_include_directories(test_overlap PRIVATE . clipper3)
target_link_libraries(test_overlap PRIVATE geometry TBB::tbb)

add_executable(test_benchmark tests/test_benchmark.cpp)
target_include_directories(test_benchmark PRIVATE . clipper3)
target_link_libraries(test_benchmark PRIVATE geometry TBB::tbb)

add_executable(test_minkowski tests/test_minkowski_gpu.cpp)
target_include_directories(test_minkowski PRIVATE . clipper3)
target_link_libraries(test_minkowski PRIVATE geometry TBB::tbb)

add_executable(test_offsets tests/test_offsets_gpu.cpp)
target_include_directories(test_offsets PRIVATE .)
if(CUDAToolkit_FOUND)
    target_link_libraries(test_offsets PRIVATE CUDA::cudart)
endif()

# Regression driver exercising GPU overlap code on JSON parts
add_executable(test_regression tests/test_regression_gpu.cpp)
target_include_directories(test_regression PRIVATE . clipper3)
target_link_libraries(test_regression PRIVATE geometry TBB::tbb)
