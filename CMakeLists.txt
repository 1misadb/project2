cmake_minimum_required(VERSION 3.14)
project(nesting)
set(CMAKE_CXX_STANDARD 17)
find_package(TBB REQUIRED)
add_library(clipper3
    clipper3/src/clipper.engine.cpp
    clipper3/src/clipper.offset.cpp
    clipper3/src/clipper.rectclip.cpp
)
target_include_directories(clipper3 PUBLIC clipper3)
add_executable(nest
    nesting_algorithm.cpp
)
target_include_directories(nest PRIVATE . clipper3)
target_link_libraries(nest PRIVATE clipper3 TBB::tbb)

add_executable(test_overlap tests/test_overlap.cpp)
target_include_directories(test_overlap PRIVATE . clipper3)
target_link_libraries(test_overlap PRIVATE clipper3 TBB::tbb)
