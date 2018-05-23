cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(model_anim)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(BENCHMARK_ENABLE_TESTING OFF)
set(BENCHMARK_ENABLE_INSTALL OFF)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
add_subdirectory(lib/benchmark)

add_library(anim src/anim_naive.cpp src/anim_smart.cpp src/anim_soa.cpp)
target_include_directories(anim PUBLIC include)
target_include_directories(anim PUBLIC SYSTEM lib/benchmark/include)
target_compile_options(anim PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(anim benchmark)

add_executable(bench test.cpp)
target_compile_options(bench PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(bench anim benchmark)

