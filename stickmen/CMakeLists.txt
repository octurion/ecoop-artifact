cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(model_anim)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(lib/benchmark)

add_library(anim src/anim/anim.cpp)
target_include_directories(anim PRIVATE include/anim INTERFACE include/)
target_compile_options(anim PRIVATE -Wall -Wextra -Wpedantic)

add_executable(bench src/test/test.cpp src/test/data.cpp)
target_include_directories(bench PRIVATE include/test)
target_compile_options(bench PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(bench anim benchmark)
