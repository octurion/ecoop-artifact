cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(model_anim)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(COMMON_FLAGS -Wall -Wextra -Wpedantic)

add_subdirectory(lib/benchmark)

add_executable(bench src/traffic.cpp)
target_compile_options(bench PRIVATE ${COMMON_FLAGS})
target_link_libraries(bench benchmark)
