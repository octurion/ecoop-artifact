cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(doors)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(lib/benchmark)

add_executable(bench src/main.cpp)
target_compile_options(bench PRIVATE -Wall -Wextra -Wpedantic)
target_link_libraries(bench m benchmark)
