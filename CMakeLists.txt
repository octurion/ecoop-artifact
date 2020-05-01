cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(ecoop_artifact)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(WARN_FLAGS -Wall -Wextra -Wpedantic)

add_subdirectory(lib/benchmark)

add_subdirectory(doors/)
add_subdirectory(forex/)
add_subdirectory(op2reimpl/)
add_subdirectory(stickmen/)
add_subdirectory(traffic/)

add_subdirectory(charts/)
