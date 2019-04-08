cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(model_anim)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(COMMON_FLAGS -Wall -Wextra -Wpedantic)

add_subdirectory(lib/benchmark)

add_executable(traffic_aos src/traffic_aos.cpp)
target_compile_options(traffic_aos PRIVATE ${COMMON_FLAGS})
target_link_libraries(traffic_aos benchmark)

add_executable(traffic_soa src/traffic_soa.cpp)
target_compile_options(traffic_soa PRIVATE ${COMMON_FLAGS})
target_link_libraries(traffic_soa benchmark)
