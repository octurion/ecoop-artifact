cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(traffic)

set(COMMON_FLAGS ${WARN_FLAGS})

add_executable(traffic src/traffic.cpp)
target_compile_options(traffic PRIVATE ${COMMON_FLAGS})
target_link_libraries(traffic benchmark)
