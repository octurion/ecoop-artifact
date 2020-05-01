cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(doors)

add_executable(doors src/main.cpp)
target_compile_options(doors PRIVATE ${WARN_FLAGS})
target_link_libraries(doors benchmark)
