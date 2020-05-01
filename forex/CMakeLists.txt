cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(forex)

add_executable(forex src/main.cpp)
target_compile_options(forex PRIVATE ${WARN_FLAGS})
target_link_libraries(forex benchmark)
