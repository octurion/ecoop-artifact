cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(stickmen)

add_library(anim src/anim/anim.cpp)
target_include_directories(anim PRIVATE include/anim INTERFACE include/)
target_compile_options(anim PRIVATE ${WARN_FLAGS})

add_executable(stickmen src/test/test.cpp)
target_include_directories(stickmen PRIVATE include/test)
target_compile_options(stickmen PRIVATE ${WARN_FLAGS})
target_link_libraries(stickmen anim benchmark)
