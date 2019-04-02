cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(op2meme)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

find_package(OpenMP REQUIRED)

set(COMMON_FLAGS -Wall -Wextra -Wpedantic ${OpenMP_CXX_FLAGS})
set(LINK_FLAGS ${OpenMP_CXX_FLAGS})

add_executable(airfoil src/airfoil.cpp)
target_compile_options(airfoil PRIVATE ${COMMON_FLAGS})
target_link_libraries(airfoil ${LINK_FLAGS})

add_executable(airfoil_aos src/airfoil_aos.cpp)
target_compile_options(airfoil_aos PRIVATE ${COMMON_FLAGS})
target_link_libraries(airfoil_aos ${LINK_FLAGS})

add_executable(airfoil_soa src/airfoil_soa.cpp)
target_compile_options(airfoil_soa PRIVATE ${COMMON_FLAGS})
target_link_libraries(airfoil_soa ${LINK_FLAGS})

add_executable(airfoil_v3 src/airfoil_v3.cpp)
target_compile_options(airfoil_v3 PRIVATE ${COMMON_FLAGS})
target_link_libraries(airfoil_v3 ${LINK_FLAGS})

add_executable(aero src/aero.cpp)
target_compile_options(aero PRIVATE ${COMMON_FLAGS})
target_link_libraries(aero ${LINK_FLAGS})

add_executable(aero_aos src/aero_aos.cpp)
target_compile_options(aero_aos PRIVATE ${COMMON_FLAGS})
target_link_libraries(aero_aos ${LINK_FLAGS})

add_executable(aero_v2 src/aero_v2.cpp)
target_compile_options(aero_v2 PRIVATE ${COMMON_FLAGS})
target_link_libraries(aero_v2 ${LINK_FLAGS})

add_executable(jac1 src/jac1.cpp)
target_compile_options(jac1 PRIVATE ${COMMON_FLAGS})
target_link_libraries(jac1 ${LINK_FLAGS})
