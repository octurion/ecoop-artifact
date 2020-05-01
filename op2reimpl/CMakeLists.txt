cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(op2reimpl)

find_package(OpenMP REQUIRED)

set(COMMON_FLAGS ${WARN_FLAGS} ${OpenMP_CXX_FLAGS})
set(LINK_FLAGS ${OpenMP_CXX_FLAGS} m)

add_executable(airfoil src/airfoil.cpp)
target_compile_options(airfoil PRIVATE ${COMMON_FLAGS})
target_link_libraries(airfoil ${LINK_FLAGS})

add_executable(airfoil_aos src/airfoil_aos.cpp)
target_compile_options(airfoil_aos PRIVATE ${COMMON_FLAGS})
target_link_libraries(airfoil_aos ${LINK_FLAGS})

add_executable(airfoil_soa src/airfoil_soa.cpp)
target_compile_options(airfoil_soa PRIVATE ${COMMON_FLAGS})
target_link_libraries(airfoil_soa ${LINK_FLAGS})

add_executable(aero src/aero.cpp)
target_compile_options(aero PRIVATE ${COMMON_FLAGS})
target_link_libraries(aero ${LINK_FLAGS})

add_executable(aero_aos src/aero_aos.cpp)
target_compile_options(aero_aos PRIVATE ${COMMON_FLAGS})
target_link_libraries(aero_aos ${LINK_FLAGS})
