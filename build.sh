#!/bin/sh
set -eu

# First build all projects but the original OP2 code
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release\
    -DBENCHMARK_ENABLE_TESTING=OFF \
    -DBENCHMARK_ENABLE_GTEST_TESTS=OFF \
    -DBENCHMARK_ENABLE_LTO=ON \
    -GNinja ..
cd ..
ninja -C build/

# Extract the OP2 input files
tar --lzma -xvf op2reimpl/op2_inputs.tar.lzma -C build/op2reimpl/

cd op2orig/op2/c/
./cmake.local -DCMAKE_BUILD_TYPE=Release
cd ../../../

mkdir op2orig/apps/c/build
cp build/op2reimpl/new_grid.dat build/op2reimpl/FE_grid.dat op2orig/apps/c/build
cd op2orig/apps/c/
./cmake.local -DCMAKE_BUILD_TYPE=Release
cd ../../../
