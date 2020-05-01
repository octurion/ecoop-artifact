#!/bin/sh

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
