#!/usr/bin/env sh

# Builds and tests the project.
# Run from project root directory
rm -rf build && mkdir build && cd build

cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..
make
ctest --output-on-failure
