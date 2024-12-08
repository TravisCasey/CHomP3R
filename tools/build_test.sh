#!/usr/bin/env sh

# Builds and tests the project.
# Run from project root directory
cmake --preset chomp-gcc-debug
make -C ./build/chomp-gcc-debug
ctest ./build/chomp-gcc-debug --preset default
