#!/usr/bin/env sh

# Builds and tests the project.
# Run from project root directory
cmake --preset chomp-gcc-x86-native-bench
make -C ./build/chomp-gcc-x86-native-bench
build/chomp-gcc-x86-native-bench/chomp/tests [!benchmark]