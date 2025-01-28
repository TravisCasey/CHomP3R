#!/usr/bin/env sh

# Builds and tests the project.
# Run from project root directory
cmake --preset chomp-gcc-x86-native-profile
make -C ./build/chomp-gcc-x86-native-profile
cd build/chomp-gcc-x86-native-profile/chomp
./tests [!benchmark] --benchmark-samples 40
gprof tests gmon.out > gprof.output
cd ../../..