# CHomP3R - Computational Homology Project 3 (Redux)
This library is designed for efficient homology computation of high-dimensional cell complexes. It is in an early state with a current focus on cubical complexes defined as a collection of top-dimensional (i.e., dimension of the embedding space) cubes.

## Usage
Due to its highly-templated nature, CHomP3R is a header-only library with no currrent dependencies except the standard library (with C++20 features). Including the `*.hpp` files from the `chomp` directory in your project (via CMake or otherwise) is all that is required for use. If you plan to modify the source code see the [LICENSE](CHomP3R/LICENSE).

### Scripts
- `build_test.sh` Compiles the tests (written in `*.test.cpp` files) using CMake and runs them using CTest and Catch2. Intended for and tested on CMake/CTest version 3.31.0 and Catch2 version 3.6.0. Currently, this is hard-coded to use my CMake preset but this will be modified in the future. This requires a C++20-compliant compiler; I am using gcc version 14.2.0.
- `benchmark.sh` Compiles and runs benchmarks using the same programs as `build_test.sh`. Currently, this is hard-coded to use my CMake preset but this will be modified in the future.
- `lint.sh` Lint a selected file using clang-tidy (version 18.1.8) with the options in [.clang-tidy](CHomP3R/.clang-tidy). Note that `build_test.sh` should be run first to generate `compile_commands.json` for the debug build.
- `profile.sh` Compile an optimized build (with debug symbols) of benchmarks and profile it using gprof (version 2.3.8).
