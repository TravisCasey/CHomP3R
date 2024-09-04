#!/usr/bin/env sh

# Runs clang-tidy on all files with .cpp and .hpp extensions in passed directory.
# Run from project root directory.
find $1 -iname '*.hpp' -o -iname '*.cpp' | xargs clang-tidy --quiet --config-file \
 .clang-tidy -p build/compile_commands.json --extra-arg='-DCHOMP_CLANG_TIDY'
