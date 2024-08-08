cmake_minimum_required(VERSION 3.22)

set(PROJECT_DESCRIPTION "Computational Homology Project 3 (Redux)")
set(ONLINE_REPOSITORY "https://github.com/TravisCasey/CHomP3R")
project(
    CHomP3R
    VERSION 0.1.0
    DESCRIPTION ${PROJECT_DESCRIPTION}
    HOMEPAGE_URL ${ONLINE_REPOSITORY}
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++23 -Wall -Wextra -O3")
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

find_package(Catch2 3 REQUIRED)
include(CTest)
include(Catch)

set(CHOMP_DIR ${PROJECT_SOURCE_DIR}/chomp)
add_subdirectory(${CHOMP_DIR})

target_include_directories(tests PRIVATE chomp)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)
