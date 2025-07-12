cmake_minimum_required(VERSION 3.31)

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

if(CATCH2_TEST)
    include(CTest)
endif()

if(CATCH2_TEST OR CATCH2_BENCH)
    find_package(Catch2 3 REQUIRED)
    include(Catch)
endif()

set(CHOMP_DIR ${PROJECT_SOURCE_DIR}/chomp)
add_subdirectory(${CHOMP_DIR})

if(CATCH2_TEST OR CATCH2_BENCH)
    target_include_directories(tests PRIVATE ${PROJECT_SOURCE_DIR})
    target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)
endif()
