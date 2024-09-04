set(TEST_SOURCES
    ${CHOMP_DIR}/chomp/algebra/algebra.test.cpp
    ${CHOMP_DIR}/chomp/algebra/cyclic.test.cpp
    ${CHOMP_DIR}/chomp/algebra/modules.test.cpp
    ${CHOMP_DIR}/chomp/util/iterators.test.cpp)

add_executable(tests ${TEST_SOURCES})
catch_discover_tests(tests)
