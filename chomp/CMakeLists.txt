set(TEST_SOURCES
    ${CHOMP_DIR}/chomp/modules/rings.test.cpp)

add_executable(tests ${TEST_SOURCES})
catch_discover_tests(tests)
