if(CATCH2_TEST OR CATCH2_BENCH)
    set(TEST_SOURCES
        ${CHOMP_DIR}/core/algebra/algebra.test.cpp
        ${CHOMP_DIR}/core/algebra/cyclic.test.cpp
        ${CHOMP_DIR}/core/algebra/modules.test.cpp
        ${CHOMP_DIR}/core/complexes/cubical.test.cpp
        ${CHOMP_DIR}/core/complexes/grading.test.cpp
        ${CHOMP_DIR}/core/complexes/morse.test.cpp
        ${CHOMP_DIR}/core/homology/cubical.test.cpp
        ${CHOMP_DIR}/core/homology/hasse.test.cpp
        ${CHOMP_DIR}/core/homology/morse.test.cpp
        ${CHOMP_DIR}/core/util/cache.test.cpp
        ${CHOMP_DIR}/core/util/concepts.test.cpp
        ${CHOMP_DIR}/core/util/iterators.test.cpp
        ${CHOMP_DIR}/core/util/morse.test.cpp
        ${CHOMP_DIR}/cycling/trajectory/spline.test.cpp)

    add_executable(tests ${TEST_SOURCES})
    catch_discover_tests(tests)
endif()
