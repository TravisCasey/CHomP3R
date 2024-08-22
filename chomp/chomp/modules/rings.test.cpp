/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <catch2/catch_test_macros.hpp>

#include <chomp/modules/concepts.hpp>
#include <chomp/modules/rings.hpp>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace chomp::modules {

TEST_CASE("BinaryRing concept on Z", "[rings]") {
    CHECK(BinaryRing<Z<unsigned int, 2>>);
    CHECK_FALSE(BinaryRing<Z<unsigned int, 3>>);
    CHECK_FALSE(BinaryRing<int>);
}

TEST_CASE("Identity functions on fundamental types output correctly,"
          "[rings]") {
    CHECK(zero<int>() == 0);
    CHECK(zero<float>() == 0.0);
    CHECK(one<int>() == 1);
    CHECK(one<float>() == 1.0);
}

TEST_CASE("Z constraints work as intended", "[rings]") {
    CHECK_NOTHROW(Z<unsigned int, 2>(3));
    CHECK_NOTHROW(Z<unsigned short, 4>(0));
    CHECK_NOTHROW(Z<unsigned short, (1 << 8) - 1>(0));
}

TEST_CASE("Z fulfills Ring requirement", "[rings]") {
    CHECK(Ring<Z<unsigned int, 2>>);
    CHECK(Ring<Z<unsigned char, 6>>);
}

TEST_CASE("Ring functions work correctly", "[rings]") {
    CHECK(zero<Z<unsigned long, 13>>().get_value() ==
          static_cast<unsigned long>(0));
    CHECK(one<Z<unsigned long long, 2>>().get_value() ==
          static_cast<unsigned long long>(1));
}

TEST_CASE("Comparison operators of Z work as intended", "[rings]") {
    CHECK(Z<unsigned int, 2>(3) == Z<unsigned int, 2>(-1));
    CHECK(Z<unsigned short, 7>(3) != Z<unsigned short, 7>(-1));
}

TEST_CASE("Arithmetic operators of Z work as intended", "[rings]") {
    CHECK(-Z<unsigned long, 14>(9) == Z<unsigned long, 14>(5));
    CHECK(Z<unsigned char, 10>(4) + Z<unsigned char, 10>(28) ==
          Z<unsigned char, 10>(4 + 28));
    CHECK(Z<unsigned int, 12>(11) - Z<unsigned int, 12>(1) ==
          Z<unsigned int, 12>(11 - 1));
    CHECK(Z<unsigned long long, 2>(1) * Z<unsigned long long, 2>(2) ==
          Z<unsigned long long, 2>(1 * 2));

    Z<unsigned short, 5> a(0);
    a += Z<unsigned short, 5>(1);
    REQUIRE(a.get_value() == 1);
    a -= Z<unsigned short, 5>(3);
    REQUIRE(a.get_value() == 3);
    a *= Z<unsigned short, 5>(2);
    REQUIRE(a.get_value() == 1);
}

} // namespace chomp::modules

#endif // DOXYGEN_SHOULD_SKIP_THIS
