/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <catch2/catch_test_macros.hpp>

#include <chomp/modules/rings.hpp>

namespace chomp::modules {

TEST_CASE("Fundamental numerical types fulfill Ring requirement", "[rings]") {
    CHECK(Ring<int>);
    CHECK(Ring<float>);
    CHECK_FALSE(Ring<void>);
    CHECK(Ring<char>);
}

TEST_CASE("Identity functions on fundamental types output correctly, [rings]") {
    CHECK(zero<int>() == 0);
    CHECK(zero<float>() == 0.0);
    CHECK(one<int>() == 1);
    CHECK(one<float>() == 1.0);
}

TEST_CASE("Z fulfills Ring requirement", "[rings]") {
    CHECK(Ring<Z<2>>);
    CHECK(Ring<Z<6>>);
}

TEST_CASE("Comparison operators of Z work correctly", "[rings]") {
    CHECK(Z<2>(3) == Z<2>(-1));
    CHECK(Z<7>(3) != Z<7>(-1));
}

TEST_CASE("Arithmetic operators of Z work correctly", "[rings]") {
    CHECK(-Z<14>(-9) == Z<14>(9));
    CHECK(Z<10>(-4) + Z<10>(28) == Z<10>(-4 + 28));
    CHECK(Z<12>(11) - Z<12>(1) == Z<12>(11 - 1));
    CHECK(Z<2>(1) * Z<2>(2) == Z<2>(1 * 2));

    Z<5> a(0);
    a += Z<5>(-1);
    REQUIRE(a == Z<5>(0 - 1));
    a -= Z<5>(-3);
    REQUIRE(a == Z<5>(0 - 1 + 3));
    a *= Z<5>(2);
    REQUIRE(a == Z<5>(2 * (0 - 1 + 3)));
}

TEST_CASE("Integral coversion of Z works correctly", "[rings]") {
    CHECK(int(Z<12>(-2)) == 10);
    CHECK(zero<Z<24>>() == Z<24>(0));
    CHECK(one<Z<3>>() == Z<3>(1));

    CHECK(int(Z<2>(-455)) == 1);
    CHECK(int(Z<2>(44)) == 0);
}

} // namespace chomp::modules
