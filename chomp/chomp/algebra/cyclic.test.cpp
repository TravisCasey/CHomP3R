/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/algebra/algebra.hpp>
#include <chomp/algebra/cyclic.hpp>

#include <catch2/catch_test_macros.hpp>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

TEST_CASE("BinaryRing concept on Z", "[algebra]") {
  CHECK(BinaryRing<Z<2>>);
  CHECK_FALSE(BinaryRing<Z<3>>);
}

TEST_CASE("Ring identify functions work correctly", "[algebra]") {
  CHECK(zero<Z<13>>().rep() == 0);
  CHECK(zero<Z<2>>().rep() == 0);
  CHECK(one<Z<22>>().rep() == 1);
  CHECK(one<Z<2>>().rep() == 1);
}

TEST_CASE("Comparison operators of Z work as intended", "[algebra]") {
  CHECK(Z<7>(-800) == Z<7>(-107));
  CHECK(Z<7>(3) != Z<7>(-1));
  CHECK(Z<2>(3) == Z<2>(-1));
  CHECK(Z<2>(-800) != Z<2>(-107));
}

TEST_CASE("Arithmetic operators of Z work as intended", "[algebra]") {
  CHECK(-Z<14>(9) == Z<14>(5));
  CHECK(-Z<2>(9) == Z<2>(9));

  CHECK(Z<10>(4) + Z<10>(28) == Z<10>(4 + 28));
  CHECK(Z<12>(11) - Z<12>(1) == Z<12>(11 - 1));
  CHECK(Z<2>(1) * Z<2>(2) == Z<2>(1 * 2));

  Z<5> a(0);
  a += Z<5>(1);
  REQUIRE(a.rep() == 1);
  a -= Z<5>(3);
  REQUIRE(a.rep() == 3);
  a *= Z<5>(2);
  REQUIRE(a.rep() == 1);
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
