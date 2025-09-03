/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/core/algebra/algebra.hpp>
#include <chomp/core/algebra/cyclic.hpp>

#include <catch2/catch_test_macros.hpp>

#include <stdexcept>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

TEST_CASE("Ring identify functions work correctly", "[algebra]") {
  CHECK(zero<Z<13>>().rep() == 0);
  CHECK(zero<Z<2>>().rep() == 0);
  CHECK(one<Z<23>>().rep() == 1);
  CHECK(one<Z<2>>().rep() == 1);
}

TEST_CASE("Comparison operators of Z work as intended", "[algebra]") {
  CHECK(Z<7>(-800) == Z<7>(-107));
  CHECK(Z<7>(3) != Z<7>(-1));
  CHECK(Z<2>(3) == Z<2>(-1));
  CHECK(Z<2>(-800) != Z<2>(-107));
}

TEST_CASE("Arithmetic operators of Z work as intended", "[algebra]") {
  CHECK(-Z<17>(9) == Z<17>(8));
  CHECK(-Z<2>(9) == Z<2>(9));

  CHECK(Z<11>(4) + Z<11>(28) == Z<11>(4 + 28));
  CHECK(Z<13>(11) - Z<13>(1) == Z<13>(11 - 1));
  CHECK(Z<2>(1) * Z<2>(2) == Z<2>(1 * 2));

  Z<5> a(0);
  a += Z<5>(1);
  REQUIRE(a.rep() == 1);
  a -= Z<5>(3);
  REQUIRE(a.rep() == 3);
  a *= Z<5>(2);
  REQUIRE(a.rep() == 1);
}

TEST_CASE("Inverses in Z", "[algebra]") {
  CHECK_FALSE(invertible(Z<2>(2)));
  CHECK_FALSE(invertible(Z<3>(3)));
  CHECK(invertible(Z<2>(3)));
  CHECK(invertible(Z<3>(2)));

  CHECK_THROWS_AS(invert(Z<2>(0)), std::domain_error);
  CHECK_THROWS_AS(invert(Z<19>(-19)), std::domain_error);

  CHECK(invert(Z<2>(3)) == Z<2>(3));
  CHECK(invert(Z<5>(1)) == Z<5>(1));
  CHECK(invert(Z<5>(2)) == Z<5>(3));
  CHECK(invert(Z<5>(3)) == Z<5>(2));
  CHECK(invert(Z<5>(4)) == Z<5>(4));
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
