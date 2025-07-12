/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/core/algebra/algebra.hpp>

#include <catch2/catch_test_macros.hpp>

#include <stdexcept>
#include <string>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

TEST_CASE("Fundamental types and algebraic concepts", "[algebra]") {
  CHECK_FALSE(Ring<std::string>);
  CHECK(Ring<int>);
  CHECK(Ring<float>);
  CHECK_FALSE(Ring<void>);
  CHECK(Ring<char>);
  CHECK_FALSE(BinaryRing<char>);
}

TEST_CASE("Identity functions on fundamental types.", "[algebra]") {
  CHECK(zero<int>() == 0);
  CHECK(zero<float>() == 0.0);
  CHECK(one<int>() == 1);
  CHECK(one<float>() == 1.0);
}

TEST_CASE("Inverse functions on fundamental types", "[algebra]") {
  CHECK_THROWS_AS(invert(0.0), std::domain_error);
  CHECK_THROWS_AS(invert(static_cast<int>(0)), std::domain_error);
  CHECK(invert(static_cast<int>(1)) == static_cast<int>(1));
  CHECK(invert(static_cast<int>(-1)) == static_cast<int>(-1));
  CHECK_THROWS_AS(invert(static_cast<unsigned int>(-1)), std::domain_error);
  CHECK(invert(static_cast<unsigned int>(1)) == static_cast<unsigned int>(1));
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
