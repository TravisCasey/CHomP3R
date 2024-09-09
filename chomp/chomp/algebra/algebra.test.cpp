/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/algebra/algebra.hpp>

#include <catch2/catch_test_macros.hpp>

#include <string>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

TEST_CASE("Fundamental types and algebraic concepts", "[algebra]") {
  CHECK_FALSE(Group<std::string>);
  CHECK(Ring<int>);
  CHECK(Ring<float>);
  CHECK_FALSE(Group<void>);
  CHECK(Ring<char>);
  CHECK_FALSE(BinaryRing<char>);
}

TEST_CASE(
    "Identity functions on fundamental types output correctly", "[algebra]"
) {
  CHECK(zero<int>() == 0);
  CHECK(zero<float>() == 0.0);
  CHECK(one<int>() == 1);
  CHECK(one<float>() == 1.0);
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
