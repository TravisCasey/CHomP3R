/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/core/algebra/cyclic.hpp>
#include <chomp/core/util/morse.hpp>

#include <catch2/catch_test_macros.hpp>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

TEST_CASE("Trichotomy initialization and comparison.", "[util]") {
  Trichotomy tri;

  CHECK(tri == tri);
  CHECK(tri == Trichotomy::ace);
  CHECK(tri != Trichotomy::queen);
  CHECK(tri != Trichotomy::king);

  CHECK(Trichotomy::ace == tri);
  CHECK(Trichotomy::ace == Trichotomy::ace);
  CHECK(Trichotomy::ace != Trichotomy::queen);
  CHECK(Trichotomy::ace != Trichotomy::king);

  CHECK(Trichotomy::queen != tri);
  CHECK(Trichotomy::queen != Trichotomy::ace);
  CHECK(Trichotomy::queen == Trichotomy::queen);
  CHECK(Trichotomy::queen != Trichotomy::king);

  CHECK(Trichotomy::king != tri);
  CHECK(Trichotomy::king != Trichotomy::ace);
  CHECK(Trichotomy::king != Trichotomy::queen);
  CHECK(Trichotomy::king == Trichotomy::king);
}

TEST_CASE("MatchResult methods behave correctly.", "[util]") {
  const MatchResult<int, Z<3>> def_MR;
  CHECK(def_MR.cell() == int());
  CHECK(def_MR.coef() == Z<3>());
  CHECK(def_MR.tri() == Trichotomy::ace);
  CHECK(def_MR.is_ace());
  CHECK_FALSE(def_MR.is_queen());
  CHECK_FALSE(def_MR.is_king());

  const MatchResult<int, Z<3>> queen_MR(3, Z<3>(1), Trichotomy::queen);
  CHECK(queen_MR.cell() == 3);
  CHECK(queen_MR.coef() == Z<3>(1));
  CHECK(queen_MR.tri() == Trichotomy::queen);
  CHECK_FALSE(queen_MR.is_ace());
  CHECK(queen_MR.is_queen());
  CHECK_FALSE(queen_MR.is_king());

  const MatchResult<int, Z<3>> king_MR(-20, Z<3>(0), Trichotomy::king);
  CHECK(king_MR.cell() == -20);
  CHECK(king_MR.coef() == Z<3>(0));
  CHECK(king_MR.tri() == Trichotomy::king);
  CHECK_FALSE(king_MR.is_ace());
  CHECK_FALSE(king_MR.is_queen());
  CHECK(king_MR.is_king());
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
