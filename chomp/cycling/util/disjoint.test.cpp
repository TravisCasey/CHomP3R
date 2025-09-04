/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/cycling/util/disjoint.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>

#ifndef CHOMP_DOXYGEN

namespace chomp::cycling::detail {

TEST_CASE("DisjointSet class", "[cycling][util]") {
  DisjointSet dset(7);

  // Currently all disjoint singleton sets
  for (std::size_t element = 0; element < 7; ++element) {
    CHECK(dset.rep(element) == element);
  }

  // Merge two sets
  dset.merge(0, 3);
  for (std::size_t element = 0; element < 7; ++element) {
    if (element != 0 && element != 3) {
      CHECK(dset.rep(element) == element);
    }
  }
  CHECK(dset.rep(0) == dset.rep(3));

  // Check for weighted union optimization
  dset.merge(5, 6);
  dset.merge(0, 5);
  // set containing 0, 3, 5, and 6 now has rank 3
  const std::size_t larger_rep = dset.rep(0);
  dset.merge(2, 4);
  // set containing 2 and 4 now has rank 2
  dset.merge(2, 0);
  // root after merging should be the root of the set of higher rank.
  CHECK(dset.rep(0) == larger_rep);
  CHECK(dset.rep(2) == larger_rep);
  CHECK(dset.rep(3) == larger_rep);
  CHECK(dset.rep(4) == larger_rep);
  CHECK(dset.rep(5) == larger_rep);
  CHECK(dset.rep(6) == larger_rep);

  // Append new disjoint singleton set and merge
  CHECK(dset.insert() == 7);
  CHECK(dset.rep(7) == 7);
  dset.merge(7, 0);
  CHECK(dset.rep(7) == larger_rep);
}

}  // namespace chomp::cycling::detail

#endif  // CHOMP_DOXYGEN
