/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/util/concepts.hpp>

#include <catch2/catch_test_macros.hpp>

#include <concepts>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

TEST_CASE("Hashable concept functions as expected", "[util]") {
  CHECK(Hashable<int>);
  CHECK(Hashable<std::vector<bool>>);
  CHECK_FALSE(Hashable<std::vector<int>>);
  CHECK(AssociativeKey<std::vector<int>>);
}

TEST_CASE("Comparable concept functions as expected", "[util]") {
  struct NonComparable {
    int val;
    explicit NonComparable(int n) : val(n) {}
    bool operator==(const NonComparable& rhs) const {
      return val == rhs.val;
    }
  };
  CHECK(Comparable<int>);
  CHECK(Comparable<std::vector<unsigned short>>);
  CHECK_FALSE(Comparable<NonComparable>);
  CHECK_FALSE(AssociativeKey<NonComparable>);
}

TEST_CASE("DefaultMap chooses correctly", "[util]") {
  // int is hashable but std::vector<int> is not
  CHECK(std::same_as<
        DefaultMap<int, std::vector<int>>,
        std::unordered_map<int, std::vector<int>>>);
  CHECK(std::same_as<
        DefaultMap<std::vector<int>, int>, std::map<std::vector<int>, int>>);
}

TEST_CASE("DefaultSet chooses correctly", "[util]") {
  // int is hashable but std::vector<int> is not
  CHECK(std::same_as<DefaultSet<int>, std::unordered_set<int>>);
  CHECK(std::same_as<DefaultSet<std::vector<int>>, std::set<std::vector<int>>>);
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
