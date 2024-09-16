/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/util/cache.hpp>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <map>
#include <tuple>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

using CacheTypes = std::tuple<LRUCache<int, int>, LRUCache<int, int, std::map>>;


TEMPLATE_LIST_TEST_CASE(
    "LRUCache access (and insertion) operations", "[util]", CacheTypes
) {
  TestType c(
      [](const int& x) {
        return 2 * x + 1;
      },
      4
  );

  REQUIRE(c.size() == 0);
  REQUIRE(c.max_size() == 4);

  REQUIRE_FALSE(c.contains(3));
  REQUIRE(c[3] == 7);
  REQUIRE(c.contains(3));
  REQUIRE(c.size() == 1);
  REQUIRE(c[3] == 7);
  REQUIRE(c.size() == 1);

  REQUIRE(c[2] == 5);
  REQUIRE(c[-1] == -1);
  REQUIRE(c[0] == 1);
  REQUIRE(c.size() == 4);

  TestType temp(c);

  REQUIRE(c.contains(3));
  REQUIRE(c.contains(2));
  REQUIRE(c[10] == 21);
  REQUIRE(c.size() == 4);
  REQUIRE_FALSE(c.contains(3));
  REQUIRE(c.contains(2));

  c = std::move(temp);

  c[3];  // 3 is not least recently used anymore
  REQUIRE(c.contains(3));
  REQUIRE(c.contains(2));
  REQUIRE(c[10] == 21);
  REQUIRE(c.size() == 4);
  REQUIRE(c.contains(3));
  REQUIRE_FALSE(c.contains(2));
}

TEST_CASE("CachedFunctionWrapper functions correctly", "[util]") {
  CachedFunctionWrapper<int, int> wrapped_func(
      [](int input) -> int {
        return input << 2;
      },
      4
  );

  REQUIRE(wrapped_func(0) == 0);
  REQUIRE(wrapped_func(0) == 0);
  REQUIRE(wrapped_func(1) == 4);
  REQUIRE(wrapped_func(2) == 8);
  REQUIRE(wrapped_func(3) == 12);
  REQUIRE(wrapped_func(4) == 16);
  REQUIRE(wrapped_func(0) == 0);
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
