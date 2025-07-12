/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/core/util/iterators.hpp>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <map>
#include <tuple>
#include <type_traits>
#include <unordered_map>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

struct ComparableCell {
  int val;
  constexpr ComparableCell(int n) noexcept : val(n) {};
  constexpr bool operator<(const ComparableCell& rhs) const noexcept {
    return val < rhs.val;
  }
  constexpr bool operator==(const ComparableCell& rhs) const noexcept {
    return val == rhs.val;
  }
};

using MapTypes = std::tuple<
    std::tuple<
        std::map<ComparableCell, short>,
        std::integral_constant<ComparableCell, ComparableCell(33)>,
        std::integral_constant<short, -20>,
        std::integral_constant<ComparableCell, ComparableCell(0)>,
        std::integral_constant<short, 1>>,

    std::tuple<
        std::unordered_map<int, int>, std::integral_constant<int, 29>,
        std::integral_constant<int, 29>, std::integral_constant<int, -29>,
        std::integral_constant<int, 0>>>;

TEMPLATE_LIST_TEST_CASE(
    "KeyIterator functions correctly as a forward iterator", "[util]", MapTypes
) {
  using Map_t = std::tuple_element_t<0, TestType>;
  using key_type = typename Map_t::key_type;
  using mapped_type = typename Map_t::mapped_type;
  const key_type key_0 = std::tuple_element_t<1, TestType>();
  mapped_type val_0 = std::tuple_element_t<2, TestType>();  // NOLINT
  key_type key_1 = std::tuple_element_t<3, TestType>();  // NOLINT
  const mapped_type val_1 = std::tuple_element_t<4, TestType>();

  Map_t m;
  m[key_0] = val_0;
  m[key_1] = val_1;

  auto it = KeyIterator(m.cbegin());
  CHECK(std::forward_iterator<decltype(it)>);

  CHECK(m.count(*it) == 1);
  CHECK(it == KeyIterator(m.cbegin()));

  auto temp = it++;
  CHECK(temp != it);
  CHECK(temp == KeyIterator(m.cbegin()));
  CHECK(m.count(*it) == 1);
  CHECK(++it == KeyIterator(m.cend()));
}

TEST_CASE("Test BFS metaprogramming routine.", "[util]") {
  SECTION("TupleConcat struct.") {
    using Tuple_1 = std::tuple<
        std::integral_constant<std::size_t, 3>,
        std::integral_constant<std::size_t, 0>>;
    using Tuple_2 = std::tuple<>;
    using Tuple_3 = std::tuple<std::integral_constant<std::size_t, 11>>;
    using ResultTuple = std::tuple<
        std::integral_constant<std::size_t, 3>,
        std::integral_constant<std::size_t, 0>,
        std::integral_constant<std::size_t, 11>>;

    CHECK(std::same_as<
          detail::TupleConcat<Tuple_1, Tuple_2, Tuple_3>, ResultTuple>);
  }

  SECTION("BitIf struct.") {
    using original = std::integral_constant<std::size_t, 0b1001>;
    const std::size_t set_bit_flag = 0b0001;
    using set_generated =
        std::tuple<std::integral_constant<std::size_t, 0b1000>>;
    using set_continued = std::tuple<original>;
    const std::size_t unset_bit_flag = 0b0100;
    using unset_generated = std::tuple<>;
    using unset_continued = std::tuple<>;

    CHECK(std::same_as<
          typename detail::BitIf<
              original, original() ^ set_bit_flag,
              bool(original() & set_bit_flag)>::generated,
          set_generated>);
    CHECK(std::same_as<
          typename detail::BitIf<
              original, original() ^ set_bit_flag,
              bool(original() & set_bit_flag)>::continued,
          set_continued>);

    CHECK(std::same_as<
          typename detail::BitIf<
              original, original() ^ unset_bit_flag,
              bool(original() & unset_bit_flag)>::generated,
          unset_generated>);
    CHECK(std::same_as<
          typename detail::BitIf<
              original, original() ^ unset_bit_flag,
              bool(original() & unset_bit_flag)>::continued,
          unset_continued>);
  }

  SECTION("BFSLoop struct.") {
    using InputTuple = std::tuple<
        std::integral_constant<std::size_t, 0b110>,
        std::integral_constant<std::size_t, 0b101>,
        std::integral_constant<std::size_t, 0b011>>;
    using ResultTuple = std::tuple<
        std::integral_constant<std::size_t, 0b100>,
        std::integral_constant<std::size_t, 0b010>,
        std::integral_constant<std::size_t, 0b001>>;

    CHECK(std::same_as<
          typename detail::BFSLoop<3, 0, InputTuple>::type, ResultTuple>);
  }

  SECTION("BFS_CUBE_ARRAY generates correctly.", "[util]") {
    std::array<std::size_t, 16> correct = {0b1111, 0b1110, 0b1101, 0b1011,
                                           0b0111, 0b1100, 0b1010, 0b0110,
                                           0b1001, 0b0101, 0b0011, 0b1000,
                                           0b0100, 0b0010, 0b0001, 0b0000};

    CHECK(BFS_CUBE_ARRAY<4> == correct);
  }
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
