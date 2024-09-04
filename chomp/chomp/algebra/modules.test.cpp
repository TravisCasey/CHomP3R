/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/algebra/algebra.hpp>
#include <chomp/algebra/cyclic.hpp>
#include <chomp/algebra/modules.hpp>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <concepts>
#include <cstddef>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

struct HashableCell {
  std::size_t val;
  explicit constexpr HashableCell(std::size_t n) noexcept : val(n) {};
  constexpr bool operator==(const HashableCell& rhs) const noexcept {
    return val == rhs.val;
  }
};

}  // namespace chomp::core

namespace std {
template <>
struct hash<chomp::core::HashableCell> {
  size_t operator()(const chomp::core::HashableCell& cell) const {
    return cell.val;
  }
};
}  // namespace std

namespace chomp::core {

struct ComparableCell {
  int val;
  explicit constexpr ComparableCell(int n) noexcept : val(n) {};
  constexpr bool operator<(const ComparableCell& rhs) const noexcept {
    return val < rhs.val;
  }
  constexpr bool operator==(const ComparableCell& rhs) const noexcept {
    return val == rhs.val;
  }
};

using ModuleTypes = std::tuple<

    std::tuple<
        UnorderedSetModule<int, Z<2>>, std::integral_constant<int, 3>,
        std::integral_constant<int, -25>>,

    std::tuple<
        SetModule<ComparableCell, Z<2>>,
        std::integral_constant<ComparableCell, ComparableCell(20)>,
        std::integral_constant<ComparableCell, ComparableCell(-3)>>,

    std::tuple<
        UnorderedMapModule<HashableCell, Z<14>>,
        std::integral_constant<HashableCell, HashableCell(2445)>,
        std::integral_constant<HashableCell, HashableCell(0)>>,

    std::tuple<
        MapModule<bool, int>, std::integral_constant<bool, true>,
        std::integral_constant<bool, false>>

    >;

TEMPLATE_LIST_TEST_CASE(
    "Module classes model Module concept", "[algebra]", ModuleTypes
) {
  using M = std::tuple_element_t<0, TestType>;
  CHECK(Module<M>);
}

TEMPLATE_LIST_TEST_CASE(
    "Module classes access, insertion, and iteration", "[algebra]", ModuleTypes
) {
  using M = std::tuple_element_t<0, TestType>;
  using T = typename M::BasisType;
  using R = typename M::RingType;
  const T cell_0 = std::tuple_element_t<1, TestType>();
  T cell_1 = std::tuple_element_t<2, TestType>();  // NOLINT
  M elem;

  REQUIRE(elem[cell_0] == zero<R>());
  REQUIRE(elem[cell_1] == zero<R>());
  REQUIRE(elem.begin() == elem.end());

  elem.insert(cell_0, zero<R>());
  elem.insert(cell_1, one<R>());

  REQUIRE(elem[cell_0] == zero<R>());
  REQUIRE(elem[cell_1] == one<R>());
  REQUIRE(elem.begin() != elem.end());

  typename M::BasisIterType it = elem.begin();
  REQUIRE((*it == cell_0 || *it == cell_1));
  it++;
  REQUIRE((it == elem.end() || *it == cell_0 || *it == cell_1));

  elem.insert(cell_0, one<R>());
  elem.insert(cell_1, -one<R>());

  REQUIRE(elem[cell_0] == one<R>());
  REQUIRE(elem[cell_1] == zero<R>());
  REQUIRE(elem.begin() != elem.end());

  it = elem.begin();
  REQUIRE((*it == cell_0 || *it == cell_1));
  it++;
  REQUIRE((it == elem.end() || *it == cell_0 || *it == cell_1));
}

// Clang-tidy false positive with LinearMap concept
#  ifndef CHOMP_CLANG_TIDY

TEMPLATE_LIST_TEST_CASE(
    "Linear function interface to modules", "[algebra]", ModuleTypes
) {
  using M = std::tuple_element_t<0, TestType>;
  using T = typename M::BasisType;
  using R = typename M::RingType;
  const T cell_0 = std::tuple_element_t<1, TestType>();
  T cell_1 = std::tuple_element_t<2, TestType>();
  M elem_0;

  std::vector<std::pair<T, R>> result;
  LinearMap<M, typename std::vector<std::pair<T, R>>::iterator> lfunc =
      [&](const T& cell) {
        result.clear();
        result.push_back(std::make_pair(cell, one<R>()));
        result.push_back(std::make_pair(cell_0, zero<R>()));
        result.push_back(std::make_pair(cell_1, one<R>()));
        return std::make_pair(result.begin(), result.end());
      };

  M elem_1 = linear_apply(elem_0, lfunc);
  REQUIRE(elem_0[cell_0] == zero<R>());
  REQUIRE(elem_0[cell_1] == zero<R>());
  REQUIRE(elem_1[cell_0] == zero<R>());
  REQUIRE(elem_1[cell_1] == zero<R>());

  elem_0.insert(cell_0, -one<R>());
  elem_1 = linear_apply(elem_0, lfunc);
  REQUIRE(elem_0[cell_0] == -one<R>());
  REQUIRE(elem_0[cell_1] == zero<R>());
  REQUIRE(elem_1[cell_0] == -one<R>());
  REQUIRE(elem_1[cell_1] == -one<R>());

  elem_0.insert(cell_0, -one<R>());
  elem_0.insert(cell_1, one<R>());
  elem_1 = linear_apply(elem_0, lfunc);
  REQUIRE(elem_0[cell_0] == -one<R>() - one<R>());
  REQUIRE(elem_0[cell_1] == one<R>());
  REQUIRE(elem_1[cell_0] == -one<R>() - one<R>());
  REQUIRE(elem_1[cell_1] == zero<R>());
}

#  endif  // CHOMP_CLANG_TIDY

TEMPLATE_LIST_TEST_CASE(
    "Comparison operators on modules", "[algebra]", ModuleTypes
) {
  using M = std::tuple_element_t<0, TestType>;
  using T = typename M::BasisType;
  using R = typename M::RingType;
  const T cell_0 = std::tuple_element_t<1, TestType>();
  M elem_0;
  M elem_1;

  REQUIRE(elem_0 == elem_1);
  REQUIRE_FALSE(elem_0 != elem_1);

  elem_0.insert(cell_0, zero<R>());
  REQUIRE(elem_0 == elem_1);
  REQUIRE_FALSE(elem_0 != elem_1);

  elem_0.insert(cell_0, one<R>());
  REQUIRE_FALSE(elem_0 == elem_1);
  REQUIRE(elem_0 != elem_1);

  elem_1.insert(cell_0, one<R>());
  REQUIRE(elem_0 == elem_1);
  REQUIRE_FALSE(elem_0 != elem_1);
}

TEMPLATE_LIST_TEST_CASE(
    "Arithmetic operators on modules", "[algebra]", ModuleTypes
) {
  using M = std::tuple_element_t<0, TestType>;
  using T = typename M::BasisType;
  using R = typename M::RingType;
  const T cell_0 = std::tuple_element_t<1, TestType>();
  T cell_1 = std::tuple_element_t<2, TestType>();  // NOLINT
  M elem_0;
  M elem_1;

  M elem_n = -elem_0;
  M elem_a = elem_0 + elem_1;
  M elem_s = elem_0 - elem_1;
  REQUIRE(elem_0[cell_0] == zero<R>());
  REQUIRE(elem_0[cell_1] == zero<R>());
  REQUIRE(elem_1[cell_0] == zero<R>());
  REQUIRE(elem_1[cell_1] == zero<R>());
  REQUIRE(elem_n[cell_0] == zero<R>());
  REQUIRE(elem_n[cell_1] == zero<R>());
  REQUIRE(elem_a[cell_0] == zero<R>());
  REQUIRE(elem_a[cell_1] == zero<R>());
  REQUIRE(elem_s[cell_0] == zero<R>());
  REQUIRE(elem_s[cell_1] == zero<R>());
  REQUIRE(elem_s == elem_0 + (-elem_1));
  REQUIRE(elem_a == -elem_n + elem_1);

  elem_0.insert(cell_0, one<R>());
  elem_1.insert(cell_1, one<R>());

  elem_n = -elem_0;
  elem_a = elem_0 + elem_1;
  elem_s = elem_0 - elem_1;
  REQUIRE(elem_0[cell_0] == one<R>());
  REQUIRE(elem_0[cell_1] == zero<R>());
  REQUIRE(elem_1[cell_0] == zero<R>());
  REQUIRE(elem_1[cell_1] == one<R>());
  REQUIRE(elem_n[cell_0] == -one<R>());
  REQUIRE(elem_n[cell_1] == zero<R>());
  REQUIRE(elem_a[cell_0] == one<R>());
  REQUIRE(elem_a[cell_1] == one<R>());
  REQUIRE(elem_s[cell_0] == one<R>());
  REQUIRE(elem_s[cell_1] == -one<R>());
  REQUIRE(elem_s == elem_0 + (-elem_1));
  REQUIRE(elem_a == -elem_n + elem_1);

  elem_0 *= one<R>() + one<R>();
  elem_1 = -one<R>() * elem_1 * -one<R>();

  elem_n = -elem_0;
  elem_a = elem_0 + elem_1;
  elem_s = elem_0 - elem_1;
  REQUIRE(elem_0[cell_0] == one<R>() + one<R>());
  REQUIRE(elem_0[cell_1] == zero<R>());
  REQUIRE(elem_1[cell_0] == zero<R>());
  REQUIRE(elem_1[cell_1] == one<R>());
  REQUIRE(elem_n[cell_0] == -(one<R>() + one<R>()));
  REQUIRE(elem_n[cell_1] == zero<R>());
  REQUIRE(elem_a[cell_0] == one<R>() + one<R>());
  REQUIRE(elem_a[cell_1] == one<R>());
  REQUIRE(elem_s[cell_0] == one<R>() + one<R>());
  REQUIRE(elem_s[cell_1] == -one<R>());
  REQUIRE(elem_s == elem_0 + (-elem_1));
  REQUIRE(elem_a == -elem_n + elem_1);

  // NOLINTNEXTLINE
  REQUIRE(elem_a + std::move(elem_s) == (one<R>() + one<R>()) * elem_0);

  elem_0 += elem_1;
  REQUIRE(elem_0[cell_0] == one<R>() + one<R>());
  REQUIRE(elem_0[cell_1] == one<R>());
  elem_1 -= std::move(elem_0);
  REQUIRE(elem_1[cell_0] == -one<R>() - one<R>());
  REQUIRE(elem_1[cell_1] == zero<R>());

  REQUIRE(elem_1 * -one<R>() == -elem_1);
}

using CellAndRingTypes = std::tuple<
    std::tuple<int, Z<2>, UnorderedSetModule<int, Z<2>>>,
    std::tuple<int, Z<3>, UnorderedMapModule<int, Z<3>>>,
    std::tuple<std::vector<short>, Z<2>, SetModule<std::vector<short>, Z<2>>>,
    std::tuple<std::vector<short>, Z<5>, MapModule<std::vector<short>, Z<5>>>>;

TEMPLATE_LIST_TEST_CASE(
    "DefaultModule chooses Module type correctly", "[algebra]", CellAndRingTypes
) {
  using T = std::tuple_element_t<0, TestType>;
  using R = std::tuple_element_t<1, TestType>;
  using M = std::tuple_element_t<2, TestType>;

  CHECK(std::same_as<DefaultModule<T, R>, M>);
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
