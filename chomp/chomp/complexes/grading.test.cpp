/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/complexes/grading.hpp>
#include <chomp/util/constants.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <utility>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

struct GradingTest {
  using InputType = std::size_t;

  GradingResultType operator()(const InputType& input) {
    if (input > 30) {
      return 4;
    }
    if (input > 10) {
      return 2;
    }
    return 12;
  }
};

struct LBGradingTest {
  using InputType = std::size_t;
  using Minimum = std::integral_constant<GradingResultType, 2>;

  GradingResultType operator()(const InputType& input) {
    if (input > 30) {
      return 4;
    }
    if (input > 10) {
      return 2;
    }
    return 12;
  }
};

struct UBGradingTest {
  using InputType = std::size_t;
  using Maximum = std::integral_constant<GradingResultType, 12>;

  GradingResultType operator()(const InputType& input) {
    if (input > 30) {
      return 4;
    }
    if (input > 10) {
      return 2;
    }
    return 12;
  }
};

struct BGradingTest {
  using InputType = std::size_t;
  using Minimum = std::integral_constant<GradingResultType, 2>;
  using Maximum = std::integral_constant<GradingResultType, 12>;

  GradingResultType operator()(const InputType& input) {
    if (input > 30) {
      return 4;
    }
    if (input > 10) {
      return 2;
    }
    return 12;
  }
};


TEST_CASE("Grading concepts differentiate correctly", "[complexes]") {
  CHECK(Grading<GradingTest>);
  CHECK(Grading<LBGradingTest>);
  CHECK(Grading<UBGradingTest>);
  CHECK(Grading<BGradingTest>);

  CHECK_FALSE(LowerBoundedGrading<GradingTest>);
  CHECK(LowerBoundedGrading<LBGradingTest>);
  CHECK_FALSE(LowerBoundedGrading<UBGradingTest>);
  CHECK(LowerBoundedGrading<BGradingTest>);

  CHECK_FALSE(UpperBoundedGrading<GradingTest>);
  CHECK_FALSE(UpperBoundedGrading<LBGradingTest>);
  CHECK(UpperBoundedGrading<UBGradingTest>);
  CHECK(UpperBoundedGrading<BGradingTest>);

  CHECK_FALSE(BoundedGrading<GradingTest>);
  CHECK_FALSE(BoundedGrading<LBGradingTest>);
  CHECK_FALSE(BoundedGrading<UBGradingTest>);
  CHECK(BoundedGrading<BGradingTest>);
}

TEST_CASE(
    "CachedGradingWrapper propagates grading concepts correctly", "[complexes]"
) {
  CHECK(Grading<CachedGradingWrapper<GradingTest>>);
  CHECK(Grading<CachedGradingWrapper<LBGradingTest>>);
  CHECK(Grading<CachedGradingWrapper<UBGradingTest>>);
  CHECK(Grading<CachedGradingWrapper<BGradingTest>>);

  CHECK_FALSE(LowerBoundedGrading<CachedGradingWrapper<GradingTest>>);
  CHECK(LowerBoundedGrading<CachedGradingWrapper<LBGradingTest>>);
  CHECK_FALSE(LowerBoundedGrading<CachedGradingWrapper<UBGradingTest>>);
  CHECK(LowerBoundedGrading<CachedGradingWrapper<BGradingTest>>);

  CHECK_FALSE(UpperBoundedGrading<CachedGradingWrapper<GradingTest>>);
  CHECK_FALSE(UpperBoundedGrading<CachedGradingWrapper<LBGradingTest>>);
  CHECK(UpperBoundedGrading<CachedGradingWrapper<UBGradingTest>>);
  CHECK(UpperBoundedGrading<CachedGradingWrapper<BGradingTest>>);

  CHECK_FALSE(BoundedGrading<CachedGradingWrapper<GradingTest>>);
  CHECK_FALSE(BoundedGrading<CachedGradingWrapper<LBGradingTest>>);
  CHECK_FALSE(BoundedGrading<CachedGradingWrapper<UBGradingTest>>);
  CHECK(BoundedGrading<CachedGradingWrapper<BGradingTest>>);
}

TEST_CASE("SetGrading functions correctly", "[complexes]") {
  const std::initializer_list<int> ilist = {-2, 5, 3, 10, 1, 0, 4};
  constexpr std::size_t MIN = 4;
  constexpr std::size_t MAX = 10;
  const SetGrading<int, MIN, MAX> grade_func(ilist);

  CHECK(Grading<decltype(grade_func)>);
  CHECK(LowerBoundedGrading<decltype(grade_func)>);
  CHECK(UpperBoundedGrading<decltype(grade_func)>);
  CHECK(BoundedGrading<decltype(grade_func)>);

  for (const int cell : ilist) {
    CHECK(grade_func(cell) == MIN);
  }

  CHECK(grade_func(-1) == MAX);
  CHECK(grade_func(-3) == MAX);
}

TEST_CASE("MapGrading functions correctly", "[complexes]") {
  const std::initializer_list<std::pair<const int, GradingResultType>> ilist = {
      std::make_pair(0, 4), std::make_pair(1, 5), std::make_pair(2, 7),
      std::make_pair(20, 10)
  };
  constexpr std::size_t MIN = 4;
  constexpr std::size_t MAX = 10;
  const MapGrading<int, MIN, MAX> grade_func(ilist);

  CHECK(Grading<decltype(grade_func)>);
  CHECK(LowerBoundedGrading<decltype(grade_func)>);
  CHECK(UpperBoundedGrading<decltype(grade_func)>);
  CHECK(BoundedGrading<decltype(grade_func)>);

  for (const std::pair<const int, GradingResultType>& cell_pair : ilist) {
    CHECK(grade_func(cell_pair.first) == cell_pair.second);
  }

  CHECK(grade_func(-1) == MAX);
  CHECK(grade_func(-3) == MAX);
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
