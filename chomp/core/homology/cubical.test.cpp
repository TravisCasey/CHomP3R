/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/core/algebra/algebra.hpp>
#include <chomp/core/algebra/cyclic.hpp>
#include <chomp/core/complexes/cubical.hpp>
#include <chomp/core/homology/cubical.hpp>
#include <chomp/core/util/concepts.hpp>
#include <chomp/core/util/morse.hpp>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <map>
#include <memory>
#include <set>
#include <tuple>
#include <type_traits>
#include <utility>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

TEST_CASE("CubicalMatching class.", "[homology]") {
  SECTION("2-dimensional S1 with Z2 coefficients") {
    using GradingType = TopCubeSetGrading<2, 0, 1>;
    using ComplexType = CubicalComplex<2, GradingType>;
    using MatchingType = CubicalMatching<ComplexType>;
    using ChainType = typename ComplexType::ChainType;

    const GradingType grading_func{
        {0, 0},
        {0, 1},
        {0, 2},
        {1, 0},
        {1, 2},
        {2, 0},
        {2, 1},
        {2, 2}
    };
    const std::shared_ptr<ComplexType> complex =
        std::make_shared<ComplexType>(Orthant<2>{3, 3}, grading_func);
    MatchingType matching(complex);

    SECTION("Morse (co)boundary operator and critical cells") {
      const std::set<Cube<2>> correct_zero_grades{
          {{1, 1}, 0b01},
          {{3, 3}, 0b00}
      };
      const std::set<Cube<2>> correct_one_grades{
          {{1, 1}, 0b11},
          {{3, 3}, 0b01}
      };
      std::set<Cube<2>> zero_grades;
      std::set<Cube<2>> one_grades;

      const auto [boundaries, coboundaries] = matching.compute_operators();
      for (const auto& [ace, boundary_chain] : boundaries) {
        if (grading_func(ace) == 0) {
          CHECK(boundary_chain == ChainType());
          zero_grades.insert(ace);
        } else {
          one_grades.insert(ace);
        }
      }

      CHECK(correct_zero_grades == zero_grades);
      CHECK(correct_one_grades == one_grades);
    }

    SECTION("Match queries") {
      MatchResult<Cube<2>, Z<2>> correct;

      correct = MatchResult<Cube<2>, Z<2>>(
          Cube<2>({0, 0}, 0b01), one<Z<2>>(), Trichotomy::queen
      );
      CHECK(matching.match(Cube<2>({0, 0}, 0b00)) == correct);

      correct = MatchResult<Cube<2>, Z<2>>(
          Cube<2>({0, 0}, 0b10), one<Z<2>>(), Trichotomy::king
      );
      CHECK(matching.match(Cube<2>({0, 0}, 0b11)) == correct);

      correct = MatchResult<Cube<2>, Z<2>>(
          Cube<2>({1, 1}, 0b01), one<Z<2>>(), Trichotomy::ace
      );
      CHECK(matching.match(Cube<2>({1, 1}, 0b01)) == correct);

      correct = MatchResult<Cube<2>, Z<2>>(
          Cube<2>({1, 1}, 0b11), one<Z<2>>(), Trichotomy::ace
      );
      CHECK(matching.match(Cube<2>({1, 1}, 0b11)) == correct);
    }
  }

  SECTION("2-dimensional figure eight and extra cubes with Z2 coefficients") {
    using GradingType = TopCubeSetGrading<2, 0, 1>;
    using ComplexType = CubicalComplex<2, GradingType>;
    using MatchingType = CubicalMatching<ComplexType>;
    using GradeLimitedMatchingType = CubicalMatching<ComplexType, 0>;
    using ChainType = typename ComplexType::ChainType;
    using RingType = typename ComplexType::RingType;

    const GradingType grading_func{
        {0, 0},
        {1, 0},
        {2, 0},
        {3, 0},
        {4, 0},
        {0, 1},
        {2, 1},
        {4, 1},
        {0, 2},
        {1, 2},
        {2, 2},
        {3, 2},
        {4, 2},
        {0, 3},
        {2, 3},
        {4, 3}
    };
    const std::shared_ptr<ComplexType> complex =
        std::make_shared<ComplexType>(Orthant<2>{7, 7}, grading_func);
    MatchingType matching(complex);
    GradeLimitedMatchingType max_grade_matching(complex);

    SECTION("Morse (co)boundary operator and critical cells") {
      const std::set<Cube<2>> correct_zero_grades{
          {{1, 1}, 0b01},
          {{3, 1}, 0b01},
          {{1, 3}, 0b01},
          {{3, 3}, 0b01},
          {{1, 4}, 0b00},
          {{3, 4}, 0b00},
          {{5, 4}, 0b00}
      };
      const std::set<Cube<2>> correct_one_grades{
          {{1, 1}, 0b11},
          {{3, 1}, 0b11},
          {{1, 3}, 0b11},
          {{3, 3}, 0b11},
          {{1, 4}, 0b01},
          {{3, 4}, 0b01},
          {{5, 4}, 0b01}
      };
      std::set<Cube<2>> zero_grades;
      std::set<Cube<2>> one_grades;

      std::map<Cube<2>, ChainType> zero_boundaries{
          {Cube<2>({1, 1}, 0b01), ChainType()},
          {Cube<2>({3, 1}, 0b01), ChainType()},
          {Cube<2>({1, 4}, 0b00), ChainType()},
          {Cube<2>({3, 4}, 0b00), ChainType()},
          {Cube<2>({5, 4}, 0b00), ChainType()}
      };

      ChainType chain;
      chain.insert(Cube<2>({1, 4}, 0b00), one<RingType>());
      chain.insert(Cube<2>({3, 4}, 0b00), one<RingType>());
      zero_boundaries[Cube<2>({1, 3}, 0b01)] = chain;

      chain = ChainType();
      chain.insert(Cube<2>({3, 4}, 0b00), one<RingType>());
      chain.insert(Cube<2>({5, 4}, 0b00), one<RingType>());
      zero_boundaries[Cube<2>({3, 3}, 0b01)] = chain;

      const auto [boundaries, coboundaries] = matching.compute_operators();
      for (const auto& [ace, boundary_chain] : boundaries) {
        if (grading_func(ace) == 0) {
          CHECK(boundary_chain == zero_boundaries[ace]);
          zero_grades.insert(ace);
        } else {
          one_grades.insert(ace);
        }
      }

      CHECK(correct_zero_grades == zero_grades);
      CHECK(correct_one_grades == one_grades);

      std::set<Cube<2>> max_grade_zero_grades;
      const auto [max_grade_boundaries, max_grade_coboundaries] =
          max_grade_matching.compute_operators();
      for (const auto& [ace, boundary_chain] : max_grade_boundaries) {
        REQUIRE(grading_func(ace) == 0);
        CHECK(boundary_chain == zero_boundaries[ace]);
        max_grade_zero_grades.insert(ace);
      }

      CHECK(correct_zero_grades == max_grade_zero_grades);
    }

    SECTION("Match queries") {
      MatchResult<Cube<2>, Z<2>> correct;

      correct = MatchResult<Cube<2>, Z<2>>(
          Cube<2>({1, 1}, 0b10), one<Z<2>>(), Trichotomy::queen
      );
      CHECK(matching.match(Cube<2>({1, 1}, 0b00)) == correct);

      correct = MatchResult<Cube<2>, Z<2>>(
          Cube<2>({1, 2}, 0b10), one<Z<2>>(), Trichotomy::king
      );
      CHECK(matching.match(Cube<2>({1, 2}, 0b11)) == correct);

      correct = MatchResult<Cube<2>, Z<2>>(
          Cube<2>({1, 3}, 0b01), one<Z<2>>(), Trichotomy::ace
      );
      CHECK(matching.match(Cube<2>({1, 3}, 0b01)) == correct);

      correct = MatchResult<Cube<2>, Z<2>>(
          Cube<2>({3, 3}, 0b11), one<Z<2>>(), Trichotomy::ace
      );
      CHECK(matching.match(Cube<2>({3, 3}, 0b11)) == correct);
    }
  }
}

template <std::size_t, typename>
struct TupleAppend;

template <
    std::size_t NewInt, template <typename...> typename TUP,
    typename... IntConsts>
struct TupleAppend<NewInt, TUP<IntConsts...>> {
  using type = TUP<std::integral_constant<std::size_t, NewInt>, IntConsts...>;
};

template <std::size_t...>
struct TupleSequence;

template <std::size_t First, std::size_t... Other>
struct TupleSequence<First, Other...> {
  using type =
      typename TupleAppend<First, typename TupleSequence<Other...>::type>::type;
};

template <>
struct TupleSequence<> {
  using type = std::tuple<>;
};

using dimension_sequence = typename TupleSequence<2, 3, 4, 5, 6, 7>::type;

TEMPLATE_LIST_TEST_CASE(
    "CubicalMatching constructions in each dimension.", "[homology]",
    dimension_sequence
) {
  constexpr std::size_t dim = TestType();

  using GradingType = TopCubeSetGrading<dim, 0, 1>;
  using ComplexType = CubicalComplex<dim, GradingType>;
  using ChainType = typename ComplexType::ChainType;

  SECTION("n-dimensional S1, Z2 coefficients") {
    using MatchingType = CubicalMatching<ComplexType, 0, 1>;

    Orthant<dim> maximum{};
    maximum[0] = 2;
    maximum[1] = 2;
    Orthant<dim> hole{};
    hole[0] = 1;
    hole[1] = 1;

    DefaultSet<Orthant<dim>> grading_set;
    LexOrthantIterator<dim> lex_it(Orthant<dim>{}, maximum);

    for (; lex_it; ++lex_it) {
      if (*lex_it != hole) {
        grading_set.insert(*lex_it);
      }
    }
    const GradingType grading_func(std::move(grading_set));

    std::set<Cube<dim>> correct_zero_grades;
    std::set<Cube<dim>> zero_grades;
    for (std::size_t axis = 0; axis != dim; ++axis) {
      ++maximum[axis];
    }
    correct_zero_grades.insert(Cube<dim>(maximum, 0));
    for (std::size_t axis = 2; axis != dim; ++axis) {
      ++hole[axis];
    }
    correct_zero_grades.insert(Cube<dim>(hole, 1));

    const std::shared_ptr<ComplexType> complex =
        std::make_shared<ComplexType>(maximum, grading_func);
    MatchingType matching(complex);
    const auto [boundaries, coboundaries] = matching.compute_operators();

    for (const auto& [ace, boundary_chain] : boundaries) {
      REQUIRE(grading_func(ace) == 0);
      zero_grades.insert(ace);
      CHECK(boundary_chain == ChainType());
    }
    CHECK(zero_grades == correct_zero_grades);
  }

  SECTION("n-dimensional Sn-1 with Z2 coefficients") {
    using MatchingType = CubicalMatching<ComplexType>;

    Orthant<dim> maximum{};
    Orthant<dim> hole{};
    for (std::size_t axis = 0; axis != dim; ++axis) {
      maximum[axis] = 2;
      hole[axis] = 1;
    }

    DefaultSet<Orthant<dim>> grading_set;
    LexOrthantIterator<dim> lex_it(Orthant<dim>{}, maximum);

    for (; lex_it; ++lex_it) {
      if (*lex_it != hole) {
        grading_set.insert(*lex_it);
      }
    }
    const GradingType grading_func(std::move(grading_set));

    std::set<Cube<dim>> correct_zero_grades;
    std::set<Cube<dim>> zero_grades;
    for (std::size_t axis = 0; axis != dim; ++axis) {
      ++maximum[axis];
    }
    correct_zero_grades.insert(Cube<dim>(maximum, 0));
    correct_zero_grades.insert(Cube<dim>(hole, (1 << (dim - 1)) - 1));

    const std::shared_ptr<ComplexType> complex =
        std::make_shared<ComplexType>(maximum, grading_func);
    MatchingType matching(complex);
    const auto [boundaries, coboundaries] = matching.compute_operators();

    for (const auto& [ace, boundary_chain] : boundaries) {
      if (grading_func(ace) == 0) {
        zero_grades.insert(ace);
        CHECK(boundary_chain == ChainType());
      }
    }
    CHECK(zero_grades == correct_zero_grades);
  }
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
