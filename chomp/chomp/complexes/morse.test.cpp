/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/algebra/algebra.hpp>
#include <chomp/algebra/cyclic.hpp>
#include <chomp/algebra/modules.hpp>
#include <chomp/complexes/cubical.hpp>
#include <chomp/complexes/morse.hpp>
#include <chomp/homology/cubical.hpp>
#include <chomp/homology/morse.hpp>
#include <chomp/util/concepts.hpp>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

template <std::size_t DIM, Ring R, Module M,
    template <typename...> typename SetType = DefaultSet>
std::shared_ptr<
    CubicalComplex<DIM, TopCubeSetGrading<DIM, 0, 1, SetType>, R, M>>
initialize_top_cubical_sphere_complex() {
  Orthant<DIM> maximum{};
  Orthant<DIM> hole{};
  for (std::size_t axis = 0; axis != DIM; ++axis) {
    maximum[axis] = 2;
    hole[axis] = 1;
  }

  SetType<Orthant<DIM>> grading_set;
  LexOrthantIterator<DIM> lex_it(Orthant<DIM>{}, maximum);
  for (; lex_it; ++lex_it) {
    if (*lex_it != hole) {
      grading_set.insert(*lex_it);
    }
  }
  const TopCubeSetGrading<DIM, 0, 1, SetType> grading_func(
      std::move(grading_set));

  for (std::size_t axis = 0; axis != DIM; ++axis) {
    ++maximum[axis];
  }
  return std::make_shared<
      CubicalComplex<DIM, TopCubeSetGrading<DIM, 0, 1, SetType>, R, M>>(maximum,
      grading_func);
}

template <std::size_t, typename>
struct TupleAppend;

template <std::size_t NewInt, template <typename...> typename TUP,
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

using dimension_sequence = typename TupleSequence<4, 5, 6, 7, 8>::type;

TEMPLATE_LIST_TEST_CASE("Top cubical sphere benchmarks", "[!benchmark]",
    dimension_sequence) {
  constexpr std::size_t DIM = TestType();
  using RingType = Z<2>;
  using GradingType = TopCubeSetGrading<DIM, 0, 1, DefaultSet>;
  using ChainType = DefaultModule<Cube<DIM>, RingType>;
  using ComplexType = CubicalComplex<DIM, GradingType, RingType, ChainType>;

  std::shared_ptr<ComplexType> complex
      = initialize_top_cubical_sphere_complex<DIM, RingType, ChainType>();

  auto compute_reduction =
      [](std::shared_ptr<ComplexType> input_complex) {
        constexpr std::size_t CACHE_SPAN = 2;
        const std::shared_ptr<PartialMatching<ComplexType, DefaultMap>>
            morse_matching = std::make_shared<
                CubicalMatching<ComplexType, DefaultMap, CACHE_SPAN>>(
                input_complex, 0);
        const std::shared_ptr<MorseComplex<ComplexType>> morse_complex
            = std::make_shared<MorseComplex<ComplexType>>(morse_matching);
        return full_reduce(morse_complex);
      };

  BENCHMARK("Top Cubical S" + std::to_string(DIM - 1) + " in "
            + std::to_string(DIM) + " dimensions") {
    return compute_reduction(complex);
  };

  SECTION("Verify totally-reduced complex") {
    const std::shared_ptr<MorseComplex<MorseComplex<ComplexType>>> reduced =
    compute_reduction(complex);
    std::set<std::size_t> correct_dimensions{0, DIM - 1};

    CHECK(reduced->size() == 2);

    std::set<std::size_t> dimensions;
    for (const Cube<DIM>& cell : *reduced) {
      dimensions.insert(cell.dimension());
      CHECK(boundary(*reduced, cell) == ChainType());
      CHECK(coboundary(*reduced, cell) == ChainType());
    }
    CHECK(dimensions == correct_dimensions);
  }
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
