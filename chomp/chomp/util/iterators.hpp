/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains custom iterators used throughout CHomP.
 */

#ifndef CHOMP_UTIL_ITERATORS_H
#define CHOMP_UTIL_ITERATORS_H

#include <array>
#include <cstddef>
#include <iterator>
#include <tuple>
#include <type_traits>

namespace chomp::core {

/**
 * @brief Constant iterator wrapper for associative containers whose iterators
 * yield key, value pairs.
 *
 * This forward iterator only fetches `const`-qualified references to keys from
 * an associative container when dereferenced rather than key, value pairs. This
 * is used for `UnorderedMapModule` and `MapModule` class templates to maintain
 * a consistent interface with the set-based module class implementations.
 *
 * @tparam I A forward iterator type returning pairs (or tuples) when
 * dereferenced.
 */
template <typename I>
requires std::forward_iterator<I> && requires(I it) { std::get<0>(*it); }
class KeyIterator {
private:
  I it;

public:
  /** @brief Difference type between iterators. */
  using difference_type = std::iter_difference_t<I>;
  /** @brief Value type when dereferenced. */
  using value_type = const std::tuple_element_t<0, std::iter_value_t<I>>;
  /** @brief Pointer type. */
  using pointer = value_type*;
  /** @brief Reference type. */
  using reference = value_type&;
  /** @brief Tag for `iterator_traits` */
  using iterator_concept = std::forward_iterator_tag;

  /**
   * @brief Default initialize a new KeyIterator object.
   *
   * Unusable in this state but can be assigned to, as normal.
   */
  KeyIterator() = default;
  /**
   * @brief Construct a new KeyIterator object using a passed iterator.
   *
   * @param it
   */
  explicit KeyIterator(I it) : it(it) {};

  /**
   * @brief Dereferencing this iterator yields a constant reference to
   * the key in the map.
   *
   * @return reference
   */
  [[nodiscard]] reference operator*() const {
    return std::get<0>(*it);
  }

  /**
   * @brief Equality operates on wrapped pointers.
   *
   * @param rhs
   * @return true
   * @return false
   */
  [[nodiscard]] bool operator==(const KeyIterator& rhs) const {
    return it == rhs.it;
  }

  /**
   * @brief Preincrememnt operates on wrapped iterator.
   *
   * @return KeyIterator&
   */
  KeyIterator& operator++() {
    ++it;
    return *this;
  }
  /**
   * @brief Postincrement operates on wrapped iterator.
   *
   */
  KeyIterator operator++(int) {
    KeyIterator temp = *this;
    ++it;
    return temp;
  }
};

#ifndef CHOMP_DOXYGEN

namespace detail {

// Concatenate (type-wise) arbitrary tuples (INPUT) together.
template <typename... INPUT>
using TupleConcat = decltype(std::tuple_cat(std::declval<INPUT>()...));

// Metaprogramming If used to determine if original (value) has a bit set;
// If it does, it returns (a tuple containing) the value without the bit set
// and the original value to continue iterating on.
template <typename original, std::size_t bit_and, bool result>
struct BitIf {
  using generated = std::tuple<std::integral_constant<std::size_t, bit_and>>;
  using continued = std::tuple<original>;
};
// Otherwise, returns empty tuples. See `BFSLoop`.
template <typename original, std::size_t bit_and>
struct BitIf<original, bit_and, false> {
  using generated = std::tuple<>;
  using continued = std::tuple<>;
};

// Empty primary template used for extracting tuple parameter pack.
template <std::size_t, std::size_t, typename>
struct BFSLoop;
// Iterate over each bit position and produce further elements in a BFS manner
// from UpperElemets.
template <
    std::size_t N, std::size_t M, template <typename...> typename UpperTuple,
    typename... UpperElements>
struct BFSLoop<N, M, UpperTuple<UpperElements...>> {
  static constexpr std::size_t bit_flag = 1 << M;
  using generated = TupleConcat<typename BitIf<
      UpperElements, UpperElements() ^ bit_flag,
      bool(UpperElements() & bit_flag)>::generated...>;
  using continued = TupleConcat<typename BitIf<
      UpperElements, UpperElements() ^ bit_flag,
      bool(UpperElements() & bit_flag)>::continued...>;
  using type =
      TupleConcat<generated, typename BFSLoop<N, M + 1, continued>::type>;
};
// Terminate when M = N.
template <
    std::size_t N, template <typename...> typename UpperTuple,
    typename... UpperElements>
struct BFSLoop<N, N, UpperTuple<UpperElements...>> {
  using type = std::tuple<>;
};

// Primary template used for extracting tuple parameter pack.
template <std::size_t, typename>
struct ArrayCreator;
// Generate an array holding `Elements`.
template <
    std::size_t N, template <typename...> typename ElementsTuple,
    typename... Elements>
struct ArrayCreator<N, ElementsTuple<Elements...>> {
  static constexpr std::array<std::size_t, (1 << N)> value = {Elements()...};
};

// The exterior loop over dimension; the interior loop over axis is run by
// `BFSLoop` struct. See `BFS_CUBE_ARRAY`.
template <std::size_t N, typename MEM, typename PREV>
consteval std::array<std::size_t, (1 << N)> generate_next() {
  if constexpr (std::tuple_size_v<PREV> == 0) {
    return ArrayCreator<N, MEM>::value;
  }
  return generate_next<
      N, TupleConcat<MEM, PREV>, typename BFSLoop<N, 0, PREV>::type>();
}


}  // namespace detail

#endif  // CHOMP_DOXYGEN

/**
 * @brief The integers from `(1 << N) - 1` to `0` ordered in a breadth-first
 * manner.
 *
 * Precisely, these integers are first ordered by popcount (number of bits set)
 * and otherwise by their value as integers. This computation is done at
 * compile time and allows for easy iteration in this order without the use of
 * runtime queues.
 *
 * @tparam N The number of bits used; and `(1 << N)` is the length of the array.
 */
template <std::size_t N>
constexpr std::array<std::size_t, (1 << N)> BFS_CUBE_ARRAY =
    detail::generate_next<
        N, std::tuple<>,
        std::tuple<std::integral_constant<std::size_t, (1 << N) - 1>>>();


}  // namespace chomp::core

#endif  // CHOMP_UTIL_ITERATORS_H
