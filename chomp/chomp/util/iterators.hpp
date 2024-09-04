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

#include <iterator>
#include <tuple>

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

}  // namespace chomp::core

#endif  // CHOMP_UTIL_ITERATORS_H
