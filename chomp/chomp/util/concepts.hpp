/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains general-purpose concepts used throughout CHomP3R.
 */

#ifndef CHOMP_UTIL_CONCEPTS_H
#define CHOMP_UTIL_CONCEPTS_H

#include <concepts>
#include <cstddef>
#include <functional>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace chomp::core {

/**
 * @brief Types `T` modeling this concept specialize `std::hash<T>` and are
 * equality comparable.
 *
 * Types modeling this concept can be stored in `std::unordered_set` and can be
 * keys in `std::unordered_map`.
 *
 * @tparam T
 */
template <typename T>
concept Hashable = requires(const T& a) {
  { std::hash<T>{}(a) } -> std::convertible_to<std::size_t>;
  std::equality_comparable<T>;
};

/**
 * @brief Types `T` modeling this concept are equality comparable and their
 * inequality operators define a strict total order.
 *
 * Types modeling this concept can be stored in `std::set` and can be keys in
 * `std::map`.
 *
 * @tparam T
 */
template <typename T>
concept Comparable = std::totally_ordered<T>;

/**
 * @brief Types `T` modeling this concept are either hashable or comparable,
 * enabling their use as keys in associative containers.
 *
 * @tparam T
 */
template <typename T>
concept AssociativeKey = Hashable<T> || Comparable<T>;


#ifndef CHOMP_DOXYGEN
namespace detail {

template <bool H, typename K, typename V>
struct MapChooser {
  using type = std::map<K, V>;
};

template <typename K, typename V>
struct MapChooser<true, K, V> {
  using type = std::unordered_map<K, V>;
};

template <AssociativeKey K, typename V>
struct DefaultMapChooser {
  using type = typename MapChooser<Hashable<K>, K, V>::type;
};

}  // namespace detail
#endif  // CHOMP_DOXYGEN

/**
 * @brief Default choice between `std::map` and `std::unordered_map` depending
 * on concepts fulfilled by `K`.
 *
 * If `K` has an `std::hash` specialization, `std::unordered_map` is preferred.
 * Otherwise, uses `std::map`.
 *
 * @tparam K Key type; must model `AssociativeKey`.
 * @tparam V Value type.
 */
template <AssociativeKey K, typename V>
using DefaultMap = typename detail::DefaultMapChooser<K, V>::type;


#ifndef CHOMP_DOXYGEN
namespace detail {

template <bool H, typename K>
struct SetChooser {
  using type = std::set<K>;
};

template <typename K>
struct SetChooser<true, K> {
  using type = std::unordered_set<K>;
};

template <AssociativeKey K>
struct DefaultSetChooser {
  using type = typename SetChooser<Hashable<K>, K>::type;
};

}  // namespace detail
#endif  // CHOMP_DOXYGEN

/**
 * @brief Default choice between `std::set` and `std::unordered_set` depending
 * on concepts fulfilled by `K`.
 *
 * If `K` has an `std::hash` specialization, `std::unordered_set` is preferred.
 * Otherwise, uses `std::set`.
 *
 * @tparam K Key type; must model `AssociativeKey`.
 */
template <AssociativeKey K>
using DefaultSet = typename detail::DefaultSetChooser<K>::type;

}  // namespace chomp::core

#endif  // CHOMP_UTIL_CONCEPTS_H
