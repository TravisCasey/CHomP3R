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

}  // namespace chomp::core

#endif  // CHOMP_UTIL_CONCEPTS_H
