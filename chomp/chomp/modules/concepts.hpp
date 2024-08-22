/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief Contains algebraic and functional concepts used for defining module
 * classes.
 */

#ifndef MODULES_CONCEPTS_H
#define MODULES_CONCEPTS_H

#include <concepts>
#include <functional>

namespace chomp::modules {

/**
 * @brief The minimal requirements for a class `G` to implement an additive
 * group.
 *
 * This defines the expected interface for group-like object. No checks are
 * made that these operations make sense or that they fulfill group axioms.
 *
 * @tparam G
 */
template <typename G>
concept Group = requires(G a, G b) {
    static_cast<G>(0); // Additive identity
    { -a } -> std::convertible_to<G>;
    { a + b } -> std::convertible_to<G>;
    { a - b } -> std::convertible_to<G>;
    { a == b } -> std::convertible_to<bool>;
    { a != b } -> std::convertible_to<bool>;
    a += b;
    a -= b;
};

/**
 * @brief The minimal requirements for a class `R` to implement a ring.
 *
 * This defines the expected interface for ring-like object. No checks are
 * made that these operations make sense or that they fulfill ring axioms.
 *
 * Importantly, we also assume (but do not explicitly check) that the ring is an
 * integral domain, i.e., it has no zero divisors.
 *
 * @tparam R
 */
template <typename R>
concept Ring = Group<R> && requires(R a, R b) {
    static_cast<R>(1); // multiplicative identity
    { a* b } -> std::convertible_to<R>;
    a *= b;
};

/**
 * @brief Classes implementing this concept represent the ring (field) with two
 * elements.
 *
 * While there is only one such ring, there may be different `Z` data
 * structures implementing it due to the `T` template argument.
 *
 * Module classes over a binary ring can be represented more efficiently by only
 * implicitly storing these coefficients.
 *
 * @tparam R The proposed class implementing a two-element ring.
 */
template <typename R>
concept BinaryRing = Ring<R> &&
    (static_cast<R>(1) + static_cast<R>(1) == static_cast<R>(0));

/**
 * @brief Types `T` modeling this concept define `std::hash<T>` as well as the
 * equality operator.
 *
 * Types modeling this concept can be stored in `std::unordered_set` and can be
 * keys in `std::unordered_map`, enabling the use of `UnorderedSetModule` and
 * `UnorderedMapModule` with cell type `T`.
 *
 * @tparam T
 */
template <typename T>
concept Hashable = requires(T a, T b) {
    { std::hash<T>{}(a) } -> std::convertible_to<std::size_t>;
    { a == b } -> std::same_as<bool>;
};

/**
 * @brief Types `T` modeling this concept have `std::less<T>` defined as well as
 * the equality operator.
 *
 * Types modeling this concept can be stored in `std::set` and can be keys in
 * `std::map`, enabling the use of `SetModule` and `MapModule` with cell type
 * `T`.
 *
 * @tparam T
 */
template <typename T>
concept Comparable = requires(T a, T b) {
    { a < b } -> std::same_as<bool>;
    { a == b } -> std::same_as<bool>;
};

/**
 * @brief Types `T` modeling this concept model either `Hashable` or
 * `Comparable`, enabling their use in one of the module classes.
 *
 * @tparam T
 */
template <typename T>
concept CellType = Hashable<T> || Comparable<T>;

} // namespace chomp::modules

#endif // MODULES_CONCEPTS_H
