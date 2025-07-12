/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains constants and type aliases used throughout
 * CHomP3R.
 */

#ifndef CHOMP_CORE_UTIL_CONSTANTS_HPP
#define CHOMP_CORE_UTIL_CONSTANTS_HPP

#include <cstddef>
#include <cstdint>
#include <limits>

namespace chomp::core {

/**
 * @brief Bit width of the largest unsigned integer representable by
 * `std::size_t` (platform-dependent).
 *
 * As `std::size_t` is used for the extent/shape parameter of `Cube` objects,
 * this is the maximal ambient dimension for `CubicalComplex` objects.
 */
constexpr std::size_t SIZE_T_BITS = std::numeric_limits<std::size_t>::digits;

/**
 * @brief Unsigned integer type used for grading function output.
 */
using GradingResultType = std::size_t;

/**
 * @brief Integer type representing coordinates for expressing the
 * location of an orthant in the hypercube as an ordered tuple.
 *
 * @sa `CubeOrthant`, `CubicalCell`
 */
using HypercubeCoordinate = std::int_fast8_t;

/**
 * @brief Prime value used in custom hash function for `CubicalCell`.
 *
 * @sa `CubicalCell`
 */
constexpr std::size_t CUBE_HASH_PRIME = 11;


}  // namespace chomp::core

#endif  // CHOMP_CORE_UTIL_CONSTANTS_HPP
