/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains a helper function for interfacing with the core
 * of CHomP3R to obtain the cohomology generators needed for cycling signature
 * computation.
 */

#ifndef CHOMP_CYCLING_UTIL_GENERATORS_HPP
#define CHOMP_CYCLING_UTIL_GENERATORS_HPP

#include <chomp/core/complexes/complexes.hpp>
#include <chomp/core/complexes/morse.hpp>
#include <chomp/core/homology/morse.hpp>
#include <chomp/core/util/concepts.hpp>

#include <concepts>
#include <memory>
#include <ranges>
#include <vector>

namespace chomp::cycling::detail {

/**
 * Compute generators of the cohomology of the chain complex `complex` subject
 * to the condition `cond`.
 *
 * @tparam PM1 The partial matching type computed for `complex`. Following
 * partial matchings (computed on the Morse complexes) are of type `PM2`.
 * @tparam PM2 The partial matching type computed after the first, on the
 * type `MorseComplex<CC>`.
 */
template <core::ChainComplex CC, typename PM1, typename PM2,
    template <typename...> typename MapType = core::DefaultMap>
requires std::derived_from<PM1, core::PartialMatching<CC, MapType>>
         && std::derived_from<PM2,
             core::PartialMatching<core::MorseComplex<CC>, MapType>>
         && std::constructible_from<PM1, std::shared_ptr<CC>>
         && std::constructible_from<PM2,
             std::shared_ptr<core::MorseComplex<CC>>>
[[nodiscard]] std::vector<typename CC::ChainType> compute_cycling_generators(
    const std::shared_ptr<CC> complex,
    const core::ConditionalType<typename CC::CellType>& cond) {
  using MC = core::MorseComplex<CC>;
  using PM1Base = core::PartialMatching<CC, MapType>;
  using PM2Base = core::PartialMatching<MC, MapType>;
  using RingType = typename CC::RingType;
  using CellType = typename CC::CellType;
  using ChainType = typename CC::ChainType;

  std::shared_ptr<PM1Base> top_matching = std::make_shared<PM1>(complex);
  std::shared_ptr<PM2Base> next_matching;
  std::vector<std::shared_ptr<PM2Base>> matchings;
  std::shared_ptr<MC> prev_complex;
  std::shared_ptr<MC> next_complex = std::make_shared<MC>(top_matching);

  while (true) {
    prev_complex = next_complex;
    next_matching = std::make_shared<PM2>(prev_complex);
    next_complex = std::make_shared<MC>(next_matching);
    if (prev_complex->size() == next_complex->size()) {
      break;
    }
    matchings.push_back(next_matching);
  }

  std::vector<ChainType> cogenerators;
  for (const CellType& cell : *prev_complex) {
    if (cond(cell)) {
      ChainType cochain;
      cochain.insert(cell, core::one<RingType>());
      cogenerators.push_back(std::move(cochain));
    }
  }

  for (std::shared_ptr<PM2Base> match : matchings | std::views::reverse) {
    for (ChainType& cogen : cogenerators) {
      cogen = match->colift(cogen);
    }
  }
  for (ChainType& cogen : cogenerators) {
    cogen = top_matching->colift(cogen);
  }

  return cogenerators;
}

}  // namespace chomp::cycling::detail

#endif  // CHOMP_CYCLING_UTIL_GENERATORS_HPP
