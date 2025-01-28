/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief The `MorseComplex` class template defined in this header is the result
 * of Morse reduction (via a `PartialMatching`) on another Chain Complex.
 */

#ifndef CHOMP_COMPLEXES_MORSE_H
#define CHOMP_COMPLEXES_MORSE_H

#include <chomp/algebra/algebra.hpp>
#include <chomp/complexes/complexes.hpp>
#include <chomp/complexes/cubical.hpp>
#include <chomp/complexes/grading.hpp>
#include <chomp/homology/cubical.hpp>
#include <chomp/homology/morse.hpp>
#include <chomp/util/concepts.hpp>
#include <chomp/util/constants.hpp>

#include <concepts>
#include <cstddef>
#include <memory>
#include <tuple>
#include <vector>

#ifndef CHOMP_DOXYGEN

namespace chomp::core::detail {

template <ChainComplex CC,
    template <typename, typename> typename MapType = DefaultMap>
struct PartialMatchingChooser {
  using type = CoreductionMatching<CC, MapType>;
};

template <std::size_t CCDIM, Grading G, Ring R, Module M,
    template <typename, typename> typename MapType>
struct PartialMatchingChooser<CubicalComplex<CCDIM, G, R, M>, MapType> {
  using type = CubicalMatching<CubicalComplex<CCDIM, G, R, M>, MapType>;
};

}  // namespace chomp::core::detail

#endif  // CHOMP_DOXYGEN

namespace chomp::core {

/**
 * @brief Compute and return a (smart pointer to) `PartialMatching` on
 * `complex`.
 *
 * Currently, this uses a specialization of `CubicalMatching` if and only if
 * `CC` is a specialization of `CubicalComplex`; otherwise, this function uses
 * the correct specialization `CoreductionMatching`.
 *
 * @tparam CC The input `ChainComplex` type.
 * @tparam MapType The associative container in which the boundary maps are
 * stored in the computed matching. Default value is `DefaultMap`.
 * @param complex A pointer to the complex on which the matching is computed.
 * @return std::shared_ptr<PartialMatching<CC, MapType>> The computed matching
 * handled through a smart pointer to its abstract parent class.
 */
template <ChainComplex CC,
    template <typename, typename> typename MapType = DefaultMap>
std::shared_ptr<PartialMatching<CC, MapType>> compute_matching(
    std::shared_ptr<CC> complex) {
  return std::make_shared<
      typename detail::PartialMatchingChooser<CC, MapType>::type>(complex);
}

/**
 * @brief A lightweight chain complex class template that is a reduced version
 * of a complex of type `CC`.
 *
 * It is assumed that this complex is small enough for its cells, their
 * boundaries and coboundaries, as well as their grades to be stored explicitly.
 * Its `RingType`, `CellType`, `ChainType`, and `GradingType` are all
 * propagated from `CC`.
 *
 * @tparam CC The original complex type modeling `ChainComplex` that has been
 * reduced to yield an instantiation of this class template.
 * @tparam MapType The associative container in which the boundary maps are
 * stored. Default value is `DefaultMap`.
 */
template <ChainComplex CC,
    template <typename, typename> typename MapType = DefaultMap>
class MorseComplex {
public:
  /** @brief Coefficient ring type propagated from `CC`. */
  using RingType = typename CC::RingType;
  /** @brief Cell type propagated from `CC`. */
  using CellType = typename CC::CellType;
  /** @brief Chain type propagated from `CC`. */
  using ChainType = typename CC::ChainType;
  /** @brief Grading type propagated from `CC`. */
  using GradingType = typename CC::GradingType;
  /** @brief Iterator type returned by `begin` and `end` member functions. */
  using CellIterType = typename std::vector<CellType>::const_iterator;

private:
  std::vector<CellType> cells;
  MapType<CellType, ChainType> boundaries;
  MapType<CellType, ChainType> coboundaries;
  MapType<CellType, GradingResultType> grades;

public:
  /**
   * @brief Construct a `MorseComplex` using the `PartialMatching` pointed to
   * by `match_ptr`.
   *
   * Explicitly passing the `PartialMatching` object to construct the
   * `MorseComplex` specialization allows greater customization of the matching
   * options.
   *
   * @tparam MC The complex type modeling `ChainComplex` that the matching in
   * `match_ptr` is computed on. The type `MC` either `CC` or
   * `MorseComplex<CC, MapType>`.
   * @param match_ptr The `PartialMatching` object used to reduce the complex
   * of type `MC` to the `MorseComplex` being constructed.
   */
  template <ChainComplex MC>
  requires std::same_as<MC, CC> || std::same_as<MC, MorseComplex<CC, MapType>>
  MorseComplex(std::shared_ptr<PartialMatching<MC, MapType>> match_ptr) {
    cells.assign(match_ptr->begin(), match_ptr->end());
    std::tie(boundaries, coboundaries) = match_ptr->compute_operators();
    for (const CellType& cell : cells) {
      grades[cell] = match_ptr->upper()->grade(cell);
    }
  }

  /**
   * @brief Construct a `MorseComplex` directly from the complex to be reduced.
   *
   * This computes a matching (via the `compute_matching` function template)
   * then constructs a `MorseComplex` instantiation from that matching. For a
   * greater degree of control over the computed matching, the other constructor
   * of this class template accepts the matching instead.
   *
   * @tparam MC The type pointed to by `complex` which models `ChainComplex`
   * @param complex The input `ChainComplex` object being reduced.
   */
  template <ChainComplex MC>
  requires std::same_as<MC, CC> || std::same_as<MC, MorseComplex<CC, MapType>>
  MorseComplex(std::shared_ptr<MC> complex) :
      MorseComplex(compute_matching<MC, MapType>(complex)) {}

  /** @brief Return an iterator to the cells of this complex. */
  [[nodiscard]] CellIterType begin() const noexcept {
    return cells.cbegin();
  }
  /** @brief Return a past-the-end iterator to the cells of this complex. */
  [[nodiscard]] CellIterType end() const noexcept {
    return cells.cend();
  }

  /** @brief The number of cells in this complex. */
  [[nodiscard]] std::size_t size() const noexcept {
    return cells.size();
  }

  /** @brief Return the grade of `cell`. */
  [[nodiscard]] GradingResultType grade(const CellType& cell) {
    return grades[cell];
  }

  /** @brief Return the boundary of `cell` subject to the condition `cond`. */
  [[nodiscard]] ChainType boundary_if(const CellType& cell,
      const ConditionalType<CellType>& cond) {
    ChainType full_boundary = boundaries[cell];
    ChainType result;
    for (const CellType& cell : full_boundary) {
      if (cond(cell)) {
        result.insert(cell, full_boundary[cell]);
      }
    }
    return result;
  }
  /** @brief Return the coboundary of `cell` subject to the condition `cond`. */
  [[nodiscard]] ChainType coboundary_if(const CellType& cell,
      const ConditionalType<CellType>& cond) {
    ChainType full_coboundary = coboundaries[cell];
    ChainType result;
    for (const CellType& cell : full_coboundary) {
      if (cond(cell)) {
        result.insert(cell, full_coboundary[cell]);
      }
    }
    return result;
  }
};

/**
 * @brief Given a (smart pointer to the) a cell complex, iterate Morse reduction
 * until it terminates. At this point, other methods must be used to compute the
 * homology.
 *
 * @tparam CC The input complex type modeling `ChainComplex`.
 * @tparam MapType The associative container in which the boundaries and
 * coboundaries are stored in both the computed matchings and the
 * `MorseComplex<CC, MapType>` instantiations.
 * @param complex A smart pointer to the input complex.
 * @return std::shared_ptr<MorseComplex<CC, MapType>> The fully-reduced Morse
 * complex of `complex`.
 */
template <ChainComplex CC,
    template <typename, typename> typename MapType = DefaultMap>
std::shared_ptr<MorseComplex<CC, MapType>> full_reduce(
    std::shared_ptr<CC> complex) {
  std::shared_ptr<MorseComplex<CC, MapType>> previous
      = std::make_shared<MorseComplex<CC, MapType>>(complex);
  std::shared_ptr<MorseComplex<CC, MapType>> next
      = std::make_shared<MorseComplex<CC, MapType>>(previous);
  while (next->size() != previous->size()) {
    previous = next;
    next = std::make_shared<MorseComplex<CC, MapType>>(previous);
  }
  return previous;
}

}  // namespace chomp::core

#endif  // CHOMP_COMPLEXES_MORSE_H
