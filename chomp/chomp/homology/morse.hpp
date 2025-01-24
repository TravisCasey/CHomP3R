/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains the partial matching (Morse function surrogate)
 * interface as well as its derived classes and helper objects.
 */

#ifndef CHOMP_HOMOLOGY_MORSE_H
#define CHOMP_HOMOLOGY_MORSE_H

#include <chomp/complexes/complexes.hpp>
#include <chomp/homology/hasse.hpp>
#include <chomp/util/concepts.hpp>
#include <chomp/util/morse.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <memory>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace chomp::core {

/**
 * @brief An abstract base class representing the mathematical partial matching,
 * which is a more computable surrogate for the discrete Morse function used to
 * reduce chain complexes.
 *
 * A partial matching uses the trichotomy induced by the Morse function to match
 * the "queen" cells (those cell with exactly one coface of lesser Morse
 * function value) to "king" cells (those with exactly one face of greater Morse
 * function value), as well as "ace" cells (neither kings nor queens, also
 * referred to as critical cells) to themselves.
 *
 * The partial matching is used to define a Morse chain complex on the critical
 * cell (aces); the `MorseComplex` class holds a pointer to `PartialMatching`,
 * which it used to determine which cells are critical as well as the boundary
 * and coboundary homomorphisms.
 *
 * For an involved treatment of this approach see: `Discrete Morse Theoretic
 * Algorithms for Computing Homology of Complexes and Maps` - Harker,
 * Mischaikow, Mrozek, Nanda.
 *
 * @tparam CC The chain complex type (models `ChainComplex`) that the partial
 * matching is computed on.
 * @tparam MapType Associative container used to store the computed boundary and
 * coboundary homomorphisms for the induced Morse complex. The cell type of `CC`
 * is the key type while the value type is the chain type of `CC`. Default value
 * is `DefaultMap`.
 */
template <ChainComplex CC, template <typename...> typename MapType = DefaultMap>
class PartialMatching {
public:
  /** @brief `ChainComplex` type this class acts upon. */
  using ComplexType = CC;
  /** @brief Coefficient ring type propagated from `CC`. */
  using RingType = typename CC::RingType;
  /** @brief Cell type propagated from `CC`. */
  using CellType = typename CC::CellType;
  /** @brief Chain type propagated from `CC`. */
  using ChainType = typename CC::ChainType;
  /** @brief Type detailing the results of a match. */
  using ResultType = MatchResult<CellType, RingType>;
  /** @brief Constant iterator type over critical cells. */
  using CellIterType = typename std::vector<CellType>::const_iterator;

protected:
  /** @brief Pointer to the `ChainComplex` this object acts upon. */
  std::shared_ptr<ComplexType> upper_complex_ptr;
  /** @brief Critical (ace) cells, i.e. the cells of the lowered complex. */
  std::vector<CellType> critical_cells;

  /** @brief Initializes base class data members. */
  PartialMatching(std::shared_ptr<CC> upper_complex_ptr) :
      upper_complex_ptr(upper_complex_ptr) {}

  /*
  Strict weak ordering on queen cells; let K_1 be the match of queen Q_1; Then,
  this method must satisfy:
    (queen_ordering(Q_1, Q_2) == true) => (Q_2 not in boundary(K_1))

  Must also satisfy `Compare` requirements for `std::priority_queue`.
  */

  /**
   * @brief Strict weak ordering on queen cells for gradient path computation.
   *
   * let K_1 be the match of queen Q_1; Then, this method must satisfy:
   * (queen_ordering(Q_1, Q_2) == true) => (Q_2 not in boundary(K_1))
   */
  [[nodiscard]] virtual bool queen_ordering(const CellType& lhs,
      const CellType& rhs) const
      = 0;

  /**
   * @brief Strict weak ordering on king cells for gradient path computation.
   *
   * let Q_1 be the match of queen K_1; Then, this method must satisfy:
   * (king_ordering(K_1, K_2) == true) => (K_2 not in coboundary(Q_1))
   */
  [[nodiscard]] virtual bool king_ordering(const CellType& lhs,
      const CellType& rhs) const
      = 0;

private:
  using MatchedPairType = std::pair<CellType, ResultType>;
  using CompareType
      = std::function<bool(const MatchedPairType&, const MatchedPairType&)>;

public:
  /** @brief Deleted; handle derived classes through pointers. */
  PartialMatching(const PartialMatching&) = delete;
  /** @brief Deleted; handle derived classes through pointers. */
  PartialMatching(PartialMatching&&) = delete;
  /** @brief Deleted; handle derived classes through pointers. */
  PartialMatching& operator=(const PartialMatching&) = delete;
  /** @brief Deleted; handle derived classes through pointers. */
  PartialMatching& operator=(PartialMatching&&) = delete;
  /** @brief Virtual default destructor for smart pointer compatability. */
  virtual ~PartialMatching() = default;

  /**
   * @brief Return a (smart) pointer to the upper chain complex on which this
   * `PartialMatching` acts.
   *
   * @return std::shared_ptr<ComplexType>
   */
  [[nodiscard]] std::shared_ptr<ComplexType> upper() const noexcept {
    return upper_complex_ptr;
  }

  /**
   * @brief Query the partial matching on `cell`. The data of this match is
   * returned in a `MatchResult` class template instance.
   *
   * Implementation depends upon the derived class.
   *
   * @param cell The cell of the upper chain complex to match.
   * @return ResultType The data of the match; specifically, the other cell that
   * `cell` is matched to, their incidence on each other, and whether `cell` is
   * a queen, ace, or king.
   */
  [[nodiscard]] virtual ResultType match(const CellType& cell) = 0;

  /**
   * @brief Constant iterator to the beginning of the critical cells (aces).
   *
   * Enables iterating over this matching (specifically, the critical cells) in
   * range-based for loops.
   *
   * @return CellIterType
   */
  [[nodiscard]] CellIterType begin() const noexcept {
    return critical_cells.cbegin();
  }
  /**
   * @brief Constant iterator to the end of the critical cells (aces).
   *
   * Enables iterating over this matching (specifically, the critical cells) in
   * range-based for loops.
   *
   * @return CellIterType
   */
  [[nodiscard]] CellIterType end() const noexcept {
    return critical_cells.cend();
  }

  /**
   * @brief The number of critical cells in the partial matching.
   *
   * Iterated partial matching terminates when the number of critical cells
   * no longer changes.
   *
   * @return std::size_t The number of critical cells in this matching.
   */
  [[nodiscard]] std::size_t size() const noexcept {
    return critical_cells.size();
  }

  /**
   * @brief Compute the boundary and coboundary homomorphisms in the Morse
   * complex associated to this partial matching.
   *
   * Used to initialize `MorseComplex` instances.
   *
   * @return std::pair<MapType<CellType, ChainType>, MapType<CellType,
   * ChainType>> A pair (boundary, coboundary) homomorphisms represented using
   * a map data structure from the basis to their result in the (co)chain
   * modules.
   */
  [[nodiscard]] virtual std::pair<MapType<CellType, ChainType>,
      MapType<CellType, ChainType>>
  compute_operators() {
    MapType<CellType, ChainType> boundaries;
    MapType<CellType, ChainType> coboundaries;

    // Initialize all coboundary chains so they can be inserted into later.
    for (const CellType& cell : critical_cells) {
      coboundaries.insert(std::make_pair(cell, ChainType()));
    }

    // Compute the lowered boundary of each cell.
    // Instead of colowering the coboundary, we use the fact that the coboundary
    // operator is dual to the boundary operator.
    for (const CellType& cell : critical_cells) {
      ChainType temp_boundary = lower(boundary(*upper_complex_ptr, cell));
      for (const CellType& boundary_cell : temp_boundary) {
        coboundaries[boundary_cell].insert(cell, temp_boundary[boundary_cell]);
      }
      boundaries.insert(std::make_pair(cell, std::move(temp_boundary)));
    }

    return std::make_pair(std::move(boundaries), std::move(coboundaries));
  }

  /**
   * @brief Lower a chain in the upper chain complex to its representative in
   * the reduced Morse complex.
   *
   * @param chain Chain in the upper complex.
   * @return ChainType Representative lowered chain in the Morse complex.
   */
  [[nodiscard]] ChainType lower(const ChainType& chain) {
    // Remaining queens and their coefficients to be eliminated.
    ChainType queen_chain;
    // Aces in this chain are the resulting lowered chain.
    ChainType ace_chain;
    // Scaled boundary of king cells to queens in queen_chain.
    ChainType boundary_chain;
    // Coefficient of queen currently being processed.
    RingType incidence;
    // Pair of cell being processed and its matching data.
    MatchedPairType match_pair;
    // Scaling coefficient of boundary chain
    RingType coef;

    // Yields maximal queens to process in the correct order.
    std::priority_queue<MatchedPairType, std::vector<MatchedPairType>,
        CompareType>
        queen_queue(
            [this](const MatchedPairType& lhs, const MatchedPairType& rhs) {
              return this->queen_ordering(lhs.first, rhs.first);
            });

    // Initialize queen_queue, queen_chain, ace_chain with elements of `chain`
    for (const CellType& cell : chain) {
      match_pair = std::make_pair(cell, match(cell));
      if (match_pair.second.is_queen()) {
        queen_queue.push(std::move(match_pair));
        queen_chain.insert(cell, chain[cell]);
      } else if (match_pair.second.is_ace()) {
        ace_chain.insert(cell, chain[cell]);
      }
    }

    // Loop until all queens are eliminated
    while (!queen_queue.empty()) {
      // Get the first (in ordering on queens) with nonzero coefficient.
      incidence = queen_chain[queen_queue.top().first];
      while (incidence == zero<RingType>()) {
        queen_queue.pop();
        if (queen_queue.empty()) {
          return ace_chain;
        }
        incidence = queen_chain[queen_queue.top().first];
      }

      // Find match of first queen and calculate its scaled boundary.
      match_pair = queen_queue.top();
      queen_queue.pop();
      coef = -incidence * invert(match_pair.second.coef());
      boundary_chain
          = coef * boundary(*upper_complex_ptr, match_pair.second.cell());

      // Add queens in boundary to queen_queue/queen_chain & aces to ace_chain.
      for (const CellType& cell : boundary_chain) {
        match_pair = std::make_pair(cell, match(cell));
        if (match_pair.second.is_queen()) {
          queen_queue.push(std::move(match_pair));
          queen_chain.insert(cell, boundary_chain[cell]);
        } else if (match_pair.second.is_ace()) {
          ace_chain.insert(cell, boundary_chain[cell]);
        }
      }
    }

    return ace_chain;
  }
  /** @overload */
  [[nodiscard]] ChainType lower(const CellType& cell) {
    ChainType cell_chain;
    cell_chain.insert(cell, one<RingType>());
    return lower(cell_chain);
  }

  /**
   * @brief Lift a chain from the Morse complex to its representative in the
   * upper chain complex.
   *
   * @param chain Chain in the Morse complex.
   * @return ChainType Representative lifted chain in the upper complex.
   */
  [[nodiscard]] ChainType lift(const ChainType& chain) {
    // Remaining queens and their coefficients to be eliminated.
    ChainType queen_chain;
    // Kings in this chain plus `chain` itself is the resuling lifted chain.
    ChainType king_chain;
    // Scaled boundary of king cells to queens in queen_chain.
    ChainType boundary_chain(boundary(*upper_complex_ptr, chain));
    // Coefficient of queen currently being processed.
    RingType incidence;
    // Pair of cell being processed and its matching data.
    MatchedPairType match_pair;
    // Scaling coefficient of boundary chain.
    RingType coef;

    // Yields maximal queens to process in the correct order.
    std::priority_queue<MatchedPairType, std::vector<MatchedPairType>,
        CompareType>
        queen_queue(
            [this](const MatchedPairType& lhs, const MatchedPairType& rhs) {
              return this->queen_ordering(lhs.first, rhs.first);
            });

    // Initialize queen_queue and queen_chain with boundary of `chain`.
    for (const CellType& cell : boundary_chain) {
      match_pair = std::make_pair(cell, match(cell));
      if (match_pair.second.is_queen()) {
        queen_queue.push(std::move(match_pair));
        queen_chain.insert(cell, boundary_chain[cell]);
      }
    }

    // Loop until all queens are eliminated
    while (!queen_queue.empty()) {
      // Get the first (in ordering on queens) with nonzero coefficient.
      incidence = queen_chain[queen_queue.top().first];
      while (incidence == zero<RingType>()) {
        queen_queue.pop();
        if (queen_queue.empty()) {
          return chain + king_chain;
        }
        incidence = queen_chain[queen_queue.top().first];
      }

      // Find match of first queen, add it to king_chain, and calculate its
      // scaled boundary.
      match_pair = queen_queue.top();
      queen_queue.pop();
      coef = -incidence * invert(match_pair.second.coef());
      boundary_chain
          = coef * boundary(*upper_complex_ptr, match_pair.second.cell());
      king_chain.insert(std::move(match_pair.second.cell()), std::move(coef));

      // Add queens in boundary to queen_queue/queen_chain.
      for (const CellType& cell : boundary_chain) {
        match_pair = std::make_pair(cell, match(cell));
        if (match_pair.second.is_queen()) {
          queen_queue.push(std::move(match_pair));
          queen_chain.insert(cell, boundary_chain[cell]);
        }
      }
    }

    return chain + king_chain;
  }
  /** @overload */
  [[nodiscard]] ChainType lift(const CellType& cell) {
    ChainType cell_chain;
    cell_chain.insert(cell, one<RingType>());
    return lift(cell_chain);
  }

  /**
   * @brief Lower a cochain in the upper chain complex to its representative in
   * the reduced Morse complex.
   *
   * @param chain Cochain in the upper complex.
   * @return ChainType Representative lowered cochain in the Morse complex.
   */
  [[nodiscard]] ChainType colower(const ChainType& chain) {
    // Remaining kings and their coefficients to be eliminated.
    ChainType king_chain;
    // Aces in this chain are the resulting (co)lowered chain.
    ChainType ace_chain;
    // Scaled coboundary of queen cells to kings in king_chain.
    ChainType coboundary_chain;
    // Coefficient of king currently being processed.
    RingType incidence;
    // Pair of cell being processed and its matching data.
    MatchedPairType match_pair;
    // Scaling coefficient of coboundary chain.
    RingType coef;

    // Yields maximal kings to process in the correct order.
    std::priority_queue<MatchedPairType, std::vector<MatchedPairType>,
        CompareType>
        king_queue(
            [this](const MatchedPairType& lhs, const MatchedPairType& rhs) {
              return this->king_ordering(lhs.first, rhs.first);
            });

    // Initialize king_queue, king_chain, ace_chain with elements of `chain`
    for (const CellType& cell : chain) {
      match_pair = std::make_pair(cell, match(cell));
      if (match_pair.second.is_king()) {
        king_queue.push(std::move(match_pair));
        king_chain.insert(cell, chain[cell]);
      } else if (match_pair.second.is_ace()) {
        ace_chain.insert(cell, chain[cell]);
      }
    }

    // Loop until all kings are eliminated
    while (!king_queue.empty()) {
      // Get the first (in ordering on kings) with nonzero coefficient.
      incidence = king_chain[king_queue.top().first];
      while (incidence == zero<RingType>()) {
        king_queue.pop();
        if (king_queue.empty()) {
          return ace_chain;
        }
        incidence = king_chain[king_queue.top().first];
      }

      // Find match of first king and calculate its scaled coboundary.
      match_pair = king_queue.top();
      king_queue.pop();
      coef = -incidence * invert(match_pair.second.coef());
      coboundary_chain
          = coef * coboundary(*upper_complex_ptr, match_pair.second.cell());

      // Add kings in coboundary to king_queue/king_chain & aces to ace_chain.
      for (const CellType& cell : coboundary_chain) {
        match_pair = std::make_pair(cell, match(cell));
        if (match_pair.second.is_king()) {
          king_queue.push(std::move(match_pair));
          king_chain.insert(cell, coboundary_chain[cell]);
        } else if (match_pair.second.is_ace()) {
          ace_chain.insert(cell, coboundary_chain[cell]);
        }
      }
    }

    return ace_chain;
  }
  /** @overload */
  [[nodiscard]] ChainType colower(const CellType& cell) {
    ChainType cell_chain;
    cell_chain.insert(cell, one<RingType>());
    return colower(cell_chain);
  }

  /**
   * @brief Lift a cochain from the Morse complex to its representative in the
   * upper chain complex.
   *
   * @param chain Cochain in the Morse complex.
   * @return ChainType Representative lifted cochain in the upper complex.
   */
  [[nodiscard]] ChainType colift(const ChainType& chain) {
    // Remaining kings and their coefficients to be eliminated.
    ChainType king_chain;
    // Queens in this chain plus `chain` itself is the resuling lifted chain.
    ChainType queen_chain;
    // Scaled coboundary of queen cells to kings in king_chain.
    ChainType coboundary_chain(coboundary(*upper_complex_ptr, chain));
    // Coefficient of king currently being processed.
    RingType incidence;
    // Pair of cell being processed and its matching data.
    MatchedPairType match_pair;
    // Scaling coefficient of coboundary chain.
    RingType coef;

    // Yields maximal kings to process in the correct order.
    std::priority_queue<MatchedPairType, std::vector<MatchedPairType>,
        CompareType>
        king_queue(
            [this](const MatchedPairType& lhs, const MatchedPairType& rhs) {
              return this->king_ordering(lhs.first, rhs.first);
            });

    // Initialize king_queue and king_chain with coboundary of `chain`.
    for (const CellType& cell : coboundary_chain) {
      match_pair = std::make_pair(cell, match(cell));
      if (match_pair.second.is_king()) {
        king_queue.push(std::move(match_pair));
        king_chain.insert(cell, coboundary_chain[cell]);
      }
    }

    // Loop until all kings are eliminated
    while (!king_queue.empty()) {
      // Get the first (in ordering on kings) with nonzero coefficient.
      incidence = king_chain[king_queue.top().first];
      while (incidence == zero<RingType>()) {
        king_queue.pop();
        if (king_queue.empty()) {
          return chain + queen_chain;
        }
        incidence = king_chain[king_queue.top().first];
      }

      // Find match of first king, add it to queen_chain, and calculate its
      // scaled coboundary.
      match_pair = king_queue.top();
      king_queue.pop();
      coef = -incidence * invert(match_pair.second.coef());
      coboundary_chain
          = coef * coboundary(*upper_complex_ptr, match_pair.second.cell());
      queen_chain.insert(std::move(match_pair.second.cell()), coef);

      // Add kings in coboundary to king_queue/king_chain.
      for (const CellType& cell : coboundary_chain) {
        match_pair = std::make_pair(cell, match(cell));
        if (match_pair.second.is_king()) {
          king_queue.push(std::move(match_pair));
          king_chain.insert(cell, coboundary_chain[cell]);
        }
      }
    }

    return chain + queen_chain;
  }
  /** @overload */
  [[nodiscard]] ChainType colift(const CellType& cell) {
    ChainType cell_chain;
    cell_chain.insert(cell, one<RingType>());
    return colift(cell_chain);
  }
};

/**
 * @brief A `PartialMatching`-derived implementation based on excising
 * coreductions in the Hasse diagram according to the face partial order.
 *
 * Applicable to any complex but may be less efficient than specialized
 * implementations. All (non-critical) matches are explicitly stored after
 * computation as this computation is difficult to query; thus it is best run on
 * smaller complexes.
 *
 * This uses a reformulation of Algorithm 3.6 in `Discrete Morse Theoretic
 * Algorithms for Computing Homology of Complexes and Maps` - Harker,
 * Mischaikow, Mrozek, Nanda.
 *
 * @tparam CC The chain complex type (models `ChainComplex`) that the partial
 * matching is computed on.
 * @tparam MapType Associative container used to store the computed boundary and
 * coboundary homomorphisms for the induced Morse complex. The cell type of `CC`
 * is the key type while the value type is the chain type of `CC`. Default value
 * is `DefaultMap`.
 */
template <ChainComplex CC, template <typename...> typename MapType = DefaultMap>
class CoreductionMatching : public PartialMatching<CC, MapType> {
public:
  /** @brief `ChainComplex` type this class acts upon. */
  using ComplexType = CC;
  /** @brief Coefficient ring type propagated from `CC`. */
  using RingType = typename CC::RingType;
  /** @brief Cell type propagated from `CC`. */
  using CellType = typename CC::CellType;
  /** @brief Type detailing the results of a match. */
  using ResultType = MatchResult<CellType, RingType>;

private:
  MapType<CellType, ResultType> matches;
  MapType<CellType, std::size_t> priority;

  using MatchIterType = decltype(matches)::const_iterator;
  using PriorityIterType = decltype(priority)::const_iterator;

  // Based on the order of matching; earlier matched (lower priority) queens are
  // maximal w.r.t. queens matched later.
  [[nodiscard]] bool queen_ordering(const CellType& lhs,
      const CellType& rhs) const override {
    if (this->upper_complex_ptr->grade(lhs)
        < this->upper_complex_ptr->grade(rhs)) {
      return true;
    }

    PriorityIterType lhs_it = priority.find(lhs);
    PriorityIterType rhs_it = priority.find(rhs);
    if (lhs_it == priority.cend() || rhs_it == priority.cend()) {
      throw std::invalid_argument(
          "Partial matching ordering function called on invalid cell; critical "
          "cells cannot be ordered.");
    }
    return lhs_it->second < rhs_it->second;
  }

  // The inverse of the queen ordering
  [[nodiscard]] bool king_ordering(const CellType& lhs,
      const CellType& rhs) const override {
    if (this->upper_complex_ptr->grade(lhs)
        > this->upper_complex_ptr->grade(rhs)) {
      return true;
    }

    PriorityIterType lhs_it = priority.find(lhs);
    PriorityIterType rhs_it = priority.find(rhs);
    if (lhs_it == priority.cend() || rhs_it == priority.cend()) {
      throw std::invalid_argument(
          "Partial matching ordering function called on invalid cell; critical "
          "cells cannot be ordered.");
    }
    return lhs_it->second > rhs_it->second;
  }

public:
  /**
   * @brief Initializes the `CoreductionMatching` object by computing the
   * partial matching on the chain complex pointed to by `complex`.
   *
   * @param complex The upper chain complex to be matched.
   */
  CoreductionMatching(std::shared_ptr<CC> complex) :
      PartialMatching<CC, MapType>(complex) {
    std::tie(matches, priority, this->critical_cells)
        = detail::HasseCoreduction<CC, MapType>::compute_matching(
            this->upper_complex_ptr);
  }

  /**
   * @brief Query the partial matching on `cell`. The data of this match is
   * returned in a `MatchResult` class template instance.
   *
   * All (non-critical) matches are stored explicitly after computation of the
   * matching. Hence, querying further matches does not require extensive
   * computation.
   *
   * @param cell The cell of the upper chain complex to match.
   * @return ResultType The data of the match; specifically, the other cell that
   * `cell` is matched to, their incidence on each other, and whether `cell` is
   * a queen, ace, or king.
   */
  [[nodiscard]] ResultType match(const CellType& cell) override {
    MatchIterType match_it = matches.find(cell);
    if (match_it != matches.cend()) {
      return match_it->second;
    }
    return ResultType(cell, one<RingType>(), Trichotomy::ace);
  }
};

}  // namespace chomp::core

#endif  // CHOMP_HOMOLOGY_MORSE_H
