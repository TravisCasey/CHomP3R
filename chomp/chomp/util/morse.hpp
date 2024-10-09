/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains the `MatchResult` class and the associated
 * `Trichotomy` class which encapsulate and conveniently expose data about
 * partial matchings used in Morse reduction routines.
 */

#ifndef CHOMP_UTIL_HASSE_H
#define CHOMP_UTIL_HASSE_H

#include <chomp/algebra/algebra.hpp>
#include <chomp/complexes/complexes.hpp>

#include <concepts>
#include <type_traits>

namespace chomp::core {

/**
 * @brief A lightweight class with three types of values corresponding to
 * queens, aces, and kings.
 *
 * These are intended to be exposed through the static inline variables `queen`,
 * `ace`, and `king`, but `ace` can also be constructed via the default
 * constructor.
 */
class Trichotomy {
private:
  signed char value;
  constexpr explicit Trichotomy(signed char value) noexcept : value(value) {}

public:
  /** @brief A value indicating the associated cell is a queen. */
  static const Trichotomy queen;
  /** @brief A value indicating the associated cell is an ace (critical). */
  static const Trichotomy ace;
  /** @brief A value indicating the associated cell is a king. */
  static const Trichotomy king;

  /** @brief Default constructor yields an equivalent object to `ace`. */
  constexpr Trichotomy() : value() {}

  /** @brief Each valid value compares equal only to themselves. */
  [[nodiscard]] constexpr bool operator==(Trichotomy rhs) const noexcept {
    return value == rhs.value;
  }
};

// Defining valid values.
inline constexpr Trichotomy Trichotomy::queen(-1);
inline constexpr Trichotomy Trichotomy::ace(0);
inline constexpr Trichotomy Trichotomy::king(1);

/**
 * @brief A lightweight class indicating the result of a partial matching
 * applied to a cell of some chain complex.
 *
 * A cell (of type `T`) is matched to another cell (type `T`), which may be
 * itself; this cell is accessible via the `cell` method.
 *
 * The type of match (ace, queen, king) is indicated by the `Trichotomy` object
 * which is accessed via the `tri` method; alternatively, the type can be
 * queried by the `is_queen`, `is_ace`, `is_king` helper methods.
 *
 * Finally, the coefficient of incidence (type `R`) of the cell on its match is
 * given by the `coef` method.
 *
 * @tparam T The cell type modeling `Cellular`.
 * @tparam R The coefficient ring type modeling `Ring`.
 */
template <Cellular T, Ring R>
class MatchResult {
private:
  T match_cell;
  R match_coef;
  Trichotomy match_tri;

public:
  /** @brief The cell type and the return type of the `cell` method. */
  using CellType = T;
  /** @brief The coefficient type and the return type of the `coef` method. */
  using RingType = R;

  /**
   * @brief Default initialize a `MatchResult` object by default initializing
   * the cell, coefficient, and trichotomy objects it wraps.
   */
  MatchResult() noexcept : match_cell(), match_coef(), match_tri() {}

  /**
   * @brief Initalize a `MatchResult` object for a given cell by passing the
   * cell it is matched to (`match_cell`), the coefficient of their incidence
   * (`match_coef`) and the `Trichotomy` object indicating the type of match.
   *
   * @tparam TFor Forwarding type for `CellType`.
   * @tparam RFor Forwarding type for `RingType`.
   * @param match_cell The matched cell; uses perfect forwarding.
   * @param match_coef The incidence coefficient of the cell being matched and
   * `match_cell`; uses perfect forwarding.
   * @param match_tri Indicator of whether the cell being matched (not
   * `match_cell`) is a queen, king, or ace.
   */
  template <typename TFor, typename RFor>
  requires std::same_as<std::remove_cvref_t<TFor>, CellType> &&
               std::same_as<std::remove_cvref_t<RFor>, RingType>
  MatchResult(
      TFor&& match_cell, RFor&& match_coef, Trichotomy match_tri
  ) noexcept :
      match_cell(std::forward<TFor>(match_cell)),
      match_coef(std::forward<RFor>(match_coef)), match_tri(match_tri) {}

  /** @brief Get the wrapped cell object. */
  [[nodiscard]] CellType cell() const noexcept {
    return match_cell;
  }
  /** @brief Get the wrapped coefficient ring object. */
  [[nodiscard]] RingType coef() const noexcept {
    return match_coef;
  }
  /** @brief Get the wrapped `Trichotomy` object. */
  [[nodiscard]] Trichotomy tri() const noexcept {
    return match_tri;
  }

  /** @brief Helper function to determine if the queried cell is a queen. */
  [[nodiscard]] bool is_queen() const noexcept {
    return match_tri == Trichotomy::queen;
  }
  /** @brief Helper function to determine if the queried cell is an ace. */
  [[nodiscard]] bool is_ace() const noexcept {
    return match_tri == Trichotomy::ace;
  }
  /** @brief Helper function to determine if the queried cell is a king. */
  [[nodiscard]] bool is_king() const noexcept {
    return match_tri == Trichotomy::king;
  }
};

}  // namespace chomp::core

#endif  // CHOMP_UTIL_HASSE_H
