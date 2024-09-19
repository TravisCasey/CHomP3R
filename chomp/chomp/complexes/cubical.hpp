/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains the implementations for cubical complexes,
 * including the `CubicalComplex` and `Cube` classes.
 */

#ifndef CHOMP_COMPLEXES_CUBICAL_H
#define CHOMP_COMPLEXES_CUBICAL_H

#include <chomp/algebra/algebra.hpp>
#include <chomp/algebra/cyclic.hpp>
#include <chomp/algebra/modules.hpp>
#include <chomp/complexes/complexes.hpp>
#include <chomp/complexes/grading.hpp>
#include <chomp/util/constants.hpp>

#include <array>
#include <compare>
#include <concepts>
#include <cstddef>
#include <functional>


namespace chomp::core {

/**
 * @brief Array type denoting the location of an orthant on the hypercubical
 * grid. Each entry is a coordinate along each axis forming an ordered tuple.
 *
 * @tparam CCDIM The ambient dimension of the hypercubical grid, i.e., the
 * number of axes and thus the number of entries in the array.
 */
template <std::size_t CCDIM>
using CubeOrthant = std::array<HypercubeCoordinate, CCDIM>;

/**
 * @brief A hypercube embedded in `CCDIM`-dimensional space. This is the cell
 * type for the `CubicalComplex` chain complex.
 *
 * Importantly, the `CCDIM` template parameter is the dimension of the ambient
 * space in which the associated cubical complex (and this cube) is embedded;
 * this is not the dimension of the cube as a cell of that complex.
 *
 * Each cube comprises an orthant (see `CubeOrthant`) and a shape integer that
 * has each bit 1 or 0 depending if the cube has extent along the corresponding
 * axis.
 *
 * @tparam CCDIM The dimension of the hypercubical grid in which the associated
 * cubical complex and this cube is embedded. Not the dimension of this cell as
 * a cell of said complex.
 *
 * @sa `CubicalComplex`, `CubeOrthant`
 */
template <std::size_t CCDIM>
class Cube {
private:
  CubeOrthant<CCDIM> cube_orthant;
  std::size_t cube_extent;

public:
  /**
   * @brief Initialize a `Cube` instance by providing an orthant and an extent
   * parameter.
   *
   * @param cube_orthant Location of the orthant this cube is in.
   * @param cube_extent Shape of this cube with each bit 1 or 0 depending on
   * if the cube has extent along the corresponding axis.
   */
  Cube(const CubeOrthant<CCDIM>& cube_orthant, std::size_t cube_extent) :
      cube_orthant(cube_orthant), cube_extent(cube_extent) {}

  /** @copydoc Cube() */
  Cube(CubeOrthant<CCDIM>&& cube_orthant, std::size_t cube_extent) :
      cube_orthant(std::move(cube_orthant)), cube_extent(cube_extent) {}

  /**
   * @brief Get the orthant of this cube.
   *
   * @return const CubeOrthant<CCDIM>&
   */
  [[nodiscard]] const CubeOrthant<CCDIM>& orthant() const noexcept {
    return cube_orthant;
  }

  /**
   * @brief Get the shape parameter of this cube.
   *
   * @return std::size_t
   */
  [[nodiscard]] std::size_t extent() const noexcept {
    return cube_extent;
  }

  /**
   * @brief Equality of `Cube` instances is based on equality of orthant and
   * extent/shape parameters.
   *
   * @param rhs
   * @return true
   * @return false
   */
  [[nodiscard]] bool operator==(const Cube& rhs) const noexcept {
    return cube_orthant == rhs.cube_orthant && cube_extent == rhs.cube_extent;
  }

  /**
   * @brief Three-way comparison operator synthesizes comparison operators
   * between `Cube` instances.
   *
   * Notably, this enables the use of `Cube` as a basis for module classes.
   *
   * The ordering scheme is based first on orthant lexographically then the
   * shape/extent parameter as an integer.
   */
  [[nodiscard]] std::strong_ordering operator<=>(const Cube& rhs
  ) const noexcept {
    if (cube_orthant < rhs.cube_orthant ||
        (cube_orthant == rhs.cube_orthant && cube_extent < rhs.cube_extent)) {
      return std::strong_ordering::less;
    }
    if (cube_orthant > rhs.cube_orthant ||
        (cube_orthant == rhs.cube_orthant && cube_extent > rhs.cube_extent)) {
      return std::strong_ordering::greater;
    }
    return std::strong_ordering::equivalent;
  }
};

}  // namespace chomp::core

namespace std {

/**
 * @brief Hash specialization for `Cube` class.
 *
 * @tparam CCDIM
 */
template <size_t CCDIM>
struct hash<chomp::core::Cube<CCDIM>> {
  /** @brief Hash the `Cube` by its orthant and its extent/shape parameter. */
  size_t operator()(const chomp::core::Cube<CCDIM>& cube) const {
    constexpr size_t PRIME = chomp::core::CUBE_HASH_PRIME;
    constexpr size_t BIT_DIFFERENCE = chomp::core::SIZE_T_BITS - CCDIM;
    size_t hash_result = 0;
    const chomp::core::CubeOrthant<CCDIM>& cube_orthant = cube.orthant();
    for (size_t idx = 0; idx < CCDIM; ++idx) {
      hash_result = PRIME * hash_result + cube_orthant[idx];
    }
    return hash_result ^ (cube.extent() << BIT_DIFFERENCE);
  }
};

}  // namespace std

namespace chomp::core {

/**
 * @brief Class implementing a cubical complex embedded in a `CCDIM`-dimensional
 * hypercubical grid.
 *
 * The corresponding cell class is `Cube`.
 *
 * The cubical complex includes all orthants in the hypercubical grid between
 * some minimum orthant (default the origin) and some user-provided maximum
 * orthant.
 *
 * @tparam CCDIM Ambient dimension of the cubical complex i.e. the maximum
 * dimension of `Cube` instances as cells and the number of axes. Notably, this
 * must be fewer than the bitwidth of `std::size_t`.
 * @tparam G The type of grading function object, which must model `Grading`.
 * @tparam R The coefficient ring type, which must model `Ring`. Default value
 * is `Z<2>`, i.e. the ring (field) with two elements.
 * @tparam M The chain type, which must model `Module`. The basis type must be
 * `Cube<CCDIM>` and the coefficient ring type must be `R`. The default type is
 * set to `DefaultModule` on these types.
 */
template <
    std::size_t CCDIM, Grading G, Ring R = Z<2>,
    Module M = DefaultModule<Cube<CCDIM>, R>>
requires requires {
  requires CCDIM <= SIZE_T_BITS;
  requires std::same_as<typename G::InputType, Cube<CCDIM>>;
  requires std::same_as<typename M::RingType, R>;
  requires std::same_as<typename M::BasisType, Cube<CCDIM>>;
}
class CubicalComplex {
private:
  CubeOrthant<CCDIM> minimum_orthant;
  CubeOrthant<CCDIM> maximum_orthant;
  G grading_function;

public:
  /** @brief Coefficient ring type for chains. */
  using RingType = R;
  /** @brief Cell type for cubical complexes. */
  using CellType = Cube<CCDIM>;
  /** @brief Chain (module) type. */
  using ChainType = M;
  /** @brief Grading function object type. */
  using GradingType = G;

  /**
   * @brief Initialize a new cubical complex with a maximum orthant and the
   * origin as minimum orthant. Also requires a grading function.
   *
   * @param maximum_orthant
   * @param grading_function
   */
  CubicalComplex(
      const CubeOrthant<CCDIM>& maximum_orthant,
      const GradingType& grading_function
  ) :
      minimum_orthant(), maximum_orthant(maximum_orthant),
      grading_function(grading_function) {}
  /** @copydoc CubicalComplex(const CubeOrthant<CCDIM>&, const GradingType&) */
  CubicalComplex(
      const CubeOrthant<CCDIM>& maximum_orthant, GradingType&& grading_function
  ) :
      minimum_orthant(), maximum_orthant(maximum_orthant),
      grading_function(std::move(grading_function)) {}
  /**
   * @brief Initialize a new cubical complex with both a maximum and minimum
   * orthant. Also requires a grading function.
   *
   * @param minimum_orthant
   * @param maximum_orthant
   * @param grading_function
   */
  CubicalComplex(
      const CubeOrthant<CCDIM>& minimum_orthant,
      const CubeOrthant<CCDIM>& maximum_orthant,
      const GradingType& grading_function
  ) :
      minimum_orthant(minimum_orthant), maximum_orthant(maximum_orthant),
      grading_function(grading_function) {}
  /**
   * @copydoc CubicalComplex(const CubeOrthant<CCDIM>&,
   * const CubeOrthant<CCDIM>&, const GradingType&)
   */
  CubicalComplex(
      const CubeOrthant<CCDIM>& minimum_orthant,
      const CubeOrthant<CCDIM>& maximum_orthant, GradingType&& grading_function
  ) :
      minimum_orthant(minimum_orthant), maximum_orthant(maximum_orthant),
      grading_function(std::move(grading_function)) {}

  /**
   * @brief Get the minimum orthant in the complex.
   *
   * @return const CubeOrthant<CCDIM>&
   */
  [[nodiscard]] const CubeOrthant<CCDIM>& minimum() const noexcept {
    return minimum_orthant;
  }
  /**
   * @brief Get the minimum coordinate in the complex along a particular axis.
   *
   * @param axis
   * @return HypercubeCoordinate
   */
  [[nodiscard]] HypercubeCoordinate minimum(std::size_t axis) const {
    return minimum_orthant.at(axis);
  }
  /**
   * @brief Get the maximum orthant in the complex.
   *
   * @return const CubeOrthant<CCDIM>&
   */
  [[nodiscard]] const CubeOrthant<CCDIM>& maximum() const noexcept {
    return maximum_orthant;
  }
  /**
   * @brief Get the maximum coordinate in the complex along a particular axis.
   *
   * @param axis
   * @return HypercubeCoordinate
   */
  [[nodiscard]] HypercubeCoordinate maximum(std::size_t axis) const {
    return maximum_orthant.at(axis);
  }

  /**
   * @brief Grade a given cell according to the complex's grading function.
   *
   * @param cell
   * @return GradingResultType
   */
  GradingResultType grade(const CellType& cell) {
    return grading_function(cell);
  }

  /**
   * @brief Get the boundary of `cell` in the complex subject to some constraint
   * `cond`.
   *
   * This method synthesizes the `boundary`, `graded_boundary`, and
   * `closure_boundary` function templates for this class along with the
   * generalizations to chain inputs.
   *
   * @param cell
   * @param cond A function taking (a constant reference to) a potential
   * boundary cell and returning a boolean value; if `true`, the cell is added
   * to the boundary.
   * @return ChainType
   */
  [[nodiscard]] ChainType
  boundary_if(const CellType& cell, const ConditionalType<CellType>& cond) {
    // Implementation follows `Computational Homology` Kaczynski et al.
    const CubeOrthant<CCDIM>& cube_orthant = cell.orthant();
    const std::size_t cube_extent = cell.extent();
    std::size_t axis_bit = 1;
    RingType coef = one<RingType>();  // axes with extent negate the coefficient
    ChainType result;

    for (std::size_t axis = 0; axis < CCDIM; ++axis, axis_bit <<= 1) {
      // cell must have extent along this axis to have a boundary
      if (cube_extent & axis_bit) {
        // The extent is the same for inner and outer cell
        std::size_t new_extent = cube_extent - axis_bit;

        // No outer cells along maximum edge of complex
        if (cube_orthant[axis] != maximum_orthant[axis]) {
          CubeOrthant<CCDIM> new_orthant = cube_orthant;
          new_orthant[axis] += 1;
          Cube<CCDIM> outer_cell(std::move(new_orthant), new_extent);
          if (cond(outer_cell)) {
            result.insert(outer_cell, coef);
          }
        }

        // Always inner cells
        Cube<CCDIM> inner_cell(cube_orthant, new_extent);
        if (cond(inner_cell)) {
          result.insert(inner_cell, -coef);
        }

        // Negate coefficient on axes with extent
        coef = -coef;
      }
    }
    return result;
  }

  /**
   * @brief Get the coboundary of `cell` in the complex subject to some
   * constraint `cond`.
   *
   * This method synthesizes the `coboundary`, `graded_coboundary`, and
   * `closure_coboundary` function templates for this class.
   *
   * @param cell
   * @param cond A function taking (a constant reference to) a potential
   * coboundary cell and returning a boolean value; if `true`, the cell is added
   * to the coboundary.
   * @return ChainType
   */
  [[nodiscard]] ChainType
  coboundary_if(const CellType& cell, const ConditionalType<CellType>& cond) {
    const CubeOrthant<CCDIM>& cube_orthant = cell.orthant();
    const std::size_t cube_extent = cell.extent();
    std::size_t axis_bit = 1;
    RingType coef = one<RingType>();
    ChainType result;

    for (std::size_t axis = 0; axis < CCDIM; ++axis, axis_bit <<= 1) {
      // cell must not have extent along this axis to have a boundary
      if (!(cube_extent & axis_bit)) {
        // The extent is the same for inner and outer cell
        std::size_t new_extent = cube_extent + axis_bit;

        // No inner cells along minimum edge of complex
        if (cube_orthant[axis] != minimum_orthant[axis]) {
          CubeOrthant<CCDIM> new_orthant = cube_orthant;
          new_orthant[axis] -= 1;
          Cube<CCDIM> inner_cell(std::move(new_orthant), new_extent);
          if (cond(inner_cell)) {
            result.insert(inner_cell, coef);
          }
        }

        // Always outer cells
        Cube<CCDIM> outer_cell(cube_orthant, new_extent);
        if (cond(outer_cell)) {
          result.insert(outer_cell, -coef);
        }

      } else {
        // Negate coefficient on axes with extent
        coef = -coef;
      }
    }
    return result;
  }
};

}  // namespace chomp::core

#endif  // CHOMP_COMPLEXES_CUBICAL_H
