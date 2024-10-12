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
#include <chomp/util/concepts.hpp>
#include <chomp/util/constants.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <compare>
#include <concepts>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <type_traits>

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
 * this is not the dimension of the cube as a cell of that complex. The
 * dimension of the cube itself is accessible via the `dimension` method.
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
  /** @brief Default initialize a `Cube` object. */
  Cube() = default;
  /**
   * @brief Initialize a `Cube` instance by providing an orthant and an extent
   * parameter.
   *
   * @param cube_orthant Location of the orthant this cube is in.
   * @param cube_extent Shape of this cube with each bit 1 or 0 depending on if
   * the cube has extent along the corresponding axis.
   */
  Cube(const CubeOrthant<CCDIM>& cube_orthant, std::size_t cube_extent) :
      cube_orthant(cube_orthant), cube_extent(cube_extent) {}
  /** @copydoc Cube(const CubeOrthant<CCDIM>&, std::size_t) */
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
   * @brief Get the dimension of the cube; this is the dimension of the chain
   * group in which it belongs, not the ambient dimension of the space.
   *
   * @return std::size_t
   */
  [[nodiscard]] std::size_t dimension() const noexcept {
    return std::popcount(cube_extent);
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
   * The ordering scheme is based first on the shape/extent parameter as an
   * integer then the orthant lexographically.
   */
  [[nodiscard]] std::strong_ordering operator<=>(const Cube& rhs
  ) const noexcept {
    if (cube_extent < rhs.cube_extent ||
        (cube_extent == rhs.cube_extent && cube_orthant < rhs.cube_orthant)) {
      return std::strong_ordering::less;
    }
    if (cube_extent > rhs.cube_extent ||
        (cube_extent == rhs.cube_extent && cube_orthant > rhs.cube_orthant)) {
      return std::strong_ordering::greater;
    }
    return std::strong_ordering::equivalent;
  }
};

}  // namespace chomp::core

namespace std {

/** @brief Hash specialization for `CubeOrthant` class template. */
template <size_t CCDIM>
struct hash<chomp::core::CubeOrthant<CCDIM>> {
  /** @brief Hash the orthant by prime combinations of its axis values. */
  size_t operator()(const chomp::core::CubeOrthant<CCDIM>& orthant) const {
    constexpr size_t PRIME = chomp::core::CUBE_HASH_PRIME;
    size_t hash_result = 0;
    for (size_t axis = 0; axis < CCDIM; ++axis) {
      hash_result = PRIME * hash_result + orthant[axis];
    }
    return hash_result;
  }
};

/** @brief Hash specialization for `Cube` class template. */
template <size_t CCDIM>
struct hash<chomp::core::Cube<CCDIM>> {
  /** @brief Hash the `Cube` by its orthant and its extent/shape parameter. */
  size_t operator()(const chomp::core::Cube<CCDIM>& cube) const {
    return hash<chomp::core::CubeOrthant<CCDIM>>{}(cube.orthant()) ^
           (cube.extent() << (chomp::core::SIZE_T_BITS - CCDIM));
  }
};

}  // namespace std

namespace chomp::core {

/**
 * @brief An iterator over all cubes in the hypercubical grid.
 *
 * This iterator functions by wrapping an iterator over the orthants, then
 * iterating over each cube in each orthant.
 *
 * This is used for cell iteration in `CubicalComplex` objects when wrapping a
 * `LexOrthantIterator`. However, note that the number of cells in a cubical
 * complex grows quickly in a high-dimensional grid, making individual iteration
 * over cells inefficient.
 *
 * @tparam CCDIM The dimension of the ambient hypercubical grid.
 * @tparam I The orthant iterator type; must be at least a forward iterator and
 * return orthants when dereferenced.
 */
template <std::size_t CCDIM, std::forward_iterator I>
requires requires(I it) {
  { *it } -> std::convertible_to<CubeOrthant<CCDIM>>;
}
class BasicCubeIterator {
private:
  I orthant_iterator;
  std::size_t current_extent;

public:
  /** @brief Difference type between iterators. */
  using difference_type = std::ptrdiff_t;
  /** @brief Value type when dereferenced. */
  using value_type = Cube<CCDIM>;
  /** @brief Pointer type. */
  using pointer = value_type*;
  /** @brief Reference type. */
  using reference = value_type&;
  /** @brief Tag for `iterator_traits`. */
  using iterator_concept = std::forward_iterator_tag;

  /**
   * @brief Default initialize a new BasicCubeIterator object.
   *
   * Unusable in this state but can be assigned to, as normal.
   */
  BasicCubeIterator() = default;
  /**
   * @brief Construct a new BasicCubeIterator by passing an orthant iterator.
   *
   * @param it A forward iterator over the orthants in the hypercubical grid.
   */
  explicit BasicCubeIterator(I it) : orthant_iterator(it), current_extent(0) {}

  /**
   * @brief Dereferencing this iterator yields a cube in the associated
   * `CubicalComplex`.
   *
   * @return Cube<CCDIM>
   */
  [[nodiscard]] Cube<CCDIM> operator*() const {
    return Cube<CCDIM>(*orthant_iterator, current_extent);
  }

  /**
   * @brief Equality operates first on the extent in this iterator and secondly
   * on the wrapped orthant iterator.
   *
   * @param rhs
   * @return true
   * @return false
   */
  [[nodiscard]] bool operator==(const BasicCubeIterator& rhs) const {
    if (current_extent != rhs.current_extent) {
      return false;
    }
    return orthant_iterator == rhs.orthant_iterator;
  }

  /**
   * @brief Preincrememnt increments extent parameter first, then increments
   * wrapped orthant iterator when extent parameter has reached its max value.
   *
   * @return BasicCubeIterator&
   */
  BasicCubeIterator& operator++() {
    ++current_extent;
    if (current_extent == 1 << CCDIM) {
      ++orthant_iterator;
      current_extent = 0;
    }
    return *this;
  }
  /**
   * @brief Postincrement operates as normal.
   */
  BasicCubeIterator operator++(int) {
    BasicCubeIterator temp(*this);
    ++*this;
    return temp;
  }
};

/**
 * @brief An iterator over the orthants of a hypercubical grid in a
 * lexicographical order.
 *
 * The grid is defined by a minimum orthant and a maximum orthant; starting at
 * any orthant between these two, the iterator moves along each axis from `0` to
 * `CCDIM - 1` and increments until it reaches the maximum orthant.
 *
 * Used for cell iteration in `CubicalComplex` objects when wrapped with a
 * `BasicCubeIterator` object.
 *
 * @tparam CCDIM The number of dimensions of the hypercubical grid.
 */
template <std::size_t CCDIM>
class LexOrthantIterator {
private:
  CubeOrthant<CCDIM> minimum_orthant;
  CubeOrthant<CCDIM> maximum_orthant;
  CubeOrthant<CCDIM> current_orthant;
  bool terminal;

public:
  /** @brief Difference type between iterators. */
  using difference_type = std::ptrdiff_t;
  /** @brief Value type when dereferenced. */
  using value_type = const CubeOrthant<CCDIM>;
  /** @brief Pointer type. */
  using pointer = value_type*;
  /** @brief Reference type. */
  using reference = value_type&;
  /** @brief Tag for `iterator_traits` */
  using iterator_concept = std::forward_iterator_tag;

  /**
   * @brief Default initialize a new LexOrthantIterator object.
   *
   * Unusable in this state but can be assigned to, as normal.
   */
  LexOrthantIterator() = default;
  /**
   * @brief Construct a new LexOrthantIterator object by passing minimal and
   * maximal orthants.
   *
   * The terminal flag can be passed as well; a value of `true` notes this as
   * a past-the-end iterator.
   *
   * @param it
   */
  explicit LexOrthantIterator(
      CubeOrthant<CCDIM> minimum_orthant, CubeOrthant<CCDIM> maximum_orthant,
      CubeOrthant<CCDIM> current_orthant, bool terminal = false
  ) :
      minimum_orthant(minimum_orthant), maximum_orthant(maximum_orthant),
      current_orthant(current_orthant), terminal(terminal) {}

  /**
   * @brief Dereferencing this iterator yields a constant reference to
   * the current orthant pointed to.
   *
   * @return reference
   */
  [[nodiscard]] reference operator*() const {
    return current_orthant;
  }

  /**
   * @brief Equality operates first on the `terminal` flag, then on the current
   * orhant.

   *
   * @param rhs
   * @return true
   * @return false
   */
  [[nodiscard]] bool operator==(const LexOrthantIterator& rhs) const {
    if (terminal != rhs.terminal) {
      return false;
    }
    if (terminal) {
      return true;
    }
    return current_orthant == rhs.current_orthant;
  }

  /**
   * @brief Preincrememnt moves the pointer to the next orthant in a
   * lexicographical order.
   *
   * If it would move past the final orthant, sets the `terminal` flag instead.
   *
   * @return LexOrthantIterator&
   */
  LexOrthantIterator& operator++() {
    for (std::size_t axis = 0; axis < CCDIM; ++axis) {
      if (current_orthant[axis] != maximum_orthant[axis]) {
        ++current_orthant[axis];
        return *this;
      }
      current_orthant[axis] = minimum_orthant[axis];
    }
    terminal = true;
    return *this;
  }
  /**
   * @brief Postincrement operates as normal.
   */
  LexOrthantIterator operator++(int) {
    LexOrthantIterator temp = *this;
    ++*this;
    return temp;
  }
};

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
  /** @brief Iterator type over cells. */
  using CellIterType = BasicCubeIterator<CCDIM, LexOrthantIterator<CCDIM>>;
  /** @brief Ambient dimension in which the complex is embedded. */
  static constexpr std::size_t dimension = CCDIM;

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
   * @brief Grade `input` according to the complex's grading function.
   *
   * @tparam T The type of the input. At minimum, this includes the complex's
   * cell type.
   * @param input
   * @return GradingResultType
   */
  template <typename T>
  requires requires(T t) {
    { grading_function(t) } -> std::convertible_to<GradingResultType>;
  }
  GradingResultType grade(const T& input) {
    return grading_function(input);
  }

  /**
   * @brief Beginning iterator for lexicographical iteration over cells. Allows
   * use of complex in range-based for loops.
   *
   * Due to the very high number of cells in a complex of higher ambient
   * dimension, it may be inadvisable to iterate over all cells separately using
   * this method.
   *
   * @return CellIterType
   */
  [[nodiscard]] CellIterType begin() const noexcept {
    return CellIterType(
        LexOrthantIterator(minimum_orthant, maximum_orthant, minimum_orthant)
    );
  }

  /**
   * @brief End sentinel for lexicographical iteration over cells.
   *
   * @return CellIterType
   */
  [[nodiscard]] CellIterType end() const noexcept {
    return CellIterType(LexOrthantIterator(
        minimum_orthant, maximum_orthant, minimum_orthant, true
    ));
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

/**
 * @brief A function object modeling `BoundedGrading` based on the cell being in
 * the closure of top dimensional (i.e. dimesnion `CCDIM`) cubes.
 *
 * @tparam CCDIM Dimension of the ambient hypercubical complex.
 * @tparam MIN Minimal value, i.e. the value returned if the queried cell is in
 * the closure of an included top-dimensional cell.
 * @tparam MAX Maximal value, i.e. the value returned if the queried cell is NOT
 * in the closure of an included top-dimensional cell.
 * @tparam SetType The underlying set data structure. Expected to be either
 * `std::set` or `std::unordered_set` but any set-like container with
 * sufficiently similar interface can work. The default type is
 * `std::unordered_set` as the stored `CubeOrthant<CCDIM>` type is hashable.
 */
template <
    std::size_t CCDIM, GradingResultType MIN = 0, GradingResultType MAX = 1,
    template <typename...> typename SetType = DefaultSet>
class TopCubeSetGrading {
public:
  /** @brief Input cell type to the grading function. */
  using InputType = Cube<CCDIM>;
  /** @brief Orthant type. */
  using OrthantType = CubeOrthant<CCDIM>;
  /** @brief Value type of the underlying orthant set. */
  using ValueType = typename SetType<OrthantType>::value_type;
  /** @brief Minimal grading value. */
  using Minimum = std::integral_constant<GradingResultType, MIN>;
  /** @brief Maximal grading value. */
  using Maximum = std::integral_constant<GradingResultType, MAX>;

private:
  SetType<OrthantType> top_cube_set;

public:
  /**
   * @brief Initialize a TopCubeSetGrading function object by supplying a set of
   * the (orthants of) top-dimensional cubes with minimal grade.
   *
   * @param top_cubes
   */
  TopCubeSetGrading(const SetType<OrthantType>& top_cubes) :
      top_cube_set(top_cubes) {}
  /** @overload */
  TopCubeSetGrading(SetType<OrthantType>&& top_cubes) :
      top_cube_set(top_cubes) {}
  /** @overload */
  TopCubeSetGrading(const std::initializer_list<ValueType>& top_cubes) :
      top_cube_set(top_cubes) {}
  /** @overload */
  TopCubeSetGrading(std::initializer_list<ValueType>&& top_cubes) :
      top_cube_set(top_cubes) {}

  /**
   * @brief Query the grade of `orthant`, i.e. the minimal value (`MIN`) if
   * `orthant` corresponds to a stored top-dimensional cube or the maximal value
   * (`MAX`) otherwise.
   *
   * @param orthant
   * @return GradingResultType
   */
  GradingResultType operator()(const OrthantType& orthant) const noexcept {
    return top_cube_set.contains(orthant) ? MIN : MAX;
  }

  /**
   * @brief Query the grade of a cell (cube) by determining if it is in the
   * closure of a top-dimensional cell of minimal grade.
   *
   * @param cube
   * @return GradingResultType
   */
  GradingResultType operator()(const InputType& cube) const noexcept {
    // Used for iterating in surrounding orthants.
    std::array<bool, CCDIM> minimal{};

    CubeOrthant<CCDIM> current_orthant = cube.orthant();
    while (true) {
      if (top_cube_set.contains(current_orthant)) {
        return MIN;
      }

      // Iterate to next orthant.
      bool broken = false;
      std::size_t axis = 0;
      std::size_t axis_flag = 1;
      for (; axis < CCDIM; ++axis, axis_flag <<= 1) {
        if (cube.extent() & axis_flag) {
          continue;
        }
        if (!minimal[axis]) {
          minimal[axis] = true;
          --(current_orthant[axis]);
          broken = true;
          break;
        }
        minimal[axis] = false;
        ++(current_orthant[axis]);
      }

      if (!broken) {
        return MAX;
      }
    }
  }
};

}  // namespace chomp::core

#endif  // CHOMP_COMPLEXES_CUBICAL_H
