/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains the implementations for cubical complexes,
 * including the `CubicalComplex` and `Cube` classes.
 */

#ifndef CHOMP_CORE_COMPLEXES_CUBICAL_HPP
#define CHOMP_CORE_COMPLEXES_CUBICAL_HPP

#include <chomp/core/algebra/algebra.hpp>
#include <chomp/core/algebra/cyclic.hpp>
#include <chomp/core/algebra/modules.hpp>
#include <chomp/core/complexes/complexes.hpp>
#include <chomp/core/complexes/grading.hpp>
#include <chomp/core/util/concepts.hpp>
#include <chomp/core/util/constants.hpp>

#include <algorithm>
#include <array>
#include <compare>
#include <concepts>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>

namespace chomp::core {

/**
 * @brief Represents an ordered tuple of points in a `CCDIM`-dimensional
 * hypercubical grid.
 *
 * A (closed, unit) orthant is the union of all cubical cells in the grid whose
 * closures share a vertex (in the grid) as their least point. This effectively
 * partitions the hypercubical grid as a cubical complex, with each equivalence
 * class corresponding to a vertex; this representation is encapsulated by the
 * `Orthant` class.
 *
 * The interface of `Orthant` follows that of `std::array`, with heap-allocated
 * members and improved move semantics.
 *
 * @tparam CCDIM The ambient dimension of the hypercubical grid, i.e., the
 * number of axes and thus the number of entries in the ordered tuple.
 *
 * @sa Cube
 */
template <std::size_t CCDIM>
class Orthant {
public:
  /** @brief Underlying array type wrapped by this class template. */
  using ArrayType = std::array<HypercubeCoordinate, CCDIM>;
  /**
   * @brief Underlying pointer type to `ArrayType` returned by reference from
   * the `data` method.
   */
  using PointerType = std::unique_ptr<ArrayType>;
  /** @brief Iterator type propagated from `ArrayType`. */
  using IterType = typename ArrayType::iterator;
  /** @brief `const` iterator type propagated from `ArrayType`. */
  using CIterType = typename ArrayType::const_iterator;
  /** @brief Reverse iterator type propagated from `ArrayType`. */
  using RIterType = typename ArrayType::reverse_iterator;
  /** @brief `const` reverse iterator type propagated from `ArrayType`. */
  using CRIterType = typename ArrayType::const_reverse_iterator;

private:
  PointerType coordinate_ptr;

public:
  /**
   * @brief Default constructor creates a new `Orthant` object with all
   * coordinates zero (i.e. the orthant at the origin).
   */
  Orthant() :
      coordinate_ptr(std::make_unique<ArrayType>()) {}
  /**
   * @brief Copy constructor heap allocates a new container with coordinates
   * identical to `other` for the newly-constructed object.
   */
  Orthant(const Orthant& other) :
      coordinate_ptr(std::make_unique<ArrayType>(*(other.coordinate_ptr))) {}
  /**
   * @brief Move constructor transfers ownership of the container to the
   * newly-constructed object without allocating new memory on the heap.
   *
   * The moved-from `Orthant` object is left in an invalid state, but may be
   * assigned-to as usual.
   */
  Orthant(Orthant&&) = default;
  /**
   * @brief Construct an `Orthant` object by populating the container with
   * the elements of `coordinate_ilist`.
   *
   * Unpopulated entries are left as zero, while entries beyond the first
   * `CCDIM` are ignored.
   */
  constexpr Orthant(
      std::initializer_list<HypercubeCoordinate> coordinate_ilist) :
      coordinate_ptr(std::make_unique<ArrayType>()) {
    std::size_t axis = 0;
    for (const HypercubeCoordinate& coord : coordinate_ilist) {
      (*coordinate_ptr)[axis++] = coord;
      if (axis == CCDIM) {
        break;
      }
    }
  }

  /**
   * @brief Copy assignment overwrites each coordinate in the container with a
   * copy of the corresponding coordinate in `other`.
   *
   * If the object is in an invalid state (namely, after being moved from), a
   * new container is allocated and copy-constructed from that of `other`.
   *
   * @return Orthant& Reference to this `Orthant` object.
   */
  Orthant& operator=(const Orthant& other) {
    if (coordinate_ptr) {
      *coordinate_ptr = *(other.coordinate_ptr);
    } else {
      coordinate_ptr = std::make_unique<ArrayType>(*(other.coordinate_ptr));
    }
    return *this;
  }
  /**
   * @brief The move assignment operator transfers ownership of the container
   * without allocating new memory om the heap.
   *
   * @return Orthant& Reference to this `Orthant` object.
   */
  Orthant& operator=(Orthant&&) = default;

  /** @brief Destruct this `Orthant` object and its heap-allocated container. */
  ~Orthant() = default;

  /**
   * @brief Get the size of the container, which is always equal to the `CCDIM`
   * template parameter.
   */
  [[nodiscard]] constexpr std::size_t size() const noexcept {
    return CCDIM;
  }

  /**
   * @brief Return by reference a pointer to the array underyling this object.
   *
   * @return PointerType& A reference to `PointerType` pointing to the
   * underlying array of coordinates.
   */
  [[nodiscard]] PointerType& data() noexcept {
    return coordinate_ptr;
  }
  /** @overload */
  [[nodiscard]] const PointerType& data() const noexcept {
    return coordinate_ptr;
  }

  /** @brief Access the coordinate at `axis` in the container. */
  [[nodiscard]] HypercubeCoordinate& operator[](const std::size_t axis) {
    return (*coordinate_ptr)[axis];
  }
  /** @overload */
  [[nodiscard]] const HypercubeCoordinate& operator[](
      const std::size_t axis) const {
    return (*coordinate_ptr)[axis];
  }

  /** @brief Return an iterator to the beginning of the container. */
  [[nodiscard]] IterType begin() noexcept {
    return coordinate_ptr->begin();
  }
  /** @overload */
  [[nodiscard]] CIterType begin() const noexcept {
    return coordinate_ptr->begin();
  }
  /** @brief Return a `const` iterator to the beginning of the container. */
  [[nodiscard]] CIterType cbegin() const noexcept {
    return coordinate_ptr->cbegin();
  }
  /** @brief Return a reverse iterator to the end of the container. */
  [[nodiscard]] RIterType rbegin() noexcept {
    return coordinate_ptr->rbegin();
  }
  /** @overload */
  [[nodiscard]] CRIterType rbegin() const noexcept {
    return coordinate_ptr->rbegin();
  }
  /** @brief Return a reverse `const` iterator to the end of the container. */
  [[nodiscard]] CRIterType crbegin() const noexcept {
    return coordinate_ptr->crbegin();
  }
  /** @brief Return an iterator to the end of the container. */
  [[nodiscard]] IterType end() noexcept {
    return coordinate_ptr->end();
  }
  /** @overload */
  [[nodiscard]] CIterType end() const noexcept {
    return coordinate_ptr->end();
  }
  /** @brief Return a `const` iterator to the end of the container. */
  [[nodiscard]] CIterType cend() const noexcept {
    return coordinate_ptr->cend();
  }
  /** @brief Return a reverse iterator to the beginning of the container. */
  [[nodiscard]] RIterType rend() noexcept {
    return coordinate_ptr->rend();
  }
  /** @overload */
  [[nodiscard]] CRIterType rend() const noexcept {
    return coordinate_ptr->rend();
  }
  /**
   * @brief Return a `const` reverse iterator to the beginning of the container.
   */
  [[nodiscard]] CRIterType crend() const noexcept {
    return coordinate_ptr->crend();
  }

  /** @brief Coordinate-wise equality comparison in this container and `rhs`. */
  [[nodiscard]] bool operator==(const Orthant& rhs) const {
    return *(coordinate_ptr) == *(rhs.coordinate_ptr);
  }
  /**
   * @brief Coordinate-wise three-way comparison on entries in this container
   * and `rhs`.
   *
   * @return std::strong_ordering The result of three-way comparison on the
   * first pair of non-equivalent entries in this object and `rhs` should such
   * entries exist. Else, returns `std::strong_ordering::equal`.
   */
  [[nodiscard]] std::strong_ordering operator<=>(const Orthant& rhs) const {
    return *(coordinate_ptr) <=> *(rhs.coordinate_ptr);
  }
};

/**
 * @brief A cubical cell residing in a `CCDIM`-dimensional hypercubical grid.
 *
 * Geometrically, a cubical cell (represented by this object) is a cartesian
 * product of `CCDIM` intervals, which are either open unit intervals or
 * signular points. It is oriented in the grid naturally, with the corners of
 * its closure each being a vertex in the hypercubical grid.
 *
 * The representation of a cubical cell in this program is a pair of orthants:
 * the base orthant and the dual orthant. Each orthant contains one
 * `CCDIM`-dimensional cubical cell (namely, the open `CCDIM`-cube). Each
 * cubical cell can be expressed as the intersection of the closure of two such
 * cubes, which correspond to the `CCDIM`-cubes in the base and dual orthants.
 * The base orthant further is the orthant in which the cubical cell resides.
 *
 * An alternative but equivalent representation offered comprises the base
 * `Orthant` and an unsigned integer (namely, `std::size_t`) representing its
 * shape that has a bit set if the cubical cell has a unit interval along the
 * corresponding axis in its cartesian product. If it instead has a singular
 * point, the bit is not set. The `extent` method calculates this integer, and a
 * `Cube` can be constructed from this paremter instead.
 *
 * Algebraically, this is the cell type of the `CubicalComplex` class. Its chain
 * group corresponds to its dimension as a cell, which is the number of unit
 * intervals in its product (accessible by the `dimension` method). It is
 * hashable and ordered by its underlying `Orthant` objects.
 *
 * @tparam CCDIM The ambient dimension of the hypercubical grid. This matches
 * the `CCDIM` template parameter for the associated `Orthant` class.
 *
 * @sa Orthant, CubicalComplex
 */
template <std::size_t CCDIM>
class Cube {
private:
  Orthant<CCDIM> dual_orthant;
  Orthant<CCDIM> base_orthant;

  [[nodiscard]] static Orthant<CCDIM> dual_from_extent(
      const Orthant<CCDIM>& base, const std::size_t extent) {
    Orthant<CCDIM> dual(base);
    std::size_t axis = 0;
    std::size_t axis_flag = 1;
    for (; axis != CCDIM; ++axis, axis_flag <<= 1) {
      if (!(extent & axis_flag)) {
        --dual[axis];
      }
    }
    return dual;
  }

public:
  /**
   * @brief Default construct a `Cube` with default-constructed `Orthant`
   * objects.
   *
   * The resulting cube is the `CCDIM`-dimensional cubical cell in the
   * orthant at the origin of the hypercubical grid.
   */
  Cube() = default;

  /**
   * @brief Construct a `Cube` by supplying a base and dual `Orthant` objects.
   *
   * See the `Cube` docstring for the geometric description of these objects.
   *
   * @param base The orthant in which the cubical cell represented by the
   * newly-constructed `Cube` object resides.
   * @param dual See the `Cube` docstring.
   */
  Cube(const Orthant<CCDIM>& base, const Orthant<CCDIM>& dual) :
      dual_orthant(dual), base_orthant(base) {}
  /** @overload */
  Cube(const Orthant<CCDIM>& base, Orthant<CCDIM>&& dual) :
      dual_orthant(std::move(dual)), base_orthant(base) {}
  /** @overload */
  Cube(Orthant<CCDIM>&& base, const Orthant<CCDIM>& dual) :
      dual_orthant(dual), base_orthant(std::move(base)) {}
  /** @overload */
  Cube(Orthant<CCDIM>&& base, Orthant<CCDIM>&& dual) :
      dual_orthant(std::move(dual)), base_orthant(std::move(base)) {}

  /**
   * @brief Construct a `Cube` by supplying a base `Orthant` object and an
   * extent paramter.
   *
   * See the `Cube` docstring for the geometric description of these objects.
   * Note that the dual `Orthant` must be constructed from `base` and `extent`,
   * which is linear in the `CCDIM` template parameter.
   *
   * @param base The orthant in which the cubical cell represented by the
   * newly-constructed `Cube` object resides.
   * @param extent See the `Cube` docstring.
   */
  Cube(const Orthant<CCDIM>& base, const std::size_t extent) :
      dual_orthant(dual_from_extent(base, extent)), base_orthant(base) {}
  /** @overload */
  Cube(Orthant<CCDIM>&& base, const std::size_t extent) :
      dual_orthant(dual_from_extent(base, extent)),
      base_orthant(std::move(base)) {}

  /** @brief Return a reference to the base `Orthant` object. */
  [[nodiscard]] Orthant<CCDIM>& base() noexcept {
    return base_orthant;
  }
  /** @overload */
  [[nodiscard]] const Orthant<CCDIM>& base() const noexcept {
    return base_orthant;
  }
  /** @brief Return a reference to the dual `Orthant` object. */
  [[nodiscard]] Orthant<CCDIM>& dual() noexcept {
    return dual_orthant;
  }
  /** @overload */
  [[nodiscard]] const Orthant<CCDIM>& dual() const noexcept {
    return dual_orthant;
  }

  /**
   * @brief Return a reference to the coordinate along `axis` in the base
   * orthant of this `Cube` object.
   */
  [[nodiscard]] HypercubeCoordinate& base(std::size_t axis) {
    return base_orthant[axis];
  }
  /** @overload */
  [[nodiscard]] const HypercubeCoordinate& base(std::size_t axis) const {
    return base_orthant[axis];
  }
  /**
   * @brief Return a reference to the coordinate along `axis` in the dual
   * orthant of this `Cube` object.
   */
  [[nodiscard]] HypercubeCoordinate& dual(std::size_t axis) {
    return dual_orthant[axis];
  }
  /** @overload */
  [[nodiscard]] const HypercubeCoordinate& dual(std::size_t axis) const {
    return dual_orthant[axis];
  }

  /**
   * @brief Compute and return the shape parameter of the cubical cell
   * represented by this `Cube` object.
   *
   * Note this parameter is not stored, and must be calculated each time; this
   * computation is linear in the `CCDIM` template parameter.
   *
   * @return std::size_t An integer that has a bit set iff the cubical cell has
   * a unit interval in its cartesian product along the corresponding axis. See
   * the `Cube` docstring for more details.
   */
  [[nodiscard]] std::size_t extent() const {
    std::size_t extent = 0;
    std::size_t axis = 0;
    std::size_t axis_flag = 1;
    for (; axis != CCDIM; ++axis, axis_flag <<= 1) {
      if (base_orthant[axis] == dual_orthant[axis]) {
        extent += axis_flag;
      }
    }
    return extent;
  }

  /**
   * @brief Compute and return the dimension of the cubical cell represented by
   * this `Cube` object.
   *
   * Note that this is the dimension of the chain group in which it resides,
   * not the ambient dimension of the space.
   *
   * Equivalently, this is the number of axes (between 0 and the template
   * parameter `CCDIM`) along which the `Cube` object has a unit interval in its
   * cartesian product.
   */
  [[nodiscard]] std::size_t dimension() const {
    std::size_t extent_count = 0;
    for (std::size_t axis = 0; axis != CCDIM; ++axis) {
      if (base_orthant[axis] == dual_orthant[axis]) {
        ++extent_count;
      }
    }
    return extent_count;
  }

  /**
   * @brief Equality of `Cube` instances is based on equality of the base and
   * dual orthants in this object and `rhs`.
   */
  [[nodiscard]] bool operator==(const Cube& rhs) const = default;
  /**
   * @brief Three-way comparison operator synthesizes comparison operators
   * between `Cube` instances.
   *
   * Notably, this enables the use of `Cube` as a basis for ordered module
   * classes.
   *
   * The ordering scheme is first three-way-comparison on the base `Orthant`
   * object. If they do not compare equal, return this result. Else, the result
   * of three-way-comparison on the dual `Orthant` object is returned.
   */
  [[nodiscard]] std::strong_ordering operator<=>(const Cube& rhs) const {
    std::strong_ordering base_comp = base_orthant <=> rhs.base_orthant;
    return base_comp != std::strong_ordering::equal
               ? base_comp
               : dual_orthant <=> rhs.dual_orthant;
  }
};

}  // namespace chomp::core

namespace std {

/** @brief `std::hash` specialization for the `Orthant` class template. */
template <size_t CCDIM>
struct hash<chomp::core::Orthant<CCDIM>> {
private:
  constexpr static size_t PRIME = 71;
  constexpr static size_t INITIAL = 15605072507422122298ULL;

public:
  /** @brief Compute and return the hash of `orthant`. */
  [[nodiscard]] size_t operator()(
      const chomp::core::Orthant<CCDIM>& orthant) const noexcept {
    size_t hash = INITIAL;
    for (chomp::core::HypercubeCoordinate coordinate : orthant) {
      hash = PRIME * hash + coordinate;
    }
    return hash;
  }
};

/** @brief `std::hash` specialization for the `Cube` class template. */
template <size_t CCDIM>
struct hash<chomp::core::Cube<CCDIM>> {
private:
  // Inspiration/Credit: Wolfgang Brehm, StackOverflow.

  static constexpr size_t ALT = 0x5555555555555555ULL;
  static constexpr size_t RAN = 17316035218449499591ULL;
  static constexpr size_t DIG = numeric_limits<size_t>::digits;

  [[nodiscard]] static size_t xorshift(size_t n, size_t shift) noexcept {
    return n ^ (n >> shift);
  }

  [[nodiscard]] static size_t distribute(size_t n) noexcept {
    return RAN * xorshift(ALT * xorshift(n, DIG / 2), DIG / 2);
  }

public:
  /** @brief Compute and return the hash of `cube`. */
  [[nodiscard]] size_t operator()(
      const chomp::core::Cube<CCDIM>& cube) const noexcept {
    hash<chomp::core::Orthant<CCDIM>> hasher;
    return hasher(cube.base()) << (DIG / 3) ^ distribute(hasher(cube.dual()));
  }
};

}  // namespace std

namespace chomp::core {

/**
 * @brief This class template is a input iterator (modeling the concept
 * `std::input_iterator`) of dual `Orthant` objects and equivalent integer shape
 * parameters associated with cubes between (in the sense of the face partial
 * order) a maximum and minimum cubical cell; all objects are assumed to reside
 * in the same base orthant.
 *
 * Combining the dereferenced dual `Orthant` object (or alternatively, the shape
 * parameter) with the base `Orthant` object in which the iteration takes place,
 * this is effectively an iterator over all cubical cells (as `Cube` objects)
 * residing in the base orthant. Under this identification, this iterator class
 * template is a breadth-first traversal of (part of) the Hasse diagram in the
 * face partial order of cubes residing in the base orthant.
 *
 * For more information on the identification of base/dual orthants, cubical
 * cells, and shape parameters, see the `Cube` class template docstring.
 *
 * @tparam CCDIM The ambient dimension of the hypercubical grid. This matches
 * the `CCDIM` template parameter of the associated `Orthant` class template.
 *
 * @sa `Cube`, `Orthant`
 */
template <std::size_t CCDIM>
class BFSDualOrthantIterator {
public:
  /** @brief Difference type between iterators. */
  using difference_type = std::ptrdiff_t;
  /**
   * @brief The iterator value type, which is a pair containing a (dual)
   * `Orthant` object and equivalent (in context of a shared base orthant) shape
   * parameter.
   *
   * When iterators of this class template are dereferenced they yield a
   * `const`-qualified references to this type.
   */
  using value_type = std::pair<Orthant<CCDIM>, std::size_t>;
  /**
   * @brief Reference type to `const` `value_type`; returned when objects of
   * this class are dereferenced.
   */
  using reference = const value_type&;
  /**
   * @brief Pointer type to `const` `value_type` returned by the `->` operator.
   */
  using pointer = const value_type*;
  /** @brief Input iterator tag for `iterator_traits` */
  using iterator_concept = std::input_iterator_tag;

private:
  /*
  Implementation Notes:

  This iterator is designed for the cubical Morse reduction-based homology
  algorithm, but may be used in other contexts. The cubical homology algorithm
  uses both the dual orthant and extent parameter interchangeably as their
  representations are both useful. This is why the iterator yields pairs of
  these objects.

  As the prior algorithm will ideally not require iterating the entire graph, we
  compute the partial Hasse diagram lazily. The graph is computed in the
  breadth-first order in the dynamic array `arr`, whose length is 2^n; n being
  the number of inequivalent axes in the supplied minimum and maximum orthants
  (which is the number of vertices in the partial Hasse diagram).

  The computation uses a slow pointer from which new values are computed and a
  fast pointer which points to the newest computed value (returned when the
  iterator is dereferenced). New values are computed from old by decrementing
  axis values along which the old value is not equal to the minimum orthant axis
  value. If the maximum and minimum orthants agree on an axis, this can never
  happen; hence `axes` holds the axes along which differences may occur. It is
  implemented as a static array of length `CCDIM` (a relatively small number),
  although the effective length (demarcated by `axes_end_it`) may be shorter.

  Iteration ends once the slow pointer passes the end of `arr`, where it
  surpasses the fast pointer. The bool conversion operator marks the end, beyond
  which incrementing the iterator is ineffective. Sentinels are not used for
  this iterator.

  Due to the pointers and dynamic arrays this class is not trivially copyable.
  It is unnecessary in the current implementation of algorithms to do so as well
  as quite inefficient as `arr` would need copying as well.

  From the lack of copying and sentinels this class only models the input
  iterator concept and not the forward iterator concept. If necessary, the class
  can be made to model the forward iterator concept but the required features
  may be significantly more inefficient or contrived.
  */

  Orthant<CCDIM> minimum_orthant;
  Orthant<CCDIM> maximum_orthant;

  std::unique_ptr<value_type[]> arr{};
  value_type* slow_it;
  value_type* fast_it;

  std::array<std::size_t, CCDIM> axes{};
  std::size_t* axes_it = axes.data();
  std::size_t* axes_end_it = axes.data();

  void iniitalize_arr(const std::size_t maximum_extent) {
    std::size_t n = 0;
    for (std::size_t axis = 0; axis != CCDIM; ++axis) {
      if (minimum_orthant[axis] != maximum_orthant[axis]) {
        *axes_end_it = axis;
        ++axes_end_it;
        ++n;
      }
    }
    arr = std::make_unique_for_overwrite<value_type[]>(1 << n);
    fast_it = arr.get();
    slow_it = arr.get();
    *fast_it = {maximum_orthant, maximum_extent};
  }

  [[nodiscard]] static Orthant<CCDIM> decrement_axes(
      const Orthant<CCDIM>& base) {
    Orthant<CCDIM> minimum(base);
    for (std::size_t axis = 0; axis != CCDIM; ++axis) {
      --minimum[axis];
    }
    return minimum;
  }

public:
  /** @brief Copy construction is not available for this class template. */
  BFSDualOrthantIterator(const BFSDualOrthantIterator&) = delete;
  /** @brief Default move construction. */
  BFSDualOrthantIterator(BFSDualOrthantIterator&&) = default;
  /** @brief Copy assignment is not available for this class template. */
  BFSDualOrthantIterator& operator=(const BFSDualOrthantIterator&) = delete;
  /** @brief Default move assignment. */
  BFSDualOrthantIterator& operator=(BFSDualOrthantIterator&&) = default;
  /** @brief Default destructor. */
  ~BFSDualOrthantIterator() = default;

  /**
   * @brief Construct a new BFSDualOrthantIterator object by supplying the
   * `minimum` and `maximum` `Orthant` objects the iterator traverses between.
   *
   * While no base `Orthant` object is explicitly supplied, all traversal must
   * occur within the same orthant. Hence, each axis of `maximum` must agree
   * with that of `minimum` or exceed it by one.
   *
   * The optional `maximum_extent` shape parameter marks the shape of the
   * maximum cell (represented by the `maximum` `Orthant` object). This
   * implicitly determines the base orthant and allows for greater control over
   * the range iterated. The default value assumes the maximum cubical cell is
   * the single `CCDIM`-dimensional cell based in each orthant.
   */
  BFSDualOrthantIterator(const Orthant<CCDIM>& minimum,
      const Orthant<CCDIM>& maximum,
      std::size_t maximum_extent = (1 << CCDIM) - 1) :
      minimum_orthant(minimum), maximum_orthant(maximum) {
    iniitalize_arr(maximum_extent);
  }
  /** @overload */
  BFSDualOrthantIterator(const Orthant<CCDIM>& minimum,
      Orthant<CCDIM>&& maximum, std::size_t maximum_extent = (1 << CCDIM) - 1) :
      minimum_orthant(minimum), maximum_orthant(std::move(maximum)) {
    iniitalize_arr(maximum_extent);
  }
  /** @overload */
  BFSDualOrthantIterator(Orthant<CCDIM>&& minimum,
      const Orthant<CCDIM>& maximum,
      std::size_t maximum_extent = (1 << CCDIM) - 1) :
      minimum_orthant(std::move(minimum)), maximum_orthant(maximum) {
    iniitalize_arr(maximum_extent);
  }
  /** @overload */
  BFSDualOrthantIterator(Orthant<CCDIM>&& minimum, Orthant<CCDIM>&& maximum,
      std::size_t maximum_extent = (1 << CCDIM) - 1) :
      minimum_orthant(std::move(minimum)), maximum_orthant(std::move(maximum)) {
    iniitalize_arr(maximum_extent);
  }
  /**
   * @brief Construct a new BFSDualOrthantIterator object iterating over every
   * cubical cell residing in the orthant represented by the `base` `Orthant`
   * object.
   */
  BFSDualOrthantIterator(const Orthant<CCDIM>& base) :
      BFSDualOrthantIterator(decrement_axes(base), base) {}

  /**
   * @brief Dereference the iterator, yielding a constant reference to the
   * currently pointed-to `Orthant` and shape parameter pair.
   *
   * Note that the pointed-to object is always valid, even if the iterator is
   * in its terminal state. Further incrementing will leave it unchanged,
   * however.
   *
   * @return reference A reference to `const` `std::pair` containing an
   * `Orthant` object and an unsigned integer; both equivalently represent a
   * cubical cell currently pointed to by this iterator.
   */
  [[nodiscard]] reference operator*() const noexcept {
    return *fast_it;
  }
  /**
   * @brief Access the pointed-to `Orthant`, shape parameter pair via a pointer.
   *
   * Note that the pointed-to object is always valid, even if the iterator is
   * in its terminal state. Further incrementing will leave it unchanged,
   * however.
   *
   * @return pointer A pointer to `const` `std::pair` containing an `Orthant`
   * object and an unsigned integer; both equivalently represent a cubical cell
   * currently pointed to by this iterator.
   */
  [[nodiscard]] pointer operator->() const noexcept {
    return fast_it;
  }

  /**
   * @brief Increment the iterator, computing the next value and moving the
   * pointer forward.
   *
   * Has no effect if the iterator is terminal, denoted by the bool conversion
   * operator returning `false`.
   *
   * @return BFSDualOrthantIterator& Reference to this object.
   */
  BFSDualOrthantIterator& operator++() {
    for (; slow_it <= fast_it; ++slow_it, axes_it = axes.data()) {
      if (axes_it != axes_end_it
          && minimum_orthant[*axes_it] != slow_it->first[*axes_it]) {
        ++fast_it;
        *fast_it = *slow_it;
        --fast_it->first[*axes_it];
        fast_it->second -= (1 << *axes_it);
        ++axes_it;
        return *this;
      }
    }
    return *this;
  }
  /**
   * @brief Postincrement functions the same as preincrement; no copies are
   * made.
   *
   * @return BFSDualOrthantIterator& Reference to this object.
   */
  BFSDualOrthantIterator& operator++(int) {
    return ++*this;
  }

  /**
   * @brief Return `true` if the iterator is not in its terminal state. Else,
   * return `false`.
   *
   * Further incrementing has no effect once this reads false, and dereferencing
   * will yield the minimum orthant only.
   */
  [[nodiscard]] explicit operator bool() const noexcept {
    return slow_it <= fast_it;
  }
};

/**
 * @brief A forward iterator class over all orthants in the hypercubical grid
 * lexicographically.
 *
 * @tparam CCDIM The ambient dimension of the hypercubical grid. This matches
 * the `CCDIM` template parameter of the associated `Orthant` class template.
 */
template <std::size_t CCDIM>
class LexOrthantIterator {
private:
  Orthant<CCDIM> minimum_orthant;
  Orthant<CCDIM> maximum_orthant;
  Orthant<CCDIM> current_orthant;

  // Set once the iterator repeats to mark the completion of one cycle.
  bool wrapped = false;

public:
  /** @brief Difference type between iterators. */
  using difference_type = std::ptrdiff_t;
  /**
   * @brief The iterator value type are `Orthant` objects. When objects of this
   * class are dereferenced they return a `const` reference to this type.
   */
  using value_type = Orthant<CCDIM>;
  /**
   * @brief Reference type to `const` `value_type`; returned when objects of
   * this class are dereferenced.
   */
  using reference = const value_type&;
  /**
   * @brief Pointer type to `const` `value_type`; returned by the `->` operator.
   */
  using pointer = const value_type*;
  /** @brief Forward iterator tag for `iterator_traits` */
  using iterator_concept = std::forward_iterator_tag;

  /**
   * @brief Default initialize a new LexOrthantIterator object.
   *
   * This is a trivial (one element only) iterator at the origin of the
   * hypercubical grid. However, it can be assigned to, as normal.
   */
  LexOrthantIterator() = default;
  /**
   * @brief Construct a new LexOrthantIterator object iterating over the
   * rectangle of orthants between the `minimum` and `maximum` `Orthant`
   * objects, beginning at `current` and proceeding lexicographically.
   *
   * @param minimum The least (in lexicographical ordering) orthant in the range
   * iterated over by the newly-constructed object. Each axis must be less than
   * or equal to that of `maximum`.
   * @param maximum The greatest (in lexicographical ordering) orthant in the
   * range iterated over by the newly-constructed object. Each axis must be
   * greater than or equal to that of `minimum`.
   * @param current The orthant at which iteration begins; expected to lie
   * within the rectangle formed by `minimum` and `maximum`.
   */
  LexOrthantIterator(const Orthant<CCDIM>& minimum,
      const Orthant<CCDIM>& maximum, const Orthant<CCDIM>& current) :
      minimum_orthant(minimum),
      maximum_orthant(maximum),
      current_orthant(current) {}
  /**
   * @brief This constructor overload uses `minimum` as the beginning orthant,
   * so that the newly-constructed object iterates over the entire rectangle
   * between `minimum` and `maximum`.
   */
  LexOrthantIterator(const Orthant<CCDIM>& minimum,
      const Orthant<CCDIM>& maximum) :
      LexOrthantIterator(minimum, maximum, minimum) {}

  /**
   * @brief Dereference the iterator, yielding a constant reference to the
   * currently pointed-to `Orthant` in the grid.
   *
   * The pointed-to object is always valid, even if the iterator has cycled.
   */
  [[nodiscard]] reference operator*() const noexcept {
    return current_orthant;
  }
  /**
   * @brief Access the pointed-to `Orthant` object via a pointer.
   *
   * The pointed-to object is always valid, even if the iterator has cycled.
   */
  [[nodiscard]] pointer operator->() const noexcept {
    return &current_orthant;
  }

  /**
   * @brief Compare two LexOrthantIterator objects based on the range they
   * iterate over (defined by minimum and maximum `Orthant` objects), the
   * currently pointed-to `Orthant` object, and whether the iterator has cycled.
   */
  [[nodiscard]] bool operator==(const LexOrthantIterator&) const = default;

  /**
   * @brief Increment the iterator, pointing to the next orthant in the grid.
   *
   * @return LexOrthantIterator& A reference to this object.
   */
  LexOrthantIterator& operator++() {
    for (std::size_t axis = 0; axis != CCDIM; ++axis) {
      if (current_orthant[axis] != maximum_orthant[axis]) {
        ++current_orthant[axis];
        return *this;
      }
      current_orthant[axis] = minimum_orthant[axis];
    }
    wrapped = true;
    return *this;
  }
  /**
   * @brief Postincrement operates as expected.
   *
   * @return LexOrthantIterator A copy of this object prior to incrementing.
   */
  LexOrthantIterator operator++(int) {
    LexOrthantIterator temp = *this;
    ++*this;
    return temp;
  }

  /**
   * @brief Determine if the iterator has reached the maximum orthant.
   *
   * This iterator will start from the minimum orthant after reaching the
   * maximum; this operator returns `true` after this iterator has been
   * incremented past the maximum orthant.
   */
  [[nodiscard]] explicit operator bool() const noexcept {
    return !wrapped;
  }
};

/**
 * @brief A forward iterator class template over all cubical cells (as `Cube`
 * objects) in a rectangle of orthants in the `CCDIM`-dimensional hypercubical
 * grid.
 *
 * The order of iteration is first lexicographical in the orthants (as iterated
 * by `LexOrthantIterator`) then lexicographical over the cubical cells in
 * each orthant (under the identification with their dual orthants).
 *
 * This is the iterator class used by the `CubicalComplex` class template.
 *
 * @tparam CCDIM The ambient dimension of the hypercubical grid. This matches
 * the `CCDIM` template parameter of the associated `Cube` class template.
 */
template <std::size_t CCDIM>
class CubeIterator {
private:
  LexOrthantIterator<CCDIM> base_it;
  Orthant<CCDIM> dual_orthant;
  bool is_sentinel = false;

  CubeIterator(bool) :
      is_sentinel(true) {}

public:
  /** @brief Difference type between iterators. */
  using difference_type = std::ptrdiff_t;
  /**
   * @brief The class type returned when iterators of this type are
   * dereferenced.
   */
  using value_type = Cube<CCDIM>;
  /**
   * @brief The class type returned when iterators of this type are
   * dereferenced. Identical to `value_type` as dereferencing returns by value.
   */
  using reference = value_type;
  /** @brief Forward iterator tag for `iterator_traits`. */
  using iterator_concept = std::forward_iterator_tag;

  /**
   * @brief Default initialize a trivial (but valid) iterator that can be
   * assigned to, as normal.
   */
  CubeIterator() = default;
  /**
   * @brief Construct a new CubeIterator object by specifying the minimum and
   * maximum orthants that will be iterated over (specifically, all cubical
   * cells in them).
   *
   * @param minimum The least `Orthant` object iterated over. All cubical cells
   * iterated over reside in orthants greater than or equal to `minimum` in all
   * axes.
   * @param maximum The greatest `Orthant` object iterated over. All cubical
   * cells iterated over reside in orthants less than or equal to `minimum` in
   * all axes.
   */
  CubeIterator(const Orthant<CCDIM>& minimum, const Orthant<CCDIM>& maximum) :
      base_it(LexOrthantIterator<CCDIM>(minimum, maximum)),
      dual_orthant(minimum) {}

  /**
   * @brief Construct and return a sentinel which evaluates equal only to other
   * sentinels and `CubeIterator` objects which have completed iteration.
   *
   * @sa operator==
   */
  [[nodiscard]] static CubeIterator sentinel() noexcept {
    return CubeIterator(true);
  }

  /**
   * @brief Dereference the iterator, returning (by value) a `Cube` object with
   * the currently pointed-to `Orthant` objects.
   */
  [[nodiscard]] value_type operator*() const {
    return Cube<CCDIM>(*base_it, dual_orthant);
  }

  /**
   * @brief Compare iterators by their minimum and maximum bounds as well as
   * the currently pointed-to cubical cell.
   *
   * If this object or `rhs` is a sentinel (see the `sentinel` static member
   * function), `true` is only returned if both are sentinels or the other is an
   * iterator that has finished iteration.
   *
   * @param rhs The `CubeIterator` instance to compare this object against.
   */
  [[nodiscard]] bool operator==(const CubeIterator& rhs) const {
    if (!is_sentinel && !rhs.is_sentinel) {
      return dual_orthant == rhs.dual_orthant && base_it == rhs.base_it;
    }
    if ((is_sentinel || !base_it) && (rhs.is_sentinel || !rhs.base_it)) {
      return true;
    }
    return false;
  }

  /**
   * @brief Increment the iterator to the next cubical cell lexicographically.
   *
   * If the iterator has iterated over all cells in the current orthant, it
   * points to the next orthant (lexicographically) and begins iterating over
   * cells (again, lexicographically) residing in that orthant.
   *
   * @return CubeIterator& A reference to this object.
   */
  CubeIterator& operator++() {
    for (std::size_t axis = 0; axis != CCDIM; ++axis) {
      if (dual_orthant[axis] == (*base_it)[axis]) {
        --dual_orthant[axis];
        return *this;
      }
      dual_orthant[axis] = (*base_it)[axis];
    }
    ++base_it;
    dual_orthant = *base_it;
    return *this;
  }
  /**
   * @brief Postincrement operator functions as expected, incrementing this
   * object but returning a copy that is not incremented.
   *
   * This involves copying several objects and may not be efficient in expensive
   * operations.
   *
   * @return CubeIterator A copy of this object, prior to this object being
   * incremented.
   */
  CubeIterator operator++(int) {
    CubeIterator temp = *this;
    ++*this;
    return temp;
  }

  /**
   * @brief Determine if the iterator has iterated over all cubical cells in
   * the rectangle.
   *
   * This iterator will start from the minimum orthant after reaching the
   * maximum; this operator returns `true` after this iterator has been
   * incremented past the maximum orthant (along with all its cells).
   */
  [[nodiscard]] explicit operator bool() const noexcept {
    return static_cast<bool>(base_it);
  }
};

/**
 * @brief A function object modeling `BoundedGrading` based on a queried cubical
 * cell residing in the closure of top-dimensional (i.e. dimension `CCDIM`).
 * cubes.
 *
 * @tparam CCDIM The ambient dimension of the hypercubical grid. This matches
 * the `CCDIM` template parameter of the associated `Cube` and `Orthant` class
 * templates.
 * @tparam MIN Minimal grading value, i.e. the value returned if the queried
 * cubical cell is in the closure of an included top-dimensional cell.
 * @tparam MAX Maximal grading value, i.e. the value returned if the queried
 * cubical cell is NOT in the closure of an included top-dimensional cell.
 * @tparam SetType The underlying set-like data structure. Expected to be either
 * `std::set` or `std::unordered_set` but any set-like container with
 * sufficiently similar interface can work. The default type is
 * `std::unordered_set` as the stored `Orthant<CCDIM>` type is hashable.
 */
template <std::size_t CCDIM, GradingResultType MIN = 0,
    GradingResultType MAX = 1,
    template <typename...> typename SetType = DefaultSet>
class TopCubeSetGrading {
public:
  /** @brief Input cubical cell type to the grading function. */
  using InputType = Cube<CCDIM>;
  /** @brief The stored `Orthant` type. `Orthant` objects of this type can also
   * be queried, equivalent to querying the top-dimensional cubical cell
   * residing in that orthant.
   */
  using OrthantType = Orthant<CCDIM>;
  /** @brief Minimal grading value. */
  using Minimum = std::integral_constant<GradingResultType, MIN>;
  /** @brief Maximal grading value. */
  using Maximum = std::integral_constant<GradingResultType, MAX>;

private:
  SetType<OrthantType> top_cube_set;

public:
  /**
   * @brief Initialize an empty TopCubeSetGrading function object. In other
   * words, all queries will return the value of the `MAX` template parameter.
   */
  TopCubeSetGrading() = default;
  /**
   * @brief Initialize a TopCubeSetGrading function object by supplying a set of
   * the `Orthant` objects corresponding to top-dimensional cubes with minimal
   * grade.
   *
   * @param top_cubes A set-like (`std::initializer_list`, or the class template
   * parameter `SetType`) collection of `Orthant` objects.
   */
  TopCubeSetGrading(const SetType<OrthantType>& top_cubes) :
      top_cube_set(top_cubes) {}
  /** @overload */
  TopCubeSetGrading(SetType<OrthantType>&& top_cubes) :
      top_cube_set(std::move(top_cubes)) {}
  /** @overload */
  TopCubeSetGrading(std::initializer_list<OrthantType> top_cubes) :
      top_cube_set(top_cubes) {}

  /**
   * @brief Query the grade of `orthant`, i.e. the minimal value (`MIN`) if
   * `orthant` corresponds to a stored top-dimensional cube or the maximal value
   * (`MAX`) otherwise.
   *
   * This specialization is used by the cubical Morse-reduction-based homology
   * algorithm to efficiently query the grades of cubical cells in each orthant.
   */
  GradingResultType operator()(const OrthantType& orthant) const noexcept {
    return top_cube_set.contains(orthant) ? MIN : MAX;
  }

  /**
   * @brief Query the grade of a cubical cell by determining if it is in the
   * closure of a top-dimensional cube of minimal grade.
   *
   * This requires querying inclusion in the underlying set a number of times
   * linear in the `CCDIM` template parameter for each cube, and can thus be
   * costly if done commonly or inefficiently.
   */
  GradingResultType operator()(const InputType& cube) const noexcept {
    for (BFSDualOrthantIterator<CCDIM> bfs_it(cube.dual(), cube.base()); bfs_it;
         ++bfs_it) {
      if (top_cube_set.contains(bfs_it->first)) {
        return MIN;
      }
    }
    return MAX;
  }
};

/**
 * @brief Class implementing a cubical complex embedded in a `CCDIM`-dimensional
 * hypercubical grid.
 *
 * The corresponding cell class is `Cube<CCDIM>`.
 *
 * The cubical complex includes all orthants (and cubical cells residing within
 * them) in the hypercubical grid between some minimum orthant (by default the
 * origin) and a specified maximum orthant.
 *
 * @tparam CCDIM Ambient dimension of the cubical complex i.e. the maximum
 * dimension of `Cube` instances as cubical cells and the number of axes. This
 * matches the `CCDIM` template parameter on the associated `Cube` and `Orthant`
 * class template specializations. Notably, this must also be less than the
 * bit-width of `std::size_t`.
 * @tparam G The type of grading function object, which must model `Grading`.
 * @tparam R The coefficient ring type, which must model `Ring`. Default value
 * is `Z<2>`, i.e. the ring (field) with two elements.
 * @tparam M The chain type, which must model `Module`. The basis type must be
 * `Cube<CCDIM>` and the coefficient ring type must be `R`. The default type is
 * set to `DefaultModule` specialized on `R` and `Cube<CCDIM>`.
 */
template <std::size_t CCDIM, Grading G, Ring R = Z<2>,
    Module M = DefaultModule<Cube<CCDIM>, R>>
requires requires {
  requires CCDIM <= SIZE_T_BITS;
  requires std::same_as<typename G::InputType, Cube<CCDIM>>;
  requires std::same_as<typename M::RingType, R>;
  requires std::same_as<typename M::BasisType, Cube<CCDIM>>;
}
class CubicalComplex {
private:
  Orthant<CCDIM> minimum_orthant;
  Orthant<CCDIM> maximum_orthant;
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
  using CellIterType = CubeIterator<CCDIM>;
  /** @brief Ambient dimension in which the complex is embedded. */
  static constexpr std::size_t dimension = CCDIM;

  /**
   * @brief Initialize a new cubical complex with a maximum orthant and the
   * origin as minimum orthant. Also requires a grading function.
   */
  CubicalComplex(const Orthant<CCDIM>& maximum_orthant,
      const GradingType& grading_function) :
      minimum_orthant(),
      maximum_orthant(maximum_orthant),
      grading_function(grading_function) {}
  /** @overload */
  CubicalComplex(const Orthant<CCDIM>& maximum_orthant,
      GradingType&& grading_function) :
      minimum_orthant(),
      maximum_orthant(maximum_orthant),
      grading_function(std::move(grading_function)) {}
  /**
   * @brief Initialize a new cubical complex with both a maximum and minimum
   * orthant. Also requires a grading function.
   */
  CubicalComplex(const Orthant<CCDIM>& minimum_orthant,
      const Orthant<CCDIM>& maximum_orthant,
      const GradingType& grading_function) :
      minimum_orthant(minimum_orthant),
      maximum_orthant(maximum_orthant),
      grading_function(grading_function) {}
  /** @overload */
  CubicalComplex(const Orthant<CCDIM>& minimum_orthant,
      const Orthant<CCDIM>& maximum_orthant, GradingType&& grading_function) :
      minimum_orthant(minimum_orthant),
      maximum_orthant(maximum_orthant),
      grading_function(std::move(grading_function)) {}

  /**
   * @brief Return (a `const` reference to) the minimum orthant in the complex.
   */
  [[nodiscard]] const Orthant<CCDIM>& minimum() const noexcept {
    return minimum_orthant;
  }
  /**
   * @brief Return the minimum coordinate in the complex along a particular
   * axis.
   */
  [[nodiscard]] HypercubeCoordinate minimum(std::size_t axis) const {
    return minimum_orthant[axis];
  }
  /**
   * @brief Return (a `const` reference to) the maximum orthant in the complex.
   */
  [[nodiscard]] const Orthant<CCDIM>& maximum() const noexcept {
    return maximum_orthant;
  }
  /**
   * @brief Return the maximum coordinate in the complex along a particular
   * axis.
   */
  [[nodiscard]] HypercubeCoordinate maximum(std::size_t axis) const {
    return maximum_orthant[axis];
  }

  /**
   * @brief Query the grading function with an object `t` of type `T`. This
   * function only participates in overload resolution if the grading function
   * can query type `T`.
   *
   * The grading function need only take `Cube<CCDIM>` objects as input, but may
   * grade more types (notably, `Orthant<CCDIM>` objects). This function
   * provides a uniform interface to the grading function object.
   */
  template <typename T>
  requires requires(T t) {
    { grading_function(t) } -> std::convertible_to<GradingResultType>;
  }
  GradingResultType grade(const T& input) {
    return grading_function(input);
  }

  /**
   * @brief Return a forward iterator to the beginning of all cells comprised
   * by this `CubicalComplex` object.
   *
   * These cells are not explicitly stored and are computed by the iterators.
   * Furthermore, `CubicalComplex` objects may comprise a huge number of cubical
   * cells, over which iteration may not be feasible.
   *
   * Enables the use of range-based for loops for objects of this type.
   */
  [[nodiscard]] CellIterType begin() const {
    return CubeIterator<CCDIM>(minimum_orthant, maximum_orthant);
  }

  /**
   * @brief Return a forward iterator past the end of all cells comprised by
   * this `CubicalComplex` object.
   *
   * Compare against the iterator returned by `begin` to demarcate the end of
   * iteration.
   *
   * Enables the use of range-based for loops for objects of this type.
   */
  [[nodiscard]] CellIterType end() const {
    return CubeIterator<CCDIM>::sentinel();
  }

  /**
   * @brief Return the boundary of `cell` in the complex subject to some
   * constraint `cond`.
   *
   * This method synthesizes the `boundary`, `graded_boundary`, and
   * `closure_boundary` function templates for this class along with the
   * generalizations to chain inputs.
   *
   * @param cell The cubical cell to take the boundary of.
   * @param cond A function taking (a `const` reference to) a potential
   * boundary cell and returning a boolean value; if `true`, the cell is added
   * to the boundary and discarded otherwise.
   * @return ChainType The boundary of `cell` subject to `cond`.
   */
  [[nodiscard]] ChainType boundary_if(const CellType& cell,
      const ConditionalType<CellType>& cond) {
    // Implementation follows `Computational Homology` Kaczynski et al.
    RingType coef = one<RingType>();  // axes with extent negate the coefficient
    ChainType result;

    for (std::size_t axis = 0; axis < CCDIM; ++axis) {
      // cell must have extent along this axis to have a boundary
      if (cell.base(axis) == cell.dual(axis)) {
        // No outer cells along maximum edge of complex
        if (cell.base(axis) != maximum_orthant[axis]) {
          Cube<CCDIM> outer_cell(cell);
          ++outer_cell.base(axis);
          if (cond(outer_cell)) {
            result.insert(std::move(outer_cell), coef);
          }
        }

        // Always inner cells
        Cube<CCDIM> inner_cell(cell);
        --inner_cell.dual(axis);
        if (cond(inner_cell)) {
          result.insert(std::move(inner_cell), -coef);
        }

        // Negate coefficient on axes with extent
        coef = -coef;
      }
    }
    return result;
  }

  /**
   * @brief Return the coboundary of `cell` in the complex subject to some
   * constraint `cond`.
   *
   * This method synthesizes the `coboundary`, `graded_coboundary`, and
   * `closure_coboundary` function templates for this class.
   *
   * @param cell The cubical cell to take the coboundary of.
   * @param cond A function taking (a constant reference to) a potential
   * coboundary cell and returning a boolean value; if `true`, the cell is added
   * to the coboundary.
   * @return ChainType The coboundary of `cell` subject to `cond`.
   */
  [[nodiscard]] ChainType coboundary_if(const CellType& cell,
      const ConditionalType<CellType>& cond) {
    RingType coef = one<RingType>();
    ChainType result;

    for (std::size_t axis = 0; axis < CCDIM; ++axis) {
      // cell must not have extent along this axis to have a boundary
      if (cell.base(axis) != cell.dual(axis)) {
        // No inner cells along minimum edge of complex
        if (cell.base(axis) != minimum_orthant[axis]) {
          Cube<CCDIM> inner_cell(cell);
          --inner_cell.base(axis);
          if (cond(inner_cell)) {
            result.insert(std::move(inner_cell), coef);
          }
        }

        // Always outer cells
        Cube<CCDIM> outer_cell(cell);
        ++outer_cell.dual(axis);
        if (cond(outer_cell)) {
          result.insert(std::move(outer_cell), -coef);
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

#endif  // CHOMP_CORE_COMPLEXES_CUBICAL_HPP
