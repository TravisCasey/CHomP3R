/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains the `CubicalMatching` class template for
 * efficiently computing a partial matching on a `CubicalComplex` object for
 * Morse reduction.
 */

#ifndef CHOMP_CORE_HOMOLOGY_CUBICAL_HPP
#define CHOMP_CORE_HOMOLOGY_CUBICAL_HPP

#include <chomp/core/algebra/modules.hpp>
#include <chomp/core/complexes/complexes.hpp>
#include <chomp/core/complexes/cubical.hpp>
#include <chomp/core/complexes/grading.hpp>
#include <chomp/core/homology/morse.hpp>
#include <chomp/core/util/concepts.hpp>
#include <chomp/core/util/constants.hpp>
#include <chomp/core/util/morse.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <concepts>
#include <cstddef>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef CHOMP_DOXYGEN

namespace chomp::core::detail {

/*
 * Simplifies later notation by checking that the type of a chain complex is a
 * specialization of CubicalComplex.
 */
template <typename CC>
concept Cubical = ChainComplex<CC>
                  && std::same_as<CC,
                      CubicalComplex<CC::dimension, typename CC::GradingType,
                          typename CC::RingType, typename CC::ChainType>>;

/*
 * Store selected entries of the grading function queried in CubicalMatching
 * in arrays for faster lookup.
 *
 * Dual orthants are associated with a number (specifically, 2^`CCDIM`, `CCDIM`
 * the dimension of the cubical complex) of base orthants, and so their grades
 * are regularly re-queried. This cache type stores the queries of all dual
 * orthants associated to a `SPAN`-dimensional subgrid of base orthants.
 * Furthermore, half of this subgrid (stored in `upper`) is relevant for the
 * subgrid of base orthants as well, potentially replacing half of the number of
 * sparse grade queries with array lookups.
 */
template <Cubical CC, std::size_t SPAN>
class SubgridCache {
private:
  static constexpr std::size_t CCDIM = CC::dimension;

  std::shared_ptr<CC> complex;
  std::array<std::size_t, SPAN> lengths;
  std::size_t CACHE_SIZE = 1;

  // Whether the cache `lower` (resp. `upper`) holds the correct grading value
  // is indicated by `stored_lower` (resp. `stored_upper`).

  std::unique_ptr<bool[]> stored_lower;
  std::unique_ptr<GradingResultType[]> lower;
  std::unique_ptr<bool[]> stored_upper;
  std::unique_ptr<GradingResultType[]> upper;

public:
  SubgridCache(std::shared_ptr<CC> complex_ptr) :
      complex(complex_ptr) {
    // The number of dual orthants; note there is a dual orthant that is less
    // than the minimum base orthant of the complex along each axis.
    for (std::size_t axis = SPAN - 1; axis < SPAN; --axis) {
      lengths[axis] = complex->maximum(axis) - complex->minimum(axis) + 2;
      CACHE_SIZE *= lengths[axis];
    }
    CACHE_SIZE *= 1 << (CCDIM - SPAN - 1);

    stored_lower = std::make_unique<bool[]>(CACHE_SIZE);
    lower = std::make_unique_for_overwrite<GradingResultType[]>(CACHE_SIZE);
    stored_upper = std::make_unique<bool[]>(CACHE_SIZE);
    upper = std::make_unique_for_overwrite<GradingResultType[]>(CACHE_SIZE);
  }

  /*
   * Query the grade of the dual orthant in the pair (possibly from the cache),
   * and update the cache.
   */
  [[nodiscard]] GradingResultType operator()(
      const std::pair<Orthant<CCDIM>, std::size_t>& cell_pair) const {
    std::size_t pos = 0;
    for (std::size_t axis = SPAN - 1; axis < SPAN; --axis) {
      pos = pos * lengths[axis] + cell_pair.first[axis] - complex->minimum(axis)
            + 1;
    }
    pos = (pos << (CCDIM - SPAN - 1)) + (cell_pair.second >> (SPAN + 1));

    constexpr std::size_t FLAG = 1 << SPAN;
    if (cell_pair.second & FLAG) {
      if (!stored_upper[pos]) {
        upper[pos] = complex->grade(cell_pair.first);
        stored_upper[pos] = true;
      }
      return upper[pos];
    }
    if (!stored_lower[pos]) {
      lower[pos] = complex->grade(cell_pair.first);
      stored_lower[pos] = true;
    }
    return lower[pos];
  }

  /* Pass half of the cache to the next subgrid when compatible. */
  void shift() {
    std::swap(stored_lower, stored_upper);
    std::swap(lower, upper);
    std::fill(stored_upper.get(), stored_upper.get() + CACHE_SIZE, false);
  }

  /* Otherwise, reset the cache completely. */
  void reset() {
    std::fill(stored_lower.get(), stored_lower.get() + CACHE_SIZE, false);
    std::fill(stored_upper.get(), stored_upper.get() + CACHE_SIZE, false);
  }
};

}  // namespace chomp::core::detail

#endif  // CHOMP_DOXYGEN

namespace chomp::core {

/**
 * @brief A `PartialMatching`-derived implementation using the structure of
 * cubical complexes to efficiently compute a partial matching for Morse
 * reduction.
 *
 * This implementation is an optimization (in the case that the grading is
 * determined by top-dimensional cubes) of the algorithm presented (roughly
 * Theorems 1.1, 1.2, 1.3) in `Morse Theoretic Templates for High Dimensional
 * Homology Computation` - Harker, Mischaikow, Spendlove.
 *
 * By considering just one orthant at a time and iterating over its cell in a
 * breadth-first manner (in the face partial order) we efficiently compute the
 * grade of each cell and match entire suborthants at a time, rather than just
 * matching two cells together. These optimizations make this algorithm
 * well-suited to high-dimensional homology computation.
 *
 * @tparam CC The chain complex type (models `ChainComplex`) that the partial
 * matching is computed on. Must be a specialization of `CubicalComplex` with
 * a `grade` member function specialization for `Orthant` objects.
 * @tparam GRADE_MAX The maximum (inclusive) grade for which critical cells are
 * computed. The (co)homology groups of the associated Morse complex will be
 * isomorphic in grades less than or equal to `GRADE_MAX`. Default value
 * includes all cell grades in the complex type `CC`.
 * @tparam DIM_MAX The maximum (inclusive) cell dimension for which critical
 * cells are computed. The (co)homology groups of the associated Morse complex
 * will be isomorphic in dimensions strictly less than `DIM_MAX`. Default value
 * includes all cell dimensions.
 * @tparam CACHE_SPAN A parameter controlling the size of the grade query cache.
 * The cache size is proportional to the product of the number of orthants along
 * the first `CACHE_SPAN` axes in the cubical complex type `CC`. Increasing it
 * reduces repeated queries at the cost of memory overhead. Default value is 0.
 * @tparam MapType Associative container template used to store the computed
 * boundary and coboundary homomorphisms for the induced Morse complex. The cell
 * type of `CC` (which is the corresponding `Cube` specialization) is the key
 * type while the value type is the chain type of `CC`. Default value is
 * `DefaultMap`.
 */
template <ChainComplex CC,
    GradingResultType GRADE_MAX
    = BoundAccess<typename CC::GradingType>::Maximum,
    std::size_t DIM_MAX = std::numeric_limits<std::size_t>::max(),
    std::size_t CACHE_SPAN = 0,
    template <typename...> typename MapType = DefaultMap>
requires detail::Cubical<CC> && requires(CC c, Orthant<CC::dimension> o) {
  { c.grade(o) } -> std::convertible_to<GradingResultType>;
}
class CubicalMatching : public PartialMatching<CC, MapType> {
public:
  /**
   * @brief `ChainComplex` type (specialization of `CubicalComplex`) that this
   * class acts upon.
   */
  using ComplexType = CC;
  /** @brief Coefficient ring type propagated from `CC`. */
  using RingType = typename CC::RingType;
  /** @brief Cell type (specialization of `Cube`) propagated from `CC`. */
  using CellType = typename CC::CellType;
  /** @brief Chain type propagated from `CC`. */
  using ChainType = typename CC::ChainType;
  /** @brief Type detailing the results of a match. */
  using ResultType = MatchResult<CellType, RingType>;

private:
  static constexpr std::size_t DIM = ComplexType::dimension;
  static constexpr GradingResultType GRADE_MIN
      = BoundAccess<typename CC::GradingType>::Minimum;

  /*
  We store the current state of gradient paths first by Orthant, then by the
  critical cell that the gradient path originates from (in `AceMapType`).
  Finally, the extents of the cells at which the gradient path is at currently
  is stored in `ExtentChainType`.

  Grade queries are cached in `CacheType` to reduce repeated queries.
  */

  using ExtentChainType = DefaultModule<std::size_t, RingType>;
  using AceMapType = std::unordered_map<std::size_t, ExtentChainType>;
  using CacheType = detail::SubgridCache<ComplexType, CACHE_SPAN>;

  /*
  A suborthant is a subset of the cubical cells in an orthant that has the same
  structure as an orthant (its Hasse diagram under the face partial order is
  isomorphic to an orthant of dimension no greater than that of the orthant it
  resides in).

  The computed matching on each `Orthant` is represented as a tree of
  `Suborthant` objects under inclusion (so that each cubical cell belongs to
  exactly one leaf `Suborthant` object). Furthermore, the leaf `Suborthant`
  objects comprise only cubical cells of the same grade, and all are matched
  along a single axis (which is `match_axis`). In the case that there is only
  one cubical cell in a `Suborthant` leaf, it is considered critical.

  Regarding cubical cells in the orthant under consideration by their extent
  shape parameter, `extent` xor `match_axis` yields the mate of the cell
  corresponding to `extent`. The `prime_extent` parameter is the first cell
  found in the top-down breadth-first search through the cubical cells whose
  grade differs from that of the top-dimensional cell in the suborthant. If
  there is no such cell, then the `Suborthant` object must be a leaf.

  Matching recurs on the suborthant between `prime_extent` and the least extent
  in the `Suborthant`; it then recurs on the remaining cubical cells in
  the `Suborthant` by partitioning it into smaller `Suborthant` objects.
  Pointers to these `Suborthants` are stored in `partition`.
  */

  struct Suborthant {
    std::size_t prime_extent;
    std::size_t match_axis;
    std::vector<std::unique_ptr<Suborthant>> partition;

    Suborthant(std::size_t maximum, std::size_t minimum) :
        prime_extent(minimum),
        // computes least bit set in maximum but not in minimum
        match_axis((maximum ^ minimum) & (minimum - maximum)) {}
  };

  /*
  LexOrthantIterator iterates lexicographically starting from axis 0. The
  `std::map` specialization `gradient` returns the next `Orthant` object (as its
  key) from its `begin` method that is least in this lexicographical order.
  This function object type implements the correct ordering.
  */

  struct ReverseLexCompare {
    bool operator()(const Orthant<DIM>& lhs, const Orthant<DIM>& rhs) const {
      return std::lexicographical_compare(lhs.crbegin(), lhs.crend(),
          rhs.crbegin(), rhs.crend());
    }
  };

  /*
  Data structure `gradient` stores the current progress of partially-computed
  gradient paths. The map orders the orthant according to `ReverseLexCompare`
  which is the order in which the gradient paths are computed. The value type is
  `{Orthant<DIM>, AceMapType}`, where `AceMapType` is an associative type with
  value type `{ace, chain}`; `chain` is the current gradient path in base
  orthant `orthant` of previously found critical cell, which is in position
  `ace` of `this->critical_cells`.

  As each orthant is iterated over, these gradient paths are computed (via the
  `flow` method) until they reach other critical cells, at which point they are
  added to `boundaries`. These are the boundary maps of the reduced Morse
  complex being computed.
  */

  CacheType cache{this->upper_complex_ptr};
  MapType<CellType, ChainType> boundaries{};
  std::map<Orthant<DIM>, AceMapType, ReverseLexCompare> gradient{};

  /*
  Adds cell `extent` with coefficient to `coef` to the extent chain in
  `gradient[base_orthant][ace_index]`. Checks to ensure these entries exist
  first, and if not, creates them.
  */

  void update_gradient(const Orthant<DIM>& base_orthant,
      const std::size_t ace_index, const std::size_t extent,
      const RingType& coef) {
    (*(*gradient.insert({base_orthant, {}}).first)
            .second.insert({ace_index, ExtentChainType()})
            .first)
        .second.insert(extent, coef);
  }

  /*
  This method computes the partial matching on each suborthant. The maximum cell
  in the suborthant is described by `maximum_orthant` (which is its dual,
  a convenient form for querying its grade)`, `maximum_extent` (convenient for
  querying grade from cache arrays and gradient path storage), and
  `maximum_grade`. The minimum cell is given by `minimum_orthant` (which is
  again its dual), `minimum_extent`, and `minimum_dim` (cell dimension). The
  `critical_extents` vector is updated with critical cells as they are computed.

  A critical cell occurs exactly when there is one cell in the suborthant, i.e.
  `maximum_extent == minimum_extent`. It is only added to `critical` if the cell
  falls below the grade and dimension cutoffs. The `Suborthant` object for a
  critical cell is just the extent itself, 0 `match_axis`, and empty `partition`
  vector.

  A key optimization occurs when `maximum_grade` is equal to the least grade
  found in the complex. As the grading function is order-preserving w.r.t. the
  face partial order, all cells in the suborthant are of the same grade and thus
  all match along the first axis in the suborthant. This is computed by the
  `Suborthant` constructor, and the returned object has no `Suborthant` objects
  in its `partition` vector.

  If those cases do not occur, we traverse the cells of the suborthant in a
  breadth-first (w.r.t. cell dimension) manner until a cell with grade different
  than `maximum_grade` is encountered. The suborthant between this cell and
  the minimum cell is called the prime suborthant; this extent is added to the
  `Suborthant` object and this method recurs on the prime suborthant and
  appended to the suborthant vector. Following that, the rest of this method
  recurs on the rest of the suborthant in a "stepping-up" manner.

  Finally, if no such prime suborthant is found, all cells must be of the same
  grade in the suborthant and so they all match as well.
  */

  [[nodiscard]] std::unique_ptr<Suborthant> match_suborthant(
      std::vector<std::size_t>& critical_extents,
      const Orthant<DIM>& maximum_dual, const std::size_t maximum_extent,
      const GradingResultType maximum_grade, const Orthant<DIM>& minimum_dual,
      const std::size_t minimum_extent, const std::size_t minimum_dim) {
    // Critical cell
    if (maximum_extent == minimum_extent) {
      if (maximum_grade <= GRADE_MAX && minimum_dim <= DIM_MAX) {
        critical_extents.push_back(minimum_extent);
      }
      return std::make_unique<Suborthant>(maximum_extent, minimum_extent);
    }

    // All cells are matched without further analysis
    if (maximum_grade == GRADE_MIN || minimum_dim > DIM_MAX) {
      return std::make_unique<Suborthant>(maximum_extent, minimum_extent);
    }

    // BFS through the cells in this suborthant, searching for grades which
    // differ from the maximum grade; the cell found is the prime extent and is
    // the maximum cell in the prime suborthant.
    BFSDualOrthantIterator bfs_it(minimum_dual, maximum_dual, maximum_extent);
    ++bfs_it;  // skip first cell as it is the maximum cell.

    for (; bfs_it; ++bfs_it) {
      GradingResultType prime_grade = cache(*bfs_it);
      if (prime_grade != maximum_grade) {
        // match the prime suborthant
        Orthant<DIM> prime_dual(bfs_it->first);
        std::size_t prime_extent = bfs_it->second;
        std::unique_ptr<Suborthant> result_ptr
            = std::make_unique<Suborthant>(prime_extent, prime_extent);
        result_ptr->partition.reserve(DIM - minimum_dim);

        result_ptr->partition.push_back(
            match_suborthant(critical_extents, prime_dual, prime_extent,
                prime_grade, minimum_dual, minimum_extent, minimum_dim));

        // match all remaining cells by "stepping up" along axes not in the
        // prime extent
        std::size_t axis = 0;
        std::size_t axis_flag = 1;
        Orthant<DIM> minimum_dual_copy(minimum_dual);
        for (; axis != DIM; ++axis, axis_flag <<= 1) {
          if (prime_dual[axis] != maximum_dual[axis]) {
            ++prime_dual[axis];
            prime_extent ^= axis_flag;
            ++minimum_dual_copy[axis];
            result_ptr->partition.push_back(match_suborthant(critical_extents,
                prime_dual, prime_extent, maximum_grade, minimum_dual_copy,
                minimum_extent ^ axis_flag, minimum_dim + 1));
            --minimum_dual_copy[axis];
          }
        }

        return result_ptr;
      }
    }

    // non-critical (checked prior, i.e. more than one element in suborthant)
    // but all elements are the same grade; no proper suborthants.
    return std::make_unique<Suborthant>(maximum_extent, minimum_extent);
  }

  /*
  Match entire Orthant using the `match_suborthant` method. Return any critical
  cells and the associated `Suborthant` object for computing gradient paths.
  */

  [[nodiscard]] std::pair<std::vector<std::size_t>, std::unique_ptr<Suborthant>>
  match_orthant(const Orthant<DIM>& base_orthant,
      const Orthant<DIM>& minimum_orthant) {
    std::vector<std::size_t> critical_extents;
    GradingResultType maximum_grade = cache({base_orthant, (1 << DIM) - 1});
    std::unique_ptr<Suborthant> result_ptr = match_suborthant(critical_extents,
        base_orthant, (1 << DIM) - 1, maximum_grade, minimum_orthant, 0, 0);
    return {std::move(critical_extents), std::move(result_ptr)};
  }

  /*
  Propagate gradient paths in `base_orthant` (as stored in `gradient`) to
  critical cells or later orthants.

  The boundary homomorphism in the reduced Morse complex is computed by
  following gradient paths, which begin from a critical cell and successively
  take the boundary homomorhpism (in the original complex) then match up along
  the computed partial matching until they reach another critical cell. King
  cells (cells that match down) are discareded in boundaries.

  The tree in `suborthant_ptr` denotes the result of a match in the
  `base_orthant` as computed by the `match_orthant` member function.
  */

  void flow(const std::size_t ace_index, const Orthant<DIM>& base_orthant,
      const std::unique_ptr<Suborthant>& suborthant_ptr,
      ExtentChainType& chain) {
    // leaf
    if (suborthant_ptr->partition.empty()) {
      const std::size_t prime_extent = suborthant_ptr->prime_extent;
      // critical cell
      if (suborthant_ptr->match_axis == 0) {
        boundaries[this->critical_cells[ace_index]].insert(
            CellType(base_orthant, prime_extent), chain[prime_extent]);
        chain.erase(prime_extent);
        return;
      }

      // matches
      // Store in temporary vector `suborthant_cells` as a workaround for
      // iterator invalidation in chain
      std::vector<std::size_t> suborthant_cells;
      for (const std::size_t extent : chain) {
        // prime extent must be face of extent
        if (!(prime_extent & ~extent)) {
          suborthant_cells.push_back(extent);
        }
      }

      Orthant<DIM> base_orthant_copy = base_orthant;
      const std::size_t match_axis = suborthant_ptr->match_axis;
      for (const std::size_t extent : suborthant_cells) {
        if (extent & match_axis) {
          chain.erase(extent);
          continue;  // king
        }

        std::size_t axis = 0;
        std::size_t axis_flag = 1;
        const std::size_t king_extent = extent ^ match_axis;
        RingType boundary_coef
            = (std::popcount(extent % match_axis) % 2 ? -one<RingType>()
                                                      : one<RingType>())
              * chain[extent];
        chain.erase(extent);
        for (; axis != DIM; ++axis, axis_flag <<= 1) {
          if (king_extent & axis_flag) {
            const std::size_t boundary_extent = king_extent ^ axis_flag;
            if (base_orthant[axis] != this->upper_complex_ptr->maximum(axis)) {
              ++base_orthant_copy[axis];
              update_gradient(base_orthant_copy, ace_index, boundary_extent,
                  boundary_coef);
              --base_orthant_copy[axis];
            }
            if (prime_extent & axis_flag) {
              chain.insert(boundary_extent, -boundary_coef);
            }
            boundary_coef = -boundary_coef;
          }
        }
      }
      return;
    }

    // Not a leaf
    // Note Reverse iteration
    for (auto it = suborthant_ptr->partition.crbegin();
         it != suborthant_ptr->partition.crend(); ++it) {
      flow(ace_index, base_orthant, *it, chain);
    }
  }

  /* Used for initialization of orthant iterators. */
  [[nodiscard]] static Orthant<DIM> decrement_axes(Orthant<DIM> orthant) {
    for (std::size_t axis = 0; axis != DIM; ++axis) {
      --orthant[axis];
    }
    return orthant;
  }


public:
  /**
   * @brief Compute an optimized partial matching on the cubical complex
   * `complex`.
   *
   * This class is typically handled through a smart pointer to the base class
   * template `PartialMatching`. It can be used to construct a reduced
   * `MorseComplex` specialization based on this partial matching.
   *
   * @param complex A (smart pointer to a) cubical complex whose grading is
   * determined by top-dimensional cubes.
   */
  CubicalMatching(std::shared_ptr<ComplexType> complex) :
      PartialMatching<CC, MapType>(complex) {
    Orthant<DIM> lower_minimum
        = decrement_axes(this->upper_complex_ptr->minimum());
    Orthant<DIM> lower_maximum
        = decrement_axes(this->upper_complex_ptr->maximum());

    LexOrthantIterator<DIM> orthant_upper_it(this->upper_complex_ptr->minimum(),
        this->upper_complex_ptr->maximum());
    LexOrthantIterator<DIM> orthant_lower_it(lower_minimum, lower_maximum);

    while (orthant_upper_it) {
      const auto [critical_extents, suborthant_ptr]
          = match_orthant(*orthant_upper_it, *orthant_lower_it);

      for (const std::size_t extent : critical_extents) {
        const CellType cell(*orthant_upper_it, extent);
        this->critical_cells.push_back(cell);
        boundaries.insert({cell, ChainType()});
        const ChainType boundary_chain
            = boundary(*(this->upper_complex_ptr), cell);
        for (const CellType& boundary_cell : boundary_chain) {
          update_gradient(boundary_cell.base(), this->critical_cells.size(),
              boundary_cell.extent(), boundary_chain[boundary_cell]);
        }
        this->critical_cells.push_back(cell);
      }
      auto outer_it = gradient.begin();
      if (outer_it != gradient.end()) {
        if ((*outer_it).first == *orthant_upper_it) {
          for (auto& [ace_index, chain] : (*outer_it).second) {
            flow(ace_index, *orthant_upper_it, suborthant_ptr, chain);
          }
          gradient.erase(gradient.cbegin());
        }
      }

      HypercubeCoordinate previous = (*orthant_upper_it)[CACHE_SPAN];
      ++orthant_upper_it;
      ++orthant_lower_it;
      if (previous != (*orthant_upper_it)[CACHE_SPAN]) {
        if (previous == this->upper_complex_ptr->maximum(CACHE_SPAN)) {
          cache.reset();
        } else {
          cache.shift();
        }
      }
    }
  }

private:
  /*
  Used for computing individual matches via the `match` interface - is not used
  during main construction of the matching. It is assumed that there are far too
  many matches to store, so they must be recomputed if they are queried
  individually. However, this method is guaranteed to yield the same matching
  computed during construction.
  */

  [[nodiscard]] ResultType match_helper(const CellType& cell,
      const std::unique_ptr<Suborthant>& suborthant_ptr,
    const std::size_t skip_axes) const {
    if (suborthant_ptr->partition.empty()) {
      if (suborthant_ptr->match_axis == 0) {
        return ResultType(cell, one<RingType>(), Trichotomy::ace);
      }
      RingType coef
          = std::popcount(cell.extent() % suborthant_ptr->match_axis) % 2
                ? -one<RingType>()
                : one<RingType>();
      CellType match_cell(cell.base(),
          cell.extent() ^ suborthant_ptr->match_axis);
      if (suborthant_ptr->match_axis & cell.extent()) {
        return ResultType(match_cell, coef, Trichotomy::king);
      }
      return ResultType(match_cell, coef, Trichotomy::queen);
    }

    std::size_t axis = 0;
    std::size_t axis_flag = 1;
    std::size_t excess_axes = cell.extent() & ~suborthant_ptr->prime_extent;
    auto it = suborthant_ptr->partition.cbegin();
    for (; axis != DIM; ++axis, axis_flag <<= 1, excess_axes >>= 1) {
      if (skip_axes & axis_flag) {
        continue;
      }
      if (excess_axes == 0) {
        return match_helper(cell, *it, ~(suborthant_ptr->prime_extent | (axis_flag - 1)));
      }
      if (!(suborthant_ptr->prime_extent & axis_flag)) {
        ++it;
      }
    }
    return match_helper(cell, *it, ~(suborthant_ptr->prime_extent | (axis_flag - 1)));
  }


public:
  /**
   * @brief A partial ordering on queen cells that is compatible with the
   * ordering on queens induced by a partial matching.
   *
   * This is used for the (co)lift and (co)lower routines.
   */
  [[nodiscard]] bool queen_ordering(const CellType& lhs,
      const CellType& rhs) const override {
    return lhs.base() > rhs.base();
  }
  /**
   * @brief A partial ordering on king cells that is compatible with the
   * ordering on kings induced by a partial matching.
   *
   * This is used for the (co)lift and (co)lower routines.
   */
  [[nodiscard]] bool king_ordering(const CellType& lhs,
      const CellType& rhs) const override {
    return lhs.base() < rhs.base();
  }

  /**
   * @brief Query the match of `cell` under this partial matching.
   *
   * Note that the results may be nonsensical if the grade or dimension of
   * `cell` exceeds the `GRADE_MAX` or `DIM_MAX` template parameters,
   * respectively.
   *
   * @param cell The cell to query the match of.
   * @return ResultType The matched cell, incidence between them, and
   * designation as queen, king, or ace.
   */
  [[nodiscard]] ResultType match(const CellType& cell) override {
    cache.reset();
    std::unique_ptr<Suborthant> suborthant_ptr
        = match_orthant(cell.base(), decrement_axes(cell.base())).second;
    return match_helper(cell, suborthant_ptr, 0);
  }

  [[nodiscard]] std::pair<MapType<CellType, ChainType>,
      MapType<CellType, ChainType>>
  compute_operators() override {
    MapType<CellType, ChainType> coboundaries;

    // Initialize all coboundary chains so they can be inserted into later.
    for (const CellType& cell : this->critical_cells) {
      coboundaries.insert(std::make_pair(cell, ChainType()));
    }

    for (const CellType& cell : this->critical_cells) {
      auto it = boundaries.find(cell);
      if (it != boundaries.cend()) {
        for (const CellType& boundary_cell : (*it).second) {
          coboundaries[boundary_cell].insert(cell, (*it).second[boundary_cell]);
        }
      }
    }

    return std::make_pair(boundaries, std::move(coboundaries));
  }
};

}  // namespace chomp::core

#endif  // CHOMP_CORE_HOMOLOGY_CUBICAL_HPP
