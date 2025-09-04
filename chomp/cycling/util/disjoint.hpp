/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains an efficient implementation of the disjoint-set
 * forest data structure (also known as the union-find or merge-find) for use
 * in cycle identification.
 */

#ifndef CHOMP_CYCLING_UTIL_DISJOINT_HPP
#define CHOMP_CYCLING_UTIL_DISJOINT_HPP

#include <concepts>
#include <cstddef>
#include <numeric>
#include <tuple>
#include <vector>

#ifndef CHOMP_DOXYGEN

namespace chomp::cycling::detail {

/**
 * @brief A disjoint-set forest implementation.
 *
 * This class template represents a collection
 * of disjoint sets of objects of type `T`. The implementation is a collection
 * of trees, each containing all the elements of a singe disjoint set. The root
 * of this tree can be obtained through the `rep` method, which enables checking
 * if two elements are in the same disjoint set. Furthermore, the `merge` method
 * replaces two (possibly distinct) disjoint sets with their union.
 *
 * This implementation uses path compression when searching and weighted unions
 * when merging to increase efficiency.
 *
 * @tparam T The type of object stored in the disjoint sets. This implementation
 * requires that `T` be implicitly convertible to `std::size_t` and that `T` may
 * be constructed from a single argument of type `std::size_t`. Default type is
 * `std::size_t`.
 */
template <typename T = std::size_t>
requires std::convertible_to<T, std::size_t>
         && std::constructible_from<T, std::size_t>
class DisjointSet {
private:
  std::vector<T> parents;
  std::vector<std::size_t> ranks;

public:
  /**
   * @brief Construct a new collection of `size` disjoint sets, each containing
   * a single element from the list `0`, `1`, ..., `size - 1`.
   *
   * @param size The initial number of disjoint singleton sets.
   */
  DisjointSet(std::size_t size = 0) :
      parents(size), ranks(size, 1) {
    std::iota(parents.begin(), parents.end(), T(0));
  }

  /**
   * @brief Return a representative of the disjoint set containing `node`.
   *
   * This same representative is a member of the same disjoint set, and is
   * returned by this method if and only if `node` is in the same disjoint set
   * as the representative. This enables checking if two objects are in the same
   * set.
   *
   * Uses path compression to reduce the path length between `node` and its
   * representative (the root of the tree containing `node`) which makes
   * related representative queries more time-efficient.
   *
   * @param node The object for which a representative is returned.
   * @exception std::out_of_range If `node` is not in any of the disjoint sets.
   * @return T The representative of the disjoint set containing `node`.
   */
  T rep(T node) {
    while (parents.at(node) != node) {
      std::tie(node, parents[node]) = {parents[node], parents[parents[node]]};
    }
    return node;
  }

  /**
   * @brief Replace the (possibly distinct) disjoint sets containing `first` and
   * `second` by their union.
   *
   * Employs a weighted union to attach the smaller disjoint set to root of the
   * larger, reducing path lengths.
   *
   * @param first A member of the first set to merge.
   * @param second A member of the second set to merge.
   * @exception std::out_of_range If either `first` or `second` are not in any
   * of the disjoint sets.
   * @return T The representative of the disjoint set which now contains `first`
   * and `second`.
   */
  T merge(T first, T second) {
    T rep_first = rep(first);
    T rep_second = rep(second);

    // Weighted union attaching the tree of lesser rank to the root of the other
    if (rep_first != rep_second) {
      if (ranks[rep_first] > ranks[rep_second]) {
        parents[rep_second] = rep_first;
        return rep_first;
      }
      if (ranks[rep_first] < ranks[rep_second]) {
        parents[rep_first] = rep_second;
        return rep_second;
      }

      // If the sets have equal rank, this increments the root rank by 1.
      parents[rep_second] = rep_first;
      ++ranks[rep_first];
      return rep_first;
    }
    return rep_first;
  }

  /**
   * @brief Create a new disjoint singleton set in this collection.
   *
   * @return T The representative (and sole element of) the newly-created
   * disjoint set.
   */
  T insert() {
    parents.push_back(T(parents.size()));
    ranks.push_back(1);
    return parents.back();
  }
};

}  // namespace chomp::cycling::detail

#endif  // CHOMP_DOXYGEN

#endif  // CHOMP_CYCLING_UTIL_DISJOINT_HPP
