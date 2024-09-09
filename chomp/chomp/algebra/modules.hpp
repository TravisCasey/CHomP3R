/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header defines module classes built on different container types
 * as well as a module type decider based on basis and ring types.
 */

#ifndef CHOMP_ALGEBRA_MODULES_H
#define CHOMP_ALGEBRA_MODULES_H

#include <chomp/algebra/algebra.hpp>
#include <chomp/util/concepts.hpp>
#include <chomp/util/iterators.hpp>

#include <concepts>
#include <map>
#include <set>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace chomp::core {

namespace detail {

/**
 * @brief Common implementation for `UnorderedMapModule` and `MapModule`.
 *
 * It is preferable to use the template aliases in the `chomp::core` namespace
 * than to directly access this.
 *
 * @tparam T Basis type.
 * @tparam R Coefficient ring type.
 * @tparam MapType Either `std::unordered_map<T, R>` or `std::map<T, R>`.
 *
 * @sa `UnorderedMapModule`, `MapModule`.
 */
template <Basis T, Ring R, typename MapType>
class AssociativeModule {
  MapType cells;
  using MapIterType = typename MapType::iterator;
  using MapCIterType = typename MapType::const_iterator;
  using NodeType = typename MapType::node_type;

public:
  /** @brief Basis element type. */
  using BasisType = T;
  /** @brief Coefficient ring type. */
  using RingType = R;
  /** @brief Iterator type over basis elements. */
  using BasisIterType = KeyIterator<MapCIterType>;

  /**
   * @brief Get the `R` coefficient of basis element `cell` in the module
   * element.
   *
   * A value of zero (in `R`) indicates that `cell` is not present in the
   * element.
   *
   * @param cell Basis element.
   * @return R Coefficient of `cell`.
   */
  [[nodiscard]] R operator[](const T& cell) const noexcept {
    const MapCIterType find_result = cells.find(cell);
    return find_result == cells.cend() ? zero<R>() : find_result->second;
  }

  /**
   * @brief Constant forward iterator to the beginning of the cells in this
   * `R`-linear combination.
   *
   * @return BasisIterType
   */
  [[nodiscard]] BasisIterType begin() const noexcept {
    return KeyIterator(cells.cbegin());
  }
  /**
   * @brief Constant forward iterator to the end of the cells in this
   * `R`-linear combination.
   *
   * @return BasisIterType
   */
  [[nodiscard]] BasisIterType end() const noexcept {
    return KeyIterator(cells.cend());
  }

  /**
   * @brief Insert a basis element `cell` with coefficient `coef` into this
   * module element.
   *
   * If `cell` is already present in the element, then `coef` is added to
   * the coefficient of `cell` in the element.
   *
   * @tparam TFor Forwarding type whose cv-unqualified value matches `T`.
   * @tparam RFor Forwarding type whose cv-unqualified value matches `R`.
   * @param cell Basis element.
   * @param coef Coefficient.
   */
  template <typename TFor, typename RFor>
  requires std::same_as<std::remove_cvref_t<TFor>, T> &&
           std::same_as<std::remove_cvref_t<RFor>, R>
  void insert(TFor&& cell, RFor&& coef) {
    if (coef == zero<R>()) {
      return;
    }

    const std::pair<MapIterType, bool> ins_result =
        cells.insert(std::make_pair(std::forward<TFor>(cell), coef));

    if (!ins_result.second) {
      (ins_result.first)->second += std::forward<RFor>(coef);
      if ((ins_result.first)->second == zero<R>()) {
        cells.erase(ins_result.first);
      }
    }
  }

  /** @brief Reset the element to default initialization state. */
  void clear() {
    cells.clear();
  }

  /**
   * @brief Compoud assignment sum operator computes the formal sum of this
   * module element and `rhs`.
   *
   * Provides move semantics to the `+=` operator implemented for all `Module`
   * types.
   *
   * @param rhs Module element to sum with this element.
   * @return AssociativeModule&
   */
  AssociativeModule& operator+=(AssociativeModule&& rhs) {
    for (MapCIterType it = rhs.cells.cbegin(); it != rhs.cells.cend();) {
      // Avoid invalidation of iterator with postincrement.
      const NodeType nh = rhs.cells.extract(it++);
      insert(std::move(nh.key()), std::move(nh.mapped()));
    }
    return *this;
  }
  /**
   * @brief Compoud assignment difference operator computes the formal
   * difference of this module element and `rhs`.
   *
   * Provides move semantics to the `-=` operator implemented for all `Module`
   * types.
   *
   * @param rhs Module element to subtract from this element.
   * @return AssociativeModule&
   */
  AssociativeModule& operator-=(AssociativeModule&& rhs) {
    for (MapCIterType it = rhs.cells.cbegin(); it != rhs.cells.cend();) {
      // Avoid invalidation of iterator with postincrement.
      const NodeType nh = rhs.cells.extract(it++);
      insert(std::move(nh.key()), -std::move(nh.mapped()));
    }
    return *this;
  }

  /**
   * @brief Equality as elements of the module.
   *
   * @param rhs
   * @return true
   * @return false
   */
  [[nodiscard]] bool operator==(const AssociativeModule& rhs) const {
    return cells == rhs.cells;
  }
};


/**
 * @brief Common implementation for `UnorderedSetModule` and `SetModule`.
 *
 * It is preferable to use the template aliases in the `chomp::core` namespace
 * than to directly access this.
 *
 * @tparam T Basis type.
 * @tparam R Coefficient ring type.
 * @tparam SetType Either `std::unordered_set<T>` or `std::set<T>`.
 *
 * @sa `UnorderedMapModule`, `MapModule`.
 */
template <Basis T, BinaryRing R, typename SetType>
class UniqueModule {
  SetType cells;
  using SetIterType = typename SetType::iterator;
  using SetCIterType = typename SetType::const_iterator;
  using NodeType = typename SetType::node_type;

public:
  /** @copydoc AssociativeModule::BasisType */
  using BasisType = T;
  /** @copydoc AssociativeModule::RingType */
  using RingType = R;
  /** @copydoc AssociativeModule::BasisIterType */
  using BasisIterType = SetCIterType;

  /** @copydoc AssociativeModule::operator[]() */
  [[nodiscard]] R operator[](const T& cell) const {
    return cells.count(cell) == 0 ? zero<R>() : one<R>();
  }

  /** @copydoc AssociativeModule::begin() */
  [[nodiscard]] BasisIterType begin() const noexcept {
    return cells.cbegin();
  }
  /** @copydoc AssociativeModule::end() */
  [[nodiscard]] BasisIterType end() const noexcept {
    return cells.cend();
  }

  /** @copydoc AssociativeModule::insert() */
  template <typename TFor>
  requires std::same_as<std::remove_cvref_t<TFor>, T>
  void insert(TFor&& cell, const R& coef) {
    if (coef == zero<R>()) {
      return;
    }
    const std::pair<SetIterType, bool> ins_result =
        cells.insert(std::forward<TFor>(cell));
    if (!ins_result.second) {
      cells.erase(ins_result.first);
    }
  }

  /** @copydoc AssociativeModule::clear() */
  void clear() {
    cells.clear();
  }

  /** @copydoc AssociativeModule::operator+=() */
  UniqueModule& operator+=(UniqueModule&& rhs) {
    for (SetCIterType it = rhs.cells.cbegin(); it != rhs.cells.cend();) {
      // Avoid invalidation of iterator with postincrement.
      const NodeType nh = rhs.cells.extract(it++);
      insert(std::move(nh.value()), one<RingType>());
    }
    return *this;
  }

  /** @copydoc AssociativeModule::operator-=() */
  UniqueModule& operator-=(UniqueModule&& rhs) {
    for (SetCIterType it = rhs.cells.cbegin(); it != rhs.cells.cend();) {
      // Avoid invalidation of iterator with postincrement.
      const NodeType nh = rhs.cells.extract(it++);
      insert(std::move(nh.value()), one<RingType>());
    }
    return *this;
  }

  /** @copydoc AssociativeModule::operator==() */
  bool operator==(const UniqueModule& rhs) const {
    return cells == rhs.cells;
  }
};


}  // namespace detail


/**
 * @brief This class template implements a free `R`-module on basis set `T`. Its
 * instantiations are elements of this `R`-module, which are formal `R`-linear
 * combinations of elements of the basis set `T`.
 *
 * The implementation uses `std::unordered_set` to store the cells. The
 * coefficient ring `R` is required to be binary-valued (`BinaryRing`). These
 * coefficients are not explicitly stored in this implementation and are instead
 * based on inclusion in the set.
 *
 * @tparam T Basis type modeling `Hashable` concept.
 * @tparam R Ring type modeling `BinaryRing` concept.
 */
template <Hashable T, BinaryRing R>
using UnorderedSetModule = detail::UniqueModule<T, R, std::unordered_set<T>>;

/**
 * @brief This class template implements a free `R`-module on basis set `T`. Its
 * instantiations are elements of this `R`-module, which are formal `R`-linear
 * combinations of elements of the basis set `T`.
 *
 * The implementation uses `std::set` to store the cells. The coefficient ring
 * `R` is required to be binary-valued (`BinaryRing`). These coefficients are
 * not explicitly stored in this implementation and are instead based on
 * inclusion in the set.
 *
 * @tparam T Basis type modeling `Comparable` concept.
 * @tparam Ring type modeling `BinaryRing` concept.
 */
template <Comparable T, BinaryRing R>
using SetModule = detail::UniqueModule<T, R, std::set<T>>;

/**
 * @brief This class template implements a free `R`-module on basis set `T`. Its
 * instantiations are elements of this `R`-module, which are formal `R`-linear
 * combinations of elements of the basis set `T`.
 *
 * The implementation uses `std::unordered_map` to store the cells.
 *
 * @tparam T Basis type modeling `Hashable` concept.
 * @tparam Ring type modeling `Ring` concept.
 */
template <Hashable T, Ring R>
using UnorderedMapModule =
    detail::AssociativeModule<T, R, std::unordered_map<T, R>>;

/**
 * @brief This class template implements a free `R`-module on basis set `T`. Its
 * instantiations are elements of this `R`-module, which are formal `R`-linear
 * combinations of elements of the basis set `T`.
 *
 * The implementation uses `std::map` to store the cells.
 *
 * @tparam T Basis type modeling `Comparable` concept.
 * @tparam Ring type modeling `Ring` concept.
 */
template <Comparable T, Ring R>
using MapModule = detail::AssociativeModule<T, R, std::map<T, R>>;


// Helper level of abstraction for DefaultModule
namespace detail {
#ifndef CHOMP_DOXYGEN

template <bool H, bool B, typename T, typename R>
struct Chooser {
  using type = MapModule<T, R>;
};

template <typename T, typename R>
struct Chooser<true, true, T, R> {
  using type = UnorderedSetModule<T, R>;
};

template <typename T, typename R>
struct Chooser<true, false, T, R> {
  using type = UnorderedMapModule<T, R>;
};

template <typename T, typename R>
struct Chooser<false, true, T, R> {
  using type = SetModule<T, R>;
};

template <Basis T, Ring R>
struct DefaultModuleChooser {
  using type = typename Chooser<Hashable<T>, BinaryRing<R>, T, R>::type;
};

#endif  // CHOMP_DOXYGEN
}  // namespace detail

/**
 * @brief Selects the optimal `Module` type at compile time based on concepts
 * modeled by `T` and `R`
 *
 * @tparam T Basis type.
 * @tparam R Coefficient ring type.
 */
template <Basis T, Ring R>
using DefaultModule = typename detail::DefaultModuleChooser<T, R>::type;

}  // namespace chomp::core

#endif  // CHOMP_ALGEBRA_MODULES_H
