/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains various concepts used for defining and
 * categorizing algebraic classes culminating in the `Module` concept.
 *
 * It also defines functions general to all classes modeling certain concepts.
 */

#ifndef CHOMP_ALGEBRA_ALGEBRA_H
#define CHOMP_ALGEBRA_ALGEBRA_H

#include <chomp/util/concepts.hpp>

#include <concepts>
#include <functional>
#include <iterator>
#include <type_traits>
#include <utility>

namespace chomp::core {

/**
 * @brief The minimal requirements for a class `Q` to implement an additive
 * quasigroup.
 *
 * This concept provides the arithmetic operations required to model a `Group`
 * without the requirement on identity. It also requires that `Q` model the
 * `std::regular` concept.
 *
 * @tparam Q
 */
template <typename Q>
concept Quasigroup = std::regular<Q> && requires(Q a, Q b) {
  { -a } -> std::convertible_to<Q>;
  { a + b } -> std::convertible_to<Q>;
  { a - b } -> std::convertible_to<Q>;
  { a += b } -> std::convertible_to<Q&>;
  { a -= b } -> std::convertible_to<Q&>;
};

/**
 * @brief Return the additive identity of the (quasi)group `G`.
 *
 * This is the primary template using `static_cast` of `0`. This works for
 * fundamental arithmetic types and can work for user classes which define an
 * `int` constructor.
 *
 * Otherwise, this template function can be specialized.
 *
 * @tparam G Models `Quasigroup` concept.
 * @return Additive identity of `G`.
 */
template <Quasigroup G>
[[nodiscard]] constexpr G zero() {
  return static_cast<G>(0);
}

/**
 * @brief The minimal requirements for a class `G` to implement an additive
 * group.
 *
 * This defines the expected interface for group-like object. No checks are
 * made that these operations make sense or that they fulfill group axioms.
 *
 * @tparam G
 */
template <typename G>
concept Group = requires(G a, G b) {
  requires Quasigroup<G>;
  { zero<G>() } -> std::convertible_to<G>;
};

/**
 * @brief The minimal requirements for a class `R` to implement a multiplicative
 * non-unital ring.
 *
 * This concept provides the arithmetic operations required to model a `Ring`
 * without the requirement on multiplicative identity.
 *
 * @tparam R
 */
template <typename R>
concept NonUnitalRing = requires(R a, R b) {
  requires Group<R>;
  { a* b } -> std::convertible_to<R>;
  { a *= b } -> std::convertible_to<R&>;
};

/**
 * @brief Return the multiplicative identity of the ring `R`.
 *
 * This is the primary template using `static_cast` of `1`. This works for
 * fundamental arithmetic types and can work for user classes which define an
 * `int` constructor. Otherwise, this template function can be specialized.
 *
 * @tparam R Models `NonUnitalRing` concept.
 * @return Multiplicative identity of `R`.
 */
template <NonUnitalRing R>
[[nodiscard]] constexpr R one() {
  return static_cast<R>(1);
}

/**
 * @brief The minimal requirements for a class `R` to implement a multiplicative
 * ring.
 *
 * This defines the expected interface for ring-like object. No checks are
 * made that these operations make sense or that they fulfill ring axioms.
 *
 * Notably, we also assume (but do not explicitly check) that the ring is an
 * integral domain, i.e., it has no zero divisors.
 *
 * @tparam R
 */
template <typename R>
concept Ring = requires(R a, R b) {
  requires NonUnitalRing<R>;
  { one<R>() } -> std::convertible_to<R>;
};

/**
 * @brief Classes implementing this concept represent the ring (field) with two
 * elements.
 *
 * While algebraically there is only one such ring, there may be different data
 * structures implementing it.
 *
 * Module classes over a class implementing `BinaryRing` can be represented
 * efficiently by only implicitly storing these coefficients.
 *
 * @tparam R
 */
template <typename R>
concept BinaryRing = requires {
  requires Ring<R>;
  requires one<R>() != zero<R>();
  requires one<R>() + one<R>() == zero<R>();
};


/**
 * @brief Types `T` modeling this concept model either `Hashable` or
 * `Comparable`, enabling their use as a basis set in one of the module classes.
 *
 * @tparam T
 */
template <typename T>
concept Basis = AssociativeKey<T>;

/**
 * @brief Requirements for an iterator type `I` to be a basis iterator for
 * module types.
 *
 * Its `value_type` (from `std::iterator_traits`) must be (convertible to) `T`,
 * which is checked to be the module type's basis type in the `Module` concept.
 *
 * @tparam I Forward iterator with value type (convertible to) `T`.
 * @tparam T
 *
 * @sa Module
 */
template <typename I, typename T>
concept ModuleIterator =
    std::forward_iterator<I> && std::convertible_to<std::iter_value_t<I>, T>;

/**
 * @brief Most non-arithmetic requirements for a class `M` to model a module.
 *
 * The class `M` must declare `BasisType` (modeling `Basis`), `RingType`
 * (modeling `Ring`), and `BasisIterType` (modeling `ModuleIterator`) public
 * member types as well as implement certain basic functions. It must also
 * model `std::regular`.
 *
 * With these requirements, the arithmetic (sum, difference, scalar product) and
 * other basic functions are implemented for all classes modeling
 * `ModulePrecursor`. These are sufficient for `M` to model `Module`, as well.
 *
 * @tparam M
 */
template <typename M>
concept ModulePrecursor = requires(M mod) {
  // Required member types
  typename M::BasisType;
  typename M::RingType;
  typename M::BasisIterType;

  // Constraints on types
  requires std::regular<M>;
  requires Ring<typename M::RingType>;
  requires Basis<typename M::BasisType>;
  requires ModuleIterator<typename M::BasisIterType, typename M::BasisType>;

  // Required methods
  {
    mod[std::declval<typename M::BasisType>()]
  } -> std::convertible_to<typename M::RingType>;
  mod.insert(
      std::declval<typename M::BasisType>(),
      std::declval<typename M::RingType>()
  );
  { mod.begin() } -> std::convertible_to<typename M::BasisIterType>;
  { mod.end() } -> std::convertible_to<typename M::BasisIterType>;
  mod.clear();
};

/**
 * @brief Compound assignment sum operator computes the formal sum of `lhs` and
 * `rhs` in the module `M`.
 *
 * Move semantics are dependent upon the module class implementation and are not
 * provided by default.
 *
 * @tparam M Module type modeling `ModulePrecursor`.
 * @param lhs
 * @param rhs
 * @return M&
 */
template <ModulePrecursor M>
M& operator+=(M& lhs, const M& rhs) {
  for (const typename M::BasisType cell : rhs) {
    lhs.insert(cell, rhs[cell]);
  }
  return lhs;
}

/**
 * @brief Compound assignment difference operator computes the formal difference
 * of `lhs` and `rhs` in the module `M`.
 *
 * Move semantics are dependent upon the module class implementation and are not
 * provided by default.
 *
 * @tparam M Module type modeling `ModulePrecursor`.
 * @param lhs
 * @param rhs
 * @return M&
 */
template <ModulePrecursor M>
M& operator-=(M& lhs, const M& rhs) {
  for (const typename M::BasisType cell : rhs) {
    lhs.insert(cell, -rhs[cell]);
  }
  return lhs;
}

/**
 * @brief Compound assignment scalar product operator computes the scalar
 * product of `lhs` and scalar `rhs` in the module `M`.
 *
 * @tparam M Module type modeling `ModulePrecursor`.
 * @param lhs
 * @param rhs
 * @return M&
 */
template <ModulePrecursor M>
M& operator*=(M& lhs, const typename M::RingType& rhs) {
  // For predefined module classes, calls to `insert` may invalidate pointers
  // if the new coefficient would be 0.
  // Presuming the coefficient ring is an integral domain, this prevents that
  // from happening.
  if (rhs == zero<typename M::RingType>()) {
    lhs.clear();
    return lhs;
  }
  for (const typename M::BasisType& cell : lhs) {
    lhs.insert(cell, (rhs - one<typename M::RingType>()) * lhs[cell]);
  }
  return lhs;
}

/**
 * @brief Negates each coefficient in the module element.
 *
 * @tparam Left Forwarding type modeling `ModulePrecursor`.
 * @param elem Module element.
 * @return New module element that is the negation of `elem`.
 */
template <typename Left>
requires ModulePrecursor<std::remove_cvref_t<Left>>
[[nodiscard]] std::remove_cvref_t<Left> operator-(Left&& elem) {
  std::remove_cvref_t<Left> result(std::forward<Left>(elem));
  result *= -one<typename std::remove_cvref_t<Left>::RingType>();
  return result;
}

/**
 * @brief Sum operator computes formal sum of `rhs` and `lhs` in their common
 * algebraic module.
 *
 * @tparam Left Forwarding type modeling `ModulePrecursor`, whose cv-unqualified
 * value matches that of `Right`.
 * @tparam Right Forwarding type modeling `ModulePrecursor`, whose
 * cv-unqualified value matches that of `Left`.
 * @param lhs Module element.
 * @param rhs Module element.
 * @return New module element that is the sum of `lhs` and `rhs`.
 */
template <typename Left, typename Right>
requires ModulePrecursor<std::remove_cvref_t<Left>> &&
         std::same_as<std::remove_cvref_t<Left>, std::remove_cvref_t<Right>>
[[nodiscard]] std::remove_cvref_t<Left> operator+(Left&& lhs, Right&& rhs) {
  std::remove_cvref_t<Left> result(std::forward<Left>(lhs));
  result += std::forward<Right>(rhs);
  return result;
}

/**
 * @brief Difference operator computes formal difference of `rhs` and `lhs` in
 * their common algebraic module.
 *
 * @tparam Left Forwarding type modeling `ModulePrecursor`, whose cv-unqualified
 * value matches that of `Right`.
 * @tparam Right Forwarding type modeling `ModulePrecursor`, whose
 * cv-unqualified value matches that of `Left`.
 * @param lhs Module element.
 * @param rhs Module element.
 * @return New module element that is the difference of `lhs` and `rhs`.
 */
template <typename Left, typename Right>
requires ModulePrecursor<std::remove_cvref_t<Left>> &&
         std::same_as<std::remove_cvref_t<Left>, std::remove_cvref_t<Right>>
[[nodiscard]] std::remove_cvref_t<Left> operator-(Left&& lhs, Right&& rhs) {
  std::remove_cvref_t<Left> result(std::forward<Left>(lhs));
  result -= std::forward<Right>(rhs);
  return result;
}

/**
 * @brief Product operator computes formal scalar product of `elem` with `coef`
 * in the algebraic module of `elem`.
 *
 * @tparam M Forwarding type modeling `ModulePrecursor`.
 * @tparam R Forwarding type that is the coefficient ring type of the module
 * (referred to by) `M`.
 * @param elem Module element.
 * @param coef Coefficient.
 * @return New module element that is the product of `elem` with `coef`.
 */
template <typename M, typename R>
requires ModulePrecursor<std::remove_cvref_t<M>> &&
         std::same_as<
             std::remove_cvref_t<R>, typename std::remove_cvref_t<M>::RingType>
[[nodiscard]] std::remove_cvref_t<M> operator*(M&& elem, R&& coef) {
  std::remove_cvref_t<M> result(std::forward<M>(elem));
  result *= coef;
  return result;
}

/**
 * @brief Product operator computes formal scalar product of `elem` with `coef`
 * in the algebraic module of `elem`.
 *
 * @tparam R Forwarding type that is the coefficient ring type of the module
 * (referred to by) `M`.
 * @tparam M Forwarding type modeling `ModulePrecursor`.
 * @param coef Coefficient.
 * @param elem Module element.
 * @return New module element that is the product of `elem` with `coef`.
 */
template <typename R, typename M>
requires ModulePrecursor<std::remove_cvref_t<M>> &&
         std::same_as<
             std::remove_cvref_t<R>, typename std::remove_cvref_t<M>::RingType>
[[nodiscard]] std::remove_cvref_t<M> operator*(R&& coef, M&& elem) {
  std::remove_cvref_t<M> result(std::forward<M>(elem));
  result *= coef;
  return result;
}

/**
 * @brief Template specialization of additive identity function for modules.
 *
 * The identity element is the empty element, which may be defined by default
 * initialization. Further specializations can be made for types where this is
 * not the case.
 *
 * @tparam M Module type
 * @return M Additive identity of `M`.
 */
template <typename M>
requires ModulePrecursor<M> && Quasigroup<M>
[[nodiscard]] constexpr M zero() {
  return M();
}

/**
 * @brief Minimal arithmetic, value, and method requirements for a class to
 * implement an algebraic module for the purposes of CHomP.
 *
 * @tparam M
 */
template <typename M>
concept Module = requires(M mod) {
  requires ModulePrecursor<M>;
  requires Group<M>;
  { std::declval<typename M::RingType>() * mod } -> std::convertible_to<M>;
  { mod* std::declval<typename M::RingType>() } -> std::convertible_to<M>;
  { mod *= std::declval<typename M::RingType>() } -> std::convertible_to<M&>;
};

/**
 * @brief Alias for the function type expected by the linear map interface to
 * module classes.
 *
 * The return value is a forward iterator that returns basis, coefficient pairs.
 * This avoids forcing allocation of these elements.
 *
 * @tparam M Module Type
 * @tparam I Forward iterator yielding basis, coefficient pairs which are the
 * result of the function.
 *
 * @sa linear_apply()
 */
template <typename M, typename I>
requires Module<M> && std::forward_iterator<I> &&
             std::convertible_to<
                 std::iter_value_t<I>,
                 std::pair<typename M::BasisType, typename M::RingType>>
using LinearMap = std::function<std::pair<I, I>(const typename M::BasisType&)>;

/**
 * @brief Apply the linear map `func` to each the module element `elem`.
 *
 * @tparam M Module Type
 * @tparam I Iterator type returned (as a begin, end pair) from `func`.
 * @param elem Input module element.
 * @param func Linear map to apply to `elem`.
 * @return M A new module element that is the result of `func` applied to
 * `elem`.
 *
 * @sa LinearMap
 */
template <Module M, std::forward_iterator I>
[[nodiscard]] M linear_apply(const M& elem, const LinearMap<M, I>& func) {
  M result;
  typename M::RingType coef;
  I begin_it, end_it;

  // Apply func to each basis element in elem; multiply the coefficients of the
  // result by the cells coefficient and insert into a new element.
  for (const typename M::BasisType& cell : elem) {
    coef = elem[cell];
    for (std::tie(begin_it, end_it) = func(cell); begin_it != end_it;
         begin_it++) {
      result.insert(
          std::move(begin_it->first), coef * std::move(begin_it->second)
      );
    }
  }

  return result;
}

}  // namespace chomp::core

#endif  // CHOMP_ALGEBRA_ALGEBRA_H
