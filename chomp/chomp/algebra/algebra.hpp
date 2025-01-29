/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains various concepts used for defining and
 * categorizing algebraic classes including the `Ring` and `Module` concepts.
 *
 * It also defines function templates specialized for all classes modeling
 * certain concepts.
 */

#ifndef CHOMP_ALGEBRA_ALGEBRA_H
#define CHOMP_ALGEBRA_ALGEBRA_H

#include <concepts>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace chomp::core {

/**
 * @brief The minimal requirements for a class `G` to model an additive group.
 *
 * All types modeling this concept have the `zero` function template specialized
 * for them; the primary template uses `static_cast<G>(0)` to yield the additive
 * identity for fundamental arithmetic types. Custom types modeling `Group` can
 * specialize the `zero` function for desired behavior.
 *
 * Note that no algebraic axioms are checked here; it is assumed that all
 * operations make sense.
 *
 * @tparam G The type modeling an addtive group.
 */
template <typename G>
concept Group = std::regular<G> && requires(G a, G b) {
  { -a } -> std::convertible_to<G>;
  { a + b } -> std::convertible_to<G>;
  { a - b } -> std::convertible_to<G>;
  { a += b } -> std::convertible_to<G&>;
  { a -= b } -> std::convertible_to<G&>;
};


/**
 * @brief The minimal requirements for a class `R` to model a multiplicative
 * ring (that is also an additive group).
 *
 * All types modeling this concept have the `zero` and `one` function templates
 * specialized for them; the primary templates use `static_cast<R>(0)` and
 * `static_cast<R>(1)`, respectively, to yield the additive and multiplicative
 * identities for fundamental arithmetic types. Custom types modeling `Ring`
 * can specialize the `one` function for desired behavior.
 *
 * Note that no algebraic axioms are checked here; it is assumed that all
 * operations make sense.
 *
 * @tparam R The type modeling a multiplicative ring.
 *
 * @sa `zero`, `one`.
 */
template <typename R>
concept Ring = std::regular<R> && requires(R a, R b) {
  { -a } -> std::convertible_to<R>;
  { a + b } -> std::convertible_to<R>;
  { a - b } -> std::convertible_to<R>;
  { a* b } -> std::convertible_to<R>;
  { a += b } -> std::convertible_to<R&>;
  { a -= b } -> std::convertible_to<R&>;
  { a *= b } -> std::convertible_to<R&>;
};

/**
 * @brief Return the additive identity on the `Ring`-modeling type `R`.
 *
 * This is the primary template using `static_cast<R>(0)`. This works for
 * fundamental arithmetic types and can work for user classes which define an
 * `int` constructor.
 *
 * Otherwise, this function template should be specialized.
 *
 * @tparam R Models `Ring` concept.
 * @return The additive identity of `R`.
 */
template <Group R>
[[nodiscard]] constexpr R zero() {
  return static_cast<R>(0);
}

/**
 * @brief Return the multiplicative identity on the `Ring`-modeling type `R`.
 *
 * This is the primary template using `static_cast<R>(1)`. This works for
 * fundamental arithmetic types and can work for user classes which define an
 * `int` constructor.
 *
 * Otherwise, this function template should be specialized.
 *
 * @tparam R Models `Ring` concept.
 * @return Multiplicative identity of `R`.
 */
template <Ring R>
[[nodiscard]] constexpr R one() {
  return static_cast<R>(1);
}

/**
 * @brief Classes implementing this concept represent the ring (in fact, field)
 * with two elements.
 *
 * While algebraically there is only one such ring, there may be different data
 * structures implementing it.
 *
 * Module classes over a class implementing `BinaryRing` can be represented
 * efficiently by only implicitly storing these coefficients.
 *
 * @tparam R `Ring`-modeling type that models the ring with two elements.
 */
template <typename R>
concept BinaryRing = requires {
  requires Ring<R>;
  requires one<R>() != zero<R>();
  requires one<R>() + one<R>() == zero<R>();
};

/**
 * @brief Return the multiplicative inverse of `value` in `R`, should it exist.
 * Throw `std::domain_error` otherwise.
 *
 * This is the primary template using division (`/` operator), which works for
 * all floating point types. For custom types this function template may need to
 * be specialized.
 *
 * The `invertible` function template is a safe way to check for invertibility
 * first.
 *
 * @tparam R Models `Ring` concept.
 * @param value The ring element to be inverted.
 * @return R The multiplicative inverse of `value`.
 */
template <Ring R>
[[nodiscard]] R invert(const R& value) {
  if (value == zero<R>()) {
    throw std::domain_error(
        "Attempted to invert a ring element that is not a unit.");
  }
  return one<R>() / value;
}
/** @brief A partial specialization of `invert` for unsigned integral types. */
template <typename R>
requires Ring<R> && std::unsigned_integral<R>
[[nodiscard]] R invert(const R& value) {
  if (value == one<R>()) {
    return one<R>();
  }
  throw std::domain_error(
      "Attempted to invert a ring element that is not a unit.");
}
/** @brief A partial specialization of `invert` for signed integral types. */
template <typename R>
requires Ring<R> && std::signed_integral<R>
[[nodiscard]] R invert(const R& value) {
  if (value == one<R>()) {
    return one<R>();
  }
  if (value == -one<R>()) {
    return -one<R>();
  }
  throw std::domain_error(
      "Attempted to invert a ring element that is not a unit.");
}

/**
 * @brief Check invertibility of `value` in `R`. If the result is `true`, the
 * `invert` function should correctly return when called on `value`.
 *
 * This is the primary template for any ring (field) `R` in which all values
 * are invertible except zero. This works for floating point types. This
 * template is specialized for unsigned integral types (at which only 1 is
 * treated as invertible) and signed integral types (at which 1 and -1 are
 * invertible).
 *
 * However, there may be inconsistencies in the machine number representation of
 * these types and it is recommended to define custom ring types of use the
 * cycling rings `Z<p>`.
 *
 * @tparam R Models `Ring` concept.
 * @param value The ring instantiation to query invertibility.
 * @return bool Whether `value` is invertible.
 */
template <Ring R>
[[nodiscard]] bool invertible(const R& value) {
  return value != zero<R>();
}
/**
 * @brief A partial specialization of `invertible` for unsigned integral types.
 */
template <typename R>
requires Ring<R> && std::unsigned_integral<R>
[[nodiscard]] bool invertible(const R& value) {
  return value == one<R>();
}
/**
 * @brief A partial specialization of `invertible` for signed integral types.
 */
template <typename R>
requires Ring<R> && std::signed_integral<R>
[[nodiscard]] bool invertible(const R& value) {
  return value == one<R>() || value == -one<R>();
}

/**
 * @brief This concept encapsulates the requirements for a class `M` to model an
 * algebraic module.
 *
 * The class `M` must behave similarly to an associative container with member
 * type `BasisType` as key and member type `RingType` as value; furthermore,
 * instantiations of the member type `BasisIterType` are returned by the `begin`
 * and `end` member functions to iterate over the keys.
 *
 * Algebraically, an object of type modeling `Module` models a linear
 * combination of `BasisType` objects with coefficients in `RingType`.
 *
 * Types modeling this concept have specializations of the sum, difference,
 * scalar product operators among other basic functions.
 *
 * @tparam M Models `Module` concept.
 */
template <typename M>
concept Module = requires(M mod) {
  // Required member types
  typename M::BasisType;
  typename M::RingType;
  typename M::BasisIterType;

  // Constraints on types
  requires std::regular<M>;
  requires Ring<typename M::RingType>;
  requires std::forward_iterator<typename M::BasisIterType>;
  requires std::convertible_to<std::iter_value_t<typename M::BasisIterType>,
      typename M::BasisType>;

  // Required methods
  {
    mod[std::declval<typename M::BasisType>()]
  } -> std::convertible_to<typename M::RingType>;
  mod.insert(std::declval<typename M::BasisType>(),
      std::declval<typename M::RingType>());
  { mod.begin() } -> std::convertible_to<typename M::BasisIterType>;
  { mod.end() } -> std::convertible_to<typename M::BasisIterType>;
  mod.clear();
  mod.erase(std::declval<typename M::BasisType>());
};

/**
 * @brief Compound assignment sum operator computes the formal sum of `lhs` and
 * `rhs` in the module `M`.
 *
 * Move semantics are dependent upon the module class implementation and are not
 * provided by default.
 *
 * @tparam M Module type modeling `Module`.
 * @return M& Reference to `lhs` following compound assignment sum with `rhs`.
 */
template <Module M>
M& operator+=(M& lhs, const M& rhs) {
  for (const typename M::BasisType& cell : rhs) {
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
 * @tparam M Module type modeling `Module`.
 * @return M& Reference to `lhs` following compound assignment difference with
 * `rhs`.
 */
template <Module M>
M& operator-=(M& lhs, const M& rhs) {
  for (const typename M::BasisType& cell : rhs) {
    lhs.insert(cell, -rhs[cell]);
  }
  return lhs;
}

/**
 * @brief Compound assignment product operator computes the scalar product of
 * `lhs` and scalar `rhs` in the module `M`.
 *
 * @tparam M Module type modeling `Module`.
 * @return M& Reference to `lhs` following compound assignment product with
 * `rhs`.
 */
template <Module M>
M& operator*=(M& lhs, const typename M::RingType& rhs) {
  // For predefined module classes, calls to `insert` may invalidate pointers
  // if the new coefficient would be 0. Presuming the coefficient ring is an
  // integral domain, this prevents that from happening.
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
 * @brief Negate each coefficient in the module object `elem`.
 *
 * @tparam Left A type whose cv-unqualified non-reference type models `Module`.
 * @return New module object that is the negation of `elem`.
 */
template <typename Left>
requires Module<std::remove_cvref_t<Left>>
[[nodiscard]] std::remove_cvref_t<Left> operator-(Left&& elem) {
  std::remove_cvref_t<Left> result(std::forward<Left>(elem));
  result *= -one<typename std::remove_cvref_t<Left>::RingType>();
  return result;
}

/**
 * @brief Sum operator computes formal sum of `rhs` and `lhs` in their common
 * algebraic module.
 *
 * @tparam Left A type whose cv-unqualified non-reference type models `Module`
 * and agrees with that of `Right`.
 * @tparam Right A type whose cv-unqualified non-reference type models `Module`
 * and agrees with that of `Left`.
 * @return New module object that is the sum of `lhs` and `rhs`.
 */
template <typename Left, typename Right>
requires Module<std::remove_cvref_t<Left>>
         && std::same_as<std::remove_cvref_t<Left>, std::remove_cvref_t<Right>>
[[nodiscard]] std::remove_cvref_t<Left> operator+(Left&& lhs, Right&& rhs) {
  std::remove_cvref_t<Left> result(std::forward<Left>(lhs));
  result += std::forward<Right>(rhs);
  return result;
}

/**
 * @brief Difference operator computes formal difference of `rhs` and `lhs` in
 * their common algebraic module.
 *
 * @tparam Left A type whose cv-unqualified non-reference type models `Module`
 * and agrees with that of `Right`.
 * @tparam Right A type whose cv-unqualified non-reference type models `Module`
 * and agrees with that of `Left`.
 * @return New module object that is the difference of `lhs` and `rhs`.
 */
template <typename Left, typename Right>
requires Module<std::remove_cvref_t<Left>>
         && std::same_as<std::remove_cvref_t<Left>, std::remove_cvref_t<Right>>
[[nodiscard]] std::remove_cvref_t<Left> operator-(Left&& lhs, Right&& rhs) {
  std::remove_cvref_t<Left> result(std::forward<Left>(lhs));
  result -= std::forward<Right>(rhs);
  return result;
}

/**
 * @brief Product operator computes formal scalar product of `elem` with `coef`
 * in the algebraic module of `elem`.
 *
 * @tparam M A type whose cv-unqualified non-reference type models `Module`.
 * @tparam R (Possibly cv-qualified reference to) the `RingType` of `M`.
 * @return New module object that is the scalar product of `elem` with `coef`.
 */
template <typename M, typename R>
requires Module<std::remove_cvref_t<M>>
         && std::same_as<std::remove_cvref_t<R>,
             typename std::remove_cvref_t<M>::RingType>
[[nodiscard]] std::remove_cvref_t<M> operator*(M&& elem, R&& coef) {
  std::remove_cvref_t<M> result(std::forward<M>(elem));
  result *= coef;
  return result;
}

/**
 * @brief Product operator computes formal scalar product of `elem` with `coef`
 * in the algebraic module of `elem`.
 *
 * @tparam R (Possibly cv-qualified reference to) the `RingType` of `M`.
 * @tparam M A type whose cv-unqualified non-reference type models `Module`.
 * @return New module object that is the scalar product of `elem` with `coef`.
 */
template <typename R, typename M>
requires Module<std::remove_cvref_t<M>>
         && std::same_as<std::remove_cvref_t<R>,
             typename std::remove_cvref_t<M>::RingType>
[[nodiscard]] std::remove_cvref_t<M> operator*(R&& coef, M&& elem) {
  std::remove_cvref_t<M> result(std::forward<M>(elem));
  result *= coef;
  return result;
}

/**
 * @brief Alias template for the function type expected by the linear map
 * interface to module classes.
 *
 * These maps take an element of the basis type of `M` as input and return a new
 * object of type `M`; these can be applied linearly over the corresponding
 * module type.
 *
 * @tparam M Type modeling `Module`.
 */
template <Module M>
using LinearMap = std::function<M(const typename M::BasisType&)>;

/**
 * @brief Apply the linear map `func` to each basis object in the module object
 * `elem`.
 *
 * @tparam M Type modeling `Module`.
 * @tparam F Function object type; must be convertible to `LinearMap<M>`.
 * @param elem Input module object.
 * @param func Linear map to apply to `elem`.
 * @return M A new module object that is the result of `func` applied to `elem`.
 *
 * @sa `LinearMap`
 */
template <Module M, typename F>
requires std::convertible_to<F, LinearMap<M>>
[[nodiscard]] M linear_apply(const M& elem, const F& func) {
  M result;
  for (const typename M::BasisType& cell : elem) {
    result += elem[cell] * func(cell);
  }
  return result;
}

}  // namespace chomp::core

#endif  // CHOMP_ALGEBRA_ALGEBRA_H
