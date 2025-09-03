/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header defines the cyclic integral ring `Z` class template.
 */

#ifndef CHOMP_CORE_ALGEBRA_CYCLIC_HPP
#define CHOMP_CORE_ALGEBRA_CYCLIC_HPP

#include <chomp/core/algebra/algebra.hpp>

#include <flint/flint.h>
#include <flint/nmod.h>

#include <concepts>
#include <stdexcept>

namespace chomp::core {

/**
 * @brief The cyclic ring of integers modulo a prime `p`.
 *
 * Instantiations of this class template are representatives of the equivalence
 * classes modulo `p`, with representative in `[0, p-1]`.
 *
 * Defines arithmetic and equality operators along with an accessor method for
 * the modulus value in the range `[0, p-1]`. Overflow is prevented by requiring
 * that the divisor value `p` satisfies `p * p <= max(int)`, the maximum value
 * that can be represented by `int`.
 *
 * @tparam p Divisor value; required p > 1 and `p * p <= max(int)`.
 */
template <unsigned int p>
requires requires { p > 1; }
class Z {
  nmod_t modulus;
  unsigned long value = 0;

  explicit Z(nmod_t modulus, unsigned long value) :
      modulus(modulus), value(value) {}

public:
  /**
   * @brief Construct a new `Z` object with value `n` modulo `p`.
   *
   * @param n
   */
  explicit Z(int n = 0) {
    nmod_init(&modulus, static_cast<unsigned long>(p));
    value += nmod_set_si(n, modulus);
  }

  /**
   * @brief Get the equivalence class representative (i.e. the value modulo
   * `p`) in range `[0, p-1]`.
   *
   * @return int
   */
  [[nodiscard]] int rep() const noexcept {
    return static_cast<int>(value);
  }
  /**
   * @brief Get the divisor `p`.
   *
   * @return int
   */
  [[nodiscard]] static uint divisor() noexcept {
    return p;
  }

  /**
   * @brief Negative operator modulo `p`.
   *
   * @return Z
   */
  [[nodiscard]] Z operator-() const noexcept {
    return Z(modulus, nmod_neg(value, modulus));
  }
  /**
   * @brief Sum operator modulo `p`.
   *
   * @param rhs
   * @return Z
   */
  [[nodiscard]] Z operator+(const Z& rhs) const noexcept {
    return Z(modulus, nmod_add(value, rhs.value, modulus));
  }
  /**
   * @brief Difference operator modulo `p`.
   *
   * @param rhs
   * @return Z
   */
  [[nodiscard]] Z operator-(const Z& rhs) const noexcept {
    return Z(modulus, nmod_sub(value, rhs.value, modulus));
  }
  /**
   * @brief Product operator modulo `p`.
   *
   * @param rhs
   * @return Z
   */
  [[nodiscard]] Z operator*(const Z& rhs) const noexcept {
    return Z(modulus, nmod_mul(value, rhs.value, modulus));
  }

  /**
   * @brief Equality operator modulo `p`.
   *
   * @param rhs
   * @return true If the `Z` instances are equivalent modulo `p`.
   * @return false Otherwise.
   */
  [[nodiscard]] bool operator==(const Z& rhs) const noexcept {
    return value == rhs.value;
  }

  /**
   * @brief Compound assignment sum operator modulo `p`.
   *
   * @param rhs
   * @return Z&
   */
  Z& operator+=(const Z& rhs) noexcept {
    value = nmod_add(value, rhs.value, modulus);
    return *this;
  }
  /**
   * @brief Compound assignment difference operator modulo `p`.
   *
   * @param rhs
   * @return Z&
   */
  Z& operator-=(const Z& rhs) noexcept {
    value = nmod_sub(value, rhs.value, modulus);
    return *this;
  }
  /**
   * @brief Compound assignment product operator modulo `p`.
   *
   * @param rhs
   * @return Z&
   */
  Z& operator*=(const Z& rhs) noexcept {
    value = nmod_mul(value, rhs.value, modulus);
    return *this;
  }

  /**
   * @brief Invert the element in `Z`. Throws `std::domain_error` if the element
   * is zero.
   *
   * @return Z
   */
  [[nodiscard]] Z invert() const {
    if (value == 0) {
      throw std::domain_error(
          "Attempted to invert a ring element that is not a unit.");
    }
    return Z(modulus, nmod_inv(value, modulus));
  }
};

/** @overload */
template <typename R>
requires Ring<R> && requires(const R& r) {
  { r.invert() } -> std::convertible_to<R>;
}
[[nodiscard]] R invert(const R& value) {
  return value.invert();
}

}  // namespace chomp::core

#endif  // CHOMP_CORE_ALGEBRA_CYCLIC_HPP
