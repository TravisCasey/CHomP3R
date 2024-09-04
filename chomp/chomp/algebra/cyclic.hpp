/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header defines the cyclic integral ring `Z` class template.
 */

#ifndef CHOMP_ALGEBRA_CYCLIC_H
#define CHOMP_ALGEBRA_CYCLIC_H

#include <limits>

namespace chomp::core {

/**
 * @brief The cyclic ring of integers modulo `p`.
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
template <int p>
requires requires {
  p > 1;
  p <= std::numeric_limits<int>::max() / p;
}
class Z {
  int value;

public:
  /**
   * @brief Construct a new `Z` object with value `n` modulo `p`.
   *
   * @param n
   */
  constexpr explicit Z(int n = 0) : value(n >= 0 ? n % p : n % p + p) {}
  // % operator truncates towards zero

  /**
   * @brief Get the equivalence class representative (i.e. the value modulo
   * `p`) in range `[0, p-1]`.
   *
   * @return int
   */
  [[nodiscard]] constexpr int rep() const noexcept {
    return value;
  }
  /**
   * @brief Get the divisor `p`.
   *
   * @return T
   */
  [[nodiscard]] static constexpr int divisor() noexcept {
    return p;
  }

  /**
   * @brief Negative operator modulo `p`.
   *
   * @return Z
   */
  [[nodiscard]] constexpr Z operator-() const noexcept {
    return Z(p - value);  // Result in [0, p-1]
  }
  /**
   * @brief Sum operator modulo `p`.
   *
   * @param rhs
   * @return Z
   */
  [[nodiscard]] constexpr Z operator+(const Z& rhs) const noexcept {
    Z result(*this);  // copy constructor avoids conditionals in constructor
    result += rhs;
    return result;
  }
  /**
   * @brief Difference operator modulo `p`.
   *
   * @param rhs
   * @return Z
   */
  [[nodiscard]] constexpr Z operator-(const Z& rhs) const noexcept {
    Z result(*this);
    result -= rhs;
    return result;
  }
  /**
   * @brief Product operator modulo `p`.
   *
   * @param rhs
   * @return Z
   */
  [[nodiscard]] constexpr Z operator*(const Z& rhs) const noexcept {
    Z result(*this);
    result *= rhs;
    return result;
  }

  /**
   * @brief Equality operator modulo `p`.
   *
   * @param rhs
   * @return true If the `Z` instances are equivalent modulo `p`.
   * @return false Otherwise.
   */
  [[nodiscard]] constexpr bool operator==(const Z& rhs) const noexcept {
    return value == rhs.value;
  }

  /**
   * @brief Compound assignment sum operator modulo `p`.
   *
   * @param rhs
   * @return Z&
   */
  constexpr Z& operator+=(const Z& rhs) noexcept {
    value += rhs.value;  // safe from overflow as p * p <= max(int)
    if (value >= p) {
      value -= p;
    }
    return *this;
  }
  /**
   * @brief Compound assignment difference operator modulo `p`.
   *
   * @param rhs
   * @return Z&
   */
  constexpr Z& operator-=(const Z& rhs) noexcept {
    value += (p - rhs.value);
    if (value >= p) {
      value -= p;
    }
    return *this;
  }
  /**
   * @brief Compound assignment product operator modulo `p`.
   *
   * @param rhs
   * @return Z&
   */
  constexpr Z& operator*=(const Z& rhs) noexcept {
    value = (value * rhs.value) % p;  // safe as p * p <= max(int)
    return *this;
  }
};

/**
 * @brief The two-element ring (field).
 *
 * Specialized for the common use case of using the two element ring `Z<2>` as
 * the coefficient ring of choice.
 */
template <>
class Z<2> {
  bool odd;

public:
  /**
   * @brief Construct a new `Z` object with value `n` modulo `2`.
   *
   * @param n
   */
  constexpr explicit Z(int n = 0) : odd(n & 1U) {}
  // Using 1U works for one's complement systems as well

  /**
   * @brief Get the equivalence class representative (i.e. the value modulo `2`)
   * in the range `{0, 1}`.
   *
   * @return int
   */
  [[nodiscard]] constexpr int rep() const noexcept {
    return odd;
  }
  /**
   * @brief Get the divisor, `2`.
   *
   * @return T
   */
  [[nodiscard]] static constexpr int divisor() noexcept {
    return 2;
  }

  /**
   * @brief Negative operator modulo `2` (no effect).
   *
   * @return Z
   */
  [[nodiscard]] constexpr Z operator-() const noexcept {
    return *this;
  }
  /**
   * @brief Sum operator modulo `2`.
   *
   * @param rhs
   * @return Z
   */
  [[nodiscard]] constexpr Z operator+(const Z& rhs) const noexcept {
    Z result(*this);  // copy constructor avoids conditionals in constructor
    result += rhs;
    return result;
  }
  /**
   * @brief Difference operator modulo `2`.
   *
   * @param rhs
   * @return Z
   */
  [[nodiscard]] constexpr Z operator-(const Z& rhs) const noexcept {
    Z result(*this);
    result -= rhs;
    return result;
  }
  /**
   * @brief Product operator modulo `2`.
   *
   * @param rhs
   * @return Z
   */
  [[nodiscard]] constexpr Z operator*(const Z& rhs) const noexcept {
    Z result(*this);
    result *= rhs;
    return result;
  }

  /**
   * @brief Equality operator modulo `2`.
   *
   * @param rhs
   * @return true If the `Z` instances are equivalent modulo `2`.
   * @return false Otherwise.
   */
  [[nodiscard]] constexpr bool operator==(const Z& rhs) const noexcept {
    return odd == rhs.odd;
  }

  /**
   * @brief Compound assignment sum operator modulo `2`.
   *
   * @param rhs
   * @return Z&
   */
  constexpr Z& operator+=(const Z& rhs) noexcept {
    odd = rhs.odd ? !odd : odd;
    return *this;
  }
  /**
   * @brief Compound assignment difference operator modulo `2`.
   *
   * @param rhs
   * @return Z&
   */
  constexpr Z& operator-=(const Z& rhs) noexcept {
    odd = rhs.odd ? !odd : odd;
    return *this;
  }
  /**
   * @brief Compound assignment product operator modulo `2`.
   *
   * @param rhs
   * @return Z&
   */
  constexpr Z& operator*=(const Z& rhs) noexcept {
    odd = rhs.odd ? odd : false;
    return *this;
  }
};

}  // namespace chomp::core

#endif  // CHOMP_ALGEBRA_CYCLIC_H
