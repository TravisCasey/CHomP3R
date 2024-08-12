/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This file handles the requirements on and definitions of the
 * coefficient rings.
 *
 */

#ifndef RINGS_H
#define RINGS_H

#include <concepts>

namespace chomp::modules {

/**
 * @brief The minimal requirements for a class to be considered a Ring.
 *
 * We require certain arithmetic operations, equality operators, and casting of
 * 0 and 1 to the additive and multiplicative identities, respectively.
 *
 * Note that the standard numerical types fulfill this requirement.
 *
 * @tparam R The proposed class implementing a ring.
 */
template <typename R>
concept Ring = requires(R a, R b) {
    static_cast<R>(0);
    static_cast<R>(1);
    { -a } -> std::same_as<R>;
    { a + b } -> std::same_as<R>;
    { a - b } -> std::same_as<R>;
    { a* b } -> std::same_as<R>;
    { a == b } -> std::same_as<bool>;
    { a != b } -> std::same_as<bool>;
    a += b;
    a -= b;
    a *= b;
};

/**
 * @brief Return the additive identity of the ring R.
 *
 * @tparam R The ring type.
 * @return R
 */
template <Ring R>
R zero() {
    return static_cast<R>(0);
}

/**
 * @brief Return the multiplicative identity of the ring R.
 *
 * @tparam R The ring type.
 * @return R
 */
template <Ring R>
R one() {
    return static_cast<R>(1);
}

/**
 * @brief The ring of integers modulo p. Defines a field if and only if p is
 * prime.
 *
 * Defines arithmetic operators in a straightforward way along with a conversion
 * to and from int. Internally, the value is stored as a unsigned int, and is
 * only taken modulo p for the purposes of equality tests and conversion to int.
 *
 * @tparam p Integer p > 1.
 */
template <unsigned int p>
class Z {
public:
    /**
     * @brief Construct a new Z object with value n modulo p.
     *
     * @param n
     */
    constexpr Z(int n) noexcept : value(n) {}

    /**
     * @brief (Explicit) conversion to int using value of n modulo p.
     *
     * Output range is [0, p-1).
     *
     * @return int
     */
    [[nodiscard]] explicit operator int() const noexcept {
        return mod(value);
    }

    /**
     * @brief Negative modulo p.
     *
     * @return Z
     */
    [[nodiscard]] constexpr Z operator-() const noexcept {
        return Z(-value);
    }
    /**
     * @brief Sum modulo p.
     *
     * @param n
     * @return Z
     */
    [[nodiscard]] constexpr Z operator+(const Z& n) const noexcept {
        return Z(value + n.value);
    }
    /**
     * @brief Difference modulo p.
     *
     * @param n
     * @return Z
     */
    [[nodiscard]] constexpr Z operator-(const Z& n) const noexcept {
        return Z(value - n.value);
    }
    /**
     * @brief Multiplication modulo p.
     *
     * @param n
     * @return Z
     */
    [[nodiscard]] constexpr Z operator*(const Z& n) const noexcept {
        return Z(value * n.value);
    }
    /**
     * @brief Equality modulo p.
     *
     * @param n
     * @return true If the Z instances are equivalent.
     * @return false Otherwise.
     */
    [[nodiscard]] constexpr bool operator==(const Z& n) const noexcept {
        return (mod(value - n.value) == 0);
    }
    /**
     * @brief Inequality modulo p.
     *
     * @param n
     * @return true If the Z instances are not equivalent.
     * @return false Otherwise.
     */
    [[nodiscard]] constexpr bool operator!=(const Z& n) const noexcept {
        return (mod(value - n.value) != 0);
    }
    /**
     * @brief Updating sum modulo p.
     *
     * @param n
     * @return Z&
     */
    constexpr Z& operator+=(const Z& n) noexcept {
        value += n.value;
        return *this;
    }
    /**
     * @brief Updating subtraction modulo p.
     *
     * @param n
     * @return Z&
     */
    constexpr Z& operator-=(const Z& n) noexcept {
        value -= n.value;
        return *this;
    }
    /**
     * @brief Updating multiplication modulo p.
     *
     * @param n
     * @return Z&
     */
    constexpr Z& operator*=(const Z& n) noexcept {
        value *= n.value;
        return *this;
    }

private:
    unsigned int value;
    [[nodiscard]] constexpr unsigned int mod(int n) const noexcept;
};

template <unsigned int p>
inline constexpr unsigned int Z<p>::mod(int n) const noexcept {
    int modulus = n % static_cast<int>(p);
    return modulus >= 0 ? modulus : modulus + p;
}

template <>
inline constexpr unsigned int Z<2u>::mod(int n) const noexcept {
    return (n & 1) == 1;
}

} // namespace chomp::modules

#endif // RINGS_H
