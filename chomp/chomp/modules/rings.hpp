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
 * to and from int.
 *
 * @tparam p Integer p > 1.
 */
template <int p>
class Z {
public:
    /**
     * @brief Construct a new Z object with value n modulo p.
     *
     * @param n
     */
    Z(int n) : value(mod(n)) {}

    /**
     * @brief (Explicit) conversion to int using value of n modulo p.
     *
     * Output range is [0, p-1).
     *
     * @return int
     */
    explicit operator int() {
        return value;
    }

    /**
     * @brief Negative modulo p.
     *
     * @return Z
     */
    Z operator-() const {
        return Z(-value);
    }
    /**
     * @brief Sum modulo p.
     *
     * @param n
     * @return Z
     */
    Z operator+(const Z& n) const {
        return Z(value + n.value);
    }
    /**
     * @brief Difference modulo p.
     *
     * @param n
     * @return Z
     */
    Z operator-(const Z& n) const {
        return Z(value - n.value);
    }
    /**
     * @brief Multiplication modulo p.
     *
     * @param n
     * @return Z
     */
    Z operator*(const Z& n) const {
        return Z(value * n.value);
    }
    /**
     * @brief Equality modulo p.
     *
     * @param n
     * @return true If the Z instances are equivalent.
     * @return false Otherwise.
     */
    bool operator==(const Z& n) const {
        return value == n.value;
    }
    /**
     * @brief Inequality modulo p.
     *
     * @param n
     * @return true If the Z instances are not equivalent.
     * @return false Otherwise.
     */
    bool operator!=(const Z& n) const {
        return value != n.value;
    }
    /**
     * @brief Updating sum modulo p.
     *
     * @param n
     * @return Z&
     */
    Z& operator+=(const Z& n) {
        value = mod(value + n.value);
        return *this;
    }
    /**
     * @brief Updating subtraction modulo p.
     *
     * @param n
     * @return Z&
     */
    Z& operator-=(const Z& n) {
        value = mod(value - n.value);
        return *this;
    }
    /**
     * @brief Updating multiplication modulo p.
     *
     * @param n
     * @return Z&
     */
    Z& operator*=(const Z& n) {
        value = mod(value * n.value);
        return *this;
    }

private:
    int value;
    int mod(int n) {
        if (n >= 0) {
            return n - (n / p) * p;
        }
        return n - (n / p) * p + p;
    }
};

} // namespace chomp::modules

#endif // RINGS_H
