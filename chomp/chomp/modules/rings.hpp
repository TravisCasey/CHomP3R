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
#include <limits>

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
    { -a } -> std::convertible_to<R>;
    { a + b } -> std::convertible_to<R>;
    { a - b } -> std::convertible_to<R>;
    { a* b } -> std::convertible_to<R>;
    { a == b } -> std::convertible_to<bool>;
    { a != b } -> std::convertible_to<bool>;
    a += b;
    a -= b;
    a *= b;
};

/**
 * @brief Classes implementing this concept represent the ring (field) with two
 * elements.
 *
 * While there is only one such ring, there may be different @c Z data
 * structures implementing it due to the @c T template argument.
 *
 * @tparam R The proposed class implementing a two-element ring.
 */
template <typename R>
concept BinaryRing = Ring<R> && R::get_p()
== 2;

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
 * @brief The ring of (unsigned) integers (of type @c T) modulo @c p.
 *
 * Defines standard arithmetic and equality operators along with an accessor
 * method for the modular value. Overflow is prevented by requiring that the
 * the modulo value @c p satisfies @c p squared is less than the maximum value
 * that can be represented by @c T.
 *
 * @tparam T Unsigned integral type
 * @tparam p Modulo value; required p > 2 and p * p < max value of T.
 */
template <typename T, T p>
requires(std::unsigned_integral<T>&& p > static_cast<T>(1) &&
         p < std::numeric_limits<T>::max() / p) class Z {
public:
    /**
     * @brief Underlying unsigned integral type.
     *
     */
    using value_type = T;

    /**
     * @brief Construct a new Z object with value n modulo p.
     *
     * @param n
     */
    constexpr Z(T n) : value(n % p) {}

    /**
     * @brief Return value in range [0, p-1]
     *
     * @return T
     */
    [[nodiscard]] constexpr T get_value() const noexcept {
        return value;
    }
    /**
     * @brief Get the modulo value p.
     *
     * @return T
     */
    [[nodiscard]] static constexpr T get_p() noexcept {
        return p;
    }

    /**
     * @brief Negative operator modulo p.
     *
     * @return Z
     */
    [[nodiscard]] constexpr Z operator-() const noexcept {
        return Z(p - value); // safe as value in [0, p-1]
    }
    /**
     * @brief Sum operator modulo p.
     *
     * @param rhs
     * @return Z
     */
    [[nodiscard]] constexpr Z operator+(const Z& rhs) const noexcept {
        T sum = value + rhs.value; // safe as p < sqrt(max value of T)
        return Z(sum >= p ? sum - p : sum);
    }
    /**
     * @brief Difference operator modulo p.
     *
     * @param rhs
     * @return Z
     */
    [[nodiscard]] constexpr Z operator-(const Z& rhs) const noexcept {
        T sum = value + (p - rhs.value); // safe; see - and + operators
        return Z(sum >= p ? sum - p : sum);
    }
    /**
     * @brief Multiplication operator modulo p.
     *
     * @param rhs
     * @return Z
     */
    [[nodiscard]] constexpr Z operator*(const Z& rhs) const noexcept {
        T prod = value * rhs.value; // safe as p < sqrt(max value of T)
        return Z(prod);
    }

    /**
     * @brief Equality operator modulo p.
     *
     * @param rhs
     * @return true If the Z instances are equivalent modulo p.
     * @return false Otherwise.
     */
    [[nodiscard]] constexpr bool operator==(const Z& rhs) const noexcept {
        T sum = value + (p - rhs.value); // safe; see - and + operators;
        return sum == p;
    }
    /**
     * @brief Inequality operator modulo p.
     *
     * @param rhs
     * @return true If the Z instances are not equivalent modulo p.
     * @return false Otherwise.
     */
    [[nodiscard]] constexpr bool operator!=(const Z& rhs) const noexcept {
        T sum = value + (p - rhs.value); // safe; see - and + operators;
        return sum != p;
    }

    /**
     * @brief Updating sum operator modulo p.
     *
     * @param rhs
     * @return Z&
     */
    constexpr Z& operator+=(const Z& rhs) noexcept {
        value += rhs.value; // safe; see + operator
        if (value >= p) {
            value -= p;
        }
        return *this;
    }
    /**
     * @brief Updating subtraction operator modulo p.
     *
     * @param rhs
     * @return Z&
     */
    constexpr Z& operator-=(const Z& rhs) noexcept {
        value += (p - rhs.value); // safe; see + and - operators
        if (value >= p) {
            value -= p;
        }
        return *this;
    }
    /**
     * @brief Updating multiplication operator modulo p.
     *
     * @param rhs
     * @return Z&
     */
    constexpr Z& operator*=(const Z& rhs) noexcept {
        value = (value * rhs.value) % p; // safe; see * operator
        return *this;
    }

private:
    T value;
};

} // namespace chomp::modules

#endif // RINGS_H
