/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This file contains custom iterators used throughout CHomP.
 */

#ifndef CHOMP_ITERATORS_H
#define CHOMP_ITERATORS_H

#include <concepts>
#include <cstddef>
#include <iterator>
#include <utility>

namespace chomp::util {

/**
 * @brief Constant iterator wrapper for associative containers whose iterators
 * are key, value pairs.
 *
 * This forward iterator only fetches `const`-qualified references to keys from
 * an associative container when dereferenced rather than key, value pairs. This
 * is used for `chomp::modules::UnorderedMapModule` and
 * `chomp::modules::MapModule` class templates to maintain a consistent
 * interface with the set-based module class implementations.
 *
 * @tparam I A forward iterator type returning pairs (or tuples) when
 * dereferenced.
 */
template <typename I>
requires std::forward_iterator<I> && requires(I it) {
    { std::get<0>(*it) } -> std::convertible_to<typename std::tuple_element_t<
            0, typename std::iterator_traits<I>::value_type>>;
}
class KeyIterator {

    I it;

public:
    /** @brief Difference type between iterators. */
    using difference_type = typename std::iterator_traits<I>::difference_type;
    /** @brief Value type when dereferenced. */
    using value_type = const std::tuple_element_t<
        0, typename std::iterator_traits<I>::value_type>;
    /** @brief Pointer type. */
    using pointer = value_type*;
    /** @brief Reference type. */
    using reference = value_type&;
    /** @brief Tag for `iterator_traits` */
    using iterator_concept = std::forward_iterator_tag;

    /**
     * @brief Default initialize a new KeyIterator object.
     *
     * Unusable in this state but can be assigned to, as normal.
     */
    KeyIterator() = default;
    /**
     * @brief Construct a new KeyIterator object using a passed iterator.
     *
     * @param it
     */
    KeyIterator(I it) : it(it) {};

    /**
     * @brief Dereferencing this iterator yields a constant reference to
     * the key in the map.
     *
     * @return reference
     */
    reference operator*() const {
        return std::get<0>(*it);
    }

    /**
     * @brief Equality operates on wrapped pointers.
     *
     * @param rhs
     * @return true
     * @return false
     */
    bool operator==(const KeyIterator& rhs) const {
        return it == rhs.it;
    }

    /**
     * @brief Preincrememnt operates on wrapped iterator.
     *
     * @return KeyIterator&
     */
    KeyIterator& operator++() {
        ++it;
        return *this;
    }
    /**
     * @brief Postincrement operates on wrapped iterator.
     *
     */
    KeyIterator operator++(int) {
        KeyIterator temp = *this;
        ++it;
        return temp;
    }
};

} // namespace chomp::util

#endif // CHOMP_ITERATORS_H
