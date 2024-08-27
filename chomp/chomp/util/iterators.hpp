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

#include <cstddef>
#include <iterator>
#include <utility>

namespace chomp::util {

/**
 * @brief Constant iterator wrapper for `std::map` and `std::unordered_map`.
 *
 * This forward iterator only fetches `const`-qualified references to keys from
 * an associative container when dereferenced rather than key, value pairs. This
 * is used for `chomp::modules::UnorderedMapModule` and
 * `chomp::modules::MapModule` class templates to maintain a consistent
 * interface with the set-based module class implementations.
 *
 * @tparam Container A specialization of `std::map` or `std::unordered_map`
 * class templates.
 */
template <typename Container>
class KeyIterator {

    typename Container::const_iterator it;

public:
    /** @brief Difference type between iterators. */
    using difference_type = std::ptrdiff_t;
    /** @brief Value type when dereferenced. */
    using value_type = const typename Container::key_type;
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
     * @brief Construct a new KeyIterator object using a constant
     * iterator.
     *
     * @param it
     */
    KeyIterator(typename Container::const_iterator it) : it(it) {};

    /**
     * @brief Dereferencing this iterator yields a constant reference to
     * the key in the map.
     *
     * @return value_type&
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
