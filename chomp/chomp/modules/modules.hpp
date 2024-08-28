/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This file contains various classes implementing elements of a free
 * module.
 */

#ifndef CHOMP_MODULES_H
#define CHOMP_MODULES_H

#include <concepts>
#include <functional>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <chomp/modules/concepts.hpp>
#include <chomp/modules/rings.hpp>
#include <chomp/util/iterators.hpp>

namespace chomp::modules {

using chomp::util::KeyIterator;

/**
 * @brief Abstract base class implementing free `R`-modules on basis set  `T`.
 *
 * Concrete classes deriving from this base classes have arithmetic and other
 * operations defined automatically. This class is not intended for
 * polymorphism; instead this is for ease of maintenance and to avoid repetition
 * of code among various similar classes.
 *
 * @tparam T Cell type that forms the basis for the free `R`-module.
 * @tparam R Ring coefficient type.
 * @tparam I Forward input iterator type over cells; implementation-dependent.
 */
template <CellType T, Ring R, std::forward_iterator I>
class AbstractModule {
public:
    /** @brief Cell iterator type */
    using cell_iter_t = I;
    /** @brief Cell type. */
    using cell_t = T;
    /** @brief Ring type. */
    using ring_t = R;

    /**
     * @brief Default Constructor.
     *
     * Initializes an empty `R`-linear combination which is the identity of the
     * module.
     *
     */
    AbstractModule() = default;
    /**
     * @brief Only used for `static_cast` of 0 to module types, satisfying the
     * `Group` concept.
     *
     * The additive identity is the default initialization. Values of `n` other
     * than zero throw `std::domain_error`.
     *
     * @param n
     */
    constexpr AbstractModule(int n) {
        if (n != 0) {
            throw std::domain_error(
                "Cannot construct Module with nonzero value.");
        }
    }

    /**
     * @brief Get the `R` coefficient of `cell` in the module element.
     *
     * A value of zero (in `R`) indicates that `cell` is not present in the
     * element.
     *
     * @param cell
     * @return R
     */
    virtual R operator[](const T& cell) const = 0;

    /**
     * @brief Constant input iterator to the beginning of the cells in this
     * element.
     *
     * @return cell_iter_t
     */
    virtual cell_iter_t cbegin() const = 0;
    /**
     * @brief Constant input iterator to the end of the cells in this element.
     *
     * @return cell_iter_t
     */
    virtual cell_iter_t cend() const = 0;

    /**
     * @brief Insert a cell, coefficient pair into the element.
     *
     * If `cell` is already present in the element, then `coef` is added to
     * the coefficient of `cell` in the element.
     *
     * @param cell
     * @param coef
     */
    virtual void insert(const T& cell, const R& coef) = 0;
    /** @overload */
    virtual void insert(T&& cell, const R& coef) = 0;

    /**
     * @brief Resets the element to the default initalization state.
     *
     */
    virtual void clear() = 0;
};

/**
 * @brief The requirements for a class `M` to model a free `R`-module on the
 * basis `T`.
 *
 * Namely, the requirements are:
 *   - Derived from AbstractModule base class with correct member types
 *   - Copy constructible
 *   - Default initializable
 *   - Equality comparable
 *
 * @sa AbstractModule
 *
 * @tparam T Cell type.
 * @tparam R Ring type.
 * @tparam M `R`-Module type with basis `T`.
 */
template <typename M>
concept Module = std::copy_constructible<M> && std::default_initializable<M> &&
    std::equality_comparable<M> &&
    std::derived_from<M, AbstractModule<typename M::cell_t, typename M::ring_t,
                                        typename M::cell_iter_t>> &&
    std::same_as<typename M::lfunc_t,
                 std::function<M(const typename M::cell_t&)>>;

/**
 * @brief Returns a (constant) iterator to the beginning of the cells in the
 * module `m`.
 *
 * This allows for range-based for loops over the cells in `m`, but does not
 * allow modification of the keys in `m`.
 *
 * @tparam M Module type.
 * @param m Module element of type `M`.
 * @return M::cell_iter_t Constant iterator to the beginning of the cells in
 * `m`.
 *
 * @sa end()
 */
template <Module M>
typename M::cell_iter_t begin(const M& m) {
    return m.cbegin();
}
/**
 * @brief Returns a (constant) iterator to the end of the cells in the module
 * `m`.
 *
 * This allows for range-based for loops over the cells in `m`, but does not
 * allow modification of the keys in `m`.
 *
 * @tparam M Module type.
 * @param m Module element of type `M`.
 * @return M::cell_iter_t Constant iterator to the end of the cells in `m`.
 *
 * @sa begin()
 */
template <Module M>
typename M::cell_iter_t end(const M& m) {
    return m.cend();
}

/**
 * @brief Apply a function `func` linearly along the cells in this element.
 *
 * @tparam M `Module` type.
 * @param elem Input module element to linear function `func`.
 * @param func Linear function object.
 * @return M A new module element that is the output of the function `func`
 * applied to `elem`.
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M linear_apply(const M& elem, const typename M::lfunc_t& func) {
    M result;
    // Apply `func` to each cell in `elem` and multiple the result by its
    // coefficient in `elem`.
    for (const typename M::cell_t& cell : elem) {
        result += elem[cell] * func(cell);
    }
    return result;
}

/**
 * @brief Updating sum operator computes the formal sum of `lhs` and `rhs` in
 * the module `M`.
 *
 * @tparam M `Module` type.
 * @param lhs
 * @param rhs
 * @return M&
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M& operator+=(M& lhs, const M& rhs) {
    for (const typename M::cell_t cell : rhs) {
        lhs.insert(cell, rhs[cell]);
    }
    return lhs;
}

/**
 * @brief Updating difference operator computes the formal difference of `lhs`
 * and `rhs` in the module `M`.
 *
 * @tparam M `Module` type.
 * @param lhs
 * @param rhs
 * @return M&
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M& operator-=(M& lhs, const M& rhs) {
    for (const typename M::cell_t cell : rhs) {
        lhs.insert(cell, -rhs[cell]);
    }
    return lhs;
}
/**
 * @copydoc operator-=()
 * @relatedalso AbstractModule
 */
template <Module M>
M& operator-=(M& lhs, M&& rhs) {
    lhs += -std::move(rhs);
    return lhs;
}

/**
 * @brief Updating scalar multiplication operator computes the scalar product of
 * `lhs` and scalar `rhs` in the module `M`.
 *
 * @tparam M `Module` type.
 * @param lhs
 * @param rhs
 * @return M&
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M& operator*=(M& lhs, const typename M::ring_t& rhs) {
    // This never adds more cells to `lhs`, but the erase call in the insert
    // methods will invalidate pointers. However, presuming the coefficient ring
    // is at least an integral domain, this can only happen if `rhs` is 0.
    // So we handle that case.
    if (rhs == zero<typename M::ring_t>()) {
        lhs.clear();
        return lhs;
    }
    for (const typename M::cell_t cell : lhs) {
        lhs.insert(cell, (rhs - one<typename M::ring_t>()) * lhs[cell]);
    }
    return lhs;
}

/**
 * @brief Negates each coefficient in the module element.
 *
 * @tparam M `Module` type.
 * @param elem
 * @return M
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M operator-(const M& elem) {
    M result(elem);
    result *= -one<typename M::ring_t>();
    return result;
}
/**
 * @brief Negates each coefficient in the module element.
 *
 * @tparam M `Module` type.
 * @param elem
 * @return M
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M operator-(M&& elem) {
    M result(std::move(elem));
    result *= -one<typename M::ring_t>();
    return result;
}

/**
 * @brief Computes formal sum of `lhs` and `rhs` in the module `M`.
 *
 * @tparam M `Module` type.
 * @param lhs
 * @param rhs
 * @return M
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M operator+(const M& lhs, const M& rhs) {
    M result(lhs);
    result += rhs;
    return result;
}
/**
 * @copydoc operator+()
 * @relatedalso AbstractModule
 */
template <Module M>
M operator+(const M& lhs, M&& rhs) {
    M result(lhs);
    result += std::move(rhs);
    return result;
}
/**
 * @copydoc operator+()
 * @relatedalso AbstractModule
 */
template <Module M>
M operator+(M&& lhs, const M& rhs) {
    M result(std::move(lhs));
    result += rhs;
    return result;
}
/**
 * @copydoc operator+()
 * @relatedalso AbstractModule
 */
template <Module M>
M operator+(M&& lhs, M&& rhs) {
    M result(std::move(lhs));
    result += std::move(rhs);
    return result;
}

/**
 * @brief Computes formal difference of `lhs` and `rhs` in the module `M`.
 *
 * @tparam M `Module` type.
 * @param lhs
 * @param rhs
 * @return M
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M operator-(const M& lhs, const M& rhs) {
    M result(lhs);
    result -= rhs;
    return result;
}
/**
 * @brief Computes formal difference of `lhs` and `rhs` in the module `M`.
 *
 * @tparam M `Module` type.
 * @param lhs
 * @param rhs
 * @return M
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M operator-(const M& lhs, M&& rhs) {
    M result(lhs);
    result -= std::move(rhs);
    return result;
}
/**
 * @brief Computes formal difference of `lhs` and `rhs` in the module `M`.
 *
 * @tparam M `Module` type.
 * @param lhs
 * @param rhs
 * @return M
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M operator-(M&& lhs, const M& rhs) {
    M result(std::move(lhs));
    result -= rhs;
    return result;
}
/**
 * @brief Computes formal difference of `lhs` and `rhs` in the module `M`.
 *
 * @tparam M `Module` type.
 * @param lhs
 * @param rhs
 * @return M
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M operator-(M&& lhs, M&& rhs) {
    M result(std::move(lhs));
    result -= std::move(rhs);
    return result;
}

/**
 * @brief Computes scalar product of `lhs` with scalar `rhs` in the module `M`.
 *
 * @tparam M `Module type`.
 * @param lhs
 * @param rhs
 * @return M
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M operator*(const M& lhs, const typename M::ring_t& rhs) {
    M result(lhs);
    result *= rhs;
    return result;
}
/**
 * @brief Computes scalar product of `lhs` with scalar `rhs` in the module `M`.
 *
 * @tparam M `Module type`.
 * @param lhs
 * @param rhs
 * @return M
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M operator*(M&& lhs, const typename M::ring_t& rhs) {
    M result(std::move(lhs));
    result *= rhs;
    return result;
}
/**
 * @brief Computes scalar product of scalar `lhs` with `rhs` in the module `M`.
 *
 * @tparam M `Module` type.
 * @param lhs
 * @param rhs
 * @return M
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M operator*(const typename M::ring_t& lhs, const M& rhs) {
    M result(rhs);
    result *= lhs;
    return result;
}
/**
 * @brief Computes scalar product of scalar `lhs` with `rhs` in the module `M`.
 *
 * @tparam M `Module` type.
 * @param lhs
 * @param rhs
 * @return M
 *
 * @relatedalso AbstractModule
 */
template <Module M>
M operator*(const typename M::ring_t& lhs, M&& rhs) {
    M result(std::move(rhs));
    result *= lhs;
    return result;
}

/**
 * @brief This class template implements a free `R`-module on basis set `T`. Its
 * instantiations are elemnts of this `R`-module, which are formal `R`-linear
 * combinations of elements of `T`.
 *
 * The implementation uses `std::unordered_set` to store the cells. The
 * coefficient ring `R` is required to be binary-valued (`BinaryRing`). These
 * coefficients are not explicitly stored in this implementation and are instead
 * based on inclusion in the set.
 *
 * @tparam T Cell type modeling `Hashable` concept.
 * @tparam R Ring type modeling `BinaryRing` concept.
 */
template <Hashable T, BinaryRing R>
class UnorderedSetModule
    : public AbstractModule<T, R,
                            typename std::unordered_set<T>::const_iterator> {
    std::unordered_set<T> cells;
    using iter_t = typename std::unordered_set<T>::iterator;
    using citer_t = typename std::unordered_set<T>::const_iterator;
    using node_t = typename std::unordered_set<T>::node_type;

public:
    /** @copydoc AbstractModule::cell_iter_t */
    using cell_iter_t = citer_t;
    /** @copydoc AbstractModule::cell_t */
    using cell_t = typename AbstractModule<T, R, cell_iter_t>::cell_t;
    /** @copydoc AbstractModule::ring_t */
    using ring_t = typename AbstractModule<T, R, cell_iter_t>::ring_t;
    /** @brief Linear function type. */
    using lfunc_t = std::function<UnorderedSetModule<T, R>(const T&)>;

    /** @copydoc AbstractModule() */
    UnorderedSetModule() = default;
    /** @copydoc AbstractModule(int) */
    constexpr UnorderedSetModule(int n)
        : AbstractModule<T, R, cell_iter_t>(n) {}

    /** @copydoc AbstractModule::operator[]() */
    R operator[](const T& cell) const {
        return cells.count(cell) == 0 ? zero<R>() : one<R>();
    }

    /** @copydoc AbstractModule::cbegin() */
    cell_iter_t cbegin() const noexcept {
        return cells.cbegin();
    }
    /** @copydoc AbstractModule::cend() */
    cell_iter_t cend() const noexcept {
        return cells.cend();
    }

    /** @copydoc AbstractModule::insert(const T&, const R&) */
    void insert(const T& cell, const R& coef) {
        if (coef != zero<R>()) {
            std::pair<iter_t, bool> ins_result = cells.insert(cell);
            if (!ins_result.second) {
                cells.erase(ins_result.first);
            }
        }
    }
    /** @overload */
    void insert(T&& cell, const R& coef) {
        if (coef != zero<R>()) {
            std::pair<iter_t, bool> ins_result = cells.insert(std::move(cell));
            if (!ins_result.second) {
                cells.erase(ins_result.first);
            }
        }
    }

    /**
     * @brief Updating sum operator computes the formal sum of `this` and `rhs`.
     *
     * @param rhs
     * @return UnorderedSetModule&
     */
    UnorderedSetModule& operator+=(UnorderedSetModule&& rhs) {
        for (citer_t it = rhs.cells.cbegin(); it != rhs.cells.cend();) {
            // Avoid invalidation of iterator with postincrement.
            node_t nh = rhs.cells.extract(it++);
            insert(std::move(nh.value()), one<ring_t>());
        }
        return *this;
    }

    /** @copydoc AbstractModule::clear() */
    void clear() {
        cells.clear();
    }

    /**
     * @brief Equality operator returns `true` if and only if the cells and
     * their coefficients are equal.
     *
     * @param rhs
     * @return bool
     */
    bool operator==(const UnorderedSetModule& rhs) const {
        return cells == rhs.cells;
    }
};

/**
 * @copybrief UnorderedSetModule
 *
 * The implementation uses `std::set` to store the cells. The
 * coefficient ring `R` is required to be binary-valued (`BinaryRing`). These
 * coefficients are not explicitly stored in this implementation and are instead
 * based on inclusion in the set.
 *
 * @tparam T A cell type with `<` implementing a total order.
 * @tparam R A two-element (`BinaryRing`) ring type.
 */
template <Comparable T, BinaryRing R>
class SetModule
    : public AbstractModule<T, R, typename std::set<T>::const_iterator> {
    std::set<T> cells;
    using iter_t = typename std::set<T>::iterator;
    using citer_t = typename std::set<T>::const_iterator;
    using node_t = typename std::set<T>::node_type;

public:
    /** @copydoc AbstractModule::cell_iter_t */
    using cell_iter_t = citer_t;
    /** @copydoc AbstractModule::cell_t */
    using cell_t = typename AbstractModule<T, R, cell_iter_t>::cell_t;
    /** @copydoc AbstractModule::ring_t */
    using ring_t = typename AbstractModule<T, R, cell_iter_t>::ring_t;
    /** @brief Linear function type. */
    using lfunc_t = std::function<SetModule<T, R>(const T&)>;

    /** @copydoc AbstractModule() */
    SetModule() = default;
    /** @copydoc AbstractModule(int) */
    constexpr SetModule(int n) : AbstractModule<T, R, cell_iter_t>(n) {}

    /** @copydoc AbstractModule::operator[]() */
    R operator[](const T& cell) const {
        return cells.count(cell) == 0 ? zero<R>() : one<R>();
    }

    /** @copydoc AbstractModule::cbegin() */
    cell_iter_t cbegin() const noexcept {
        return cells.cbegin();
    }
    /** @copydoc AbstractModule::cend() */
    cell_iter_t cend() const noexcept {
        return cells.cend();
    }

    /** @copydoc AbstractModule::insert(const T&, const R&) */
    void insert(const T& cell, const R& coef) {
        if (coef != zero<R>()) {
            std::pair<iter_t, bool> ins_result = cells.insert(cell);
            if (!ins_result.second) {
                cells.erase(ins_result.first);
            }
        }
    }
    /** @overload */
    void insert(T&& cell, const R& coef) {
        if (coef != zero<R>()) {
            std::pair<iter_t, bool> ins_result = cells.insert(std::move(cell));
            if (!ins_result.second) {
                cells.erase(ins_result.first);
            }
        }
    }

    /**
     * @brief Updating sum operator computes the formal sum of `this` and `rhs`.
     *
     * @param rhs
     * @return SetModule&
     */
    SetModule& operator+=(SetModule&& rhs) {
        for (citer_t it = rhs.cells.cbegin(); it != rhs.cells.cend();) {
            // Avoid invalidation of iterator with postincrement.
            node_t nh = rhs.cells.extract(it++);
            insert(std::move(nh.value()), one<ring_t>());
        }
        return *this;
    }

    /** @copydoc AbstractModule::clear() */
    void clear() {
        cells.clear();
    }

    /** @copydoc UnorderedSetModule::operator==() */
    bool operator==(const SetModule& rhs) const {
        return cells == rhs.cells;
    }
};

/**
 * @copybrief UnorderedSetModule
 *
 * The implementation uses `std::unordered_map` to store the cells (as keys) as
 * well as their corresponding coefficients.
 *
 * @tparam T Cell type modeling `Hashable` concept.
 * @tparam R Ring type.
 */
template <Hashable T, Ring R>
class UnorderedMapModule
    : public AbstractModule<
          T, R,
          KeyIterator<typename std::unordered_map<T, R>::const_iterator>> {

    std::unordered_map<T, R> cells;
    using iter_t = typename std::unordered_map<T, R>::iterator;
    using citer_t = typename std::unordered_map<T, R>::const_iterator;
    using node_t = typename std::unordered_map<T, R>::node_type;

public:
    /** @copydoc AbstractModule::cell_iter_t */
    using cell_iter_t =
        KeyIterator<typename std::unordered_map<T, R>::const_iterator>;
    /** @copydoc AbstractModule::cell_t */
    using cell_t = typename AbstractModule<T, R, cell_iter_t>::cell_t;
    /** @copydoc AbstractModule::ring_t */
    using ring_t = typename AbstractModule<T, R, cell_iter_t>::ring_t;
    /** @brief Linear function type. */
    using lfunc_t = std::function<UnorderedMapModule<T, R>(const T&)>;

    /** @copydoc AbstractModule() */
    UnorderedMapModule() = default;
    /** @copydoc AbstractModule(int) */
    constexpr UnorderedMapModule(int n)
        : AbstractModule<T, R, cell_iter_t>(n) {}

    /** @copydoc AbstractModule::operator[]() */
    R operator[](const T& cell) const noexcept {
        citer_t find_result = cells.find(cell);
        if (find_result == cells.cend()) {
            return zero<R>();
        } else {
            return find_result->second;
        }
    }

    /** @copydoc AbstractModule::cbegin() */
    cell_iter_t cbegin() const noexcept {
        return KeyIterator(cells.cbegin());
    }
    /** @copydoc AbstractModule::cend() */
    cell_iter_t cend() const noexcept {
        return KeyIterator(cells.cend());
    }

    /** @copydoc AbstractModule::insert(const T&, const R&) */
    void insert(const T& cell, const R& coef) {
        if (coef != zero<R>()) {
            std::pair<iter_t, bool> ins_result =
                cells.insert(std::make_pair(cell, coef));
            if (!ins_result.second) {
                (ins_result.first)->second += coef;
                if ((ins_result.first)->second == zero<R>()) {
                    cells.erase(ins_result.first);
                }
            }
        }
    }
    /** @overload */
    void insert(T&& cell, const R& coef) {
        if (coef != zero<R>()) {
            std::pair<iter_t, bool> ins_result =
                cells.insert(std::make_pair(std::move(cell), coef));
            if (!ins_result.second) {
                (ins_result.first)->second += coef;
                if ((ins_result.first)->second == zero<R>()) {
                    cells.erase(ins_result.first);
                }
            }
        }
    }

    /**
     * @brief Updating sum operator computes the formal sum of `this` and `rhs`.
     *
     * @param rhs
     * @return UnorderedMapModule&
     */
    UnorderedMapModule& operator+=(UnorderedMapModule&& rhs) {
        for (citer_t it = rhs.cells.cbegin(); it != rhs.cells.cend();) {
            // Avoid invalidation of iterator with postincrement.
            node_t nh = rhs.cells.extract(it++);
            insert(std::move(nh.key()), std::move(nh.mapped()));
        }
        return *this;
    }

    /** @copydoc AbstractModule::clear() */
    void clear() {
        cells.clear();
    }

    /** @copydoc UnorderedSetModule::operator==() */
    bool operator==(const UnorderedMapModule& rhs) const {
        return cells == rhs.cells;
    }
};

/**
 * @copybrief UnorderedSetModule
 *
 * The implementation uses `std::map` to store the cells (as keys) as
 * well as their corresponding coefficients.
 *
 * @tparam T A cell type with `<` implementing a total order.
 * @tparam R Ring type.
 */
template <Comparable T, Ring R>
class MapModule
    : public AbstractModule<
          T, R, KeyIterator<typename std::map<T, R>::const_iterator>> {

    std::map<T, R> cells;
    using iter_t = typename std::map<T, R>::iterator;
    using citer_t = typename std::map<T, R>::const_iterator;
    using node_t = typename std::map<T, R>::node_type;

public:
    /** @copydoc AbstractModule::cell_iter_t */
    using cell_iter_t = KeyIterator<typename std::map<T, R>::const_iterator>;
    /** @copydoc AbstractModule::cell_t */
    using cell_t = typename AbstractModule<T, R, cell_iter_t>::cell_t;
    /** @copydoc AbstractModule::ring_t */
    using ring_t = typename AbstractModule<T, R, cell_iter_t>::ring_t;
    /** @brief Linear function type. */
    using lfunc_t = std::function<MapModule<T, R>(const T&)>;

    /** @copydoc AbstractModule() */
    MapModule() = default;
    /** @copydoc AbstractModule(int) */
    constexpr MapModule(int n) : AbstractModule<T, R, cell_iter_t>(n) {}

    /** @copydoc AbstractModule::operator[]() */
    R operator[](const T& cell) const noexcept {
        citer_t find_result = cells.find(cell);
        if (find_result == cells.cend()) {
            return zero<R>();
        } else {
            return find_result->second;
        }
    }

    /** @copydoc AbstractModule::cbegin() */
    cell_iter_t cbegin() const noexcept {
        return KeyIterator(cells.cbegin());
    }
    /** @copydoc AbstractModule::cend() */
    cell_iter_t cend() const noexcept {
        return KeyIterator(cells.cend());
    }

    /** @copydoc AbstractModule::insert(const T&, const R&) */
    void insert(const T& cell, const R& coef) {
        if (coef != zero<R>()) {
            std::pair<iter_t, bool> ins_result =
                cells.insert(std::make_pair(cell, coef));
            if (!ins_result.second) {
                (ins_result.first)->second += coef;
                if ((ins_result.first)->second == zero<R>()) {
                    cells.erase(ins_result.first);
                }
            }
        }
    }
    /** @overload */
    void insert(T&& cell, const R& coef) {
        if (coef != zero<R>()) {
            std::pair<iter_t, bool> ins_result =
                cells.insert(std::make_pair(std::move(cell), coef));
            if (!ins_result.second) {
                (ins_result.first)->second += coef;
                if ((ins_result.first)->second == zero<R>()) {
                    cells.erase(ins_result.first);
                }
            }
        }
    }

    /**
     * @brief Updating sum operator computes the formal sum of `this` and `rhs`.
     *
     * @param rhs
     * @return MapModule&
     */
    MapModule& operator+=(MapModule&& rhs) {
        for (citer_t it = rhs.cells.cbegin(); it != rhs.cells.cend();) {
            // Avoid invalidation of iterator with postincrement.
            node_t nh = rhs.cells.extract(it++);
            insert(std::move(nh.key()), std::move(nh.mapped()));
        }
        return *this;
    }

    /** @copydoc AbstractModule::clear() */
    void clear() {
        cells.clear();
    }

    /** @copydoc UnorderedSetModule::operator==() */
    bool operator==(const MapModule& rhs) const {
        return cells == rhs.cells;
    }
};

// Implementation detail does not need to be documented.
#ifndef DOXYGEN_SHOULD_SKIP_THIS

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

#endif // DOXYGEN_SHOULD_SKIP_THIS

/**
 * @brief Selects best `Module` type at compile time based on concepts
 * fulfilled by `T` and `R`
 *
 * @tparam T Cell type.
 * @tparam R Ring type.
 */
template <CellType T, Ring R>
struct DefaultModule {
    /** @brief Chosen Module type. */
    using type = typename Chooser<Hashable<T>, BinaryRing<R>, T, R>::type;
};

/**
 * @brief Alias for `DefaultModule` chosen type.
 *
 * @tparam T Cell type.
 * @tparam R Ring type.
 */
template <CellType T, Ring R>
using DefaultModule_t = typename DefaultModule<T, R>::type;

} // namespace chomp::modules

#endif // CHOMP_MODULES_H
