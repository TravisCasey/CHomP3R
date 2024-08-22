/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This file contains various classes implementing elements of a free
 * module.
 */

#ifndef MODULES_H
#define MODULES_H

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

#include <chomp/modules/rings.hpp>

namespace chomp::modules {

/**
 * @brief Types `T` modeling this concept have `std::hash<T>` defined
 * as well as the equality operator.
 *
 * Types modeling this concept can be stored in `std::unordered_set` and can be
 * keys in `std::unordered_map`, enabling the use of `UnorderedSetModule` and
 * `UnorderedMapModule` with cell type `T`.
 *
 * @tparam T
 */
template <typename T>
concept Hashable = requires(T a, T b) {
    { std::hash<T>{}(a) } -> std::convertible_to<std::size_t>;
    { a == b } -> std::same_as<bool>;
};

/**
 * @brief Types `T` modeling this concept have `std::less<T>` defined as well as
 * the equality operator.
 *
 * Types modeling this concept can be stored in `std::set` and can be keys in
 * `std::map`, enabling the use of `SetModule` and `MapModule` with cell type
 * `T`.
 *
 * @tparam T
 */
template <typename T>
concept Comparable = requires(T a, T b) {
    { a < b } -> std::same_as<bool>;
    { a == b } -> std::same_as<bool>;
};

/**
 * @brief Types `T` modeling this concept either model `Hashable` or
 * `Comparable`, enabling their use in one of the module types.
 *
 * @tparam T
 */
template <typename T>
concept CellType = Hashable<T> || Comparable<T>;

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
 * @tparam I Input iterator type over cells; implementation-dependent.
 */
template <CellType T, Ring R, std::input_iterator I>
class AbstractModule {
public:
    /** @brief Cell iterator type */
    using cell_iter_t = I;
    /** @brief Cell type. */
    using cell_t = T;
    /** @brief Ring type. */
    using ring_t = R;
    /** @brief Linear function type. */
    using lfunc_t = std::function<std::vector<std::pair<const T, R>>(const T&)>;

    /**
     * @brief Default Constructor.
     *
     * Initializes an empty `R`-linear combination which is the identity of the
     * module.
     *
     */
    AbstractModule() = default;
    /**
     * @brief Only used for `static_cast` of 0 to module types,
     * satisfying the `Group` concept.
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
    virtual cell_iter_t cell_cbegin() const = 0;
    /**
     * @brief Constant input iterator to the end of the cells in this element.
     *
     * @return cell_iter_t
     */
    virtual cell_iter_t cell_cend() const = 0;

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
 * @sa AbstractModule
 *
 * @tparam T Cell type.
 * @tparam R Ring type.
 * @tparam M `R`-Module type with basis `T`.
 */
template <typename M>
concept Module = Ring<typename M::ring_t> && CellType<typename M::cell_t> &&
    std::derived_from<M, AbstractModule<typename M::cell_t, typename M::ring_t,
                                        typename M::cell_iter_t>> &&
    std::same_as<typename M::lfunc_t,
                 typename AbstractModule<typename M::cell_t, typename M::ring_t,
                                         typename M::cell_iter_t>::lfunc_t> &&
    std::copy_constructible<M> && std::default_initializable<M> &&
    std::equality_comparable<M>;

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
    for (typename M::cell_iter_t it = elem.cell_cbegin();
         it != elem.cell_cend(); it++) {
        std::vector<std::pair<const typename M::cell_t, typename M::ring_t>>
            func_result = func(*it);
        for (auto func_it = std::make_move_iterator(func_result.begin());
             func_it != std::make_move_iterator(func_result.end()); func_it++) {
            func_it->second *= elem[*it];
            result.insert(std::get<0>(*func_it), std::get<1>(*func_it));
        }
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
    for (typename M::cell_iter_t it = rhs.cell_cbegin(); it != rhs.cell_cend();
         it++) {
        lhs.insert(*it, rhs[*it]);
    }
    return lhs;
}
/**
 * @copydoc operator+=()
 * @relatedalso AbstractModule
 */
template <Module M>
M& operator+=(M& lhs, M&& rhs) {
    for (std::move_iterator<typename M::cell_iter_t> it =
             std::make_move_iterator(rhs.cell_cbegin());
         it != std::make_move_iterator(rhs.cell_cend()); it++) {
        lhs.insert(*it, rhs[*it]);
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
    for (typename M::cell_iter_t it = rhs.cell_cbegin(); it != rhs.cell_cend();
         it++) {
        lhs.insert(*it, -rhs[*it]);
    }
    return lhs;
}
/**
 * @copydoc operator-=()
 * @relatedalso AbstractModule
 */
template <Module M>
M& operator-=(M& lhs, M&& rhs) {
    for (std::move_iterator<typename M::cell_iter_t> it =
             std::make_move_iterator(rhs.cell_cbegin());
         it != std::make_move_iterator(rhs.cell_cend()); it++) {
        lhs.insert(*it, -rhs[*it]);
    }
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
    for (typename M::cell_iter_t it = lhs.cell_cbegin(); it != lhs.cell_cend();
         it++) {
        lhs.insert(*it, (rhs - one<typename M::ring_t>()) * lhs[*it]);
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

public:
    /** @copydoc AbstractModule::cell_iter_t */
    using cell_iter_t = citer_t;
    /** @copydoc AbstractModule::cell_t */
    using cell_t = typename AbstractModule<T, R, cell_iter_t>::cell_t;
    /** @copydoc AbstractModule::ring_t */
    using ring_t = typename AbstractModule<T, R, cell_iter_t>::ring_t;
    /** @copydoc AbstractModule::lfunc_t */
    using lfunc_t = typename AbstractModule<T, R, cell_iter_t>::lfunc_t;

    /** @copydoc AbstractModule() */
    UnorderedSetModule() = default;
    /** @copydoc AbstractModule(int) */
    constexpr UnorderedSetModule(int n)
        : AbstractModule<T, R, cell_iter_t>(n) {}

    /** @copydoc AbstractModule::operator[]() */
    R operator[](const T& cell) const {
        return cells.count(cell) == 0 ? zero<R>() : one<R>();
    }

    /** @copydoc AbstractModule::cell_cbegin() */
    cell_iter_t cell_cbegin() const noexcept {
        return cells.cbegin();
    }
    /** @copydoc AbstractModule::cell_cend() */
    cell_iter_t cell_cend() const noexcept {
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

public:
    /** @copydoc AbstractModule::cell_iter_t */
    using cell_iter_t = citer_t;
    /** @copydoc AbstractModule::cell_t */
    using cell_t = typename AbstractModule<T, R, cell_iter_t>::cell_t;
    /** @copydoc AbstractModule::ring_t */
    using ring_t = typename AbstractModule<T, R, cell_iter_t>::ring_t;
    /** @copydoc AbstractModule::lfunc_t */
    using lfunc_t = typename AbstractModule<T, R, cell_iter_t>::lfunc_t;

    /** @copydoc AbstractModule() */
    SetModule() = default;
    /** @copydoc AbstractModule(int) */
    constexpr SetModule(int n) : AbstractModule<T, R, cell_iter_t>(n) {}

    /** @copydoc AbstractModule::operator[]() */
    R operator[](const T& cell) const {
        return cells.count(cell) == 0 ? zero<R>() : one<R>();
    }

    /** @copydoc AbstractModule::cell_cbegin() */
    cell_iter_t cell_cbegin() const noexcept {
        return cells.cbegin();
    }
    /** @copydoc AbstractModule::cell_cend() */
    cell_iter_t cell_cend() const noexcept {
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
 * @brief Constant iterator wrapper for `std::map` and `std::unordered_map`.
 *
 * This iterator only fetches keys when dereferenced rather than key, value
 * pairs.
 *
 * @tparam Map_t Either `std::map` class or `std::unordered_map` class.
 */
template <typename Map_t>
class MapCellIterator {
    typename Map_t::const_iterator iterator;

public:
    /** @brief Constant iterator type. */
    using citer_t = typename Map_t::const_iterator;
    /** @brief Difference type between iterators. */
    using difference_type = typename Map_t::difference_type;
    /** @brief Value type when dereferenced. */
    using value_type = typename Map_t::value_type::first_type;

    /**
     * @brief Construct a new Map Cell Iterator object using a constant
     * iterator.
     *
     * @param it
     */
    MapCellIterator(citer_t it) : iterator(it){};

    /**
     * @brief Dereferencing this iterator yields only the key in the map.
     *
     * @return value_type&
     */
    value_type& operator*() const noexcept {
        return std::get<0>(*iterator);
    }

    /**
     * @brief Equality is determined by equality on wrapped pointers.
     *
     * @param rhs
     * @return true
     * @return false
     */
    bool operator==(const MapCellIterator& rhs) const {
        return iterator == rhs.iterator;
    }

    /**
     * @brief Preincrememnt has standard behavior.
     *
     * @return MapCellIterator&
     */
    MapCellIterator& operator++() {
        iterator++;
        return *this;
    }
    /**
     * @brief Postincrement has standard behavior
     *
     */
    void operator++(int) {
        ++*this;
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
    : public AbstractModule<T, R, MapCellIterator<std::unordered_map<T, R>>> {
    std::unordered_map<T, R> cells;
    using iter_t = typename std::unordered_map<T, R>::iterator;
    using citer_t = typename std::unordered_map<T, R>::const_iterator;

public:
    /** @copydoc AbstractModule::cell_iter_t */
    using cell_iter_t = MapCellIterator<std::unordered_map<T, R>>;
    /** @copydoc AbstractModule::cell_t */
    using cell_t = typename AbstractModule<T, R, cell_iter_t>::cell_t;
    /** @copydoc AbstractModule::ring_t */
    using ring_t = typename AbstractModule<T, R, cell_iter_t>::ring_t;
    /** @copydoc AbstractModule::lfunc_t */
    using lfunc_t = typename AbstractModule<T, R, cell_iter_t>::lfunc_t;

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

    /** @copydoc AbstractModule::cell_cbegin() */
    cell_iter_t cell_cbegin() const noexcept {
        return MapCellIterator<std::unordered_map<T, R>>(cells.cbegin());
    }
    /** @copydoc AbstractModule::cell_cend() */
    cell_iter_t cell_cend() const noexcept {
        return MapCellIterator<std::unordered_map<T, R>>(cells.cend());
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
class MapModule : public AbstractModule<T, R, MapCellIterator<std::map<T, R>>> {
    std::map<T, R> cells;
    using iter_t = typename std::map<T, R>::iterator;
    using citer_t = typename std::map<T, R>::const_iterator;

public:
    /** @copydoc AbstractModule::cell_iter_t */
    using cell_iter_t = MapCellIterator<std::map<T, R>>;
    /** @copydoc AbstractModule::cell_t */
    using cell_t = typename AbstractModule<T, R, cell_iter_t>::cell_t;
    /** @copydoc AbstractModule::ring_t */
    using ring_t = typename AbstractModule<T, R, cell_iter_t>::ring_t;
    /** @copydoc AbstractModule::lfunc_t */
    using lfunc_t = typename AbstractModule<T, R, cell_iter_t>::lfunc_t;

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

    /** @copydoc AbstractModule::cell_cbegin() */
    cell_iter_t cell_cbegin() const noexcept {
        return MapCellIterator<std::map<T, R>>(cells.cbegin());
    }
    /** @copydoc AbstractModule::cell_cend() */
    cell_iter_t cell_cend() const noexcept {
        return MapCellIterator<std::map<T, R>>(cells.cend());
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

} // namespace chomp::modules

#endif // CHAIN_H
