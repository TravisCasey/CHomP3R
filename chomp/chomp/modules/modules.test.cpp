/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <cstddef>
#include <exception>
#include <iterator>
#include <map>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <chomp/modules/modules.hpp>
#include <chomp/modules/rings.hpp>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace chomp::modules {

struct HashableCell {
    std::size_t val;
    constexpr HashableCell(std::size_t n) noexcept : val(n){};
    constexpr bool operator==(const HashableCell& rhs) const noexcept {
        return val == rhs.val;
    }
};

} // namespace chomp::modules

namespace std {
template <>
struct hash<chomp::modules::HashableCell> {
    size_t operator()(const chomp::modules::HashableCell& cell) const {
        return cell.val;
    }
};
} // namespace std

namespace chomp::modules {

struct ComparableCell {
    int val;
    constexpr ComparableCell(int n) noexcept : val(n){};
    constexpr bool operator<(const ComparableCell& rhs) const noexcept {
        return val < rhs.val;
    }
    constexpr bool operator==(const ComparableCell& rhs) const noexcept {
        return val == rhs.val;
    }
};

using ModuleTypes = std::tuple<

    std::tuple<UnorderedSetModule<int, Z<unsigned int, 2>>,
               std::integral_constant<int, 3>,
               std::integral_constant<int, -25>>,

    std::tuple<SetModule<ComparableCell, Z<unsigned long long, 2>>,
               std::integral_constant<ComparableCell, ComparableCell(20)>,
               std::integral_constant<ComparableCell, ComparableCell(-3)>>,

    std::tuple<UnorderedMapModule<HashableCell, Z<unsigned short, 14>>,
               std::integral_constant<HashableCell, HashableCell(2445)>,
               std::integral_constant<HashableCell, HashableCell(0)>>,

    std::tuple<MapModule<bool, int>,
               std::integral_constant<bool, true>,
               std::integral_constant<bool, false>>

>;

TEMPLATE_LIST_TEST_CASE("Module classes model Module concept", "[modules]",
                        ModuleTypes) {
    using M = std::tuple_element_t<0, TestType>;
    CHECK(Module<M>);
}

TEMPLATE_LIST_TEST_CASE("Non-default constructor functions correctly",
                        "[modules]", ModuleTypes) {
    using M = std::tuple_element_t<0, TestType>;
    M elem_0;
    M elem_1(0);

    CHECK_THROWS_AS(M(2), std::domain_error);
    CHECK(elem_0 == elem_1);
}

TEMPLATE_LIST_TEST_CASE("Module classes access, insertion, and iteration",
                        "[modules]", ModuleTypes) {
    using M = std::tuple_element_t<0, TestType>;
    using T = typename M::cell_t;
    using R = typename M::ring_t;
    const T cell_0 = std::tuple_element_t<1, TestType>();
    T cell_1 = std::tuple_element_t<2, TestType>();
    M elem;

    REQUIRE(elem[cell_0] == zero<R>());
    REQUIRE(elem[cell_1] == zero<R>());
    REQUIRE(elem.cell_cbegin() == elem.cell_cend());

    elem.insert(cell_0, zero<R>());
    elem.insert(cell_1, one<R>());

    REQUIRE(elem[cell_0] == zero<R>());
    REQUIRE(elem[cell_1] == one<R>());
    REQUIRE(elem.cell_cbegin() != elem.cell_cend());

    typename M::cell_iter_t it = elem.cell_cbegin();
    REQUIRE((*it == cell_0 || *it == cell_1));
    it++;
    REQUIRE((it == elem.cell_cend() || *it == cell_0 || *it == cell_1));

    elem.insert(cell_0, one<R>());
    elem.insert(cell_1, -one<R>());

    REQUIRE(elem[cell_0] == one<R>());
    REQUIRE(elem[cell_1] == zero<R>());
    REQUIRE(elem.cell_cbegin() != elem.cell_cend());

    it = elem.cell_cbegin();
    REQUIRE((*it == cell_0 || *it == cell_1));
    it++;
    REQUIRE((it == elem.cell_cend() || *it == cell_0 || *it == cell_1));
}

TEMPLATE_LIST_TEST_CASE("Linear function interface to modules", "[modules]",
                        ModuleTypes) {
    using M = std::tuple_element_t<0, TestType>;
    using T = typename M::cell_t;
    using R = typename M::ring_t;
    const T cell_0 = std::tuple_element_t<1, TestType>();
    T cell_1 = std::tuple_element_t<2, TestType>();
    M elem_0;

    typename M::lfunc_t lfunc = [&](const T& cell) {
        M result;
        result.insert(cell, one<R>());
        result.insert(cell_0, zero<R>());
        result.insert(cell_1, one<R>());
        return result;
    };

    M elem_1 = linear_apply(elem_0, lfunc);
    REQUIRE(elem_0[cell_0] == zero<R>());
    REQUIRE(elem_0[cell_1] == zero<R>());
    REQUIRE(elem_1[cell_0] == zero<R>());
    REQUIRE(elem_1[cell_1] == zero<R>());

    elem_0.insert(cell_0, -one<R>());
    elem_1 = linear_apply(elem_0, lfunc);
    REQUIRE(elem_0[cell_0] == -one<R>());
    REQUIRE(elem_0[cell_1] == zero<R>());
    REQUIRE(elem_1[cell_0] == -one<R>());
    REQUIRE(elem_1[cell_1] == -one<R>());

    elem_0.insert(cell_0, -one<R>());
    elem_0.insert(cell_1, one<R>());
    elem_1 = linear_apply(elem_0, lfunc);
    REQUIRE(elem_0[cell_0] == -one<R>() - one<R>());
    REQUIRE(elem_0[cell_1] == one<R>());
    REQUIRE(elem_1[cell_0] == -one<R>() - one<R>());
    REQUIRE(elem_1[cell_1] == zero<R>());
}

TEMPLATE_LIST_TEST_CASE("Comparison operators on modules", "[modules]",
                        ModuleTypes) {
    using M = std::tuple_element_t<0, TestType>;
    using T = typename M::cell_t;
    using R = typename M::ring_t;
    const T cell_0 = std::tuple_element_t<1, TestType>();
    M elem_0;
    M elem_1;

    REQUIRE(elem_0 == elem_1);
    REQUIRE_FALSE(elem_0 != elem_1);

    elem_0.insert(cell_0, zero<R>());
    REQUIRE(elem_0 == elem_1);
    REQUIRE_FALSE(elem_0 != elem_1);

    elem_0.insert(cell_0, one<R>());
    REQUIRE_FALSE(elem_0 == elem_1);
    REQUIRE(elem_0 != elem_1);

    elem_1.insert(cell_0, one<R>());
    REQUIRE(elem_0 == elem_1);
    REQUIRE_FALSE(elem_0 != elem_1);
}

TEMPLATE_LIST_TEST_CASE("Arithmetic operators on modules", "[modules]",
                        ModuleTypes) {
    using M = std::tuple_element_t<0, TestType>;
    using T = typename M::cell_t;
    using R = typename M::ring_t;
    const T cell_0 = std::tuple_element_t<1, TestType>();
    T cell_1 = std::tuple_element_t<2, TestType>();
    M elem_0;
    M elem_1;

    M elem_n = -elem_0;
    M elem_a = elem_0 + elem_1;
    M elem_s = elem_0 - elem_1;
    REQUIRE(elem_0[cell_0] == zero<R>());
    REQUIRE(elem_0[cell_1] == zero<R>());
    REQUIRE(elem_1[cell_0] == zero<R>());
    REQUIRE(elem_1[cell_1] == zero<R>());
    REQUIRE(elem_n[cell_0] == zero<R>());
    REQUIRE(elem_n[cell_1] == zero<R>());
    REQUIRE(elem_a[cell_0] == zero<R>());
    REQUIRE(elem_a[cell_1] == zero<R>());
    REQUIRE(elem_s[cell_0] == zero<R>());
    REQUIRE(elem_s[cell_1] == zero<R>());
    REQUIRE(elem_s == elem_0 + (-elem_1));
    REQUIRE(elem_a == -elem_n + elem_1);

    elem_0.insert(cell_0, one<R>());
    elem_1.insert(cell_1, one<R>());

    elem_n = -elem_0;
    elem_a = elem_0 + elem_1;
    elem_s = elem_0 - elem_1;
    REQUIRE(elem_0[cell_0] == one<R>());
    REQUIRE(elem_0[cell_1] == zero<R>());
    REQUIRE(elem_1[cell_0] == zero<R>());
    REQUIRE(elem_1[cell_1] == one<R>());
    REQUIRE(elem_n[cell_0] == -one<R>());
    REQUIRE(elem_n[cell_1] == zero<R>());
    REQUIRE(elem_a[cell_0] == one<R>());
    REQUIRE(elem_a[cell_1] == one<R>());
    REQUIRE(elem_s[cell_0] == one<R>());
    REQUIRE(elem_s[cell_1] == -one<R>());
    REQUIRE(elem_s == elem_0 + (-elem_1));
    REQUIRE(elem_a == -elem_n + elem_1);

    elem_0 *= one<R>() + one<R>();
    elem_1 = -one<R>() * elem_1 * -one<R>();

    elem_n = -elem_0;
    elem_a = elem_0 + elem_1;
    elem_s = elem_0 - elem_1;
    REQUIRE(elem_0[cell_0] == one<R>() + one<R>());
    REQUIRE(elem_0[cell_1] == zero<R>());
    REQUIRE(elem_1[cell_0] == zero<R>());
    REQUIRE(elem_1[cell_1] == one<R>());
    REQUIRE(elem_n[cell_0] == -(one<R>() + one<R>()));
    REQUIRE(elem_n[cell_1] == zero<R>());
    REQUIRE(elem_a[cell_0] == one<R>() + one<R>());
    REQUIRE(elem_a[cell_1] == one<R>());
    REQUIRE(elem_s[cell_0] == one<R>() + one<R>());
    REQUIRE(elem_s[cell_1] == -one<R>());
    REQUIRE(elem_s == elem_0 + (-elem_1));
    REQUIRE(elem_a == -elem_n + elem_1);

    REQUIRE(elem_a + std::move(elem_s) == (one<R>() + one<R>()) * elem_0);

    elem_0 += elem_1;
    REQUIRE(elem_0[cell_0] == one<R>() + one<R>());
    REQUIRE(elem_0[cell_1] == one<R>());
    elem_1 -= std::move(elem_0);
    REQUIRE(elem_1[cell_0] == -one<R>() - one<R>());
    REQUIRE(elem_1[cell_1] == zero<R>());

    REQUIRE(elem_1 * -one<R>() == -elem_1);
}

using CellAndRingTypes = std::tuple<

    std::tuple<int, Z<unsigned int, 2>,
               UnorderedSetModule<int, Z<unsigned int, 2>>>,

    std::tuple<int, Z<unsigned int, 3>,
               UnorderedMapModule<int, Z<unsigned int, 3>>>,

    std::tuple<std::vector<short>, Z<unsigned long long, 2>,
               SetModule<std::vector<short>, Z<unsigned long long, 2>>>,

    std::tuple<std::vector<short>, Z<unsigned long long, 5>,
               MapModule<std::vector<short>, Z<unsigned long long, 5>>>

>;

TEMPLATE_LIST_TEST_CASE("DefaultModule chooses Module type correctly",
                        "[modules]", CellAndRingTypes) {
    using T = std::tuple_element_t<0, TestType>;
    using R = std::tuple_element_t<1, TestType>;
    using M = std::tuple_element_t<2, TestType>;

    CHECK(std::same_as<typename DefaultModule<T, R>::type, M>);
}

using MapTypes = std::tuple<
    std::tuple<std::map<ComparableCell, short>,
        std::integral_constant<ComparableCell, ComparableCell(33)>,
        std::integral_constant<short, -20>,
        std::integral_constant<ComparableCell, ComparableCell(0)>,
        std::integral_constant<short, 1>>,

    std::tuple<std::unordered_map<int, int>,
        std::integral_constant<int, 29>,
        std::integral_constant<int, 29>,
        std::integral_constant<int, -29>,
        std::integral_constant<int, 0>>
>;

TEMPLATE_LIST_TEST_CASE("KeyIterator functions correctly as a forward iterator",
                        "[modules]", MapTypes) {
    using Map_t = std::tuple_element_t<0, TestType>;
    using key_type = typename Map_t::key_type;
    using mapped_type = typename Map_t::mapped_type;
    key_type key_0 = std::tuple_element_t<1, TestType>();
    mapped_type val_0 = std::tuple_element_t<2, TestType>();
    key_type key_1 = std::tuple_element_t<3, TestType>();
    mapped_type val_1 = std::tuple_element_t<4, TestType>();

    using I = KeyIterator<Map_t>;
    CHECK(std::forward_iterator<I>);

    Map_t m;
    m[key_0] = val_0;
    m[key_1] = val_1;

    I it(m.cbegin());
    CHECK(m.count(*it) == 1);
    CHECK(it == m.cbegin());

    I temp = it++;
    CHECK(temp != it);
    CHECK(temp == m.cbegin());
    CHECK(m.count(*it) == 1);
    CHECK(++it == m.cend());
}

} // namespace chomp::modules

#endif // DOXYGEN_SHOULD_SKIP_THIS
