/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <iterator>
#include <map>
#include <tuple>
#include <type_traits>
#include <unordered_map>

#include <catch2/catch_template_test_macros.hpp>

#include <chomp/util/iterators.hpp>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace chomp::util {

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

    Map_t m;
    m[key_0] = val_0;
    m[key_1] = val_1;

    auto it = KeyIterator(m.cbegin());
    CHECK(std::forward_iterator<decltype(it)>);

    CHECK(m.count(*it) == 1);
    CHECK(it == m.cbegin());

    auto temp = it++;
    CHECK(temp != it);
    CHECK(temp == m.cbegin());
    CHECK(m.count(*it) == 1);
    CHECK(++it == m.cend());
}

} // namespace chomp::util

#endif // DOXYGEN_SHOULD_SKIP_THIS