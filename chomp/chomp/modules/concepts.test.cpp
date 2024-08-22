/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <chomp/modules/concepts.hpp>

#ifndef DOXYGEN_SHOULD_SKIP_THIS

namespace chomp::modules {

TEST_CASE("Fundamental types and Group/Ring concept", "[rings]") {
    CHECK_FALSE(Group<std::string>);
    CHECK(Ring<int>);
    CHECK(Ring<float>);
    CHECK_FALSE(Group<void>);
    CHECK(Ring<char>);
}

TEST_CASE("Hashable concept functions as expected", "[modules]") {
    CHECK(Hashable<int>);
    CHECK(Hashable<std::vector<bool>>);
    CHECK_FALSE(Hashable<std::vector<int>>);
}

TEST_CASE("Comparable concept functions as expected", "[modules]") {
    struct NonComparable {
        int val;
        NonComparable(int n) : val(n) {}
        bool operator==(const NonComparable& rhs) {
            return val == rhs.val;
        }
    };
    CHECK(Comparable<int>);
    CHECK(Comparable<std::vector<unsigned short>>);
    CHECK_FALSE(Comparable<NonComparable>);
}

} // namespace chomp::modules

#endif // DOXYGEN_SHOULD_SKIP_THIS
