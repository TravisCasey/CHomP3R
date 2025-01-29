/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/algebra/algebra.hpp>
#include <chomp/algebra/cyclic.hpp>
#include <chomp/algebra/modules.hpp>
#include <chomp/complexes/complexes.hpp>
#include <chomp/complexes/cubical.hpp>
#include <chomp/complexes/grading.hpp>
#include <chomp/util/concepts.hpp>
#include <chomp/util/constants.hpp>

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <set>
#include <unordered_set>
#include <utility>
#include <vector>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

TEST_CASE("Orthant class", "[complexes]") {
  CHECK(std::regular<Orthant<5>>);
  CHECK(Hashable<Orthant<20>>);

  SECTION("Default construction") {
    Orthant<3> def;
    CHECK(*def.data() == std::array<HypercubeCoordinate, 3>{0, 0, 0});
  }

  SECTION("Initialization list construction") {
    const Orthant<5> orth{5, -6};
    CHECK(*orth.data() == std::array<HypercubeCoordinate, 5>{5, -6, 0, 0, 0});
  }

  SECTION("Copy construction and copy assignment") {
    Orthant<4> original{0, 1, 2, 3};
    Orthant<4> copy(original);
    Orthant<4> copy_assignment;
    copy_assignment = original;

    CHECK(*original.data() == *copy.data());
    CHECK(*original.data() == *copy_assignment.data());
  }

  SECTION("Access and iteration methods") {
    std::array<HypercubeCoordinate, 3> arr{-1, 2, 9};
    const Orthant<3> const_orth{-1, 2, 9};
    Orthant<3> plain_orth{-1, 2, 9};

    std::size_t axis = 0;
    typename std::array<HypercubeCoordinate, 3>::iterator arr_it = arr.begin();
    typename std::array<HypercubeCoordinate, 3>::reverse_iterator rarr_it
        = arr.rbegin();
    typename Orthant<3>::CIterType const_it = const_orth.begin();
    typename Orthant<3>::IterType plain_it = plain_orth.begin();
    typename Orthant<3>::CRIterType rconst_it = const_orth.rbegin();
    typename Orthant<3>::RIterType rplain_it = plain_orth.rbegin();

    CHECK(arr[axis] == const_orth[axis]);
    CHECK(arr[axis] == plain_orth[axis]);
    CHECK(*arr_it == *const_it);
    CHECK(*arr_it == *plain_it);
    CHECK(*rarr_it == *rconst_it);
    CHECK(*rarr_it == *rplain_it);

    ++axis;
    ++arr_it;
    ++rarr_it;
    ++const_it;
    ++plain_it;
    ++rconst_it;
    ++rplain_it;

    CHECK(arr[axis] == const_orth[axis]);
    CHECK(arr[axis] == plain_orth[axis]);
    CHECK(*arr_it == *const_it);
    CHECK(*arr_it == *plain_it);
    CHECK(*rarr_it == *rconst_it);
    CHECK(*rarr_it == *rplain_it);

    ++axis;
    ++arr_it;
    ++rarr_it;
    ++const_it;
    ++plain_it;
    ++rconst_it;
    ++rplain_it;

    CHECK(arr[axis] == const_orth[axis]);
    CHECK(arr[axis] == plain_orth[axis]);
    CHECK(*arr_it == *const_it);
    CHECK(*arr_it == *plain_it);
    CHECK(*rarr_it == *rconst_it);
    CHECK(*rarr_it == *rplain_it);

    ++axis;
    ++const_it;
    ++plain_it;
    ++rconst_it;
    ++rplain_it;

    CHECK(axis == const_orth.size());
    CHECK(const_it == const_orth.end());
    CHECK(plain_it == plain_orth.end());
    CHECK(rconst_it == const_orth.rend());
    CHECK(rplain_it == plain_orth.rend());
  }

  SECTION("Comparison operators") {
    Orthant<5> first{0, 1, 2, 3, 4};
    Orthant<5> second = first;
    CHECK(first == second);
    CHECK((first <=> second) == std::strong_ordering::equal);

    second[2] = 5;
    CHECK(!(first == second));
    CHECK((first <=> second) == std::strong_ordering::less);

    first[0] = 2;
    CHECK(!(first == second));
    CHECK((first <=> second) == std::strong_ordering::greater);
  }
}

TEST_CASE("Cube class", "[complexes]") {
  CHECK(Hashable<Cube<5>>);
  CHECK(std::regular<Cube<3>>);

  SECTION("Default construction") {
    Cube<5> def;
    CHECK(def.base() == Orthant<5>{0, 0, 0, 0, 0});
    CHECK(def.dual() == Orthant<5>{0, 0, 0, 0, 0});
    CHECK(def.extent() == 0b11111);
    CHECK(def.dimension() == 5);
  }

  SECTION("Construction and access methods") {
    Orthant<3> base{2, 1, -4};
    Orthant<3> dual{2, 0, -4};
    Orthant<3> base_copy(base);
    Orthant<3> dual_copy(dual);
    const std::size_t extent = 0b101;
    const std::size_t dim = 2;
    std::vector<Cube<3>> cube_vec;

    cube_vec.emplace_back(base, dual);
    cube_vec.emplace_back(base, std::move(dual_copy));
    dual_copy = dual;
    cube_vec.emplace_back(std::move(base_copy), dual);
    base_copy = base;
    cube_vec.emplace_back(std::move(base_copy), std::move(dual_copy));

    base_copy = base;
    cube_vec.emplace_back(base, extent);
    cube_vec.emplace_back(std::move(base_copy), extent);

    for (const Cube<3>& current_cube : cube_vec) {
      CHECK(current_cube.base() == base);
      CHECK(current_cube.dual() == dual);
      CHECK(current_cube.extent() == extent);
      CHECK(current_cube.dimension() == dim);

      CHECK(current_cube.base(0) == 2);
      CHECK(current_cube.base(1) == 1);
      CHECK(current_cube.base(2) == -4);
      CHECK(current_cube.dual(0) == 2);
      CHECK(current_cube.dual(1) == 0);
      CHECK(current_cube.dual(2) == -4);
    }

    for (Cube<3>& current_cube : cube_vec) {
      CHECK(current_cube.base() == base);
      CHECK(current_cube.dual() == dual);
      CHECK(current_cube.extent() == extent);
      CHECK(current_cube.dimension() == dim);

      CHECK(current_cube.base(0) == 2);
      CHECK(current_cube.base(1) == 1);
      CHECK(current_cube.base(2) == -4);
      CHECK(current_cube.dual(0) == 2);
      CHECK(current_cube.dual(1) == 0);
      CHECK(current_cube.dual(2) == -4);
    }
  }

  SECTION("Comparison operators") {
    Cube<4> first({0, 1, 0, -1}, 0b0111);
    Cube<4> second({0, 1, 0, -1}, {0, 1, 0, -2});

    CHECK(first == second);
    CHECK((first <=> second) == std::strong_ordering::equal);

    --first.dual()[2];
    CHECK(first != second);
    CHECK((first <=> second) == std::strong_ordering::less);

    --second.base()[2];
    CHECK(first != second);
    CHECK((first <=> second) == std::strong_ordering::greater);
  }
}

TEST_CASE("BFSDualOrthantIterator class", "[complexes]") {
  CHECK(std::input_iterator<BFSDualOrthantIterator<10>>);

  SECTION("BFS iteration over full orthant") {
    const Orthant<3> base{2, 1, 2};
    std::vector<std::pair<Orthant<3>, std::size_t>> correct{{{2, 1, 2}, 0b111},
        {{1, 1, 2}, 0b110}, {{2, 0, 2}, 0b101}, {{2, 1, 1}, 0b011},
        {{1, 0, 2}, 0b100}, {{1, 1, 1}, 0b010}, {{2, 0, 1}, 0b001},
        {{1, 0, 1}, 0b000}};
    std::vector<std::pair<Orthant<3>, std::size_t>> it_result;

    for (BFSDualOrthantIterator<3> bfs_it(base); bfs_it; ++bfs_it) {
      it_result.push_back(*bfs_it);
    }
    CHECK(correct == it_result);
  }

  SECTION("BFS iteration with maximum and minimum orthants") {
    Orthant<5> minimum{-1, 3, 8, -11, 10};
    Orthant<5> maximum{-1, 3, 9, -10, 10};
    const std::size_t maximum_extent = 0b11110;
    BFSDualOrthantIterator<5> bfs_it(std::move(minimum), std::move(maximum),
        maximum_extent);

    std::vector<std::pair<Orthant<5>, std::size_t>> correct{
        {{-1, 3, 9, -10, 10}, 0b11110}, {{-1, 3, 8, -10, 10}, 0b11010},
        {{-1, 3, 9, -11, 10}, 0b10110}, {{-1, 3, 8, -11, 10}, 0b10010}};
    std::vector<std::pair<Orthant<5>, std::size_t>> it_result;

    for (; bfs_it; bfs_it++) {
      it_result.push_back(*bfs_it);
    }
    CHECK(correct == it_result);
  }
}

TEST_CASE("LexOrthantIterator class.", "[complexes]") {
  CHECK(std::forward_iterator<LexOrthantIterator<10>>);

  SECTION("Iteration over partial grid") {
    LexOrthantIterator<4> it({0, -1, -1, -2}, {2, 1, 4, -2}, {1, 1, 1, -2});

    // Bool conversion, increment, dereference, member access via pointer
    REQUIRE(it);
    CHECK(*it == Orthant<4>{1, 1, 1, -2});
    ++it;
    CHECK(*it == Orthant<4>{2, 1, 1, -2});
    CHECK(*(it++) == Orthant<4>{2, 1, 1, -2});
    CHECK(*it == Orthant<4>{0, -1, 2, -2});
    CHECK(it->size() == 4);

    // Bool conversion returns false after correct number of steps
    for (int i = 0; i < 26; ++i) {
      ++it;
      REQUIRE(it);
    }
    ++it;
    REQUIRE(!it);

    // Iterator loops back to the minimum orthant
    CHECK(*it == Orthant<4>{0, -1, -1, -2});
  }

  SECTION("Iteration over full grid") {
    LexOrthantIterator<2> it({0, 0}, {1, 3});

    // Bool conversion, increment, dereference, member access via pointer
    REQUIRE(it);
    CHECK(*it == Orthant<2>{0, 0});
    ++it;
    CHECK(*it == Orthant<2>{1, 0});
    CHECK(*(it++) == Orthant<2>{1, 0});
    CHECK(*it == Orthant<2>{0, 1});
    CHECK(it->size() == 2);

    // Bool conversion returns false after correct number of steps
    for (int i = 0; i < 5; ++i) {
      ++it;
      REQUIRE(it);
    }
    ++it;
    REQUIRE(!it);

    // Iterator loops back to the minimum orthant
    CHECK(*it == Orthant<2>{0, 0});
  }

  SECTION("Equality operator") {
    LexOrthantIterator<3> first({2, 5, -2}, {14, 10, 0});
    LexOrthantIterator<3> second({2, 5, -2}, {14, 10, 0}, {3, 5, -2});
    LexOrthantIterator<3> different({2, 5, -2}, {3, 7, 3});

    // Equal bounds, pointed-to orthants are not equal
    CHECK(*first != *second);
    CHECK(first != second);

    // Equal pointed-to orthants, inequal bounds
    CHECK(*first == *different);
    CHECK(first != different);

    // Equal following increment, and equality preservation under increment
    CHECK(++first == second);
    CHECK(++first == ++second);
  }
}

TEST_CASE("CubeIterator class", "[complexes]") {
  CHECK(std::forward_iterator<CubeIterator<8>>);

  SECTION("Iteration, dereferencing, and bool conversion") {
    CubeIterator<2> it({-2, 1}, {0, 2});

    REQUIRE(it);
    CHECK(*(it++) == Cube<2>({-2, 1}, 0b11));
    CHECK(*it == Cube<2>({-2, 1}, 0b10));

    std::vector<Cube<2>> correct{Cube<2>({-2, 1}, 0b10), Cube<2>({-2, 1}, 0b01),
        Cube<2>({-2, 1}, 0b00), Cube<2>({-1, 1}, 0b11), Cube<2>({-1, 1}, 0b10),
        Cube<2>({-1, 1}, 0b01), Cube<2>({-1, 1}, 0b00), Cube<2>({0, 1}, 0b11),
        Cube<2>({0, 1}, 0b10), Cube<2>({0, 1}, 0b01), Cube<2>({0, 1}, 0b00),
        Cube<2>({-2, 2}, 0b11), Cube<2>({-2, 2}, 0b10), Cube<2>({-2, 2}, 0b01),
        Cube<2>({-2, 2}, 0b00), Cube<2>({-1, 2}, 0b11), Cube<2>({-1, 2}, 0b10),
        Cube<2>({-1, 2}, 0b01), Cube<2>({-1, 2}, 0b00), Cube<2>({0, 2}, 0b11),
        Cube<2>({0, 2}, 0b10), Cube<2>({0, 2}, 0b01), Cube<2>({0, 2}, 0b00)};

    std::vector<Cube<2>> it_result;
    for (; it; ++it) {
      it_result.push_back(*it);
    }

    CHECK(it_result == correct);
  }

  SECTION("Equality operator and sentinels") {
    CubeIterator<3> first({2, 3, 2}, {2, 3, 2});
    CubeIterator<3> second = first;
    CubeIterator<3> different({2, 3, 2}, {2, 3, 3});
    CubeIterator<3> sent = CubeIterator<3>::sentinel();

    // Equal pointed-to cube, equal bounds (copy)
    CHECK(*first == *second);
    CHECK(first == second);

    // Equal pointed-to cube, inequal bounds
    CHECK(*first == *different);
    CHECK(first != different);

    // Inequal pointed-to-cube, equal bounds
    ++second;
    CHECK(*first != *second);
    CHECK(first != second);

    // Equality preservation by increment
    ++first;
    CHECK(*first == *second);
    CHECK(first == second);

    // Sentinel and sentinel
    CHECK(sent == CubeIterator<3>::sentinel());

    // Sentinel and nonterminal iterators
    REQUIRE(sent != first);
    REQUIRE(sent != second);
    REQUIRE(sent != different);

    // Sentinel and non-terminal iterators
    for (; first; ++first) {
      REQUIRE(first != sent);
    }
    for (; second; ++second) {
      REQUIRE(second != sent);
    }
    for (; different; ++different) {
      REQUIRE(different != sent);
    }

    // Sentinel and terminal iterators
    REQUIRE(sent == first);
    REQUIRE(sent == second);
    REQUIRE(sent == different);
  }
}

TEST_CASE("TopCubeSetGrading class.", "[complexes]") {
  CHECK(BoundedGrading<TopCubeSetGrading<8, 4, 10, std::set>>);
  CHECK(std::semiregular<TopCubeSetGrading<2>>);

  SECTION("Bound access") {
    using BA_Default = BoundAccess<TopCubeSetGrading<10>>;
    CHECK(BA_Default::Minimum == 0);
    CHECK(BA_Default::Maximum == 1);

    using BA = BoundAccess<TopCubeSetGrading<4, 10, 19>>;
    CHECK(BA::Minimum == 10);
    CHECK(BA::Maximum == 19);
  }

  SECTION("`set` copy construction and querying") {
    const std::set<Orthant<2>> top_cube_set{{0, 1}};
    const TopCubeSetGrading<2, 2, 4, std::set> grading_func(top_cube_set);

    // Cube querying
    const std::vector<Cube<2>> correct_min_cubes{{{0, 1}, 0b11}, {{0, 1}, 0b10},
        {{0, 1}, 0b01}, {{0, 1}, 0b00}, {{1, 1}, 0b10}, {{1, 1}, 0b00},
        {{0, 2}, 0b01}, {{0, 2}, 0b00}, {{1, 2}, 0b00}};
    std::vector<Cube<2>> min_cubes;
    for (CubeIterator<2> it({0, 0}, {3, 3}); it; ++it) {
      if (grading_func(*it) == decltype(grading_func)::Minimum()) {
        min_cubes.push_back(*it);
      } else {
        CHECK(grading_func(*it) == decltype(grading_func)::Maximum());
      }
    }
    CHECK(correct_min_cubes == min_cubes);

    // Orthant querying
    const std::vector<Orthant<2>> correct_min_orthants{{0, 1}};
    std::vector<Orthant<2>> min_orthants;
    for (LexOrthantIterator<2> it({0, 0}, {3, 3}); it; ++it) {
      if (grading_func(*it) == decltype(grading_func)::Minimum()) {
        min_orthants.push_back(*it);
      } else {
        CHECK(grading_func(*it) == decltype(grading_func)::Maximum());
      }
    }
    CHECK(correct_min_orthants == min_orthants);
  }

  SECTION("`unordered_set` move construction and querying.") {
    std::unordered_set<Orthant<4>> top_cube_set{{1, 1, 1, 1}, {1, 1, 1, 0},
        {0, 0, 0, 0}};
    const TopCubeSetGrading<4> grading_func(std::move(top_cube_set));

    // Cube Querying
    for (std::size_t extent = 0; extent < (1 << 4); ++extent) {
      CHECK(grading_func(Cube<4>({1, 1, 1, 1}, extent)) == 0);
    }
    CHECK(grading_func(Cube<4>({2, 2, 2, 2}, 0)) == 0);
    for (std::size_t extent = 1; extent < (1 << 4); ++extent) {
      CHECK(grading_func(Cube<4>({2, 2, 2, 2}, extent)) == 1);
    }
    CHECK(grading_func(Cube<4>({1, 0, 0, 0}, 14)) == 0);
    CHECK(grading_func(Cube<4>({1, 0, 0, 0}, 1)) == 1);

    // Orthant querying
    const std::vector<Orthant<4>> correct_min_orthants{{0, 0, 0, 0},
        {1, 1, 1, 0}, {1, 1, 1, 1}};
    std::vector<Orthant<4>> min_orthants;
    for (LexOrthantIterator<4> it({0, 0, 0, 0}, {2, 2, 2, 2}); it; ++it) {
      if (grading_func(*it) == 0) {
        min_orthants.push_back(*it);
      } else {
        CHECK(grading_func(*it) == 1);
      }
    }
    CHECK(correct_min_orthants == min_orthants);
  }

  SECTION("`initializer_list` construction and querying") {
    const TopCubeSetGrading<2> grading_func{{0, 0}, {1, 0}, {2, 0}, {0, 1},
        {2, 1}, {0, 2}, {1, 2}, {2, 2}};

    // Cube Querying
    const std::vector<Cube<2>> correct_max_cubes{{{1, 1}, 0b11}};
    std::vector<Cube<2>> max_cubes;
    for (CubeIterator<2> it({0, 0}, {2, 2}); it; ++it) {
      if (grading_func(*it) == 1) {
        max_cubes.push_back(*it);
      } else {
        CHECK(grading_func(*it) == 0);
      }
    }
    CHECK(correct_max_cubes == max_cubes);

    // Orthant Querying
    const std::vector<Orthant<2>> correct_max_orthants{{1, 1}};
    std::vector<Orthant<2>> max_orthants;
    for (LexOrthantIterator<2> it({0, 0}, {2, 2}); it; ++it) {
      if (grading_func(*it) == 1) {
        max_orthants.push_back(*it);
      } else {
        CHECK(grading_func(*it) == 0);
      }
    }
    CHECK(correct_max_orthants == max_orthants);
  }
}


TEST_CASE("CubicalComplex class", "[complexes]") {
  CHECK(ChainComplex<CubicalComplex<5, TopCubeSetGrading<5>>>);
  CHECK(std::copyable<CubicalComplex<3, MapGrading<Cube<3>, 0, 45>, Z<7>>>);

  SECTION("Orthant and dimension access") {
    Orthant<5> minimum{-2, -2, -2, 0, 1};
    Orthant<5> maximum{0, 2, 2, 2, 5};
    const TrivialGrading<Cube<5>, 0> grading_func;
    const CubicalComplex<5, decltype(grading_func)> complex(minimum, maximum,
        grading_func);

    CHECK(decltype(complex)::dimension == 5);
    CHECK(complex.minimum() == minimum);
    CHECK(complex.maximum() == maximum);

    for (std::size_t axis = 0; axis < 5; ++axis) {
      CHECK(complex.minimum(axis) == minimum[axis]);
      CHECK(complex.maximum(axis) == maximum[axis]);
    }
  }

  SECTION("Grading, boundary and coboundary operators") {
    // Forms a closed 2-cell with least vertex at (0, 0, 0) along axes 1 and 3
    const std::initializer_list<Cube<3>> zero_cube_ilist
        = {{{0, 0, 0}, 0b000}, {{0, 0, 0}, 0b001}, {{0, 0, 0}, 0b100},
            {{0, 0, 0}, 0b101}, {{1, 0, 0}, 0b000}, {{1, 0, 0}, 0b100},
            {{0, 0, 1}, 0b000}, {{0, 0, 1}, 0b001}, {{1, 0, 1}, 0b000}};
    const std::vector<Cube<3>> zero_cube_vec(zero_cube_ilist);
    SetGrading<Cube<3>, 0, 1> grading_func(zero_cube_ilist);
    CubicalComplex<3, decltype(grading_func), Z<5>> complex({2, 4, 5},
        std::move(grading_func));
    using ChainType = typename decltype(complex)::ChainType;

    std::vector<Cube<3>>::const_iterator it = zero_cube_vec.cbegin();
    CHECK(complex.grade(*it) == 0);

    // *it = Cube<3>({0, 0, 0}, 0b000)
    // Boundary of vertex is empty
    ChainType boundary_result;
    CHECK(boundary(complex, *it) == boundary_result);
    CHECK(graded_boundary(complex, *it) == boundary_result);
    CHECK(closure_boundary(complex, *it) == boundary_result);

    ++it;
    CHECK(complex.grade(*it) == 0);

    // *it = Cube<3>({0, 0, 0}, 0b001)
    // Boundary of 1-cube is outer vertex - inner vertex
    boundary_result.insert(Cube<3>({1, 0, 0}, 0b000), one<Z<5>>());
    boundary_result.insert(Cube<3>({0, 0, 0}, 0b000), -one<Z<5>>());
    CHECK(boundary(complex, *it) == boundary_result);
    CHECK(graded_boundary(complex, *it) == boundary_result);
    CHECK(closure_boundary(complex, *it) == boundary_result);

    ++it;
    CHECK(complex.grade(*it) == 0);
    ++it;
    CHECK(complex.grade(*it) == 0);

    // *it = Cube<3>({0, 0, 0}, 0b101)
    // Boundary of this 2-cube is:
    // outer 0b100 - inner 0b100 - (inner 0b001 - outer 0b001)
    boundary_result.clear();
    boundary_result.insert(Cube<3>({0, 0, 1}, 0b001), -one<Z<5>>());
    boundary_result.insert(Cube<3>({1, 0, 0}, 0b100), one<Z<5>>());
    boundary_result.insert(Cube<3>({0, 0, 0}, 0b001), one<Z<5>>());
    boundary_result.insert(Cube<3>({0, 0, 0}, 0b100), -one<Z<5>>());
    CHECK(boundary(complex, *it) == boundary_result);
    CHECK(graded_boundary(complex, *it) == boundary_result);
    CHECK(closure_boundary(complex, *it) == boundary_result);

    // Check boundary of boundary is empty
    CHECK(boundary(complex, boundary_result) == ChainType());

    // 3-cube is not in grading set
    CHECK(complex.grade(Cube<3>({0, 0, 0}, 0b111)) == 1);

    // Complex is open along the maximum orthants; hence boundary of this 1-cube
    // is only - inner vertex
    boundary_result.clear();
    boundary_result.insert(Cube<3>({2, 4, 5}, 0b000), -one<Z<5>>());
    CHECK(boundary(complex, {{2, 4, 5}, 0b001}) == boundary_result);
    CHECK(graded_boundary(complex, {{2, 4, 5}, 0b001}) == boundary_result);
    CHECK(closure_boundary(complex, {{2, 4, 5}, 0b001}) == boundary_result);

    // Cube<3>({0, 0, 0}, 0b010) is not in grading set, but its inner vertex is;
    // Both are in its closure, however.
    boundary_result.clear();
    boundary_result.insert(Cube<3>({0, 1, 0}, 0b000), one<Z<5>>());
    CHECK(graded_boundary(complex, {{0, 0, 0}, 0b010}) == boundary_result);
    boundary_result.insert(Cube<3>({0, 0, 0}, 0b000), -one<Z<5>>());
    CHECK(boundary(complex, {{0, 0, 0}, 0b010}) == boundary_result);
    CHECK(closure_boundary(complex, {{0, 0, 0}, 0b010}) == boundary_result);

    // Coboundary of 3-cube is empty (in 3 dimensional complex).
    ChainType coboundary_result;
    CHECK(coboundary(complex, {{0, 0, 0}, 0b111}) == coboundary_result);
    CHECK(graded_coboundary(complex, {{0, 0, 0}, 0b111}) == coboundary_result);
    CHECK(closure_coboundary(complex, {{0, 0, 0}, 0b111}) == coboundary_result);

    // Complex is closed along minimum orthant; coboundary of least vertex has
    // just outer 1-cubes along each axis.
    // Only axis 1 and 3 are in grading set, though.
    coboundary_result.insert(Cube<3>({0, 0, 0}, 0b001), -one<Z<5>>());
    coboundary_result.insert(Cube<3>({0, 0, 0}, 0b100), -one<Z<5>>());
    CHECK(graded_coboundary(complex, {{0, 0, 0}, 0b000}) == coboundary_result);
    CHECK(closure_coboundary(complex, {{0, 0, 0}, 0b000}) == coboundary_result);
    coboundary_result.insert(Cube<3>({0, 0, 0}, 0b010), -one<Z<5>>());
    CHECK(coboundary(complex, {{0, 0, 0}, 0b000}) == coboundary_result);

    // Test linear application of boundary and coboundary operators on a chain.
    // This is a cycle of 1-cubes and hence has boundary zero but nonzero
    // coboundary.
    ChainType chain;
    chain.insert(Cube<3>({0, 0, 0}, 0b001), one<Z<5>>());
    chain.insert(Cube<3>({1, 0, 0}, 0b100), one<Z<5>>());
    chain.insert(Cube<3>({0, 0, 1}, 0b001), -one<Z<5>>());
    chain.insert(Cube<3>({0, 0, 0}, 0b100), -one<Z<5>>());

    // Cycle has 0 boundary
    boundary_result.clear();
    CHECK(boundary(complex, chain) == boundary_result);
    CHECK(graded_boundary(complex, chain) == boundary_result);
    CHECK(closure_boundary(complex, chain) == boundary_result);

    // Coboundary is more complex; in grade 0, each cell contributes a copy of
    // Cube<3>({0, 0, 0}, 0b101)
    // In the full complex, many more cells are present.
    coboundary_result.clear();
    coboundary_result.insert(Cube<3>({0, 0, 0}, 0b101), Z<5>(4));
    CHECK(graded_coboundary(complex, chain) == coboundary_result);
    CHECK(closure_coboundary(complex, chain) == coboundary_result);
    coboundary_result.insert(Cube<3>({1, 0, 0}, 0b101), -one<Z<5>>());
    coboundary_result.insert(Cube<3>({0, 0, 1}, 0b101), -one<Z<5>>());
    coboundary_result.insert(Cube<3>({0, 0, 0}, 0b011), one<Z<5>>());
    coboundary_result.insert(Cube<3>({0, 0, 0}, 0b110), one<Z<5>>());
    coboundary_result.insert(Cube<3>({1, 0, 0}, 0b110), -one<Z<5>>());
    coboundary_result.insert(Cube<3>({0, 0, 1}, 0b011), -one<Z<5>>());
    CHECK(coboundary(complex, chain) == coboundary_result);

    // Check that coboundary of coboundary is empty
    CHECK(coboundary(complex, coboundary_result) == ChainType());
  }
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
