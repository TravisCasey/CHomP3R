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

#include <catch2/catch_test_macros.hpp>

#include <concepts>
#include <initializer_list>
#include <utility>
#include <vector>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

TEST_CASE(
    "CubicalComplex boundary and coboundary operators function correctly",
    "[complexes]"
) {
  // Create grading: these cells are 0, all else are 1
  // Form a complete 2-cell with least vertex at (0, 0, 0) along axes 1 and 3
  const std::initializer_list<Cube<3>> zero_cube_ilist = {
      Cube<3>({0, 0, 0}, 0b000), Cube<3>({0, 0, 0}, 0b001),
      Cube<3>({0, 0, 0}, 0b100), Cube<3>({0, 0, 0}, 0b101),
      Cube<3>({1, 0, 0}, 0b000), Cube<3>({1, 0, 0}, 0b100),
      Cube<3>({0, 0, 1}, 0b000), Cube<3>({0, 0, 1}, 0b001),
      Cube<3>({1, 0, 1}, 0b000)
  };
  const std::vector<Cube<3>> zero_cube_vec(zero_cube_ilist);
  CachedGradingWrapper grading_func(
      SetGrading<Cube<3>, 0, 1>(zero_cube_ilist), 4
  );

  // Create cubical complex and ensure initialization and type correctness
  CubicalComplex<3, decltype(grading_func), Z<5>> complex(
      CubeOrthant<3>{2, 4, 5}, std::move(grading_func)
  );
  using TestChainType = typename decltype(complex)::ChainType;
  REQUIRE(std::same_as<TestChainType, UnorderedMapModule<Cube<3>, Z<5>>>);
  REQUIRE(complex.minimum() == CubeOrthant<3>({0, 0, 0}));
  REQUIRE(complex.maximum() == CubeOrthant<3>({2, 4, 5}));

  std::vector<Cube<3>>::const_iterator it = zero_cube_vec.cbegin();
  REQUIRE(complex.grade(*it) == 0);

  // *it = Cube<3>({0, 0, 0}, 0b000)
  // Boundary of vertex is empty
  TestChainType boundary_result;
  REQUIRE(boundary(complex, *it) == boundary_result);
  REQUIRE(graded_boundary(complex, *it) == boundary_result);
  REQUIRE(closure_boundary(complex, *it) == boundary_result);

  ++it;
  REQUIRE(complex.grade(*it) == 0);

  // *it = Cube<3>({0, 0, 0}, 0b001)
  // Boundary of 1-cube is outer vertex - inner vertex
  boundary_result.insert(Cube<3>({1, 0, 0}, 0b000), one<Z<5>>());
  boundary_result.insert(Cube<3>({0, 0, 0}, 0b000), -one<Z<5>>());
  REQUIRE(boundary(complex, *it) == boundary_result);
  REQUIRE(graded_boundary(complex, *it) == boundary_result);
  REQUIRE(closure_boundary(complex, *it) == boundary_result);

  ++it;
  REQUIRE(complex.grade(*it) == 0);
  ++it;
  REQUIRE(complex.grade(*it) == 0);

  // *it = Cube<3>({0, 0, 0}, 0b101)
  // Boundary of this 2-cube is:
  // outer 0b100 - inner 0b100 - (inner 0b001 - outer 0b001)
  boundary_result.clear();
  boundary_result.insert(Cube<3>({0, 0, 1}, 0b001), -one<Z<5>>());
  boundary_result.insert(Cube<3>({1, 0, 0}, 0b100), one<Z<5>>());
  boundary_result.insert(Cube<3>({0, 0, 0}, 0b001), one<Z<5>>());
  boundary_result.insert(Cube<3>({0, 0, 0}, 0b100), -one<Z<5>>());
  REQUIRE(boundary(complex, *it) == boundary_result);
  REQUIRE(graded_boundary(complex, *it) == boundary_result);
  REQUIRE(closure_boundary(complex, *it) == boundary_result);

  // Check boundary of boundary is empty
  REQUIRE(boundary(complex, boundary_result) == TestChainType());

  // 3-cube is not in grading set
  REQUIRE(complex.grade(Cube<3>({0, 0, 0}, 0b111)) == 1);

  // Complex is open along the maximum orthants; hence boundary of this 1-cube
  // is only - inner vertex
  boundary_result.clear();
  boundary_result.insert(Cube<3>({2, 4, 5}, 0b000), -one<Z<5>>());
  REQUIRE(boundary(complex, Cube<3>({2, 4, 5}, 0b001)) == boundary_result);
  REQUIRE(
      graded_boundary(complex, Cube<3>({2, 4, 5}, 0b001)) == boundary_result
  );
  REQUIRE(
      closure_boundary(complex, Cube<3>({2, 4, 5}, 0b001)) == boundary_result
  );

  // Cube<3>({0, 0, 0}, 0b010) is not in grading set, but its inner vertex is;
  // Both are in its closure, however.
  boundary_result.clear();
  boundary_result.insert(Cube<3>({0, 1, 0}, 0b000), one<Z<5>>());
  REQUIRE(
      graded_boundary(complex, Cube<3>({0, 0, 0}, 0b010)) == boundary_result
  );
  boundary_result.insert(Cube<3>({0, 0, 0}, 0b000), -one<Z<5>>());
  REQUIRE(boundary(complex, Cube<3>({0, 0, 0}, 0b010)) == boundary_result);
  REQUIRE(
      closure_boundary(complex, Cube<3>({0, 0, 0}, 0b010)) == boundary_result
  );

  // Coboundary of 3-cube is empty (in 3 dimensional complex).
  TestChainType coboundary_result;
  REQUIRE(coboundary(complex, Cube<3>({0, 0, 0}, 0b111)) == coboundary_result);
  REQUIRE(
      graded_coboundary(complex, Cube<3>({0, 0, 0}, 0b111)) == coboundary_result
  );
  REQUIRE(
      closure_coboundary(complex, Cube<3>({0, 0, 0}, 0b111)) ==
      coboundary_result
  );

  // Complex is closed along minimum orthant; coboundary of least vertex has
  // just outer 1-cubes along each axis.
  // Only axis 1 and 3 are in grading set, though.
  coboundary_result.insert(Cube<3>({0, 0, 0}, 0b001), -one<Z<5>>());
  coboundary_result.insert(Cube<3>({0, 0, 0}, 0b100), -one<Z<5>>());
  REQUIRE(
      graded_coboundary(complex, Cube<3>({0, 0, 0}, 0b000)) == coboundary_result
  );
  REQUIRE(
      closure_coboundary(complex, Cube<3>({0, 0, 0}, 0b000)) ==
      coboundary_result
  );
  coboundary_result.insert(Cube<3>({0, 0, 0}, 0b010), -one<Z<5>>());
  REQUIRE(coboundary(complex, Cube<3>({0, 0, 0}, 0b000)) == coboundary_result);

  // Test linear application of boundary and coboundary operators on a chain.
  // This is a cycle of 1-cubes and hence has boundary zero but nonzero
  // coboundary.
  TestChainType chain;
  chain.insert(Cube<3>({0, 0, 0}, 0b001), one<Z<5>>());
  chain.insert(Cube<3>({1, 0, 0}, 0b100), one<Z<5>>());
  chain.insert(Cube<3>({0, 0, 1}, 0b001), -one<Z<5>>());
  chain.insert(Cube<3>({0, 0, 0}, 0b100), -one<Z<5>>());

  // Cycle has 0 boundary
  boundary_result.clear();
  REQUIRE(boundary(complex, chain) == boundary_result);
  REQUIRE(graded_boundary(complex, chain) == boundary_result);
  REQUIRE(closure_boundary(complex, chain) == boundary_result);

  // Coboundary is more complex; in grade 0, each cell contributes a copy of
  // Cube<3>({0, 0, 0}, 0b101)
  // In the full complex, many more cells are present.
  coboundary_result.clear();
  coboundary_result.insert(Cube<3>({0, 0, 0}, 0b101), Z<5>(4));
  REQUIRE(graded_coboundary(complex, chain) == coboundary_result);
  REQUIRE(closure_coboundary(complex, chain) == coboundary_result);
  coboundary_result.insert(Cube<3>({1, 0, 0}, 0b101), -one<Z<5>>());
  coboundary_result.insert(Cube<3>({0, 0, 1}, 0b101), -one<Z<5>>());
  coboundary_result.insert(Cube<3>({0, 0, 0}, 0b011), one<Z<5>>());
  coboundary_result.insert(Cube<3>({0, 0, 0}, 0b110), one<Z<5>>());
  coboundary_result.insert(Cube<3>({1, 0, 0}, 0b110), -one<Z<5>>());
  coboundary_result.insert(Cube<3>({0, 0, 1}, 0b011), -one<Z<5>>());
  REQUIRE(coboundary(complex, chain) == coboundary_result);

  // Check that coboundary of coboundary is empty
  REQUIRE(coboundary(complex, coboundary_result) == TestChainType());
}


}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
