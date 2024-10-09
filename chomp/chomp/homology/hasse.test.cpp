/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/complexes/complexes.hpp>
#include <chomp/complexes/cubical.hpp>
#include <chomp/complexes/grading.hpp>
#include <chomp/homology/hasse.hpp>
#include <chomp/util/morse.hpp>

#include <catch2/catch_test_macros.hpp>

#include <memory>
#include <set>
#include <unordered_map>
#include <vector>

#ifndef CHOMP_DOXYGEN

namespace chomp::core::detail {

TEST_CASE("Hasse diagram coreductions on cubical complexes", "[homology]") {
  // Declare grading function (forms an open square in orthant (0, 0))
  // and 2D cubical complex.
  const SetGrading<Cube<2>, 0, 1> grading_func(
      {Cube<2>({0, 0}, 0b00), Cube<2>({0, 0}, 0b01), Cube<2>({0, 0}, 0b10),
       Cube<2>({1, 0}, 0b00), Cube<2>({1, 0}, 0b10), Cube<2>({0, 1}, 0b00),
       Cube<2>({0, 1}, 0b01), Cube<2>({1, 1}, 0b00)}
  );
  const std::shared_ptr<CubicalComplex<2, decltype(grading_func), int>>
      complex =
          std::make_shared<CubicalComplex<2, decltype(grading_func), int>>(
              CubeOrthant<2>({1, 1}), grading_func
          );
  using CellType = decltype(complex)::element_type::CellType;
  using ChainType = decltype(complex)::element_type::ChainType;

  // Declare hasse diagram coreduction handler but do not match yet.
  HasseCoreduction<
      CubicalComplex<2, decltype(grading_func), int>, std::unordered_map>
      hasse(complex);

  // Check the leaves after construction are correct.
  std::set<CellType> leaf_cells;
  std::set<CellType> correct(
      {Cube<2>({0, 0}, 0), Cube<2>({1, 0}, 0), Cube<2>({0, 1}, 0),
       Cube<2>({1, 1}, 0), Cube<2>({0, 0}, 3), Cube<2>({1, 0}, 1),
       Cube<2>({0, 1}, 2), Cube<2>({1, 1}, 1), Cube<2>({1, 1}, 2)}
  );
  for (const auto& node : hasse.get_diagram()) {
    leaf_cells.insert(node->cell);
  }
  REQUIRE(leaf_cells == correct);

  // Compute matching and check that critical cells have expected dimensions.
  hasse.match();
  REQUIRE(hasse.get_diagram().empty());

  std::vector<int> dimension_counts(3);
  for (const CellType& critical_cell : hasse.get_critical_cells()) {
    REQUIRE(!hasse.get_matches().contains(critical_cell));
    ++dimension_counts[critical_cell.dimension()];
  }
  REQUIRE(dimension_counts == std::vector<int>({1, 2, 1}));

  auto matches = hasse.get_matches();
  auto priority = hasse.get_priority();


  for (const auto& match_pair : matches) {
    // Check that queens satisfy
    // (notation: match(queen_i) = king_i); assume queen_1 != queen_2
    // queen_2 in boundary of king_1 -> priority(queen_2) < priority(queen_1)
    if (match_pair.second.is_queen()) {
      const ChainType neighbors =
          graded_boundary(*complex, match_pair.second.cell());
      for (const CellType& cell : neighbors) {
        if (cell == match_pair.first) {
          continue;
        }
        auto match_it = matches.find(cell);
        if (match_it != matches.end()) {
          if (match_it->second.is_queen()) {
            REQUIRE(priority[cell] > priority[match_pair.first]);
          }
        }
      }
      // Check that kigns satisfy the dual relation
      // (notation: match(king_i) = queen_i); assume king_1 != king_2
      // king_2 in coboundary of king_1 -> priority(king_2) > priority(king_1)
    } else {
      const ChainType neighbors =
          graded_coboundary(*complex, match_pair.second.cell());
      for (const CellType& cell : neighbors) {
        if (cell == match_pair.first) {
          continue;
        }
        auto match_it = matches.find(cell);
        if (match_it != matches.end()) {
          if (match_it->second.is_king()) {
            REQUIRE(priority[cell] < priority[match_pair.first]);
          }
        }
      }
    }
  }
}

}  // namespace chomp::core::detail

#endif  // CHOMP_DOXYGEN
