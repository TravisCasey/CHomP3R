/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/algebra/algebra.hpp>
#include <chomp/algebra/cyclic.hpp>
#include <chomp/complexes/complexes.hpp>
#include <chomp/complexes/cubical.hpp>
#include <chomp/complexes/grading.hpp>
#include <chomp/homology/morse.hpp>

#include <catch2/catch_test_macros.hpp>

#include <memory>
#include <utility>

#ifndef CHOMP_DOXYGEN

namespace chomp::core {

TEST_CASE("Simple CubicalComplex with CoreductionMatching", "[homology]") {
  const SetGrading<Cube<2>, 0, 1> grading_func(
      {{Cube<2>({0, 0}, 0b00)},
       {Cube<2>({0, 0}, 0b01)},
       {Cube<2>({0, 0}, 0b10)},
       {Cube<2>({1, 0}, 0b00)},
       {Cube<2>({1, 0}, 0b10)},
       {Cube<2>({0, 1}, 0b00)},
       {Cube<2>({0, 1}, 0b01)},
       {Cube<2>({1, 1}, 0b00)}}
  );
  const auto complex =
      std::make_shared<CubicalComplex<2, decltype(grading_func), Z<2>>>(
          CubeOrthant<2>({1, 1}), grading_func
      );

  using CC = decltype(complex)::element_type;
  using RingType = typename CC::RingType;
  using CellType = typename CC::CellType;
  using ChainType = typename CC::ChainType;

  const std::shared_ptr<PartialMatching<CC>> matching =
      std::make_shared<CoreductionMatching<CC>>(complex);

  auto operator_pair = matching->compute_operators();
  auto morse_boundary = operator_pair.first;
  auto morse_coboundary = operator_pair.second;

  SECTION("Computed Morse (co)boundary operators compose to zero.") {
    // Check that the morse (co)boundary operators satisfy requirements for a
    // chain complex; namely, the image lies in the kernel of the next.
    for (const CellType& cell : *matching) {
      ChainType second_boundary;
      ChainType second_coboundary;
      for (const CellType& boundary_cell : morse_boundary[cell]) {
        CHECK(matching->match(boundary_cell).is_ace());
        second_boundary +=
            morse_boundary[cell][boundary_cell] * morse_boundary[boundary_cell];
      }
      for (const CellType& coboundary_cell : morse_coboundary[cell]) {
        CHECK(matching->match(coboundary_cell).is_ace());
        second_coboundary += morse_coboundary[cell][coboundary_cell] *
                             morse_coboundary[coboundary_cell];
      }
      CHECK(second_boundary == ChainType());
      CHECK(second_coboundary == ChainType());
    }
  }

  SECTION("Lower is a chain map.") {
    // That is, it commutes with the boundary operators.
    for (const CellType& cell : *complex) {
      const ChainType lowered = matching->lower(cell);
      ChainType lowered_boundary;
      for (const CellType& lowered_cell : lowered) {
        lowered_boundary +=
            lowered[lowered_cell] * morse_boundary[lowered_cell];
      }
      CHECK(lowered_boundary == matching->lower(boundary(*complex, cell)));
    }
  }
  SECTION("Colower is a chain map.") {
    // That is, it commutes with the boundary operators.
    for (const CellType& cell : *complex) {
      const ChainType colowered = matching->colower(cell);
      ChainType colowered_coboundary;
      for (const CellType& colowered_cell : colowered) {
        colowered_coboundary +=
            colowered[colowered_cell] * morse_coboundary[colowered_cell];
      }
      CHECK(
          colowered_coboundary == matching->colower(coboundary(*complex, cell))
      );
    }
  }
  SECTION("Lift is a chain map.") {
    // That is, it commutes with the boundary operators.
    for (const CellType& cell : *matching) {
      CHECK(
          boundary(*complex, matching->lift(cell)) ==
          matching->lift(morse_boundary[cell])
      );
    }
  }
  SECTION("Colift is a chain map.") {
    // That is, it commutes with the boundary operators.
    for (const CellType& cell : *matching) {
      CHECK(
          coboundary(*complex, matching->colift(cell)) ==
          matching->colift(morse_coboundary[cell])
      );
    }
  }

  SECTION("(Co)lowering a (co)lifted chain is the identity map.") {
    for (const CellType& cell : *matching) {
      ChainType cell_chain;
      cell_chain.insert(cell, one<RingType>());
      CHECK(matching->lower(matching->lift(cell)) == cell_chain);
      CHECK(matching->colower(matching->colift(cell)) == cell_chain);
    }
  }

  SECTION("Computed (co)boundary operators align with technical definition.") {
    // Check computed boundary and coboundary operators in the morse complex
    // align with the definition in:
    // Discrete Morse Theoretic Algorithms for Computing Homology of Complexes
    // and Maps - Harker, Mischaikow, Mrozek, Nanda

    for (const CellType& cell : *matching) {
      CHECK(
          matching->lower(boundary(*complex, matching->lift(cell))) ==
          morse_boundary[cell]
      );
      CHECK(
          matching->colower(coboundary(*complex, matching->colift(cell))) ==
          morse_coboundary[cell]
      );
    }
  }
}

TEST_CASE(
    "More complex CubicalComplex with CoreductionMatching", "[homology]"
) {
  const MapGrading<Cube<3>, 1, 3> grading_func({
      {Cube<3>({0, 0, 0}, 0b000), 1},
      {Cube<3>({0, 0, 0}, 0b010), 1},
      {Cube<3>({0, 0, 0}, 0b100), 1},
      {Cube<3>({0, 1, 0}, 0b000), 1},
      {Cube<3>({0, 1, 0}, 0b010), 1},
      {Cube<3>({0, 1, 0}, 0b100), 1},
      {Cube<3>({0, 2, 0}, 0b000), 1},
      {Cube<3>({0, 2, 0}, 0b100), 1},
      {Cube<3>({0, 0, 1}, 0b000), 1},
      {Cube<3>({0, 0, 1}, 0b010), 1},
      {Cube<3>({0, 1, 1}, 0b000), 1},
      {Cube<3>({0, 1, 1}, 0b010), 1},
      {Cube<3>({0, 2, 1}, 0b000), 1},
      {Cube<3>({0, 0, 0}, 0b001), 1},
      {Cube<3>({0, 1, 0}, 0b001), 1},
      {Cube<3>({0, 2, 0}, 0b001), 1},
      {Cube<3>({0, 0, 1}, 0b001), 1},
      {Cube<3>({0, 1, 1}, 0b001), 1},
      {Cube<3>({0, 2, 1}, 0b001), 1},
      {Cube<3>({1, 0, 0}, 0b000), 1},
      {Cube<3>({1, 1, 0}, 0b000), 1},
      {Cube<3>({1, 2, 0}, 0b000), 1},
      {Cube<3>({1, 0, 1}, 0b000), 1},
      {Cube<3>({1, 1, 1}, 0b000), 1},
      {Cube<3>({1, 2, 1}, 0b000), 1},
      {Cube<3>({0, 0, 0}, 0b110), 2},
      {Cube<3>({0, 1, 0}, 0b110), 2}
  });
  const auto complex =
      std::make_shared<CubicalComplex<3, decltype(grading_func), Z<5>>>(
          CubeOrthant<3>({-1, -1, -1}), CubeOrthant<3>({4, 4, 4}), grading_func
      );

  using CC = decltype(complex)::element_type;
  using RingType = typename CC::RingType;
  using CellType = typename CC::CellType;
  using ChainType = typename CC::ChainType;

  const std::shared_ptr<PartialMatching<CC>> matching =
      std::make_shared<CoreductionMatching<CC>>(complex);

  auto operator_pair = matching->compute_operators();
  auto morse_boundary = operator_pair.first;
  auto morse_coboundary = operator_pair.second;

  SECTION("Computed Morse (co)boundary operators compose to zero.") {
    // Check that the morse (co)boundary operators satisfy requirements for a
    // chain complex; namely, the image lies in the kernel of the next.
    for (const CellType& cell : *matching) {
      ChainType second_boundary;
      ChainType second_coboundary;
      for (const CellType& boundary_cell : morse_boundary[cell]) {
        CHECK(matching->match(boundary_cell).is_ace());
        second_boundary +=
            morse_boundary[cell][boundary_cell] * morse_boundary[boundary_cell];
      }
      for (const CellType& coboundary_cell : morse_coboundary[cell]) {
        CHECK(matching->match(coboundary_cell).is_ace());
        second_coboundary += morse_coboundary[cell][coboundary_cell] *
                             morse_coboundary[coboundary_cell];
      }
      CHECK(second_boundary == ChainType());
      CHECK(second_coboundary == ChainType());
    }
  }

  SECTION("Lower is a chain map.") {
    // That is, it commutes with the boundary operators.
    for (const CellType& cell : *complex) {
      const ChainType lowered = matching->lower(cell);
      ChainType lowered_boundary;
      for (const CellType& lowered_cell : lowered) {
        lowered_boundary +=
            lowered[lowered_cell] * morse_boundary[lowered_cell];
      }
      CHECK(lowered_boundary == matching->lower(boundary(*complex, cell)));
    }
  }
  SECTION("Colower is a chain map.") {
    // That is, it commutes with the boundary operators.
    for (const CellType& cell : *complex) {
      const ChainType colowered = matching->colower(cell);
      ChainType colowered_coboundary;
      for (const CellType& colowered_cell : colowered) {
        colowered_coboundary +=
            colowered[colowered_cell] * morse_coboundary[colowered_cell];
      }
      CHECK(
          colowered_coboundary == matching->colower(coboundary(*complex, cell))
      );
    }
  }
  SECTION("Lift is a chain map.") {
    // That is, it commutes with the boundary operators.
    for (const CellType& cell : *matching) {
      CHECK(
          boundary(*complex, matching->lift(cell)) ==
          matching->lift(morse_boundary[cell])
      );
    }
  }
  SECTION("Colift is a chain map.") {
    // That is, it commutes with the boundary operators.
    for (const CellType& cell : *matching) {
      CHECK(
          coboundary(*complex, matching->colift(cell)) ==
          matching->colift(morse_coboundary[cell])
      );
    }
  }

  SECTION("(Co)lowering a (co)lifted chain is the identity map.") {
    for (const CellType& cell : *matching) {
      ChainType cell_chain;
      cell_chain.insert(cell, one<RingType>());
      CHECK(matching->lower(matching->lift(cell)) == cell_chain);
      CHECK(matching->colower(matching->colift(cell)) == cell_chain);
    }
  }

  SECTION("Computed (co)boundary operators align with technical definition.") {
    // Check computed boundary and coboundary operators in the morse complex
    // align with the definition in:
    // Discrete Morse Theoretic Algorithms for Computing Homology of Complexes
    // and Maps - Harker, Mischaikow, Mrozek, Nanda

    for (const CellType& cell : *matching) {
      CHECK(
          matching->lower(boundary(*complex, matching->lift(cell))) ==
          morse_boundary[cell]
      );
      CHECK(
          matching->colower(coboundary(*complex, matching->colift(cell))) ==
          morse_coboundary[cell]
      );
    }
  }
}

}  // namespace chomp::core

#endif  // CHOMP_DOXYGEN
