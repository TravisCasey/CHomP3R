/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <chomp/cycling/trajectory/spline.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cstddef>
#include <stdexcept>
#include <vector>

#ifndef CHOMP_DOXYGEN

namespace chomp::cycling {

TEST_CASE("CubicSpline class.", "[cycling][trajectory]") {
  std::vector<double> t{0.9, 1.3, 1.9, 2.1};
  std::vector<double> y{1.3, -2.0, 12.0, 1.5, 21.4, 2.0, 1.85, 13.2, 3.4, 2.1,
      -2.1, -20.0};

  std::vector<double> sample_t{1.1, 1.4, 1.95, 2.09};
  std::vector<std::vector<double>> sample_out{{1.4056, 11.343, 4.5634},
      {1.5406, 24.235, 3.7498}, {1.9066, 9.7557, -1.2712},
      {2.0866, -1.2771, -18.651}};

  const CubicSpline<double, 3> spline(t.cbegin(), t.cend(), y.cbegin(),
      y.cend());

  for (std::size_t i = 0; i < t.size(); ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      CHECK_THAT(spline(t[i])(j),
          Catch::Matchers::WithinAbs(y[i * 3 + j], 0.0000001));
    }
  }

  for (std::size_t i = 0; i < sample_t.size(); ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      CHECK_THAT(spline(sample_t[i])(j),
          Catch::Matchers::WithinRel(sample_out[i][j], 0.001));
    }
  }

  CHECK_THROWS_AS(spline(2.2), std::domain_error);
  using ErrorSplineType = CubicSpline<double, 4>;  // Catch2 workaround
  CHECK_THROWS_AS(ErrorSplineType(t.cbegin(), t.cend(), y.cbegin(), y.cend()),
      std::invalid_argument);
}

}  // namespace chomp::cycling

#endif  // CHOMP_DOXYGEN
