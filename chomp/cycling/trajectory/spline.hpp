/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains the `CubicSpline` class template, implementing
 * cubic spline interpolation through points in a variable number of dimensions.
 * This class template is used to ensure connectedness of the cubical
 * trajectory in the computation of cycling signatures.
 */

#ifndef CHOMP_CYCLING_TRAJECTORY_SPLINE_HPP
#define CHOMP_CYCLING_TRAJECTORY_SPLINE_HPP

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <eigen3/Eigen/Core>
#include <format>
#include <iterator>
#include <numeric>
#include <stdexcept>

namespace chomp::cycling {

/**
 * @brief An implementation of cubic spline interpolation for numerical types
 * `T` in `DIM` dimensions.
 *
 * Uses dense Eigen arrays for numerical implementation.
 *
 * @tparam T The numerical type of the data.
 * @tparam DIM The dimension of the data; the implementation creates a spline
 * curve for each dimension.
 */
template <typename T, std::size_t DIM>
class CubicSpline {
private:
  Eigen::Array<T, 1, Eigen::Dynamic> inputs;
  Eigen::Array<T, DIM, Eigen::Dynamic> prev_cube_coeffs;
  Eigen::Array<T, DIM, Eigen::Dynamic> next_cube_coeffs;
  Eigen::Array<T, DIM, Eigen::Dynamic> prev_coeffs;
  Eigen::Array<T, DIM, Eigen::Dynamic> next_coeffs;

  void compute_spline(const Eigen::Array<T, DIM, Eigen::Dynamic>& outputs) {
    using Eigen::Array, Eigen::Dynamic, Eigen::Index, Eigen::all, Eigen::last,
        Eigen::seq;

    // Linear-time and space implementation of natural cubic splines. See the
    // following for inspiration, though the implementation is not identical.
    // https://cse.engineering.nyu.edu/~mleung/CS3734/s03/ch07/cubicSpline.pdf

    const Array<T, 1, Dynamic> h
        = inputs(seq(1, last)) - inputs(seq(0, last - 1));
    Array<T, 1, Dynamic> v = 2 * (h(seq(1, last)) + h(seq(0, last - 1)));
    const Array<T, DIM, Dynamic> b
        = (outputs(all, seq(1, last)) - outputs(all, seq(0, last - 1)))
              .rowwise()
          / h;
    Array<T, DIM, Dynamic> u
        = 6 * (b(all, seq(1, last)) - b(all, seq(0, last - 1)));
    Array<T, DIM, Dynamic> z = Array<T, DIM, Dynamic>::Zero(DIM, inputs.cols());

    for (Index i = 0; i < inputs.cols() - 3; ++i) {
      v(i + 1) -= h(i + 1) * h(i + 1) / v(i);
      u(all, i + 1) -= u(all, i) * h(i + 1) / v(i);
    }

    for (Index i = inputs.cols() - 2; i > 0; --i) {
      z(all, i) = (u(all, i - 1) - h(i) * z(all, i + 1)) / v(i - 1);
    }

    prev_cube_coeffs = (z(all, seq(1, last)).rowwise() / h) / 6;
    next_cube_coeffs = (z(all, seq(0, last - 1)).rowwise() / h) / 6;
    prev_coeffs = outputs(all, seq(1, last)).rowwise() / h
                  - z(all, seq(1, last)).rowwise() * h / 6;
    next_coeffs = outputs(all, seq(0, last - 1)).rowwise() / h
                  - z(all, seq(0, last - 1)).rowwise() * h / 6;
  }

public:
  /**
   * @brief Compute natural cubic spline interpolation in `DIM` dimensions.
   *
   * The construction of the spline is linear in the number of input points in
   * both time and memory.
   *
   * The computed spline can be queried at intermediate points using the
   * `operator()` method, at various orders.
   *
   * @param input_begin Beginning of the iterator yielding input ("time")
   * values. The yielded inputs are expected to be monotonically increasing.
   * There should be exactly `DIM` output values for each input value.
   * @param input_end The end sentinel for `input_begin`.
   *
   * @param output_begin Beginning of the iterator yielding output values,
   * with each `DIM` values corresponding to one input value from `input_begin`.
   * @param output_end The end sentinal for `output_end`.
   *
   * @exception std::invalid_argument If the iterator distance is negative or if
   * there are not exactly `DIM` output values for each input value.
   */
  template <std::forward_iterator IN, std::forward_iterator OUT>
  requires std::convertible_to<typename std::iterator_traits<IN>::value_type, T>
           && std::convertible_to<
               typename std::iterator_traits<OUT>::value_type, T>
  CubicSpline(IN input_begin, IN input_end, OUT output_begin, OUT output_end) {
    using Eigen::Array, Eigen::Dynamic;
    using InDiffType = typename std::iterator_traits<IN>::difference_type;
    using OutDiffType = typename std::iterator_traits<OUT>::difference_type;

    InDiffType input_distance = std::distance(input_begin, input_end);
    OutDiffType output_distance = std::distance(output_begin, output_end);
    if (output_distance <= 0 || input_distance <= 0) {
      throw std::invalid_argument(
          std::format("The length of the input range ({}) and the length of "
                      "the output range ({}) must be positive.",
              input_distance, output_distance));
    }
    std::size_t input_count = input_distance;
    std::size_t output_count = output_distance;
    if (output_count != input_count * DIM) {
      throw std::invalid_argument(std::format(
          "The length of the output range ({}) must equal the product of the "
          "DIM template parameter ({}) and the length of the input range ({}).",
          output_count, DIM, input_count));
    }

    inputs.resize(1, input_count);
    std::copy(input_begin, input_end, inputs.begin());
    Array<T, DIM, Dynamic> outputs(DIM, input_count);
    std::copy(output_begin, output_end, outputs.reshaped().begin());
    compute_spline(outputs);
  }

  /**
   * @brief Compute natural cubic spline interpolation in `DIM` dimensions with
   * implicit inputs (0, 1, 2, ...) and provided output values.
   *
   * The construction of the spline is linear in the number of input points in
   * both time and memory.
   *
   * The computed spline can be queried at intermediate points using the
   * `operator()` method, at various orders.
   *
   * @param output_begin Beginning of the iterator yielding output values.
   * @param output_end The end sentinal for `output_end`.
   *
   * @exception std::invalid_argument if the iterator distance is not a positive
   * integer multiple of `DIM`.
   */
  template <std::forward_iterator OUT>
  requires std::convertible_to<typename std::iterator_traits<OUT>::value_type,
      T>
  CubicSpline(OUT output_begin, OUT output_end) {
    using Eigen::Array, Eigen::Dynamic;
    using OutDiffType = typename std::iterator_traits<OUT>::difference_type;

    OutDiffType output_distance = std::distance(output_begin, output_end);
    if (output_distance <= 0) {
      throw std::invalid_argument(
          std::format("The length of the output range ({}) must be positive.",
              output_distance));
    }
    std::size_t output_count = output_distance;
    std::size_t input_count = output_count / DIM;
    if (output_count % DIM != 0) {
      throw std::invalid_argument(
          std::format("The length of the output range ({}) must be a positive "
                      "integer multiple of the DIM template parameter ({}).",
              output_count, DIM));
    }

    inputs.resize(1, input_count);
    std::iota(inputs.begin(), inputs.end(), static_cast<T>(0));
    Array<T, DIM, Dynamic> outputs(DIM, input_count);
    std::copy(output_begin, output_end, outputs.reshaped().begin());
    compute_spline(outputs);
  }

  /**
   * @brief Construct the cubic spline with the output values in an Eigen array.
   *
   * The provided values in the other constructors are used to populate an
   * Eigen array anyway; if your data is already in this form, this convenience
   * constructor skips that step and accepts the data in an array.
   */
  CubicSpline(const Eigen::Array<T, DIM, Eigen::Dynamic>& outputs) {
    if (outputs.cols() == 0) {
      throw std::invalid_argument("The output array may not have 0 columns.");
    }
    inputs.resize(1, outputs.cols());
    std::iota(inputs.begin(), inputs.end(), static_cast<T>(0));
    compute_spline(outputs);
  }

  /**
   * @brief Evaluate the interpolating spline curve at a given input value.
   *
   * @param input_value The value at which to evaluate the spline.
   * @param order The order of the derivative to evaluate. Default is 0, which
   * evaluates the spline itself. 1, 2, and 3 evaluate the first, second, and
   * third derivatives, respectively.
   *
   * @exception std::domain_error If `input_value` is not between the least and
   * greatest input values provided in construction of this class template.
   * @return Eigen::Array<T, 1, DIM> The point on the interpolating curve
   * corresponding to `input_value` given as a dense `Eigen` vector of
   * dimension `DIM`.
   */
  [[nodiscard]] Eigen::Array<T, 1, DIM> operator()(const T& input_value,
      const std::size_t order = 0) const {
    using Eigen::Array, Eigen::Dynamic, Eigen::all, Eigen::last;
    typename Array<T, 1, Dynamic>::const_iterator next_it;
    if (input_value == inputs(0)) {
      next_it = inputs.cbegin() + 1;
    } else {
      next_it = std::lower_bound(inputs.cbegin(), inputs.cend(), input_value);
    }

    if (next_it == inputs.cbegin() || next_it == inputs.cend()) {
      throw std::domain_error(std::format(
          "Input value {} is not in the cubic spline domain [{}, {}].",
          input_value, inputs(0), inputs(last)));
    }

    const T prev_diff = input_value - *(next_it - 1);
    const T next_diff = *next_it - input_value;
    const std::size_t col = std::distance(inputs.cbegin(), next_it - 1);

    if (order == 0) {
      return prev_cube_coeffs(all, col) * std::pow(prev_diff, 3)
             + next_cube_coeffs(all, col) * std::pow(next_diff, 3)
             + prev_coeffs(all, col) * prev_diff
             + next_coeffs(all, col) * next_diff;
    }
    if (order == 1) {
      return 3 * prev_cube_coeffs(all, col) * std::pow(prev_diff, 2)
             - 3 * next_cube_coeffs(all, col) * std::pow(next_diff, 2)
             + prev_coeffs(all, col) - next_coeffs(all, col);
    }
    if (order == 2) {
      return 6 * prev_cube_coeffs(all, col) * prev_diff
             + 6 * next_cube_coeffs(all, col) * next_diff;
    }
    if (order == 3) {
      return 6 * prev_cube_coeffs(all, col) - 6 * next_cube_coeffs(all, col);
    }
    return Eigen::Array<T, 1, DIM>::Zero(1, DIM);
  }
};

/**
 * @brief Specialized cubic spline interpolation for points in `DIM`-dimensional
 * sphere bundle.
 *
 * The provided points are the spatial coordinates only, being `DIM / 2`-
 * dimensional. Querying this spline (using the `operator()` method) yields the
 * typical cubic spline point for the first `DIM / 2` dimensions, and a
 * normalized sphere bundle-compliant derivative of the cubic spline for the
 * remaining dimensions.
 */
template <typename T, std::size_t DIM>
requires(DIM % 2 == 0)
class SphereBundleCubicSpline {
private:
  CubicSpline<T, DIM / 2> spline;
  double radius;

public:
  /**
   * @brief Construct the cubic spline. The spatial cubic spline is created as
   * in the related `CubicSpline` constructor with `DIM / 2`.
   *
   * @param sphere_radius Radius of the sphere bundle, and the value to which
   * the maximum norm of the derivative is normalized to when queried via the
   * `operator()` method.
   */
  template <std::forward_iterator IN, std::forward_iterator OUT>
  requires std::convertible_to<typename std::iterator_traits<IN>::value_type, T>
               && std::convertible_to<
                   typename std::iterator_traits<OUT>::value_type, T>
  SphereBundleCubicSpline(IN input_begin, IN input_end, OUT output_begin,
      OUT output_end, const double sphere_radius) :
      spline(input_begin, input_end, output_begin, output_end),
      radius(sphere_radius) {}

  /**
   * @brief Construct the cubic spline. The spatial cubic spline is created as
   * in the related `CubicSpline` constructor with `DIM / 2`.\
   *
   * @param sphere_radius Radius of the sphere bundle, and the value to which
   * the maximum norm of the derivative is normalized to when queried via the
   * `operator()` method.
   */
  template <std::forward_iterator OUT>
  requires std::convertible_to<typename std::iterator_traits<OUT>::value_type,
               T>
  SphereBundleCubicSpline(OUT output_begin, OUT output_end,
      const double sphere_radius) :
      spline(output_begin, output_end), radius(sphere_radius) {}

  /**
   * @brief Construct the cubic spline. The spatial cubic spline is created as
   * in the related `CubicSpline` constructor with `DIM / 2`.
   *
   * @param sphere_radius Radius of the sphere bundle, and the value to which
   * the maximum norm of the derivative is normalized to when queried via the
   * `operator()` method.
   */
  SphereBundleCubicSpline(
      const Eigen::Array<T, DIM / 2, Eigen::Dynamic>& outputs,
      const double sphere_radius) :
      spline(outputs), radius(sphere_radius) {}

  /**
   * @brief Evaluate the interpolating sphere bundle spline curve at a given
   * input value.
   *
   * The first `DIM / 2` values returned are the spatial interpolation, while
   * the remaining values are the first derivative of the spline curve at the
   * input value with maximum norm normalized to the sphere radius provided at
   * construction to correctly lie within the sphere bundle.
   *
   * @param input_value The value at which to evaluate the spline.
   * @param order The order of the derivative to evaluate. Default is 0, which
   * evaluates the spline itself. 1, 2, and 3 evaluate the first, second, and
   * third derivatives, respectively.
   *
   * @exception std::domain_error If `input_value` is not between the least and
   * greatest input values provided in construction of this class template.
   * @return Eigen::Array<T, 1, DIM> The point on the interpolating curve
   * through the sphere bundle corresponding to `input_value` given as a dense
   * `Eigen` vector of dimension `DIM`.
   */
  [[nodiscard]] Eigen::Array<T, 1, DIM> operator()(const T& input_value) const {
    Eigen::Array<T, 1, DIM> result;
    result.template head<DIM / 2>() = spline(input_value, 0);
    Eigen::Array<T, 1, DIM / 2> bundle = spline(input_value, 1);
    result.template tail<DIM / 2>()
        = radius * (bundle / bundle.abs().maxCoeff());
    return result;
  }
};


}  // namespace chomp::cycling

#endif  // CHOMP_CYCLING_TRAJECTORY_SPLINE_HPP
