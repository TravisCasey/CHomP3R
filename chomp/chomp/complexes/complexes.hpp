/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains the `ChainComplex` concept defining the minimal
 * requirements on a class to implement a chain complex. It also includes
 * shared functionality on grading and boundary/coboundary operators by
 * free function templates general to all types modeling `ChainComplex`.
 */


#ifndef CHOMP_COMPLEXES_COMPLEXES_H
#define CHOMP_COMPLEXES_COMPLEXES_H

#include <chomp/algebra/algebra.hpp>
#include <chomp/complexes/grading.hpp>
#include <chomp/util/constants.hpp>

#include <concepts>
#include <functional>

namespace chomp::core {

/**
 * @brief Alias for the function object type used as test conditions in the
 * `boundary_if` and `coboundary_if` free functions and methods in complex
 * classes.
 *
 * @tparam C Cell type of the associated chain complex class.
 */
template <typename C>
using ConditionalType = std::function<bool(const C&)>;

/**
 * @brief Requirements on a class to implement a general chain complex.
 *
 * The `boundary`, `graded_boundary`, `closure_boundary`, `coboundary`
 * `graded_coboundary`, and `closure_coboundary` free functions are synthesized
 * from the required `boundary_if` and `coboundary_if` methods.
 *
 * @tparam CC
 */
template <typename CC>
concept ChainComplex = requires(
    CC complex, const typename CC::CellType cell,
    const ConditionalType<typename CC::CellType> cond
) {
  Ring<typename CC::RingType>;
  Basis<typename CC::CellType>;
  Module<typename CC::ChainType>;
  Grading<typename CC::GradingType>;

  { complex.grade(cell) } -> std::convertible_to<GradingResultType>;
  {
    complex.boundary_if(cell, cond)
  } -> std::convertible_to<typename CC::ChainType>;
  {
    complex.coboundary_if(cell, cond)
  } -> std::convertible_to<typename CC::ChainType>;
};

/**
 * @brief Free function calling the `grade` method of classes modeling
 * `ChainComplex` similar to the `boundary` and `coboundary` free functions.
 *
 * @tparam CC Class modeling `ChainComlex`.
 * @param complex The chain complex to which `cell` belongs and in which the
 * grading is done.
 * @param cell The cell to be graded.
 * @return GradingResultType The grade of `cell` in `complex`.
 */
template <ChainComplex CC>
[[nodiscard]] inline GradingResultType
grade(CC& complex, const typename CC::CellType& cell) {
  return complex.grade(cell);
}

/**
 * @brief The boundary of `cell` in the chain complex `complex` subject to the
 * condition `cond`.
 *
 * @tparam CC Class modeling `ChainComplex`.
 * @param complex The chain complex to which `cell` belongs and in which the
 * boundary is taken.
 * @param cell The cell for which the boundary is computed.
 * @param cond Conditional statement on proposed boundary cells; the cell is
 * in the resulting boundary if and only if this condition returns `true`.
 * @return CC::ChainType The boundary chain of `cell` in `complex`.
 */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType boundary_if(
    CC& complex, const typename CC::CellType& cell,
    const ConditionalType<typename CC::CellType>& cond
) {
  return complex.boundary_if(cell, cond);
}
/** @brief Overload of `boundary_if` linearly applied on `chain`. */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType boundary_if(
    CC& complex, const typename CC::ChainType& chain,
    const ConditionalType<typename CC::CellType>& cond
) {
  return linear_apply(
      chain,
      [&complex, &cond](const typename CC::CellType& cell) {
        return boundary_if(complex, cell, cond);
      }
  );
}

/**
 * @brief The boundary of `cell` in the chain complex `complex`.
 *
 * @tparam CC Class modeling `ChainComplex`.
 * @param complex The chain complex to which `cell` belongs and in which the
 * boundary is taken.
 * @param cell The cell for which the boundary is computed.
 * @return CC::ChainType The boundary chain of `cell` in `complex`.
 */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
boundary(CC& complex, const typename CC::CellType& cell) {
  return complex.boundary_if(
      cell,
      []([[maybe_unused]] const typename CC::CellType boundary_cell) {
        return true;
      }
  );
}
/** @brief Overload of `boundary` linearly applied on `chain`. */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
boundary(CC& complex, const typename CC::ChainType& chain) {
  return linear_apply(chain, [&complex](const typename CC::CellType cell) {
    return boundary(complex, cell);
  });
}

/**
 * @brief The boundary of `cell` in the chain complex `complex` subject to the
 * constraint that the boundary cells must have the same grade as `cell`.
 *
 * @tparam CC Class modeling `ChainComplex`.
 * @param complex The chain complex to which `cell` belongs and in which the
 * boundary is taken.
 * @param cell The cell for which the boundary is computed.
 * @return CC::ChainType The boundary chain of `cell` in `complex`; all cells
 * present must have the same grade as `cell`.
 */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
graded_boundary(CC& complex, const typename CC::CellType& cell) {
  GradingResultType current_grade = complex.grade(cell);
  return complex.boundary_if(
      cell,
      [current_grade, &complex](const typename CC::CellType boundary_cell) {
        return complex.grade(boundary_cell) == current_grade;
      }
  );
}
/** @brief Specialization for when `cell` has minimal grade. */
template <typename CC>
requires ChainComplex<CC> && LowerBoundedGrading<typename CC::GradingType>
[[nodiscard]] inline typename CC::ChainType
graded_boundary(CC& complex, const typename CC::CellType& cell) {
  GradingResultType current_grade = complex.grade(cell);
  if (current_grade == CC::GradingType::Minimum::value) {
    return boundary(complex, cell);
  }
  return complex.boundary_if(
      cell,
      [current_grade, &complex](const typename CC::CellType boundary_cell) {
        return complex.grade(boundary_cell) == current_grade;
      }
  );
}
/** @brief Overload of `graded_boundary` linearly applied on `chain`. */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
graded_boundary(CC& complex, const typename CC::ChainType& chain) {
  return linear_apply(chain, [&complex](const typename CC::CellType cell) {
    return graded_boundary(complex, cell);
  });
}

/**
 * @brief The boundary of `cell` in the chain complex `complex` subject to the
 * constraint that the boundary cells must have grade at most that of `cell`.
 *
 * Due to closure requirements on grading functions, this is equivalent to
 * the `boundary` function.
 *
 * @tparam CC Class modeling `ChainComplex`.
 * @param complex The chain complex to which `cell` belongs and in which the
 * boundary is taken.
 * @param cell The cell for which the boundary is computed.
 * @return CC::ChainType The boundary chain of `cell` in `complex`; all cells
 * present must have grade less than or equal to as `cell`.
 */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
closure_boundary(CC& complex, const typename CC::CellType& cell) {
  return boundary(complex, cell);
}
/** @brief Overload of `closure_boundary` linearly applied on `chain`. */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
closure_boundary(CC& complex, const typename CC::ChainType& chain) {
  return linear_apply(chain, [&complex](const typename CC::CellType& cell) {
    return closure_boundary(complex, cell);
  });
}

/**
 * @brief The coboundary of `cell` in the chain complex `complex` subject to the
 * condition `cond`.
 *
 * @tparam CC Class modeling `ChainComplex`.
 * @param complex The chain complex to which `cell` belongs and in which the
 * coboundary is taken.
 * @param cell The cell for which the coboundary is computed.
 * @param cond Conditional statement on proposed coboundary cells; the cell is
 * in the resulting coboundary if and only if this condition returns `true`.
 * @return CC::ChainType The coboundary chain of `cell` in `complex`.
 */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType coboundary_if(
    CC& complex, const typename CC::CellType& cell,
    const ConditionalType<typename CC::CellType>& cond
) {
  return complex.coboundary_if(cell, cond);
}
/** @brief Overload of `coboundary_if` linearly applied on `chain`. */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType coboundary_if(
    CC& complex, const typename CC::ChainType& chain,
    const ConditionalType<typename CC::CellType>& cond
) {
  return linear_apply(
      chain,
      [&complex, &cond](const typename CC::CellType& cell) {
        return coboundary_if(complex, cell, cond);
      }
  );
}

/**
 * @brief The coboundary of `cell` in the chain complex `complex`.
 *
 * @tparam CC Class modeling `ChainComplex`.
 * @param complex The chain complex to which `cell` belongs and in which the
 * coboundary is taken.
 * @param cell The cell for which the coboundary is computed.
 * @return CC::ChainType The coboundary chain of `cell` in `complex`.
 */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
coboundary(CC& complex, const typename CC::CellType& cell) {
  return complex.coboundary_if(
      cell,
      []([[maybe_unused]] const typename CC::CellType coboundary_cell) {
        return true;
      }
  );
}
/** @brief Overload of `coboundary` linearly applied on `chain`. */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
coboundary(CC& complex, const typename CC::ChainType& chain) {
  return linear_apply(chain, [&complex](const typename CC::CellType cell) {
    return coboundary(complex, cell);
  });
}

/**
 * @brief The coboundary of `cell` in the chain complex `complex` subject to the
 * constraint that the coboundary cells must have the same grade as `cell`.
 *
 * @tparam CC Class modeling `ChainComplex`.
 * @param complex The chain complex to which `cell` belongs and in which the
 * coboundary is taken.
 * @param cell The cell for which the coboundary is computed.
 * @return CC::ChainType The coboundary chain of `cell` in `complex`; all cells
 * present must have the same grade as `cell`.
 */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
graded_coboundary(CC& complex, const typename CC::CellType& cell) {
  GradingResultType current_grade = complex.grade(cell);
  return complex.coboundary_if(
      cell,
      [current_grade, &complex](const typename CC::CellType coboundary_cell) {
        return complex.grade(coboundary_cell) == current_grade;
      }
  );
}
/** @brief Specialization for when `cell` has maximal grade. */
template <typename CC>
requires ChainComplex<CC> && UpperBoundedGrading<typename CC::GradingType>
[[nodiscard]] inline typename CC::ChainType
graded_coboundary(CC& complex, const typename CC::CellType& cell) {
  GradingResultType current_grade = complex.grade(cell);
  if (current_grade == CC::GradingType::Maximum::value) {
    return coboundary(complex, cell);
  }
  return complex.coboundary_if(
      cell,
      [current_grade, &complex](const typename CC::CellType coboundary_cell) {
        return complex.grade(coboundary_cell) == current_grade;
      }
  );
}
/** @brief Overload of `graded_coboundary` linearly applied on `chain`. */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
graded_coboundary(CC& complex, const typename CC::ChainType& chain) {
  return linear_apply(chain, [&complex](const typename CC::CellType cell) {
    return graded_coboundary(complex, cell);
  });
}

/**
 * @brief The coboundary of `cell` in the chain complex `complex` subject to the
 * constraint that the coboundary cells must have grade at most that of `cell`.
 *
 * @tparam CC Class modeling `ChainComplex`.
 * @param complex The chain complex to which `cell` belongs and in which the
 * coboundary is taken.
 * @param cell The cell for which the coboundary is computed.
 * @return CC::ChainType The coboundary chain of `cell` in `complex`; all cells
 * present must have grade less than or equal to as `cell`.
 */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
closure_coboundary(CC& complex, const typename CC::CellType& cell) {
  GradingResultType current_grade = complex.grade(cell);
  return complex.coboundary_if(
      cell,
      [current_grade, &complex](const typename CC::CellType coboundary_cell) {
        return complex.grade(coboundary_cell) == current_grade;
      }
  );
}
/** @brief Specialization for when `cell` has maximal grade. */
template <typename CC>
requires ChainComplex<CC> && UpperBoundedGrading<typename CC::GradingType>
[[nodiscard]] inline typename CC::ChainType
closure_coboundary(CC& complex, const typename CC::CellType& cell) {
  GradingResultType current_grade = complex.grade(cell);
  if (current_grade == CC::GradingType::Maximum::value) {
    return coboundary(complex, cell);
  }
  return complex.coboundary_if(
      cell,
      [current_grade, &complex](const typename CC::CellType coboundary_cell) {
        return complex.grade(coboundary_cell) == current_grade;
      }
  );
}
/** @brief Overload of `closure_coboundary` linearly applied on `chain`. */
template <ChainComplex CC>
[[nodiscard]] inline typename CC::ChainType
closure_coboundary(CC& complex, const typename CC::ChainType& chain) {
  return linear_apply(chain, [&complex](const typename CC::CellType& cell) {
    return closure_coboundary(complex, cell);
  });
}

}  // namespace chomp::core

#endif  // CHOMP_COMPLEXES_COMPLEXES_H
