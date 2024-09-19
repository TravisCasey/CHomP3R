/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains the concepts related to grading functions as well
 * as class templates for common grading function objects.
 */

#ifndef CHOMP_COMPLEXES_GRADING_H
#define CHOMP_COMPLEXES_GRADING_H

#include <chomp/util/cache.hpp>
#include <chomp/util/concepts.hpp>
#include <chomp/util/constants.hpp>

#include <concepts>
#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <utility>

namespace chomp::core {

/**
 * @brief The requirements on grading function objects, namely that they take
 * some const-qualified cell type (which must be the member `G::InputType`) and
 * return `GradingResultType`.
 *
 * An important requirement not checked is that the grading functions for chain
 * complexes must respect the face poset of its corresponding complex. In other
 * words, if cell `a` is a face of cell `b`, then the grading function `g` must
 * satisfy `g(a) <= g(b)`.
 *
 * @tparam G Grading function object.
 */
template <typename G>
concept Grading = requires(G g, const typename G::InputType a) {
  typename G::InputType;
  { g(a) } -> std::convertible_to<GradingResultType>;
};

/**
 * @brief A grading function object that has a lower bound on its output values.
 *
 * This is implemented as `G::Minimum` being `std::integral_constant` of type
 * `GradingResultType` with the minimal value.
 *
 * Complexes with grading functions modeling `LowerBoundedGrading` admit some
 * specialized optimizations.
 *
 * @tparam G Grading function object.
 */
template <typename G>
concept LowerBoundedGrading =
    Grading<G> &&
    std::same_as<
        typename G::Minimum,
        std::integral_constant<GradingResultType, G::Minimum::value>>;
/**
 * @brief A grading function object that has an upper bound on its output
 * values.
 *
 * This is implemented as `G::Maximum` being `std::integral_constant` of type
 * `GradingResultType` with the maximal value.
 *
 * Complexes with grading functions modeling `UpperBoundedGrading` admit some
 * specialized optimizations.
 *
 * @tparam G Grading function object.
 */
template <typename G>
concept UpperBoundedGrading =
    Grading<G> &&
    std::same_as<
        typename G::Maximum,
        std::integral_constant<GradingResultType, G::Maximum::value>>;
/**
 * @brief A grading function object that has both a lower and an upper bound on
 * its output values.
 *
 * This is implemented as `G::Minimum` and `G::Maximum` being
 * `std::integral_constant` of type `GradingResultType` with the corresponding
 * minimal and maximal values.
 *
 * Complexes with grading functions modeling `BoundedGrading` admit some
 * specialized optimizations.
 *
 * @tparam G Grading function object.
 */
template <typename G>
concept BoundedGrading =
    Grading<G> && LowerBoundedGrading<G> && UpperBoundedGrading<G>;

/**
 * @brief A function object modeling `BoundedGrading` based on a map from the
 * input type `T` and the output type `GradingResultType`.
 *
 * @tparam T Input type to the function object. Must either have comparison
 * operators implemented or `std::hash` specialized so that it can be used as a
 * key type in the associative containers `std::unordered_map` and `std::map`.
 * @tparam MIN The minimal value in the output of the map.
 * @tparam MAX The maximal value in the output of the map. If the input to the
 * function is not in the map, returns this value.
 * @tparam MapType The underlying map data structure. Expected to be either
 * `std::map` or `std::unordered_map`, but any associative container with
 * sufficiently similar interface can work. The default type is
 * `std::unordered_map` if `T` is hashable and `std::map` otherwise.
 */
template <
    AssociativeKey T, GradingResultType MIN, GradingResultType MAX,
    template <typename...> typename MapType = DefaultMap>
class MapGrading {
private:
  MapType<T, GradingResultType> grading_map;

public:
  /** @brief The type expected as input to the call operator. */
  using InputType = T;
  /** @brief The minimal output value in the map. */
  using Minimum = std::integral_constant<GradingResultType, MIN>;
  /** @brief The maximal output value in the map. */
  using Maximum = std::integral_constant<GradingResultType, MAX>;

  /** @brief Initialize the grading by explicitly providing the map. */
  MapGrading(MapType<T, GradingResultType> grading_map) :
      grading_map(grading_map) {}
  /** @brief Initialize the grading by providing an initializer list. */
  MapGrading(std::initializer_list<std::pair<T, GradingResultType>> grading_list
  ) : grading_map(grading_list) {}

  /**
   * @brief Call the function object with `input` and return the grade.
   *
   * If `input` is not in the map, returns the maximal value `MAX`; else,
   * returns the mapped value.
   *
   * @param input
   * @return GradingResultType
   */
  GradingResultType operator()(const InputType& input) const {
    typename MapType<T, GradingResultType>::const_iterator it =
        grading_map.find(input);
    return it != grading_map.cend() ? it->second : MAX;
  }
};

/**
 * @brief A function object modeling `BoundedGrading` based on inclusion of
 * objects of the input type `T` in a set`.
 *
 * @tparam T Input type to the function object. Must either have comparison
 * operators implemented or `std::hash` specialized so that it can be used as a
 * key type in the associative containers `std::unordered_set` and `std::set`.
 * @tparam MIN The output value if the input object is in the set. Should
 * satisfy `MIN <= MAX`.
 * @tparam MAX The output value if the input object is not in the set. Should
 * satisfy `MIN <= MAX`.
 * @tparam SetType The underlying set data structure. Expected to be either
 * `std::set` or `std::unordered_set` but any set-like container with
 * sufficiently similar interface can work. The default type is
 * `std::unordered_set` if `T` is hashable and `std::set` otherwise.
 */
template <
    AssociativeKey T, GradingResultType MIN, GradingResultType MAX,
    template <typename...> typename SetType = DefaultSet>
class SetGrading {
private:
  SetType<T> grading_set;

public:
  /** @brief The type expected as input to the call operator. */
  using InputType = T;
  /** @brief The output value if the input object is in the set. */
  using Minimum = std::integral_constant<GradingResultType, MIN>;
  /** @brief The output value if the input object is not in the set. */
  using Maximum = std::integral_constant<GradingResultType, MAX>;

  /** @brief Initialize the grading by explicitly providing the set. */
  SetGrading(SetType<T> grading_set) : grading_set(grading_set) {}
  /** @brief Initialize the grading by providing an initializer list. */
  SetGrading(std::initializer_list<T> grading_list) :
      grading_set(grading_list) {}

  /**
   * @brief Call the function object with `input` and return the grade.
   *
   * If `input` is not in the set, returns the maximal value `MAX`; else,
   * returns `MIN`.
   *
   * @param input
   * @return GradingResultType
   */
  GradingResultType operator()(const InputType& input) const {
    return grading_set.contains(input) ? MIN : MAX;
  }
};

/**
 * @brief A wrapper for a grading function object that caches the results with
 * the `LRUCache` container.
 *
 * @tparam G Function object type modeling `Grading`.
 * @tparam MapType Type of map used for the `LRUCache`.
 *
 * @sa `LRUCache`.
 */
template <Grading G, template <typename...> typename MapType = DefaultMap>
class CachedGradingWrapper {
private:
  LRUCache<typename G::InputType, GradingResultType, MapType> cache;

public:
  /** @brief Wrapped function object type. */
  using FunctionType = G;
  /** @brief The type expected as input to the call operator. */
  using InputType = typename FunctionType::InputType;

  /**
   * @brief Initialize the wrapper with the function object to wrap and a
   * maximum number of elements to store in the cache.
   *
   * @param wrapped_function
   * @param cache_max_size
   */
  CachedGradingWrapper(G wrapped_function, std::size_t cache_max_size) :
      cache(wrapped_function, cache_max_size) {}

  /**
   * @brief Call the function object with `input` and return the grade.
   *
   * The result may be from the cache if present. Otherwise, it is added to the
   * cache.
   *
   * @tparam InFor Forwarding type.
   * @param input The input object to the grading call.
   * @return GradingResultType The grade of `input`.
   */
  template <typename InFor>
  requires std::same_as<std::remove_cvref_t<InFor>, InputType>
  GradingResultType operator()(InFor&& input) {
    return cache[std::forward<InFor>(input)];
  }
};

/** @brief Specialization for G modeling `LowerBoundedGrading`. */
template <LowerBoundedGrading G, template <typename...> typename MapType>
class CachedGradingWrapper<G, MapType> {
private:
  LRUCache<typename G::InputType, GradingResultType, MapType> cache;

public:
  /** @brief Wrapped function object type. */
  using FunctionType = G;
  /** @brief The type expected as input to the call operator. */
  using InputType = typename FunctionType::InputType;
  /** @brief The minimal value output by the wrapped function. */
  using Minimum = typename FunctionType::Minimum;

  /**
   * @brief Initialize the wrapper with the function object to wrap and a
   * maximum number of elements to store in the cache.
   *
   * @param wrapped_function
   * @param cache_max_size
   */
  CachedGradingWrapper(G wrapped_function, std::size_t cache_max_size) :
      cache(wrapped_function, cache_max_size) {}

  /**
   * @brief Call the function object with `input` and return the grade.
   *
   * The result may be from the cache if present. Otherwise, it is added to the
   * cache.
   *
   * @tparam InFor Forwarding type.
   * @param input The input object to the grading call.
   * @return GradingResultType The grade of `input`.
   */
  template <typename InFor>
  requires std::same_as<std::remove_cvref_t<InFor>, InputType>
  GradingResultType operator()(InFor&& input) {
    return cache[std::forward<InFor>(input)];
  }
};

/** @brief Specialization for G modeling `UpperBoundedGrading`. */
template <UpperBoundedGrading G, template <typename...> typename MapType>
class CachedGradingWrapper<G, MapType> {
private:
  LRUCache<typename G::InputType, GradingResultType, MapType> cache;

public:
  /** @brief Wrapped function object type. */
  using FunctionType = G;
  /** @brief The type expected as input to the call operator. */
  using InputType = typename FunctionType::InputType;
  /** @brief The maximal value output by the wrapped function. */
  using Maximum = typename FunctionType::Maximum;

  /**
   * @brief Initialize the wrapper with the function object to wrap and a
   * maximum number of elements to store in the cache.
   *
   * @param wrapped_function
   * @param cache_max_size
   */
  CachedGradingWrapper(G wrapped_function, std::size_t cache_max_size) :
      cache(wrapped_function, cache_max_size) {}

  /**
   * @brief Call the function object with `input` and return the grade.
   *
   * The result may be from the cache if present. Otherwise, it is added to the
   * cache.
   *
   * @tparam InFor Forwarding type.
   * @param input The input object to the grading call.
   * @return GradingResultType The grade of `input`.
   */
  template <typename InFor>
  requires std::same_as<std::remove_cvref_t<InFor>, InputType>
  GradingResultType operator()(InFor&& input) {
    return cache[std::forward<InFor>(input)];
  }
};

/** @brief Specialization for G modeling `BoundedGrading`. */
template <BoundedGrading G, template <typename...> typename MapType>
class CachedGradingWrapper<G, MapType> {
private:
  LRUCache<typename G::InputType, GradingResultType, MapType> cache;

public:
  /** @brief Wrapped function object type. */
  using FunctionType = G;
  /** @brief The type expected as input to the call operator. */
  using InputType = typename FunctionType::InputType;
  /** @brief The minimal value output by the wrapped function. */
  using Minimum = typename FunctionType::Minimum;
  /** @brief The maximal value output by the wrapped function. */
  using Maximum = typename FunctionType::Maximum;

  /**
   * @brief Initialize the wrapper with the function object to wrap and a
   * maximum number of elements to store in the cache.
   *
   * @param wrapped_function
   * @param cache_max_size
   */
  CachedGradingWrapper(G wrapped_function, std::size_t cache_max_size) :
      cache(wrapped_function, cache_max_size) {}

  /**
   * @brief Call the function object with `input` and return the grade.
   *
   * The result may be from the cache if present. Otherwise, it is added to the
   * cache.
   *
   * @tparam InFor Forwarding type.
   * @param input The input object to the grading call.
   * @return GradingResultType The grade of `input`.
   */
  template <typename InFor>
  requires std::same_as<std::remove_cvref_t<InFor>, InputType>
  GradingResultType operator()(InFor&& input) {
    return cache[std::forward<InFor>(input)];
  }
};


}  // namespace chomp::core

#endif  // CHOMP_COMPLEXES_GRADING_H
