/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header contains a LRU cache implementation specialized for memory
 * efficiency in CHomP3R data structures.
 */

#ifndef CHOMP_UTIL_CACHE_H
#define CHOMP_UTIL_CACHE_H

#include <chomp/util/concepts.hpp>

#include <concepts>
#include <cstddef>
#include <functional>
#include <list>
#include <map>
#include <type_traits>
#include <unordered_map>
#include <utility>

namespace chomp::core {

#ifndef CHOMP_DOXYGEN
namespace detail {

template <bool H, typename K, typename V>
struct CacheMapChooser {
  using type = std::map<K, V>;
};

template <typename K, typename V>
struct CacheMapChooser<true, K, V> {
  using type = std::unordered_map<K, V>;
};

template <AssociativeKey K, typename V>
struct DefaultCacheMapChooser {
  using type = typename CacheMapChooser<Hashable<K>, K, V>::type;
};

}  // namespace detail
#endif  // CHOMP_DOXYGEN

/**
 * @brief Default choice between `std::map` and `std::unordered_map` depending
 * on concepts fulfilled by `K`.
 *
 * If `K` has a `std::hash` specialization, `std::unordered_map` is preferred.
 * Otherwise, uses `std::map`.
 *
 * @tparam K Key type.
 * @tparam V Value type.
 */
template <AssociativeKey K, typename V>
using DefaultCacheMap = typename detail::DefaultCacheMapChooser<K, V>::type;


/**
 * @brief A specialized LRU cache data structure designed for memory efficiency
 * of commonly referenced objects.
 *
 * This LRU cache implementation has no `insertion` method, only a get method
 * (through `operator[]`) that will insert the item automatically if it does
 * not exist and return it.
 *
 * @tparam K Key type.
 * @tparam V Value type.
 * @tparam MapType The type of map, which is expected to be either `std::map`,
 * `std::unordered_map`, or another type structurally equivalent.
 */
template <
    AssociativeKey K, typename V,
    template <typename...> typename MapType = DefaultCacheMap>
class LRUCache {
private:
  std::list<std::pair<K, V>> cache_list;

  using ListIterType = typename std::list<std::pair<K, V>>::const_iterator;
  MapType<K, ListIterType> cache_map;
  using MapIterType = typename MapType<K, ListIterType>::iterator;

  std::function<V(const K&)> construct_value;
  std::size_t cache_max_size;

  void recompute_iterators() {
    // When copy constructed/assigned, the iterators in the map point to the
    // old list. This recomputes those iterators to the new list.
    for (ListIterType it = cache_list.cbegin(); it != cache_list.cend(); ++it) {
      cache_map[it->first] = it;
    }
  }

public:
  /** @brief Key type. */
  using KeyType = K;
  /** @brief Value type. */
  using ValueType = V;

  /**
   * @brief Construct the associative container `LRUCache` with a maximum size
   * and function constructing an object of the value type `V` from the key type
   * `K`.
   *
   * The specialization of this `LRUCache` implementation requires the ability
   * to construct new objects from the key when requested.
   *
   * @param construct_value `std::function` object taking a constant reference
   * to the key type and returning a value type object.
   * @param max_size The maximum number of entries in the cache before the least
   * recently accessed entries begin to be removed.
   */
  LRUCache(
      std::function<V(const K&)> construct_value, std::size_t max_size
  ) noexcept : construct_value(construct_value), cache_max_size(max_size) {}

  /**
   * @brief Copy constructor.
   *
   * Internal iterators must be recomputed when copy constructed or copy
   * assigned. For this reason, it's best to use move semantics when possible.
   *
   * @param other
   */
  LRUCache(const LRUCache& other) noexcept :
      cache_list(other.cache_list), cache_map(other.cache_map),
      construct_value(other.construct_value),
      cache_max_size(other.cache_max_size) {
    recompute_iterators();
  }

  /**
   * @brief Move constructor.
   *
   * @param other
   */
  LRUCache(LRUCache&& other) noexcept = default;

  /**
   * @brief Copy assignment operator.
   *
   * Internal iterators must be recomputed when copy constructed or copy
   * assigned. For this reason, it's best to use move semantics when possible.
   *
   * @param other
   * @return LRUCache&
   */
  LRUCache& operator=(const LRUCache& other) {
    cache_list = other.cache_list;
    cache_map = other.cache_map;
    construct_value = other.construct_value;
    cache_max_size = other.cache_max_size;
    recompute_iterators();
    return *this;
  }

  /**
   * @brief Move assignment operator.
   *
   * @param other
   * @return LRUCache&
   */
  LRUCache& operator=(LRUCache&& other) = default;

  /**
   * @brief Deconstructor deconstructs the stored objects and deallocates
   * storage.
   */
  ~LRUCache() = default;

  /**
   * @brief Return a constant reference to the value corresponding to `key` by
   * either accessing it in the cache or by constructing the object and adding
   * it to the cache first.
   *
   * Uses forwarding semantics for handling reference types of `K` as the key
   * parameter.
   *
   * @tparam TFor Forwarding type.
   * @param key Possibly `const`-qualified lvalue or rvalue reference to an
   * object of key type `K`.
   * @return const V& Constant reference to the corresponding value type object.
   */
  template <typename TFor>
  requires std::same_as<std::remove_cvref_t<TFor>, K>
  const V& operator[](TFor&& key) {
    MapIterType find_result = cache_map.find(key);

    // Value is not present in the cache.
    // Construct it, add to the cache, and return (a reference to) the object.
    if (find_result == cache_map.end()) {
      cache_list.push_front(std::make_pair(key, construct_value(key)));
      cache_map[std::forward<TFor>(key)] = cache_list.cbegin();
      // Remove from end of cache if necessary.
      if (cache_map.size() > cache_max_size) {
        ListIterType it = cache_list.cend();
        it--;
        cache_map.erase(it->first);
        cache_list.pop_back();
      }
      return cache_list.front().second;
    }

    // Value is present in the cache.
    // Move the list node corresponding to this key to the front and return.
    cache_list.splice(cache_list.cbegin(), cache_list, find_result->second);
    return find_result->second->second;
  }

  /**
   * @brief Query whether the cache contains `key`.
   *
   * @param key
   * @return true
   * @return false
   */
  [[nodiscard]] bool contains(const K& key) const {
    return cache_map.contains(key);
  }

  /**
   * @brief Get the current size of the container.
   *
   * @return std::size_t Current number of elements in the cache.
   */
  [[nodiscard]] std::size_t size() const noexcept {
    return cache_map.size();
  }
  /**
   * @brief Get the maximum size of the container.
   *
   * @return std::size_t Maximum number of elements in the cache.
   */
  [[nodiscard]] std::size_t max_size() const noexcept {
    return cache_max_size;
  }
};

}  // namespace chomp::core

#endif  // CHOMP_UTIL_CACHE_H
