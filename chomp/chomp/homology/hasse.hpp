/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

/** @file
 * @brief This header is not a part of the public API but serves as the matching
 * implementation for `CoreductionMatching`.
 */

#ifndef CHOMP_HOMOLOGY_HASSE_H
#define CHOMP_HOMOLOGY_HASSE_H

#include <chomp/algebra/algebra.hpp>
#include <chomp/complexes/complexes.hpp>

#include <compare>
#include <concepts>
#include <cstddef>
#include <list>
#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#ifndef CHOMP_DOXYGEN

namespace chomp::core::detail {

/*
 * Forms the Hasse diagram of the face partial order (with respect to the
 * grading) and computes coreductions with invertible coefficients.
 *
 * The graph is formed with `Node` instances containing cells of the complex and
 * `Edge` instances containing the boundary/coboundary coefficients connecting
 * them. Smart pointers enable traversal of the graph and manage object
 * lifetimes.
 *
 * The algorithm is to remove all coreductions (a cell which has exactly one
 * cell in its boundary; furthermore we require the coefficient of the boundary
 * be invertible in the coefficient ring), then remove a leaf (no boundary),
 * store it as a critical cell, and repeat until all cells are matched or
 * critical.
 *
 * This is a reformulation of Algorithm 3.6 in `Discrete Morse Theoretic
 * Algorithms for Computing Homology of Complexes and Maps` - Harker,
 * Mischaikow, Mrozek, Nanda.
 */
template <ChainComplex CC, template <typename...> typename MapType>
class HasseCoreduction {
public:
  using RingType = typename CC::RingType;
  using CellType = typename CC::CellType;
  using ChainType = typename CC::ChainType;

private:
  // Forward declarations as these classes refer to each other.
  struct Node;
  struct Edge;

  using NodeListIter =
      typename std::list<std::shared_ptr<Node>>::const_iterator;
  using EdgeListIter =
      typename std::list<std::shared_ptr<Edge>>::const_iterator;

  struct Node {
    std::list<std::shared_ptr<Edge>> coboundary;
    std::list<std::shared_ptr<Edge>> boundary;
    CellType cell;
    bool leaf{false};  // Specifically whether it is already in `leaves`.
    NodeListIter leaves_it;  // Iterator to leaves, if present.

    template <typename TFor>
    requires std::same_as<std::remove_cvref_t<TFor>, CellType>
    Node(TFor&& cell) : cell(std::forward<TFor>(cell)) {}
  };

  struct Edge {
    std::shared_ptr<Node> parent;
    EdgeListIter parent_it;  // List iterator in parent->boundary.
    std::shared_ptr<Node> child;
    EdgeListIter child_it;  // List iterator in parent->coboundary.
    RingType coef;

    template <typename RFor>
    requires std::same_as<std::remove_cvref_t<RFor>, RingType>
    Edge(
        std::shared_ptr<Node> parent, std::shared_ptr<Node> child, RFor&& coef
    ) : parent(parent), child(child), coef(std::forward<RFor>(coef)) {}
  };

  // Pointer to complex being matched. Used for boundary computation.
  std::shared_ptr<CC> upper_complex_ptr;
  // The result of matching; Cell matched to a cell, with a coefficient, and an
  // indicator of whether it matches up or down (in cell dimension).
  MapType<CellType, std::tuple<CellType, RingType, std::strong_ordering>>
      matches;
  // Unmatched cells, i.e. aces.
  std::vector<CellType> critical_cells;
  // Defines the partial order on kings and queens based on order of matching.
  MapType<CellType, std::size_t> priority;
  // Collection of pointers to nodes used when constructing the graph.
  MapType<CellType, std::shared_ptr<Node>> node_pointers;
  // Collection of leaves that can be excised and made critical.
  std::list<std::shared_ptr<Node>> leaves;

  using NodePtrIt = typename decltype(node_pointers)::iterator;

  // Find a (pointer to) the node corresponding to `cell` or create it if it
  // does not exist. If it is created, return `true` as well, else false. Note
  // that the returned pointer may be null if the node has been excised.
  template <typename TFor>
  requires std::same_as<std::remove_cvref_t<TFor>, CellType>
  [[nodiscard]] std::pair<std::shared_ptr<Node>, bool>
  construct_node(TFor&& cell) {
    NodePtrIt cell_it = node_pointers.find(cell);
    // found
    if (cell_it != node_pointers.end()) {
      return std::make_pair(cell_it->second, false);
    }

    // not found; must construct
    std::shared_ptr<Node> new_node =
        std::make_shared<Node>(std::forward<TFor>(cell));
    node_pointers[new_node->cell] = new_node;
    return std::make_pair(new_node, true);
  }

  // Create an edge linking parent and child nodes with given coefficient coef.
  // Required linkings are made within the edge and nodes.
  template <typename RFor>
  requires std::same_as<std::remove_cvref_t<RFor>, RingType>
  void construct_edge(
      std::shared_ptr<Node> parent, std::shared_ptr<Node> child, RFor&& coef
  ) const {
    std::shared_ptr<Edge> new_edge =
        std::make_shared<Edge>(parent, child, std::forward<RFor>(coef));
    parent->boundary.push_front(new_edge);
    new_edge->parent_it = parent->boundary.cbegin();
    child->coboundary.push_front(new_edge);
    new_edge->child_it = child->coboundary.cbegin();
  }

  // Match the node `upper` and its singular boundary element, if the
  // coefficient of the boundary map between them is invertible. Return true,
  // excise the nodes, and match them if it is invertible, else return false
  // and do nothing.
  bool match_pair(std::shared_ptr<Node> upper) {
    // Require invertible coefficient
    RingType coef = upper->boundary.front()->coef;
    if (!invertible(coef)) {
      return false;
    }

    std::shared_ptr<Node> lower = upper->boundary.front()->child;

    if (lower->leaf) {
      leaves.erase(lower->leaves_it);
    }
    excise_node(upper);
    excise_node(lower);

    priority[upper->cell] = priority.size();
    priority[lower->cell] = priority.size();

    CellType upper_cell(std::move(upper->cell));
    CellType lower_cell(std::move(lower->cell));
    matches.insert(std::make_pair(
        upper_cell,
        std::make_tuple(lower_cell, coef, std::strong_ordering::greater)
    ));
    matches.insert(std::make_pair(
        std::move(lower_cell),
        std::make_tuple(std::move(upper_cell), coef, std::strong_ordering::less)
    ));

    return true;
  }

  // Remove `node` from the graph along with all its connections.
  // Any nodes reduced to 0 boundary by this excision are added to `leaves`
  // Any nodes reduced to 1 boundary by this excision are (attempted) matched
  void excise_node(std::shared_ptr<Node> node) {
    std::shared_ptr<Node> node_parent;
    std::shared_ptr<Node> node_child;

    // Any subsequent nodes attemtping to refer to this node will know it has
    // been excised.
    if (!node_pointers.empty()) {
      node_pointers[node->cell] = nullptr;
    }

    // Update child counts; this only affects coboundary so no other actions
    // are required.
    for (std::shared_ptr<Edge> boundary_edge : node->boundary) {
      node_child = boundary_edge->child;
      node_child->coboundary.erase(boundary_edge->child_it);
    }

    // Update parent counts first
    for (std::shared_ptr<Edge> coboundary_edge : node->coboundary) {
      node_parent = coboundary_edge->parent;
      node_parent->boundary.erase(coboundary_edge->parent_it);
      if (node_parent->boundary.empty()) {
        leaves.push_front(node_parent);
        node_parent->leaf = true;
        node_parent->leaves_it = leaves.cbegin();
      }
    }

    // Then attempt to make matches
    for (std::shared_ptr<Edge> coboundary_edge : node->coboundary) {
      node_parent = coboundary_edge->parent;
      if (node_parent->boundary.size() == 1) {
        match_pair(node_parent);
      }
    }
  }

  // Connect a node to the graph by its boundary. May generate boundary nodes if
  // they are not yet present (or excised).
  void connect_node(std::shared_ptr<Node> parent) {
    // Connect to all boundary nodes.
    // Boundary nodes are potentially not fully initialized yet other than their
    // connection to `parent` and their corresponding cell.
    ChainType boundary_chain =
        graded_boundary(*upper_complex_ptr, parent->cell);
    std::shared_ptr<Node> child;
    bool constructed = false;  // whether each node is newly created
    std::list<std::shared_ptr<Node>> additions;  // stores newly created nodes
    for (const CellType& boundary_cell : boundary_chain) {
      std::tie(child, constructed) = construct_node(boundary_cell);
      if (child) {
        construct_edge(parent, child, boundary_chain[boundary_cell]);
        if (constructed) {
          additions.push_front(child);
        }
      }
    }

    // Add to leaves
    if (parent->boundary.empty()) {
      leaves.push_front(parent);
      parent->leaf = true;
      parent->leaves_it = leaves.cbegin();
      return;
    }
    // Attempt to match; if match is successful the node and its child are
    // excised so the function returns.
    if (parent->boundary.size() == 1) {
      if (match_pair(parent)) {
        return;
      }
    }

    // Otherwise, connect all newly created child nodes.
    for (std::shared_ptr<Node> new_node : additions) {
      connect_node(new_node);
    }
  }

public:
  // Create diagram; matches are made in the first pass (i.e. during creation)
  // but are not completed until `match` is called.
  // Note this means that the graph has no available coreductions after
  // construction and a leaf must be excised to continue.
  HasseCoreduction(std::shared_ptr<CC> upper_complex_ptr) :
      upper_complex_ptr(upper_complex_ptr) {
    std::shared_ptr<Node> node;
    bool constructed = false;
    for (const CellType& cell : *upper_complex_ptr) {
      std::tie(node, constructed) = construct_node(cell);
      if (constructed) {
        connect_node(node);
      }
    }
    node_pointers.clear();  // Marks construction as finished.
  }

  // Compute the entirety of the matching.
  void match() {
    std::shared_ptr<Node> leaf;
    while (!leaves.empty()) {
      leaf = leaves.front();
      leaves.pop_front();
      excise_node(leaf);
      critical_cells.push_back(std::move(leaf->cell));
    }
  }

  // Simple exterior interface to get the matching. Returns the matches,
  // priority, and critical cells in a tuple.
  [[nodiscard]] static std::tuple<
      MapType<CellType, std::tuple<CellType, RingType, std::strong_ordering>>,
      MapType<CellType, std::size_t>, std::vector<CellType>>
  compute_matching(std::shared_ptr<CC> complex) {
    HasseCoreduction hasse(complex);
    hasse.match();
    return std::make_tuple(
        std::move(hasse.matches), std::move(hasse.priority),
        std::move(hasse.critical_cells)
    );
  }

  // Used primarily for testing; leaves will be empty after match has been
  // called.
  [[nodiscard]] std::list<std::shared_ptr<Node>> get_diagram() const noexcept {
    return leaves;
  }

  // For simple interfacing prefer using `compute_matching`.
  [[nodiscard]] std::vector<CellType> get_critical_cells() const noexcept {
    return critical_cells;
  }

  // For simple interfacing prefer using `compute_matching`.
  [[nodiscard]] MapType<
      CellType, std::tuple<CellType, RingType, std::strong_ordering>>
  get_matches() const noexcept {
    return matches;
  }

  // For simple interfacing prefer using `compute_matching`.
  [[nodiscard]] MapType<CellType, std::size_t> get_priority() const noexcept {
    return priority;
  }
};

}  // namespace chomp::core::detail

#endif  // CHOMP_DOXYGEN

#endif  // CHOMP_HOMOLOGY_HASSE_H
