// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A node in a directed acyclic graph representing a collection of subsplits with their
// corresponding parent-child relationships.
//
// Each node represents a sorted osubsplit, which is stored as a bitset in `subsplit_`.
// The leafward edges are divided into two groups based on if they split apart the
// right clade (in which case they are called sorted children) or if they split apart
// the left clade (in which case they are called rotated children).
// Similarly, the rootward edges are divided into two groups based on if the child of
// the edge splits apart the right clade of the parent (in which case they are called
// sorted parents) or if they split apart the left clade of the parent (in which case
// they are called rotated parents).

#ifndef SRC_SUBSPLIT_DAG_NODE_HPP_
#define SRC_SUBSPLIT_DAG_NODE_HPP_

#include "bitset.hpp"
#include "reindexer.hpp"
#include "sugar.hpp"

class SubsplitDAGNode {
 public:
  SubsplitDAGNode(size_t id, Bitset subsplit)
      : id_(id), subsplit_(std::move(subsplit)) {
    // std::cout << "MAKE_NEW_DAG_NODE => " << id_ << " : " <<
    // subsplit_.SubsplitToString() << std::endl;
  }

  // Compare SubsplitDAGNode's by their ids.
  static int Compare(const SubsplitDAGNode &node_a, const SubsplitDAGNode &node_b);
  // Compare SubsplitDAGNode's by their subsplit representations.
  static int CompareBySubsplit(const SubsplitDAGNode &node_a,
                               const SubsplitDAGNode &node_b);

  // Get Node Id.
  size_t Id() const { return id_; }
  // Set Subsplit representation of node.
  void SetBitset(const Bitset subsplit) { subsplit_ = subsplit; }
  // Get Subsplit representation of node.
  const Bitset &GetBitset() const { return subsplit_; }
  const Bitset GetBitset(bool rotated) const {
    return rotated ? subsplit_.SubsplitRotate() : subsplit_;
  }
  // Check if node is the DAG root.
  bool IsDAGRootNode() const {
    return (rootward_rightside_.empty() && rootward_leftside_.empty());
  }
  // Check if node is a DAG rootsplit.
  bool IsRootsplit() const { return subsplit_.SubsplitIsRootsplit(); }
  // Check if node is a DAG leaf.
  bool IsLeaf() const {
    return leafward_leftside_.empty() && leafward_rightside_.empty();
  }

  // Add edge from given node to adjacent node with specified relationship.
  // is_leafward vs is_rootward, is_leftside vs is_rootward
  void AddEdge(size_t adjacent_node_id, bool is_leafward, bool is_leftside) {
    if (is_leafward) {
      is_leftside ? AddLeafwardLeftside(adjacent_node_id)
                  : AddLeafwardRightside(adjacent_node_id);
    }
    is_leftside ? AddRootwardLeftside(adjacent_node_id)
                : AddRootwardRightside(adjacent_node_id);
  }
  void AddLeafwardLeftside(size_t node_id) { leafward_leftside_.push_back(node_id); }
  void AddLeafwardRightside(size_t node_id) { leafward_rightside_.push_back(node_id); }
  void AddRootwardLeftside(size_t node_id) { rootward_leftside_.push_back(node_id); }
  void AddRootwardRightside(size_t node_id) { rootward_rightside_.push_back(node_id); }

  // Remove edge from given node to adjacent node with specified relationship.
  // is_leafward vs is_rootward, is_leftside vs is_rootward
  void RemoveEdge(size_t adjacent_node_id, bool is_leafward, bool is_leftside) {
    if (is_leafward) {
      is_leftside ? RemoveLeafwardLeftside(adjacent_node_id)
                  : RemoveLeafwardRightside(adjacent_node_id);
    }
    is_leftside ? RemoveRootwardLeftside(adjacent_node_id)
                : RemoveRootwardRightside(adjacent_node_id);
  }
  void RemoveLeafwardLeftside(size_t node_id) { leafward_leftside_.push_back(node_id); }
  void RemoveLeafwardRightside(size_t node_id) {
    leafward_rightside_.push_back(node_id);
  }
  void RemoveRootwardLeftside(size_t node_id) { rootward_leftside_.push_back(node_id); }
  void RemoveRootwardRightside(size_t node_id) {
    rootward_rightside_.push_back(node_id);
  }

  // #350 use enumerated types for rotated?
  // Get vector of all adjacent node vectors along the specified direction.
  const SizeVector GetEdge(bool is_leafward, bool is_leftside) const {
    if (is_leafward) {
      return (is_leftside ? GetLeafwardLeftside() : GetLeafwardRightside());
    } else {
      return (is_leftside ? GetRootwardLeftside() : GetRootwardRightside());
    }
  }
  const SizeVector &GetLeafwardOrRootward(bool is_leafward, bool is_leftside) const {
    return is_leafward ? GetLeafward(is_leftside) : GetRootward(is_leftside);
  };
  const SizeVector &GetLeafwardLeftside() const { return leafward_leftside_; }
  const SizeVector &GetLeafwardRightside() const { return leafward_rightside_; }
  const SizeVector &GetLeafward(bool is_leftside) const {
    return is_leftside ? GetLeafwardLeftside() : GetLeafwardRightside();
  }
  const SizeVector &GetRootwardLeftside() const { return rootward_leftside_; }
  const SizeVector &GetRootwardRightside() const { return rootward_rightside_; }
  const SizeVector &GetRootward(bool is_leftside) const {
    return is_leftside ? GetRootwardLeftside() : GetRootwardRightside();
  }
  void RemapNodeIds(const SizeVector node_reindexer) {
    id_ = node_reindexer.at(id_);
    Reindexer::RemapIdVector(leafward_leftside_, node_reindexer);
    Reindexer::RemapIdVector(leafward_rightside_, node_reindexer);
    Reindexer::RemapIdVector(rootward_leftside_, node_reindexer);
    Reindexer::RemapIdVector(rootward_rightside_, node_reindexer);
  }

  // TODO:
  // Check if given DAGNode is connected to this DAGNode (in any direction).
  bool IsNodeConnected(const size_t node_id) const;
  // Check if node is connected, with relation specified.
  bool IsNodeConnected(const size_t node_id, const bool is_leafward,
                       const bool is_leftside) const;

  // Check if node is in a valid state for the SubsplitDAG.
  // To be valid, a node must have at least one parent, one left child, and one right
  // child.
  bool IsValid() const;
  // Output node to string description.
  // Lists id, subpslit, and all neighboring node ids.
  std::string ToString() const;

 private:
  // Node unique identifier.  Corresponds to the positional index in SubsplitDAG's
  // dag_nodes_.
  size_t id_;
  // Node bitset subsplit clades.
  Bitset subsplit_;

  // List of adjacent nodes in all directions.
  SizeVector leafward_leftside_;
  SizeVector leafward_rightside_;
  SizeVector rootward_leftside_;
  SizeVector rootward_rightside_;
};

#endif  // SRC_SUBSPLIT_DAG_NODE_HPP_
