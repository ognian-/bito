// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag_node.hpp"

int SubsplitDAGNode::Compare(const SubsplitDAGNode &node_a,
                             const SubsplitDAGNode &node_b) {
  return node_a.Id() - node_b.Id();
}

int CompareBySubsplit(const SubsplitDAGNode &node_a, const SubsplitDAGNode &node_b) {
  return Bitset::SubsplitCompare(node_a.GetBitset(), node_b.GetBitset());
}

std::string GetNeighborString(SizeVector neighbors) {
  std::string str;
  for (size_t i : neighbors) {
    str += std::to_string(i) + " ";
  }
  return str;
}

bool SubsplitDAGNode::IsNodeConnected(const size_t node_id) const {
  for (const bool is_leafward : {true, false}) {
    for (const bool is_leftside : {true, false}) {
      if (IsNodeConnected(node_id, is_leafward, is_leftside)) {
        return true;
      }
    }
  }
  return false;
};

bool SubsplitDAGNode::IsNodeConnected(const size_t node_id, const bool is_leafward,
                                      const bool is_leftside) const {
  for (const auto &adj_node_id : GetEdge(is_leafward, is_leftside)) {
    if (adj_node_id == node_id) {
      return true;
    }
  }
  return false;
};

bool SubsplitDAGNode::IsValid() const {
  // If node is a leaf, then a valid node should have no parents.
  if (IsLeaf()) {
    return (GetLeafwardRightside().size() + GetLeafwardLeftside().size() == 0);
  }
  // If node is a root, then a valid node should have no children.
  else if (IsDAGRootNode()) {
    return (GetRootwardRightside().size() + GetRootwardLeftside().size() == 0);
  }
  // If neither, then node should either have:
  // (1) Zero parents and zero children.
  // (2) 1+ parents, 1+ sorted children, and 1+ rotated children.
  size_t parent_node_count =
      GetRootwardRightside().size() + GetRootwardLeftside().size();
  if (parent_node_count > 0) {
    if (GetLeafwardRightside().size() == 0 || GetLeafwardRightside().size() == 0) {
      return false;
    }
  } else {
    if (GetLeafwardRightside().size() > 0 || GetLeafwardRightside().size() > 0) {
      return false;
    }
  }
  return true;
}

std::string SubsplitDAGNode::ToString() const {
  std::string str = std::to_string(id_) + ": " + GetBitset().SubsplitToString() + "\n";
  str += "Rootward Rightside: " + GetNeighborString(rootward_rightside_) + "\n";
  str += "Rootward Leftside: " + GetNeighborString(rootward_leftside_) + "\n";
  str += "Leafward Rightside: " + GetNeighborString(leafward_rightside_) + "\n";
  str += "Leafward Leftside: " + GetNeighborString(leafward_leftside_) + "\n";
  return str;
}
