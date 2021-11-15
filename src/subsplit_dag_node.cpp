// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag_node.hpp"

std::string GetNeighborString(SizeVector neighbors) {
  std::string str;
  for (size_t i : neighbors) {
    str += std::to_string(i) + " ";
  }
  return str;
}

// A node
bool SubsplitDAGNode::IsValid() const {
  // If node is a leaf, then a valid node should have no parents.
  if (IsLeaf()) {
    return (GetLeafwardSorted().size() + GetLeafwardRotated().size() == 0);
  }
  // If node is a root, then a valid node should have no children.
  else if (IsDAGRootNode()) {
    return (GetRootwardSorted().size() + GetRootwardRotated().size() == 0);
  }
  // If neither, then node should either have:
  // (1) Zero parents and zero children.
  // (2) 1+ parents, 1+ sorted children, and 1+ rotated children.
  size_t parent_node_count = GetRootwardSorted().size() + GetRootwardRotated().size();
  if (parent_node_count > 0) {
    if (GetLeafwardSorted().size() == 0 || GetLeafwardSorted().size() == 0) {
      return false;
    }
  } else {
    if (GetLeafwardSorted().size() > 0 || GetLeafwardSorted().size() > 0) {
      return false;
    }
  }
  return true;
}

std::string SubsplitDAGNode::ToString() const {
  std::string str = std::to_string(id_) + ": " + GetBitset().SubsplitToString() + "\n";
  str += "Rootward Sorted: " + GetNeighborString(rootward_sorted_) + "\n";
  str += "Rootward Rotated: " + GetNeighborString(rootward_rotated_) + "\n";
  str += "Leafward Sorted: " + GetNeighborString(leafward_sorted_) + "\n";
  str += "Leafward Rotated: " + GetNeighborString(leafward_rotated_) + "\n";
  return str;
}
