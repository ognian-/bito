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

bool SubsplitDAGNode::IsValid() const {
  // If node is a leaf, then a valid node should have no parents.
  if (IsLeaf()) {
    return (GetLeafwardRightward().size() + GetLeafwardLeftward().size() == 0);
  }
  // If node is a root, then a valid node should have no children.
  else if (IsDAGRootNode()) {
    return (GetRootwardRightward().size() + GetRootwardLeftward().size() == 0);
  }
  // If neither, then node should either have:
  // (1) Zero parents and zero children.
  // (2) 1+ parents, 1+ sorted children, and 1+ rotated children.
  size_t parent_node_count =
      GetRootwardRightward().size() + GetRootwardLeftward().size();
  if (parent_node_count > 0) {
    if (GetLeafwardRightward().size() == 0 || GetLeafwardRightward().size() == 0) {
      return false;
    }
  } else {
    if (GetLeafwardRightward().size() > 0 || GetLeafwardRightward().size() > 0) {
      return false;
    }
  }
  return true;
}

std::string SubsplitDAGNode::ToString() const {
  std::string str = std::to_string(id_) + ": " + GetBitset().SubsplitToString() + "\n";
  str += "Rootward Rightward: " + GetNeighborString(rootward_rightward_) + "\n";
  str += "Rootward Leftward: " + GetNeighborString(rootward_leftward_) + "\n";
  str += "Leafward Rightward: " + GetNeighborString(leafward_rightward_) + "\n";
  str += "Leafward Leftward: " + GetNeighborString(leafward_leftward_) + "\n";
  return str;
}
