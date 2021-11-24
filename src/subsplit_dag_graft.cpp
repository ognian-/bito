// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag_graft.hpp"

#include "subsplit_dag.hpp"

// ** Constructors:

// ** Modify graft methods:

SubsplitDAGGraft::SubsplitDAGGraft(SubsplitDAG &dag) : host_dag_(&dag) {

}

void SubsplitDAGGraft::AddNodePair(Bitset &parent_subsplit, Bitset &child_subsplit) {
  // Iterate through dag nodes to find adjacent nodes.
  host_dag_.IterateOverRealNodes([&parent_subsplit, &child_subsplit](SubsplitDAGNode *node) {

  });
}

void SubsplitDAGGraft::RemoveNodePair(Bitset &parent_subsplit, Bitset &child_subsplit) {
  
}

size_t SubsplitDAGGraft::NodeCount() const {
  return HostNodeCount() + GraftNodeCount();
}

size_t SubsplitDAGGraft::GraftNodeCount() const { return graft_nodes_.size(); }

size_t SubsplitDAGGraft::HostNodeCount() const { return host_dag_->NodeCount(); }

size_t SubsplitDAGGraft::GPCSPCount() const {
  return HostGPCSPCount() + GraftGPCSPCount();
}

size_t SubsplitDAGGraft::GraftGPCSPCount() const {
  return HostGPCSPCount() + GraftGPCSPCount();
}

size_t SubsplitDAGGraft::HostGPCSPCount() const {
  return HostGPCSPCount() + GraftGPCSPCount();
}
