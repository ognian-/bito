// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag_graft.hpp"

#include "subsplit_dag.hpp"

// ** Constructors

SubsplitDAGGraft::SubsplitDAGGraft(SubsplitDAG &dag) : host_dag_(&dag) {}

// ** Comparators

int SubsplitDAGGraft::Compare(SubsplitDAGGraft &dag_a, SubsplitDAGGraft &dag_b) {}

int SubsplitDAGGraft::Compare(SubsplitDAGGraft &dag_a, SubsplitDAG &dag_b) {}

// ** Modify DAG

void SubsplitDAGGraft::CreateAndInsertNode(const Bitset &node_subsplit) {}

void SubsplitDAGGraft::CreateAndInsertEdge(const Bitset &parent_subsplit,
                                           const Bitset &child_subsplit) {}

void SubsplitDAGGraft::AddNodePair(const Bitset &parent_subsplit,
                                   const Bitset &child_subsplit) {
  // Assert that parent and child are a valid node pair to add to DAG.
  // host_dag_->IsValidAddNodePair(parent_subsplit, child_subsplit);
  // Iterate through dag nodes to find adjacent nodes.
  // host_dag_->IterateOverRealNodes([](SubsplitDAGNode* node){});
}

void SubsplitDAGGraft::RemoveNodePair(const Bitset &parent_subsplit,
                                      const Bitset &child_subsplit) {}

void SubsplitDAGGraft::Clear() {
  graft_nodes_.clear();
  graft_edges_.clear();
}

// ** Getters

size_t SubsplitDAGGraft::GetDAGNodeId(const Bitset &node_subsplit) const {
  // Check if subsplit is in the host DAG.
  // Check if subsplit is in the graft.
}

// ** Counts

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
