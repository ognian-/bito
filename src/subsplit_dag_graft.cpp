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

void SubsplitDAGGraft::CreateGraftNode(const Bitset &node_subsplit) {
  // Add node to node list.
  size_t node_id = NodeCount();
  graft_nodes_.push_back(std::make_unique<SubsplitDAGNode>(node_id, subsplit));
  SafeInsert(subsplit_)
}

void SubsplitDAGGraft::CreateGraftEdge(const Bitset &parent_subsplit,
                                      const Bitset &child_subsplit) {
  // Add edge to edge list.

}

void SubsplitDAGGraft::DestroyGraftNode(const Bitset &node_subsplit) {
  // Remove node from node list.

}

void SubsplitDAGGraft::DestroyGraftEdge(const Bitset &parent_subsplit,
                                                const Bitset &child_subsplit) {
  // Remove edge from edge list.

}

void SubsplitDAGGraft::AddGraftNodePair(const Bitset &parent_subsplit,
                                   const Bitset &child_subsplit) {
  // Assert that parent and child are a valid node pair to add to DAG.
  if (host_dag_->IsValidAddNodePair(parent_subsplit, child_subsplit) == false) {
    
  }
  // Iterate through dag nodes to find adjacent nodes.
  // host_dag_->IterateOverRealNodes([](SubsplitDAGNode* node){});
}

void SubsplitDAGGraft::RemoveGraftNodePair(const Bitset &parent_subsplit,
                                           const Bitset &child_subsplit) {
  // Find nodes in node list.
  // Remove graft edges
}

void SubsplitDAGGraft::RemoveAllGrafts() {
  graft_nodes_.clear();
  graft_edges_.clear();
}

// ** Getters

SubsplitDAGNode* SubsplitDAGGraft::GetNode(const size_t node_id) const {
  // Check if subsplit is in the host DAG.
  // Check if subsplit is in the graft.
}

size_t SubsplitDAGGraft::GetNodeId(const Bitset &node_subsplit) const {
  // Check if subsplit is in the host DAG.
  if (host_dag_->ContainsNode(node_subsplit)) {
    return host_dag_->GetNodeId(node_subsplit);
  }
  // Check if subsplit is in the graft.
  
}

size_t SubsplitDAGGraft::GetRootNodeId() const {
  return host_dag_->GetRootNodeId();
}

// ** Counts

size_t SubsplitDAGGraft::NodeCount() const {
  return HostNodeCount() + GraftNodeCount();
}

size_t SubsplitDAGGraft::GraftNodeCount() const { return graft_nodes_.size(); }

size_t SubsplitDAGGraft::HostNodeCount() const { return host_dag_->NodeCount(); }

size_t SubsplitDAGGraft::EdgeCount() const {
  return HostEdgeCount() + GraftEdgeCount();
}

size_t SubsplitDAGGraft::GraftEdgeCount() const {
  return graft_nodes_.size();
}

size_t SubsplitDAGGraft::HostEdgeCount() const {
  return host_dag_->EdgeCount();
}

// ** Contains

bool SubsplitDAGGraft::ContainsNode(const Bitset subsplit) const {

}

bool SubsplitDAGGraft::ContainsEdge(const Bitset edge_pair) const {

}

// ** Traversals



// ** Validation Tests

bool SubsplitDAGGraft::IsValid() const {
  return true;
}

bool SubsplitDAGGraft::IsValidAddNodePair(const Bitset parent_subsplit, const Bitset child_subsplit) const {
  return true;
}

bool SubsplitDAGGraft::IsValidRemoveNodePair(const Bitset parent_subsplit, const Bitset child_subsplit) const {
  return true;
}
