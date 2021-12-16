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
  size_t node_id = GraftNodeCount();
  graft_nodes_.push_back(std::make_unique<SubsplitDAGNode>(node_id, node_subsplit));
  SafeInsert(subsplit_to_id_, node_subsplit, node_id);
}

void SubsplitDAGGraft::CreateGraftEdge(const size_t &parent_id,
                                      const size_t &child_id) {
  // Add edge to edge list.
  Assert(ContainsGraftNode(parent_id), "Node with the given parent_id does not exist.");
  Assert(ContainsGraftNode(child_id), "Node with the given child_id does not exist.");
  // Insert edge between parent and child.
  ConnectGivenNodes(parent_id, child_id, rotated);
  SafeInsert(dag_edges_, {parent_id, child_id}, EdgeCountWithLeafSubsplits());
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
  // Assert that at least one of the nodes are not in the host DAG.
  const bool parent_is_new = !host_dag_->ContainsNode(parent_subsplit);
  const bool child_is_new = !host_dag_->ContainsNode(child_subsplit);
  // Assert that parent and child are a valid node pair to add to DAG.
  if (host_dag_->IsValidAddNodePair(parent_subsplit, child_subsplit) == false) {
    
  }
  // Assert that parent and child don't already exist in the DAG.
  if (!parent_is_new && !child_is_new) {
    return;
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
  if (node_id < host_dag_->NodeCount()) {
    return host_dag_->GetDAGNode(node_id);
  }
  // Check if subsplit is in the graft.
  return graft_nodes_.at(node_id - host_dag_->NodeCount()).get();
}

size_t SubsplitDAGGraft::GetNodeId(const Bitset &node_subsplit) const {
  // Check if subsplit is in the host DAG.
  if (host_dag_->ContainsNode(node_subsplit)) {
    return host_dag_->GetNodeId(node_subsplit);
  }
  // Check if subsplit is in the graft.
  return subsplit_to_id_.at(node_subsplit);
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
  return graft_edges_.size();
}

size_t SubsplitDAGGraft::HostEdgeCount() const {
  return host_dag_->EdgeCount();
}

// ** Contains

bool SubsplitDAGGraft::ContainsGraftNode(const Bitset node_subsplit) const {
  return subsplit_to_id_.find(node_subsplit) != subsplit_to_id_.end();
}

bool SubsplitDAGGraft::ContainsGraftNode(const size_t node_id) const {
  return node_id >= HostNodeCount() && node_id < NodeCount();
}

bool SubsplitDAGGraft::ContainsGraftEdge(const size_t parent_id, const size_t child_id) const {
  return graft_edges_.find({parent_id, child_id}) != graft_edges_.end();
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
