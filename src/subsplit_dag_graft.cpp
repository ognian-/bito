// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag_graft.hpp"

#include "subsplit_dag.hpp"
#include "gp_dag.hpp"

// ** Constructors

SubsplitDAGGraft::SubsplitDAGGraft(SubsplitDAG &dag) : host_dag_(&dag) {}

// ** Comparators

int SubsplitDAGGraft::Compare(SubsplitDAGGraft &dag_a, SubsplitDAGGraft &dag_b) {
  SubsplitDAG &host_a = *(dag_a.host_dag_);
  SubsplitDAG &host_b = *(dag_b.host_dag_);
  // (1) Compare host DAGs.
  int host_diff = SubsplitDAG::Compare(host_a, host_b);
  if (host_diff != 0) {
    return host_diff;
  }
  // (2) Compare graft nodes.

  // (3) Compare graft edges.
}


int SubsplitDAGGraft::CompareToSubsplitDAG(SubsplitDAGGraft &dag_a, SubsplitDAG &dag_b) {
  // auto host_a = dag_a.GetHostDAG();
  // Compare taxon counts.
  // if (host_a->TaxonCount() - dag_b.TaxonCount() != 0) {
  //   return host_a->TaxonCount() - dag_b.TaxonCount();
  // }
  // Compare nodes.
  auto nodes_a = dag_a.GetSortedVectorOfNodeBitsets();
  auto nodes_b = dag_b.GetSortedVectorOfNodeBitsets();
  // TODO: Add node translation.
  for (int i = 0; i < std::min(nodes_a.size(), nodes_b.size()); i++) {
    if (Bitset::Compare(nodes_a[i], nodes_b[i]) != 0) {
      return Bitset::Compare(nodes_a[i], nodes_b[i]);
    }
  }
  if (nodes_a.size() - nodes_b.size() != 0) {
    return nodes_a.size() - nodes_b.size();
  }
  // Compare edges.
  auto edges_a = dag_a.GetSortedVectorOfEdgeBitsets();
  auto edges_b = dag_b.GetSortedVectorOfEdgeBitsets();
  // TODO: Add node translation.
  for (int i = 0; i < std::min(edges_a.size(), edges_b.size()); i++) {
    if (Bitset::Compare(edges_a[i], edges_b[i]) != 0) {
      return Bitset::Compare(edges_a[i], edges_b[i]);
    }
  }
  if (edges_a.size() - edges_b.size() != 0) {
    return edges_a.size() - edges_b.size();
  }
}

// ** Modify DAG

void SubsplitDAGGraft::CreateGraftNode(const Bitset &node_subsplit) {
  // Add node to node list.
  size_t node_id = GraftNodeCount();
  graft_nodes_.push_back(std::make_unique<SubsplitDAGNode>(node_id, node_subsplit));
  SafeInsert(subsplit_to_id_, node_subsplit, node_id);
}

void SubsplitDAGGraft::CreateGraftEdge(const size_t parent_id, const size_t child_id) {
  // find relationship between nodes.
}

void SubsplitDAGGraft::CreateGraftEdge(const size_t parent_id, const size_t child_id,
                                       const bool rotated) {
  // Add edge to edge list.
  Assert(ContainsGraftNode(parent_id), "Node with the given parent_id does not exist.");
  Assert(ContainsGraftNode(child_id), "Node with the given child_id does not exist.");
  auto parent_node = GetDAGNode(parent_id);
  auto child_node = GetDAGNode(child_id);
  // Add edge to SubsplitDAGNodes
  if (rotated) {
    parent_node->AddLeafwardLeftside(child_node->Id());
    child_node->AddRootwardLeftside(parent_node->Id());
  } else {
    parent_node->AddLeafwardRightside(child_node->Id());
    child_node->AddRootwardRightside(parent_node->Id());
  }

  // Add edge to grafted edges.
  SafeInsert(graft_edges_, {parent_id, child_id}, EdgeCount());
}

void SubsplitDAGGraft::DestroyGraftNode(const Bitset &node_subsplit) {
  // Remove node from node list.
}

void SubsplitDAGGraft::DestroyGraftEdge(const Bitset &parent_subsplit,
                                        const Bitset &child_subsplit) {
  // Assert nodes and edge exists.

  // Remove edge from edge list.
}

void SubsplitDAGGraft::AddGraftNodePair(const Bitset &parent_subsplit,
                                        const Bitset &child_subsplit) {
  // Assert that at least one of the nodes are not in the host DAG.
  const bool parent_is_new = !host_dag_->ContainsNode(parent_subsplit);
  const bool child_is_new = !host_dag_->ContainsNode(child_subsplit);
  // Assert that parent and child are a valid node pair to add to DAG.
  if (host_dag_->IsValidAddNodePair(parent_subsplit, child_subsplit) == false) {
    Failwith("AddGraftNodePair(): Is not valid node pair.");
  }
  // Assert that parent and child don't already exist in the DAG.
  if (!parent_is_new && !child_is_new) {
    return;
  }
  // Find relationship between nodes.
  bool is_leftside = Bitset::SubsplitIsWhichChildOf(parent_subsplit, child_subsplit);
  // Add nodes.
  CreateGraftNode(parent_subsplit);
  const size_t parent_id = GetNodeId(parent_subsplit);
  CreateGraftNode(child_subsplit);
  const size_t child_id = GetNodeId(child_subsplit);
  // Find adjacent nodes and add edges.
  const auto [left_parents_of_parent, right_parents_of_parent] =
      host_dag_->BuildParentIdVector(parent_subsplit);
  ConnectNodeToAdjacentNodes(parent_id, left_parents_of_parent, false, true, child_id,
                             true);
  ConnectNodeToAdjacentNodes(parent_id, right_parents_of_parent, false, false, child_id,
                             true);
  const auto [left_children_of_parent, right_children_of_parent] =
      host_dag_->BuildChildIdVector(parent_subsplit);
  ConnectNodeToAdjacentNodes(parent_id, left_children_of_parent, true, true, child_id,
                             true);
  ConnectNodeToAdjacentNodes(parent_id, left_children_of_parent, true, false, child_id,
                             true);
  const auto [left_parents_of_child, right_parents_of_child] =
      host_dag_->BuildParentIdVector(parent_subsplit);
  ConnectNodeToAdjacentNodes(child_id, left_parents_of_child, false, true, parent_id,
                             true);
  ConnectNodeToAdjacentNodes(child_id, right_parents_of_child, false, false, parent_id,
                             true);
  const auto [left_children_of_child, right_children_of_child] =
      host_dag_->BuildChildIdVector(parent_subsplit);
  ConnectNodeToAdjacentNodes(child_id, left_parents_of_child, true, true, parent_id,
                             true);
  ConnectNodeToAdjacentNodes(child_id, right_parents_of_child, true, false, parent_id,
                             true);
  // Connect parent to child.
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

void SubsplitDAGGraft::SortGraftNodes() {}

void SubsplitDAGGraft::ConnectNodeToAdjacentNodes(const size_t main_node_id,
                                                  const SizeVector adjacent_node_ids,
                                                  const bool is_main_node_parent,
                                                  const bool is_left_child,
                                                  const size_t ignored_node_id,
                                                  const bool is_node_ignored) {
  for (const auto adjacent_node_id : adjacent_node_ids) {
    if (is_node_ignored && adjacent_node_id == ignored_node_id) {
      continue;
    }
    const size_t parent_id = (is_main_node_parent ? main_node_id : adjacent_node_id);
    const size_t child_id = (is_main_node_parent ? adjacent_node_id : main_node_id);
    CreateGraftEdge(parent_id, child_id, is_left_child);
  }
}

// ** Getters

// SubsplitDAG *SubsplitDAGGraft::GetHostDAG() {
//   return host_dag_;
// }

SubsplitDAGNode *SubsplitDAGGraft::GetDAGNode(const size_t node_id) const {
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

size_t SubsplitDAGGraft::GetRootNodeId() const { return host_dag_->GetRootNodeId(); }

BitsetVector SubsplitDAGGraft::GetSortedVectorOfNodeBitsets(bool include_host) {
  std::vector<Bitset> nodes;
  // Get graft node bitsets.
  for (const auto &subsplit_id_pair : subsplit_to_id_) {
    nodes.push_back(subsplit_id_pair.first);
  }
  // Get host nodes bitsets.
  if (include_host) {
    std::vector<Bitset> host_nodes = host_dag_->GetSortedVectorOfNodeBitsets();
    nodes.insert(nodes.end(), host_nodes.begin(), host_nodes.end());
  }
  std::sort(nodes.begin(), nodes.end());
  return nodes;
}

BitsetVector SubsplitDAGGraft::GetSortedVectorOfEdgeBitsets(bool include_host) {
  std::vector<Bitset> edges;
  // Get graft edge bitsets.
  for (const auto &idpair_idx_pair : graft_edges_) {
    auto id_pair = idpair_idx_pair.first;
    auto parent_bitset = GetDAGNode(id_pair.first)->GetBitset();
    auto child_bitset = GetDAGNode(id_pair.second)->GetBitset();
    Bitset edge_bitset = Bitset::Edge(parent_bitset, child_bitset);
    edges.push_back(edge_bitset);
  }
  // Get host edge bitsets.
  if (include_host) {
    std::vector<Bitset> host_edges = host_dag_->GetSortedVectorOfEdgeBitsets();
    edges.insert(edges.end(), host_edges.begin(), host_edges.end());
  }
  std::sort(edges.begin(), edges.end());
  return edges;
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

size_t SubsplitDAGGraft::GraftEdgeCount() const { return graft_edges_.size(); }

size_t SubsplitDAGGraft::HostEdgeCount() const {
  return host_dag_->EdgeCountWithLeafSubsplits();
}

// ** Contains

bool SubsplitDAGGraft::ContainsGraftNode(const Bitset node_subsplit) const {
  return subsplit_to_id_.find(node_subsplit) != subsplit_to_id_.end();
}

bool SubsplitDAGGraft::ContainsGraftNode(const size_t node_id) const {
  return node_id >= HostNodeCount() && node_id < NodeCount();
}

bool SubsplitDAGGraft::ContainsGraftEdge(const size_t parent_id,
                                         const size_t child_id) const {
  return graft_edges_.find({parent_id, child_id}) != graft_edges_.end();
}

// ** Traversals

// ** Validation Tests

bool SubsplitDAGGraft::IsValid() const { return true; }

bool SubsplitDAGGraft::IsValidAddNodePair(const Bitset parent_subsplit,
                                          const Bitset child_subsplit) const {
  return true;
}

bool SubsplitDAGGraft::IsValidRemoveNodePair(const Bitset parent_subsplit,
                                             const Bitset child_subsplit) const {
  return true;
}
