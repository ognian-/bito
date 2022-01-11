// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag_graft.hpp"

#include "gp_dag.hpp"
#include "reindexer.hpp"
#include "subsplit_dag.hpp"

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

std::pair<BitsetVector, BitsetVector> BitsetVectorDiff(BitsetVector &vector_a,
                                                       BitsetVector &vector_b) {
  auto get_vec_a_not_b = [](BitsetVector &vec_a, BitsetVector &vec_b) {
    BitsetVector a_not_b;
    for (size_t i = 0; i < vec_a.size(); i++) {
      bool is_found = false;
      for (size_t j = 0; j < vec_b.size(); j++) {
        if (vec_a[i] == vec_b[j]) {
          is_found = true;
          break;
        }
      }
      if (is_found == false) {
        a_not_b.push_back(vec_a[i]);
      }
    }
    return a_not_b;
  };

  BitsetVector a_not_b = get_vec_a_not_b(vector_a, vector_b);
  BitsetVector b_not_a = get_vec_a_not_b(vector_b, vector_a);
  return {a_not_b, b_not_a};
}

int SubsplitDAGGraft::CompareToDAG(SubsplitDAGGraft &dag_a, GPDAG &dag_b) {
  // Compare taxon counts.
  if (dag_a.TaxonCount() - dag_b.TaxonCount() != 0) {
    return dag_a.TaxonCount() - dag_b.TaxonCount();
  }
  // Compare nodes.
  BitsetVector nodes_a = dag_a.GetSortedVectorOfNodeBitsets();
  BitsetVector nodes_b = dag_b.GetSortedVectorOfNodeBitsets();
  // TODO: Add taxon translation.
  for (size_t i = 0; i < std::min(nodes_a.size(), nodes_b.size()); i++) {
    if (Bitset::Compare(nodes_a[i], nodes_b[i]) != 0) {
      return Bitset::Compare(nodes_a[i], nodes_b[i]);
    }
  }
  if (nodes_a.size() - nodes_b.size() != 0) {
    return nodes_a.size() - nodes_b.size();
  }
  // Compare edges.
  BitsetVector edges_a = dag_a.GetSortedVectorOfEdgeBitsets();
  BitsetVector edges_b = dag_b.GetSortedVectorOfEdgeBitsets();
  // TODO: Add taxon translation.
  for (int i = 0; i < std::min(edges_a.size(), edges_b.size()); i++) {
    if (Bitset::Compare(edges_a[i], edges_b[i]) != 0) {
      return Bitset::Compare(edges_a[i], edges_b[i]);
    }
  }
  if (edges_a.size() - edges_b.size() != 0) {
    return edges_a.size() - edges_b.size();
  }
  return 0;
}

// ** Modify DAG

void SubsplitDAGGraft::CreateGraftNode(const Bitset &node_subsplit) {
  // Do not add node to DAG if already exists in DAG.
  bool node_exists_in_dag = host_dag_->ContainsNode(node_subsplit);
  bool node_exists_in_graft = ContainsGraftNode(node_subsplit);
  if (node_exists_in_dag || node_exists_in_graft) {
    return;
  }
  // else, add node to graft node list.
  size_t node_id = NodeCount();
  graft_nodes_.push_back(std::make_unique<SubsplitDAGNode>(node_id, node_subsplit));
  GetDAGNode(node_id)->SetBitset(node_subsplit);
  // insert node into subsplit map.
  SafeInsert(subsplit_to_id_, node_subsplit, node_id);
  // insert node's clades into clade map.
  AddNodeToClades(node_id, node_subsplit);
}

void SubsplitDAGGraft::CreateGraftEdge(const size_t parent_id, const size_t child_id) {
  // find relationship between nodes.
  Bitset parent_subsplit = GetDAGNode(parent_id)->GetBitset();
  Bitset child_subsplit = GetDAGNode(child_id)->GetBitset();
  bool rotated = Bitset::SubsplitIsWhichChildOf(parent_subsplit, child_subsplit);
  CreateGraftEdge(parent_id, child_id, rotated);
}

void SubsplitDAGGraft::CreateGraftEdge(const size_t parent_id, const size_t child_id,
                                       const bool rotated) {
  // Add edge to edge list.
  Assert(ContainsNode(parent_id), "Node with the given parent_id does not exist.");
  Assert(ContainsNode(child_id), "Node with the given child_id does not exist.");
  auto parent_node = GetDAGNode(parent_id);
  auto child_node = GetDAGNode(child_id);
  //
  Bitset parent_bitset = parent_node->GetBitset();
  Bitset child_bitset = child_node->GetBitset();
  Bitset edge_bitset = Bitset::Edge(parent_bitset, child_bitset);
  // Add edge to grafted edges.
  bool edge_does_not_exist = !ContainsEdge(parent_id, child_id);
  if (edge_does_not_exist) {
    Bitset edge_bitset = Bitset::Edge(parent_bitset, child_bitset);
    graft_edges_.insert(
        std::make_pair<SizePair, size_t>({parent_id, child_id}, EdgeCount()));
    // Add edge to individual SubsplitDAGNodes (only adds them to graft nodes).
    bool parent_node_is_in_dag = parent_id >= HostNodeCount();
    bool child_node_is_in_dag = child_id >= HostNodeCount();
    if (rotated) {
      if (!parent_node_is_in_dag) {
        parent_node->AddLeafwardLeftside(child_node->Id());
      }
      if (!child_node_is_in_dag) {
        child_node->AddRootwardLeftside(parent_node->Id());
      }
    } else {
      if (!parent_node_is_in_dag) {
        parent_node->AddLeafwardRightside(child_node->Id());
      }
      if (!child_node_is_in_dag) {
        child_node->AddRootwardRightside(parent_node->Id());
      }
    }
  }
}

void SubsplitDAGGraft::DestroyGraftNode(const Bitset &node_subsplit) {
  // Remove node from node list.
}

void SubsplitDAGGraft::DestroyGraftEdge(const Bitset &parent_subsplit,
                                        const Bitset &child_subsplit) {
  // Assert nodes and edge exists.

  // Remove edge from edge list.
}

SubsplitDAG::ModificationResult SubsplitDAGGraft::AddGraftNodePair(
    const Bitset &parent_subsplit, const Bitset &child_subsplit) {
  std::cout << "Begin AddGraft" << std::endl;
  // Initialize results.
  SubsplitDAG::ModificationResult results;
  // Assert that at least one of the nodes are not in the host DAG or graft DAG.
  const bool parent_is_new =
      !(host_dag_->ContainsNode(parent_subsplit) || ContainsGraftNode(parent_subsplit));
  const bool child_is_new =
      !(host_dag_->ContainsNode(child_subsplit) || ContainsGraftNode(child_subsplit));
  // Assert that parent and child are a valid node pair to add to DAG.
  if (host_dag_->IsValidAddNodePair(parent_subsplit, child_subsplit) == false) {
    Failwith("AddGraftNodePair(): Is not valid node pair.");
  }
  // Assert that parent and child don't already exist in the DAG.
  if (!parent_is_new && !child_is_new) {
    return results;
  }
  // Find relationship between nodes.
  bool is_leftside = Bitset::SubsplitIsWhichChildOf(parent_subsplit, child_subsplit);
  // Add nodes.
  if (parent_is_new) {
    CreateGraftNode(parent_subsplit);
  }
  const size_t parent_id = GetNodeId(parent_subsplit);
  if (child_is_new) {
    CreateGraftNode(child_subsplit);
  }
  const size_t child_id = GetNodeId(child_subsplit);
  // Find all adjacent host nodes and add edges to fully connect to new graft nodes.
  if (parent_is_new) {
    const auto [left_parents_of_parent, right_parents_of_parent] =
        host_dag_->BuildParentIdVector(parent_subsplit);
    ConnectNodeToAdjacentHostNodes(parent_id, left_parents_of_parent, false, true,
                                   child_id);
    ConnectNodeToAdjacentHostNodes(parent_id, right_parents_of_parent, false, false,
                                   child_id);
    const auto [left_children_of_parent, right_children_of_parent] =
        host_dag_->BuildChildIdVector(parent_subsplit);
    ConnectNodeToAdjacentHostNodes(parent_id, left_children_of_parent, true, true,
                                   child_id);
    ConnectNodeToAdjacentHostNodes(parent_id, right_children_of_parent, true, false,
                                   child_id);
  }
  if (child_is_new) {
    const auto [left_parents_of_child, right_parents_of_child] =
        host_dag_->BuildParentIdVector(child_subsplit);
    ConnectNodeToAdjacentHostNodes(child_id, left_parents_of_child, false, true,
                                   parent_id);
    ConnectNodeToAdjacentHostNodes(child_id, right_parents_of_child, false, false,
                                   parent_id);
    const auto [left_children_of_child, right_children_of_child] =
        host_dag_->BuildChildIdVector(child_subsplit);
    ConnectNodeToAdjacentHostNodes(child_id, left_children_of_child, true, true,
                                   parent_id);
    ConnectNodeToAdjacentHostNodes(child_id, right_children_of_child, true, false,
                                   parent_id);
  }
  // Connect parent to child.
  CreateGraftEdge(parent_id, child_id);

  std::cout << "End AddGraft" << std::endl;

  return results;
}

// TODO:
SubsplitDAG::ModificationResult SubsplitDAGGraft::RemoveGraftNodePair(
    const Bitset &parent_subsplit, const Bitset &child_subsplit) {
  // Initialize results.
  SubsplitDAG::ModificationResult results;

  // Find nodes in node list.
  // Remove graft edges

  return results;
}

void SubsplitDAGGraft::RemoveAllGrafts() {
  graft_nodes_.clear();
  graft_edges_.clear();
  subsplit_to_id_.clear();
  clade_to_ids_.clear();
}

SizeVector SubsplitDAGGraft::SortGraftNodes() {
  SizeVector node_reindexer = Reindexer::IdentityReindexer(GraftNodeCount());

  // sort graft nodes using SubsplitCompare.

  return node_reindexer;
}

void SubsplitDAGGraft::ConnectNodeToAdjacentHostNodes(
    const size_t main_node_id, const SizeVector adjacent_node_ids,
    const bool is_main_node_parent, const bool is_left_child,
    std::optional<size_t> ignored_node_id_opt) {
  for (const auto adjacent_node_id : adjacent_node_ids) {
    if (ignored_node_id_opt.has_value() && (adjacent_node_id == *ignored_node_id_opt)) {
      continue;
    }
    const size_t parent_id = (is_main_node_parent ? main_node_id : adjacent_node_id);
    const size_t child_id = (is_main_node_parent ? adjacent_node_id : main_node_id);
    CreateGraftEdge(parent_id, child_id, is_left_child);
  }
}

// TODO:
// ** Clades

void SubsplitDAGGraft::InitClades() {
  clade_to_ids_.clear();
  // Add host clades.
  for (const auto &subsplit_id_pair : host_dag_->GetSubsplitToIdMap()) {
    size_t node_id = subsplit_id_pair.second;
    Bitset node_subsplit = subsplit_id_pair.first;
    AddNodeToClades(node_id, node_subsplit);
  }
  // Add graft clades.
  for (const auto &subsplit_id_pair : subsplit_to_id_) {
    size_t node_id = subsplit_id_pair.second;
    Bitset node_subsplit = subsplit_id_pair.first;
    AddNodeToClades(node_id, node_subsplit);
  }
}

void SubsplitDAGGraft::AddNodeToClades(const size_t node_id,
                                       const Bitset &node_subsplit,
                                       const bool sort_ids) {
  // Add leftside and rightside clade of node.
  Bitset left_clade = node_subsplit.SubsplitGetClade(0);
  Bitset right_clade = node_subsplit.SubsplitGetClade(1);
  for (const auto clade : {left_clade, right_clade}) {
    // If this is the first occurance, add clade and id to map.
    if (clade_to_ids_.find(clade) == clade_to_ids_.end()) {
      SizeVector ids = {node_id};
      clade_to_ids_.insert(std::make_pair(clade, ids));
    }
    // If clade is already in map, add id to list.
    else {
      auto &ids = clade_to_ids_.find(clade)->second;
      ids.push_back(node_id);
      std::sort(ids.begin(), ids.end());
    }
  }
}

void SubsplitDAGGraft::RemoveNodeFromClades(const size_t node_id,
                                            const Bitset &node_subsplit,
                                            const bool sort_ids) {
  // Add leftside and rightside clade of node.
  Bitset left_clade = node_subsplit.SubsplitGetClade(0);
  Bitset right_clade = node_subsplit.SubsplitGetClade(1);
  for (const auto clade : {left_clade, right_clade}) {
    auto &ids = clade_to_ids_.find(clade)->second;
    for (size_t i = 0; i < ids.size(); i++) {
      if (ids[i] == node_id) {
        ids.erase(ids.begin() + i);
      }
    }
  }
}

SizeVector SubsplitDAGGraft::GetAllChildrenOfNode(const Bitset &node_subsplit,
                                                  const bool which_child) const {}

SizeVector SubsplitDAGGraft::GetAllParentsOfNode(const Bitset &node_subsplit,
                                                 const bool which_child) const {}

// ** Getters

SubsplitDAG *SubsplitDAGGraft::GetHostDAG() { return host_dag_; }

SubsplitDAGNode *SubsplitDAGGraft::GetDAGNode(const size_t node_id) const {
  Assert(node_id < NodeCount(), "Node Id is out of valid range.");
  // Check if node is in the host DAG.
  if (node_id < HostNodeCount()) {
    return host_dag_->GetDAGNode(node_id);
  }
  // else, node is in the graft.
  return graft_nodes_.at(node_id - host_dag_->NodeCount()).get();
}

size_t SubsplitDAGGraft::GetNodeId(const Bitset &node_subsplit) const {
  // Check if subsplit is in the host DAG.
  if (host_dag_->ContainsNode(node_subsplit)) {
    return host_dag_->GetNodeId(node_subsplit);
  }
  // Check if subsplit is in the graft.
  Assert(subsplit_to_id_.find(node_subsplit) != subsplit_to_id_.end(),
         "GetNodeId(): Subsplit does not correspond to a node in GraftDAG.");
  return subsplit_to_id_.at(node_subsplit);
}

size_t SubsplitDAGGraft::GetRootNodeId() const { return host_dag_->GetRootNodeId(); }

size_t SubsplitDAGGraft::GetEdgeIdx(const Bitset &parent_subsplit,
                                    const Bitset &child_subsplit) const {
  size_t parent_id = GetNodeId(parent_subsplit);
  size_t child_id = GetNodeId(child_subsplit);
}

size_t SubsplitDAGGraft::GetEdgeIdx(const size_t parent_id,
                                    const size_t child_id) const {
  // Check for edge in graft.
  if (graft_edges_.find(std::make_pair(parent_id, child_id)) != graft_edges_.end()) {
    return graft_edges_.at(std::make_pair(parent_id, child_id));
  }
  // Check for edge in host.
  Assert(host_dag_->ContainsEdge(parent_id, child_id),
         "GetEdgeIdx: There is no edge corresponding to given parent/child pair in "
         "GraftDAG.");
  return host_dag_->GetEdgeIdx(parent_id, child_id);
}

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
  BitsetVector edges;
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
    BitsetVector host_edges = host_dag_->GetSortedVectorOfEdgeBitsets();
    edges.insert(edges.end(), host_edges.begin(), host_edges.end());
  }
  std::sort(edges.begin(), edges.end());

  return edges;
}

// ** Counts

size_t SubsplitDAGGraft::TaxonCount() const { return host_dag_->TaxonCount(); }

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

//    These two data structures maintain sorted vectors. There one difference is that
//    SortedVector allows for duplicate data, SortedUniqueVector does not. SortedVectors
//    are an alternative to std::set. std::set is useful when you have a large dataset
//    that you are frequently inserting into.  But vectors have superior cache
//    performance during access than the Red/Black Trees of std::set.// ** Contains

bool SubsplitDAGGraft::ContainsGraftNode(const Bitset node_subsplit) const {
  return subsplit_to_id_.find(node_subsplit) != subsplit_to_id_.end();
}

bool SubsplitDAGGraft::ContainsGraftNode(const size_t node_id) const {
  return (node_id >= HostNodeCount()) && (node_id < NodeCount());
}

bool SubsplitDAGGraft::ContainsNode(const Bitset node_subsplit) const {
  return host_dag_->ContainsNode(node_subsplit) || ContainsGraftNode(node_subsplit);
}

bool SubsplitDAGGraft::ContainsNode(const size_t node_id) const {
  return host_dag_->ContainsNode(node_id) || ContainsGraftNode(node_id);
}

bool SubsplitDAGGraft::ContainsGraftEdge(const size_t parent_id,
                                         const size_t child_id) const {
  return graft_edges_.find({parent_id, child_id}) != graft_edges_.end();
}

bool SubsplitDAGGraft::ContainsEdge(const size_t parent_id,
                                    const size_t child_id) const {
  return host_dag_->ContainsEdge(parent_id, child_id) ||
         ContainsGraftEdge(parent_id, child_id);
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
