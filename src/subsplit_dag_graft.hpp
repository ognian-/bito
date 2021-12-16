// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The SubsplitDAGGraft is a proposed addition (graft) to SubsplitDAG (host), which we
// can perform computations on without the need for adding nodes and edges and
// reindexing the full DAG.
// NOTE: This may add overhead to DAG traversals.

#include "subsplit_dag.hpp"

#ifndef SRC_SUBSPLIT_DAG_GRAFT_HPP
#define SRC_SUBSPLIT_DAG_GRAFT_HPP

class SubsplitDAGGraft {
 private: 
  // DAG that the graft is proposed to be connected to.
  SubsplitDAG *host_dag_;
  // Nodes in the graft.
  std::vector<std::unique_ptr<SubsplitDAGNode>> graft_nodes_;
  // Edges in the graft.
  std::map<SizePair, size_t> graft_edges_;
  // Nodes that are adjacent to a graft node (connect by edge).
  std::vector<size_t> graft_adjacent_dag_nodes_;
  // A map from a clade to the vector of node ids containing that clade.
  std::map<Bitset, SizeVector> clade_to_ids_;
  
 public:
  // ** Constructors:

  // Initialize empty DAGGraft.
  SubsplitDAGGraft(SubsplitDAG &dag);
  // TODO:
  // Initialize DAGGraft with an initial graft.
  SubsplitDAGGraft(SubsplitDAG &dag, NNIOperation &nni);
  SubsplitDAGGraft(SubsplitDAG &dag, Bitset &parent_subsplit, Bitset &child_subsplit);

  // ** Comparators
  // Uses same method of comparison as SubsplitDAG (node and edge sets).
  static int Compare(SubsplitDAGGraft &dag_a, SubsplitDAGGraft &dag_b);
  // Treats SubsplitDAGGraft as completed SubsplitDAG to compare against normal
  // SubsplitDAG.
  static int Compare(SubsplitDAGGraft &dag_a, SubsplitDAG &dag_b);

  // ** DAG Traversals
  // These traversals cover the DAG as though graft were formally added to host into a
  // single SubsplitDAG.

  void IterateOverRealNodes(const SubsplitDAG::NodeLambda &f) const;

  // ** Modify DAG

  // Add node pair to graft.
  void AddGraftNodePair(const Bitset &parent_subsplit, const Bitset &child_subsplit);
  // TODO:
  // Remove node pair from graft.
  void RemoveGraftNodePair(const Bitset &parent_subsplit, const Bitset &child_subsplit);
  // TODO:
  // Clear all nodes and edges from graft for reuse.
  void RemoveAllGrafts();

  // ** Getters

  // Get node based on node id.
  SubsplitDAGNode *GetDAGNode(const size_t node_id) const;
  // Get the node id based on the subsplit bitset.
  size_t GetNodeId(const Bitset &subsplit) const;
  // Gets the node id of the DAG root.
  size_t GetRootNodeId() const;
  // Get the node based on the nodes id.
  SubsplitDAGNode* GetNode(const size_t node_id) const;
  // Return the node ids corresponding to the rootsplits.
  const SizeVector &GetRootsplitIds() const;
  // Get the GPCSP/edge index by its parent-child pair subsplits from the DAG nodes.
  size_t GetGPCSPEdgeIdx(const Bitset &parent_subsplit,
                         const Bitset &child_subsplit) const;
  // Get the GPCSP edge index by its parent-child pair id from the DAG nodes.
  size_t GetGPCSPEdgeIdx(const size_t parent_id, const size_t child_id) const;
  // Get a sorted vector of all node subsplit's bitset representation. Optionally only
  // graft nodes, or graft and host nodes.
  BitsetVector GetSortedVectorOfNodeBitsets(bool include_host = true);
  // Get a sorted vector of all edge PCSP's bitset representation. Optionally only graft
  // edges, or graft and host edges.
  BitsetVector GetSortedVectorOfEdgeBitsets(bool include_host = true);

  // ** Counts

  // Total number of nodes in full proposed DAG.
  size_t NodeCount() const;
  // Total number of nodes in graft only.
  size_t GraftNodeCount() const;
  // Total number of nodes in host DAG only.
  size_t HostNodeCount() const;
  // Total number of edges in full proposed DAG.
  size_t EdgeCount() const;
  // Total number of edges in graft only.
  size_t GraftEdgeCount() const;
  // Total number of edges in host DAG only.
  size_t HostEdgeCount() const;

  // ** Contains

  // 
  bool ContainsNode(const Bitset subsplit) const;
  // 
  bool ContainsEdge(const Bitset edge) const;

  // ** Traversal

  // 
  SizeVector LeafwardPostorderTraversalTrace() const;
  // 
  SizeVector LeafwardTopologicalTraversalTrace() const;

  // ** Validity Tests

  // Checks if host and graft are in a valid state.
  bool IsValid() const;
  // Checks if host and graft are in a valid state to add given node pair.
  bool IsValidAddNodePair(const Bitset parent_subsplit, const Bitset child_subsplit) const;
  // Checks if host and graft are in a valid state to remove given node pair.
  bool IsValidRemoveNodePair(const Bitset parent_subsplit, const Bitset child_subsplit) const;

 private:
  // ** Modify DAG

  // Add node to the graft.
  void CreateGraftNode(const Bitset &node_subsplit);
  // Add edge to the graft.
  void CreateGraftEdge(const Bitset &parent_subsplit, const Bitset &child_subsplit);
  // Add node to the graft.
  void DestroyGraftNode(const Bitset &node_subsplit);
  void DestroyGraftNode(const size_t node_id);
  // Add edge to the graft.
  void DestroyGraftEdge(const Bitset &parent_subsplit, const Bitset &child_node_subsplit);
  void DestroyGraftEdge(const size_t parent_id, const size_t child_id);
  void DestroyGraftEdge(const size_t edge_idx);
};

#endif  // SRC_SUBSPLIT_DAG_GRAFT_HPP
