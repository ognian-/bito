// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The SubsplitDAGGraft is a proposed addition (graft) to SubsplitDAG (host), which we
// can perform computations on without the need for adding nodes and edges and
// reindexing the full DAG.
// NOTE: This may add overhead to DAG traversals.

#include "gp_dag.hpp"
#include "subsplit_dag.hpp"

#ifndef SRC_SUBSPLIT_DAG_GRAFT_HPP
#define SRC_SUBSPLIT_DAG_GRAFT_HPP

class GPDAG;

class SubsplitDAGGraft {
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
  static int CompareToDAG(SubsplitDAGGraft &dag_a, GPDAG &dag_b);

  // ** DAG Traversals
  // These traversals cover the DAG as though graft were formally added to host into a
  // single SubsplitDAG.

  void IterateOverRealNodes(const SubsplitDAG::NodeLambda &f) const;

  // ** Modify DAG

  // Add node pair to graft.
  SubsplitDAG::ModificationResult AddGraftNodePair(const Bitset &parent_subsplit,
                                                   const Bitset &child_subsplit);
  // TODO:
  // Remove node pair from graft.
  SubsplitDAG::ModificationResult RemoveGraftNodePair(const Bitset &parent_subsplit,
                                                      const Bitset &child_subsplit);
  // TODO:
  // Clear all nodes and edges from graft for reuse.
  void RemoveAllGrafts();
  // TODO:
  // Create sorted list of all node in graft. Returns node reindexer.
  SizeVector SortGraftNodes();

  // TODO:
  // ** Clades

  // Clear clade map and add all host and graft nodes to clade map.
  void InitClades();
  // Add both of node's clades to the clade map.
  void AddNodeToClades(const size_t node_id, const Bitset &node_subsplit,
                       const bool sort_ids = true);
  // Remove both of node's clades from the clade map.
  void RemoveNodeFromClades(const size_t node_id, const Bitset &node_subsplit,
                            const bool sort_ids = true);
  // Get all the child nodes of given subsplit (specify left or right child).
  SizeVector GetAllChildrenOfNode(const Bitset &node_subsplit,
                                  const bool which_child) const;
  // Get all parent nodes of a given subsplit (specify left or right child).
  SizeVector GetAllParentsOfNode(const Bitset &node_subsplit,
                                 const bool which_child) const;

  // ** Build Indexers/Vectors

  // Get the rotated and sorted parents of the node with the given subsplit.
  std::pair<SizeVector, SizeVector> BuildParentIdVector(const Bitset &subsplit) const;
  // Get the rotated and sorted children of the node with the given subsplit.
  std::pair<SizeVector, SizeVector> BuildChildIdVector(const Bitset &subsplit) const;

  // ** Getters

  // Get pointer to the host DAG.
  SubsplitDAG *GetHostDAG();
  // Get node based on node id.
  SubsplitDAGNode *GetDAGNode(const size_t node_id) const;
  // Get the node id based on the subsplit bitset.
  size_t GetNodeId(const Bitset &subsplit) const;
  // Gets the node id of the DAG root.
  size_t GetRootNodeId() const;
  // Get the node based on the nodes id.
  SubsplitDAGNode *GetNode(const size_t node_id) const;
  // Return the node ids corresponding to the rootsplits.
  const SizeVector &GetRootsplitIds() const;
  // Get the GPCSP/edge index by its parent/child pair.
  size_t GetEdgeIdx(const Bitset &parent_subsplit, const Bitset &child_subsplit) const;
  size_t GetEdgeIdx(const size_t parent_id, const size_t child_id) const;
  //
  SizeVector GetAdjacentEdges(const size_t node_id, const bool is_leafward,
                              const bool is_leftward) const;
  // Get a sorted vector of all node subsplit's bitset representation. Optionally only
  // graft nodes, or graft and host nodes.
  BitsetVector GetSortedVectorOfNodeBitsets(bool include_host = true);
  // Get a sorted vector of all edge PCSP's bitset representation. Optionally only graft
  // edges, or graft and host edges.
  BitsetVector GetSortedVectorOfEdgeBitsets(bool include_host = true);

  // ** Counts

  // Total number of taxa in DAG (same as Subsplit length).
  size_t TaxonCount() const;
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

  // Checks whether the node is in the host DAG or the graft.
  bool ContainsNode(const Bitset node_subsplit) const;
  bool ContainsNode(const size_t node_id) const;
  // Checks whether the node is in the graft only.
  bool ContainsGraftNode(const Bitset node_subsplit) const;
  bool ContainsGraftNode(const size_t node_id) const;
  // Checks whether the edge is in the host only.
  bool ContainsHostNode(const Bitset node_subsplit) const;
  bool ContainsHostNode(const size_t node_id) const;
  // Checks whether the edge is in the host DAG or the graft.
  bool ContainsEdge(const size_t parent_id, const size_t child_id) const;
  // Checks whether the edge is in the graft only.
  bool ContainsGraftEdge(const size_t parent_id, const size_t child_id) const;
  // Checks whether the edge is in the host only.
  bool ContainsHostEdge(const size_t parent_id, const size_t child_id) const;

  // ** Traversal

  //
  SizeVector LeafwardPostorderTraversalTrace() const;
  //
  SizeVector LeafwardTopologicalTraversalTrace() const;

  // ** Validity Tests

  // Checks if host and graft are in a valid state.
  bool IsValid() const;
  // Checks if host and graft are in a valid state to add given node pair.
  bool IsValidAddNodePair(const Bitset parent_subsplit,
                          const Bitset child_subsplit) const;
  // Checks if host and graft are in a valid state to remove given node pair.
  bool IsValidRemoveNodePair(const Bitset parent_subsplit,
                             const Bitset child_subsplit) const;

 protected:
  // ** Modify DAG
  // These modification do not ensure a valid, consistent state for DAG.

  // Add node to the graft.
  void CreateGraftNode(const Bitset &node_subsplit);
  // Add edge to the graft.
  void CreateGraftEdge(const size_t parent_id, const size_t child_id);
  // Add edge to the graft (Overload for when edge relation is known).
  void CreateGraftEdge(const size_t parent_id, const size_t child_id,
                       const bool rotated);
  // Add node to the graft.
  void DestroyGraftNode(const Bitset &node_subsplit);
  void DestroyGraftNode(const size_t node_id);
  // Add edge to the graft.
  void DestroyGraftEdge(const Bitset &parent_subsplit,
                        const Bitset &child_node_subsplit);
  void DestroyGraftEdge(const size_t parent_id, const size_t child_id);
  void DestroyGraftEdge(const size_t edge_idx);
  // Connect main node to all adjacent nodes in vector.
  void ConnectNodeToAdjacentHostNodes(
      const size_t main_node_id, const SizeVector adjacent_node_ids,
      const bool is_main_node_parent, const bool is_left_child,
      std::optional<size_t> ignored_node_id_opt = std::nullopt);

 protected:
  // DAG that the graft is proposed to be connected to.
  SubsplitDAG *host_dag_;
  // Nodes in the graft.
  std::vector<std::unique_ptr<SubsplitDAGNode>> graft_nodes_;
  // Edges in the graft.
  std::map<SizePair, size_t> graft_edges_;
  // Map of all DAG Nodes:
  //   - [ Node Subsplit (Bitset) => Node Id ]
  BitsetSizeMap subsplit_to_id_;
  // A map from a clade to the vector of node ids containing that clade.
  std::map<Bitset, SizeVector> clade_to_ids_;
};

#endif  // SRC_SUBSPLIT_DAG_GRAFT_HPP
