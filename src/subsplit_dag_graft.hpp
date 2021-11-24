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
 public:
  // ** Constructors:

  // Initialize empty graft.
  SubsplitDAGGraft(SubsplitDAG &dag);
  // TODO:
  SubsplitDAGGraft(SubsplitDAG &dag, NNIOperation &nni);
  // TODO:
  SubsplitDAGGraft(SubsplitDAG &dag, Bitset &parent_subsplit, Bitset &child_subsplit);

  // ** DAG Traversals methods:

  // ** Modify graft methods:

  // Add node pair to graft.
  void AddNodePair(Bitset &parent_subsplit, Bitset &child_subsplit);
  // Remove node pair from graft.
  void RemoveNodePair(Bitset &parent_subsplit, Bitset &child_subsplit);
  // Clear all nodes and edges from graft for reuse.
  void Clear();

  // ** Counts

  // Total number of nodes in full proposed DAG.
  size_t NodeCount() const;
  // Total number of nodes in graft only.
  size_t GraftNodeCount() const;
  // Total number of nodes in host DAG only.
  size_t HostNodeCount() const;
  // Total number of edges in full proposed DAG.
  size_t GPCSPCount() const;
  // Total number of edges in graft only.
  size_t GraftGPCSPCount() const;
  // Total number of edges in host DAG only.
  size_t HostGPCSPCount() const;

 private:
  // DAG that the graft is proposed to be connected to.
  SubsplitDAG *host_dag_;
  // Nodes in the graft.
  std::vector<std::unique_ptr<SubsplitDAGNode>> graft_nodes_;
  // Edges in the graft.
  std::map<SizePair, size_t> graft_edges_;
  // Nodes that are adjacent to a graft node (connect by edge).
  std::vector<size_t> graft_adjacent_dag_nodes_;
};

#endif  // SRC_SUBSPLIT_DAG_GRAFT_HPP
