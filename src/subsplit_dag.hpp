#ifndef __SUBSPLIT_DAG_H__
#define __SUBSPLIT_DAG_H__

// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The purpose of this class is to hold a DAG built from the parent-child relationships
// of the subsplits. We wish to have information associated with both the nodes and
// edges of the DAG. Our strategy for doing that is via non-negative integer indices:
// nodes have a unique size_t `Id`, and each edge has a unique `edge_idx`. We can then
// store arbitrary associated information in other data structures associated with these
// indices.
//
// The DAG has a well-defined notion of rootward and leafward.
//
// The data structure for this DAG is as follows:
// - The nodes of the DAG are vectors of indices representing their edges to other
// nodes. These include edges to the children (leafward edges) and also the edges that
// are connecting to that node (rootward edges). Furthermore, these sets of edges are
// separated into two subsets: those that descend from the left component of the
// subsplit
// ("rotated" edges) and those that descend from the right component of the subsplit
// ("sorted" edges). The nodes of the DAG also include bitsets describing the taxa they
// contain.
// - The edges of the DAG are indexed separately, and there is a map (`dag_edges_`)
// which maps from pairs of node Ids to its edge index. These DAG edges are indexed
// such that all of the sorted edges descending from a given node have a contiguous set
// of indices, as do all of the rotated indices. The range of indices for such a set of
// edges is given by the `parent_to_child_range_` map.
// There is also a `subsplit_to_id_` map that maps from the subsplit bitset of the DAG
// node (and its rotated version) to the node Id.
// - To clarify terminology, the "DAG root node" refers to the universal ancestor and is
// representated by a bitset where the first half is all ones and the second half is all
// zeros (e.g. for taxon_count_ = 4, 1111|0000). Note that this implies that the DAG
// root node only has rotated children. Children of the DAG root node are called
// "rootsplits" and partition the whole taxon set.

#ifndef SRC_SUBSPLIT_DAG_HPP_
#define SRC_SUBSPLIT_DAG_HPP_

#include "reindexer.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "subsplit_dag_action.hpp"
#include "subsplit_dag_nni.hpp"
#include "subsplit_dag_node.hpp"

class SubsplitDAG {
 public:
  // ** Constructor

  // Build empty Subsplit DAG with no topologies and no taxa.
  SubsplitDAG();
  // Build a Subsplit DAG expressing all tree topologies from tree_collection.
  explicit SubsplitDAG(const RootedTreeCollection &tree_collection);

  // ** Comparator

  // This compare ensures that both DAGs have the same topology according to their
  // set of node and edge bitsets.  However, it does ensure that DAGs have the same id
  // and idxs for their respective nodes and edges, only that they contain the
  // same set of nodes and edges (as long as taxon positions in the
  // clades have the same mapping). This compare ensures a total ordering.
  static int Compare(const SubsplitDAG &lhs, const SubsplitDAG &rhs);
  friend bool operator==(const SubsplitDAG &lhs, const SubsplitDAG &rhs);
  friend bool operator!=(const SubsplitDAG &lhs, const SubsplitDAG &rhs);

  // ** Counts

  // The total number of individual taxa in the DAG.
  size_t TaxonCount() const;
  // The total number of nodes in the DAG (including the root and leaves).
  size_t NodeCount() const;
  // The total number of nodes in the DAG (excluding the root, but including the
  // leaves).
  size_t NodeCountWithoutDAGRoot() const;
  // The total number of rootsplits in DAG. These count all direct descendents of the
  // root (also, the union of each rootsplits clades cover the set of all taxa in the
  // DAG).
  size_t RootsplitCount() const;
  // The total number of edges in the DAG (excluding edges which terminate at a root or
  // leaf node).
  size_t EdgeCount() const;
  // The total number of edges in the DAG (including edges which terminat at a root of
  // leaf node).
  size_t EdgeCountWithLeafSubsplits() const;
  // The total number of tree topologies expressable by the DAG.
  double TopologyCount() const;

  // ** Print

  // Print all nodes:
  // - One line is given for (node_id | subsplit_bitset).
  // - One line is given for each combination of leafward/rootward, leftside/rightside,
  // for all adjacent node_ids.
  void Print() const;
  // Print all nodes as (node_id | node_bitsets) pairs, one-per-line.
  void PrintNodes() const;
  // Print all edges/PCSP, as (pcsp_bitset | edge_idx) pairs, one-per-line.
  void PrintEdgeIndexer() const;
  // Print all edges/PCSP, as (parent_node_id | child_node_id) pairs, one-per-line.
  void PrintDAGEdges() const;
  // Print all nodes, as (bitset | range_begin | range_end)
  void PrintParentToRange() const;
  // Output DOT format graph of DAG to file.
  void ToDot(const std::string file_path, bool show_index_labels = true) const;
  // Output DOT format graph of DAG to a string.
  std::string ToDot(bool show_index_labels = true) const;

  // ** Build Indexers/Vectors

  // Create a EdgeIndexer representing the DAG.
  // The EdgeIndexer is a map (edge/PCSP bitset -> edge/PCSP index).
  // The edge/PCSP indexer is "expanded", meaning it contains leafs and rootsplits.
  BitsetSizeMap BuildEdgeIndexer() const;
  // Get the leftside and rightside nodes ids of all nodes adjacent to given node
  // (rootward or leafward as specified).
  std::pair<SizeVector, SizeVector> BuildAdjacentNodeIdVector(
      const Bitset &subsplit) const;
  // Get the rotated and sorted parents of the node with the given subsplit.
  std::pair<SizeVector, SizeVector> BuildParentIdVector(const Bitset &subsplit) const;
  // Get the rotated and sorted children of the node with the given subsplit.
  std::pair<SizeVector, SizeVector> BuildChildIdVector(const Bitset &subsplit) const;
  // Output RootedIndexerRepresentation of DAG (from RootedSBNMaps).
  // RootedIndexerRepresentation is a vector of edge idxs in topological preorder.
  RootedIndexerRepresentation IndexerRepresentationOf(const BitsetSizeMap &indexer,
                                                      const Node::NodePtr &topology,
                                                      size_t default_index) const;
  // Each node in a topology is constructed with SubsplitDAGNode ID as Node ID.
  Node::NodePtrVec GenerateAllTopologies() const;

  // ** Getter

  // Get Taxon's bitset clade positional id.
  size_t GetTaxonId(const std::string &name) const;
  // Get node based on node id.
  SubsplitDAGNode *GetDAGNode(const size_t node_id) const;
  // Get the Node ID based on the subsplit bitset.
  size_t GetNodeId(const Bitset &subsplit) const;
  // Gets the Node ID of the DAG root.
  size_t GetRootNodeId() const;
  // Return the Node ids corresponding to the rootsplits.
  const SizeVector &GetRootsplitIds() const;
  // Get the Edge index by its parent-child pair subsplits from the DAG nodes.
  size_t GetEdgeIdx(const Bitset &parent_subsplit, const Bitset &child_subsplit) const;
  // Get the Edge index by its parent-child pair id from the DAG nodes.
  size_t GetEdgeIdx(size_t parent_id, size_t child_id) const;
  // Get the range of outgoing Edge indexes from the given clade of a subsplit.
  SizePair GetEdgeRange(const Bitset &subsplit, const bool rotated) const;
  // Get a vector of all node ids which are adjacent to the given bitset.
  SizeVector GetNodesAdjacentTo(const Bitset &node_subsplit) const;
  SizeVector GetNodesAdjacentTo(const size_t node_id) const;
  // Get set of all taxon names.
  std::vector<std::string> GetSortedVectorOfTaxonNames() const;
  // Get set of all node Subsplit bitsets.
  std::vector<Bitset> GetSortedVectorOfNodeBitsets() const;
  // Get set of all edge PCSP bitsets.
  std::vector<Bitset> GetSortedVectorOfEdgeBitsets() const;
  // Get reference to taxon map.
  std::map<std::string, size_t> &GetTaxonMap();
  // Get subsplit to id map.
  BitsetSizeMap &GetSubsplitToIdMap();

  // ** DAG Lambda Iterator
  // These methods iterate over the nodes and take lambda functions with arguments
  // relative to current node.

  // NodeLambda is for iterating over nodes.
  using NodeLambda = std::function<void(const SubsplitDAGNode *)>;
  // EdgeDestinationLambda takes in a rotation status (true is rotated, false is not)
  // and a "destination" node. For iterating over DAG edges with a rotation status.
  using EdgeDestinationLambda = std::function<void(bool, const SubsplitDAGNode *)>;
  // EdgeAndNodeLambda takes a PCSP index of an edge, its rotation status, and an index
  // of the node on the other side of the edge.
  using EdgeAndNodeLambda = std::function<void(const size_t, const bool, const size_t)>;
  // ParentEdgeChildLambda takes: the parent id in the DAG, the rotation status of the
  // edge, the child id, and the GCPSP index of the edge.
  using ParentRotationChildEdgeLambda =
      std::function<void(const size_t, const bool, const size_t, const size_t)>;
  //
  // Iterate over the "real" nodes, i.e. those that do not correspond to
  // leaf subsplits or the DAG root node.
  void IterateOverRealNodes(const NodeLambda &f) const;
  // Iterate over the all leafward edges, rotated and sorted, of node using an
  // EdgeDestinationLambda.
  void IterateOverLeafwardEdges(const SubsplitDAGNode *node,
                                const EdgeDestinationLambda &f) const;
  // Iterate over only the rotated/sorted leafward edges of node using a NodeLambda.
  void IterateOverLeafwardEdges(const SubsplitDAGNode *node, bool rotated,
                                const NodeLambda &f) const;
  // Iterate over the leafward edges, supplying both the a PCSP index of the edge and
  // the SubsplitDAGNode of the corresponding child.
  void IterateOverLeafwardEdgesAndChildren(const SubsplitDAGNode *node,
                                           const EdgeAndNodeLambda &f) const;
  // Iterate over the rootward edges of node using an EdgeDestinationLambda.
  // Excludes edges to the DAG root node.
  void IterateOverRootwardEdges(const SubsplitDAGNode *node,
                                const EdgeDestinationLambda &f) const;
  // Iterate over the rootward edges, supplying both the a PCSP index of the edge and
  // the SubsplitDAGNode of the corresponding child.
  void IterateOverRootwardEdgesAndParents(const SubsplitDAGNode *node,
                                          const EdgeAndNodeLambda &f) const;
  // Iterate over the leafward edges, supplying the parent node id, child node id,
  // rotation of child, and the PCSP index of the rootward edge connecting the two.
  void IterateOverParentAndChildAndLeafwardEdges(
      const SubsplitDAGNode *node, const ParentRotationChildEdgeLambda &f) const;

  // ** DAG Traversals
  // These perform a traversal of the DAG and takes an argument TraversalActionT.
  // TraversalActionT passes four functions of the following form and order:
  // (1)  void BeforeNode(size_t node_id)
  // (2)  void AfterNode(size_t node_id)
  // (3)  void BeforeNodeClade(size_t node_id, bool rotated)
  // (4)  void VisitEdge(size_t node_id, size_t child_id, bool rotated)
  // See below for order of usage.

  // Apply a TraversalAction via a depth first traversal. Do not visit leaf nodes.
  // Applied to a given node, we:
  // - Apply BeforeNode()
  // - For each of the clades of the node, we:
  //     - Apply BeforeNodeClade()
  //     - For each edge descending from that clade, we:
  //         - Recur into the child node of the clade if it is not a leaf
  //         - Apply VisitEdge() to the edge
  // - Apply AfterNode()
  template <typename TraversalActionT>
  void DepthFirstWithAction(const SizeVector &starting_nodes,
                            const TraversalActionT &action) const {
    std::unordered_set<size_t> visited_nodes;
    for (const auto &node_id : starting_nodes) {
      DepthFirstWithActionForNode(action, node_id, visited_nodes);
    }
  };
  // The portion of the traversal that is below a given node.
  template <typename TraversalActionT>
  void DepthFirstWithActionForNode(const TraversalActionT &action, size_t node_id,
                                   std::unordered_set<size_t> &visited_nodes) const {
    action.BeforeNode(node_id);
    DepthFirstWithActionForNodeClade(action, node_id, false, visited_nodes);
    DepthFirstWithActionForNodeClade(action, node_id, true, visited_nodes);
    action.AfterNode(node_id);
  };
  // The portion of the traversal that is below a given clade of a given node.
  // Does not recurse into leaf nodes.
  template <typename TraversalActionT>
  void DepthFirstWithActionForNodeClade(
      const TraversalActionT &action, size_t node_id, bool rotated,
      std::unordered_set<size_t> &visited_nodes) const {
    action.BeforeNodeClade(node_id, rotated);
    const auto node = GetDAGNode(node_id);
    for (const size_t child_id : node->GetLeafward(rotated)) {
      if (visited_nodes.count(child_id) == 0) {
        visited_nodes.insert(child_id);
        if (!GetDAGNode(child_id)->IsLeaf()) {
          DepthFirstWithActionForNode(action, child_id, visited_nodes);
        }
      }
      action.VisitEdge(node_id, child_id, rotated);
    }
  };

  // ** Edge Index Traversal Traces
  // #288: consider renaming these re leafward and rootward.
  // Also, note that they could be called "TraversalTrace" to signify that they are
  // recording the trace of a traversal, not doing actual work.
  [[nodiscard]] SizeVector LeafwardPassTraversal(bool include_dag_root_node) const;
  [[nodiscard]] SizeVector RootwardPassTraversal(bool include_dag_root_node) const;
  // Note: This is leafward topological ordering.
  [[nodiscard]] SizeVector ReversePostorderTraversal() const;
  // Do a reverse postorder traversal on the edges of the DAG, including edges from the
  // DAG root node to the rootsplits, supplying the relevant indices to a lambda.
  void ReversePostorderIndexTraversal(ParentRotationChildEdgeLambda f) const;

  // ** Statistical

  // Discrete uniform distribution over each subsplit.
  [[nodiscard]] EigenVectorXd BuildUniformQ() const;
  // Uniform prior over all topologies in the subsplit support.
  [[nodiscard]] EigenVectorXd BuildUniformOnTopologicalSupportPrior() const;
  // Uniform prior over all topologies, whether or not they are in the support.
  // Thus, this will only be a normalized probability distribution for each subsplit if
  // all topologies are in the support.
  [[nodiscard]] EigenVectorXd BuildUniformOnAllTopologiesPrior() const;
  // Get a vector from each DAG node index to the probability of sampling that DAG node
  // with the supplied SBN parameters. These SBN parameters must be indexed in a
  // compatible way as the dag_edges_ of the subsplit DAG.
  EigenVectorXd UnconditionalNodeProbabilities(
      EigenConstVectorXdRef normalized_sbn_parameters) const;
  // Get a map from each non-leaf subsplit to the probability of observing that
  // subsplit with the supplied SBN parameters. See
  // UnconditionalSubsplitProbabilityVector for notes.
  BitsetDoubleMap UnconditionalSubsplitProbabilities(
      EigenConstVectorXdRef normalized_sbn_parameters) const;
  // Get a vector from each PCSP index to the Bayes-inverted probability of sampling
  // the parent given the child.
  EigenVectorXd InvertedGPCSPProbabilities(
      EigenConstVectorXdRef normalized_sbn_parameters,
      EigenConstVectorXdRef node_probabilities) const;

  // ** Contains

  // Does a taxon with the given name exist?
  bool ContainsTaxon(const std::string &name) const;
  // Does a node with the given subsplit exist?
  bool ContainsNode(const Bitset &subsplit) const;
  // Does a node with the given id exist?
  bool ContainsNode(const size_t node_id) const;
  // Does an edge that connects the two nodes exist?
  bool ContainsEdge(const size_t parent_id, const size_t child_id) const;

  // ** Modify DAG
  // These methods are for directly modifying the DAG by adding or removing nodes and
  // edges. All these methods promise that all data owned by the SubsplitDAG will be
  // valid and up-to-date at their conclusion. This includes maintaining the following
  // DAG invariants:
  // - Tips have ids 0 to taxon_count_.
  // - Parents have higher ids than their children.
  // - The DAG root node has highest id.
  // - Edges descending from the same node clade have contiguous idxs.
  // - There are no gaps in node ids or edge idxs.
  // Returned is the ModificationResult, which will specify changes and reordering to
  // data that occurred during modification.

  // ModificationResult is the return value of Modification methods.
  // Contains the output needed to update all related data to reflect modifications to
  // DAG.
  struct ModificationResult {
    // Nodes that were added or removed by modification.
    SizeVector added_node_ids;
    // Edges that were added or removed by modification.
    SizeVector added_edge_idxs;
    // New ordering of node ids relative to their ordering before DAG modification.
    SizeVector node_reindexer;
    // New ordering of edge idxs relative to their ordering before DAG modification.
    SizeVector edge_reindexer;
  };

  // Add an adjacent node pair to the DAG.
  ModificationResult AddNodePair(const Bitset &parent_subsplit,
                                 const Bitset &child_subsplit);
  // Remove an adjacent node pair to the DAG.
  ModificationResult RemoveNodePair(const Bitset &parent_supsplit,
                                    const Bitset &child_subsplit);
  // Add all tree topologies in topology_counter to DAG.
  std::tuple<BitsetSizeMap, SizeBitsetMap, BitsetVector> ProcessTopologyCounter(
      const Node::TopologyCounter &topology_counter);

  // ** Reindexer
  // These methods are for building and applying reindexers.
  // A reindexer describes a remapping/transform of identifiers (i.e. node ids, edge
  // idxs) before and after a modification to the DAG. The index is before the
  // modification, the value is after the modification. Reindexers are used to perform a
  // transform on data vectors,

  // Uses a depth first DAG traversal to build a reindexer for node_ids.
  SizeVector BuildNodeReindexer(const size_t prev_node_count);
  // Get the reindexer for edges idxs.
  SizeVector BuildEdgeReindexer(const size_t prev_edge_count);
  // Remap all node ids according to the node_reindexer.
  void RemapNodeIds(const SizeVector &node_reindexer);
  // Remap all edge idxs according to the edge_reindexer.
  void RemapEdgeIdxs(const SizeVector &edge_reindexer);

  // ** Validation Tests
  // These methods are used to assert that a DAG is in a valid state or that given
  // operation will result in a valid DAG.

  // Checks if SubsplitDAG's corresponding data is consistent and up-to-date.
  // Specifically, checks that:
  // - subsplit_to_id_ map consistent with nodes in dag_nodes_.
  // - parent_to_child_range_ map consistent with nodes in dag_nodes_'s (adjacent.
  // - parent_to_child_range_ map matches each edge in dag_edges_.
  bool IsConsistent() const;
  // Checks if SubsplitDAG is in a valid state (assumes that DAG is consistent).
  // Specifically, checks that:
  // - Each node is valid. That is, either:
  //    - Node has zero parents and zero children.
  //    - Node has 1+ parents, 1+ sorted children, and 1+ rotated children.
  // -
  bool IsValid() const;
  // Check if it is valid to add given node pair. Specifically, check that:
  // - The nodes are adjacent.
  // - The nodes do not add/remove any taxa.
  // - The parent node has at least one parent.
  // - Including the child node, each clade of the parent node has at least one child.
  // - Each clade of the child node has at least 1 child.
  bool IsValidAddNodePair(const Bitset &parent_subsplit,
                          const Bitset &child_subsplit) const;
  // Check if it is valid to remove given node pair. Specifically, check that:
  // - The nodes are adjacent.
  // - The nodes both exist in the DAG.
  // - For parents and children of all removed nodes contain alternate parents and
  // children (so, more than one edge of that type).
  //   (e.g., if node X is being removed from DAG and is a sorted child of Y, then Y
  //   needs to have more than one outgoing sorted edge, otherwise the removal of X
  //   will make Y an invalid node.)
  bool IsValidRemoveNodePair(const Bitset &parent_subsplit,
                             const Bitset &child_subsplit) const;

  // ** Miscellaneous methods:

  // Rotates Node Subsplit only if it is out of sorted order (rotated).
  Bitset SubsplitToSortedOrder(const Bitset &subsplit, bool rotated) const;
  // Fixes Edge PCSP only if child is out of sorted order (rotated).
  Bitset EdgeToSortedOrder(const Bitset &pcsp) const;
  // Build vector between from SubsplitDAGs dag_a to dag_b corresponding to their taxon
  // ids. Can be treated as a "map" with indices representing keys. Requires that both
  // SubsplitDAGs use the same taxon set.
  // - Map: [ index: dag_a's taxon_id => value: dag_b's taxon_id ]
  static SizeVector BuildTaxonTranslationMap(const SubsplitDAG &dag_a,
                                             const SubsplitDAG &dag_b);
  // Compare two bitsets (subsplits or PCSPs) from two different DAGs using a taxon
  // map.
  static int BitsetCompareViaTaxonTranslationMap(const Bitset &bitset_a,
                                                 const Bitset &bitset_b,
                                                 const SizeVector &taxon_map);
  // Translate bitset using the taxon translation map output positions.
  // NOTE:  forward_translate uses index in map as input bit locations, values in map as
  // output bit locations.
  //        If false, map is used vice versa.
  static Bitset BitsetTranslateViaTaxonTranslationMap(
      const Bitset &bitset, const SizeVector &taxon_map,
      const bool forward_translate = true);

 protected:
  // ** Constructor methods:

  // Build a Subsplit DAG on given number of taxa, expressing all tree topologies from
  // tree_collection, with trees on the given taxa names/labels.
  SubsplitDAG(size_t taxon_count, const Node::TopologyCounter &topology_counter,
              const TagStringMap &tag_taxon_map);

  // Builds a vector of subsplits of all children, optionally including leaf nodes.
  std::vector<Bitset> GetChildSubsplits(const SizeBitsetMap &index_to_child,
                                        const Bitset &subsplit,
                                        bool include_leaf_subsplits = false);

  // ** DAG Traversal methods:
  // NOTES: - visit_order is the output vector of node ids in specified traversal order.
  //        - visited_nodes is a set of all node ids reached by traversal.

  // Creates vector of node ids in leafward depth first post-order traversal of DAG.
  void LeafwardDepthFirst(size_t node_id, SizeVector &visit_order,
                          std::unordered_set<size_t> &visited_nodes) const;
  // Create vector of node ids in rootward depth first post-order traversal of DAG.
  void RootwardDepthFirst(size_t node_id, SizeVector &visit_order,
                          std::unordered_set<size_t> &visited_nodes) const;

  // ** Modify DAG Helper methods:
  // These methods help calling functions to modify the DAG, but do NOT ensure a valid
  // state at the end of the function. Do not necessarily handle remapping ids and idxs,
  // removing old references to deleted objects, etc.

  // Create Node and insert it into the DAG.
  // WARNING: Does not check for redundant nodes.
  void CreateAndInsertNode(const Bitset &subsplit);
  // Create Edge between given nodes and insert it into the DAG.
  // WARNING: Does not check for redundant edges.
  void CreateAndInsertEdge(const size_t parent_id, const size_t child_id, bool rotated);
  // TODO:
  // Delete Node from the DAG.
  void DeleteNode(const Bitset &subsplit);
  // Delete Edge from the DAG.
  void DeleteEdge(const size_t parent_id, const size_t child_id, bool rotated);

  // Reorder nodes so that DAG nodes sorted in a topological ordering. 
  // Returns a node reindexer. If optional node_reindexer is passed, then 
  // that is used for the initial base reindexer, otherwise uses identity.
  SizeVector SortNodes(std::optional<SizeVector> node_reindexer = std::nullopt);

  // Add edge relationship to each node's respective SubsplitDAGNode.
  void ConnectGivenNodes(const size_t parent_id, const size_t child_id, bool rotated);
  // Add edge between node_id and nodes according to index to child.
  void ConnectNodes(const SizeBitsetMap &index_to_child, size_t node_id, bool rotated);
  //
  void BuildNodes(const SizeBitsetMap &index_to_child, const BitsetVector &rootsplits);
  //
  void BuildNodesDepthFirst(const SizeBitsetMap &index_to_child, const Bitset &subsplit,
                            std::unordered_set<Bitset> &visited_subsplits);
  //
  void BuildEdges(const SizeBitsetMap &index_to_child);
  //
  void BuildDAGEdgesFromEdgeIndexer(BitsetSizeMap &edge_indexer);
  // Adds edges to DAG and to DAGNodes, connecting focal node to every node in
  // adjacent_node_ids, with the option to ignore one node id (will not create connect
  // edge).
  void ConnectNodeToAdjacentNodes(
      const size_t focal_node_id, const SizeVector adjacent_node_ids,
      const bool is_main_node_parent, const bool is_left_child,
      std::optional<SizeVector> added_edge_idxs = std::nullopt,
      const std::optional<size_t> ignored_node_id_opt = std::nullopt);
  // Add edges to DAG and to DAGNodes, connecting focal node to all parent nodes or all
  // child nodes.
  void ConnectNodeToChildOrParentNodes(
      const size_t focal_node_id, const bool is_focal_parent,
      std::optional<SizeVector> added_edge_idxs = std::nullopt,
      const std::optional<size_t> ignored_node_id_opt = std::nullopt);
  // Connect the child to all of its children. Push all new edges to
  // added_edge_idxs.
  void ConnectChildToAllChildren(const Bitset &child_subsplit,
                                 SizeVector &added_edge_idxs);
  // Connect the parent to all of its children except for the given child node. Insert
  // all new edges to added_edge_idxs vector.
  void ConnectParentToAllChildrenExcept(const Bitset &parent_subsplit,
                                        const Bitset &child_subsplit,
                                        SizeVector &added_edge_idxs);
  // Connect the child to all of its parents except for the given parent node. Insert
  // all new edge to added_edge_idxs vector.
  void ConnectChildToAllParentsExcept(const Bitset &parent_subsplit,
                                      const Bitset &child_subsplit,
                                      SizeVector &added_edge_idxs);
  // Connect the parent to all of its parents. Insert all new edges to
  // added_edge_idxs vector.
  void ConnectParentToAllParents(const Bitset &parent_subsplit,
                                 SizeVector &added_edge_idxs);
  // Expand dag_edges_ and parent_to_child_range_ with leaf subsplits at the end.
  void AddLeafSubsplitsToDAGEdgesAndParentToRange();
  // To be called after modifying the DAG.
  // Traverses the DAG and refreshes count of the number of topologies contained in the
  // DAG. Updates topology_count_ and topology_count_below_.
  void CountTopologies();

  // ** Clades

  // Initialize clades
  void InitClades();
  // Add both of node's clades to the clade map.
  void AddNodeToClades(const size_t node_id, const Bitset &node_subsplit);
  // Remove both of node's clades from the clade map.
  void RemoveNodeFromClades(const size_t node_id, const Bitset &node_subsplit);
  // Get all nodes containing clade.
  SizeVector GetAllNodesContainingClade(const Bitset &subsplit_clade) const;
  // Get all the child nodes of given subsplit (specify left or right children).
  SizeVector GetAllChildrenOfNode(const Bitset &node_subsplit,
                                  const bool which_child) const;
  // Get all parent nodes of a given subsplit.
  SizeVector GetAllParentsOfNode(const Bitset &node_subsplit) const;

 protected:
  // NOTE: When using unique identifiers, for DAG nodes (aka Subsplits) we use the term
  // "ids", and for edges (aka PCSPs) we use the term index or "idx", to more easily
  // distinguish the two. This corresponds to the analogous concept for topologies.

  // - Map of Taxon Names
  //    - [ Taxon Name => Taxon Id (Bitset Clade Position) ]
  std::map<std::string, size_t> dag_taxons_;
  // - Vector of DAG Nodes: each Node index in the vector corresponds to their Id.
  // - NOTE: For all parent/child ids, child_id < parent_id.
  // - This can be viewed as a map:
  //    - [ Node Id => Pointer to SubsplitDAGNode ]
  std::vector<std::unique_ptr<SubsplitDAGNode>> dag_nodes_;
  // TODO: Rename to edges_idpair_to_idx_?
  // - Map of all DAG Edges:
  //    - [ (parent_id, child_id) Node Id Pairs => Edge/PCSP Idxs. ]
  std::map<SizePair, size_t> dag_edges_;
  // TODO: FIX THIS DESCRIPTION
  // - Map of all DAG Nodes:
  //    - [ Node Subsplit (Bitset) => Node Id ]
  // A node's id is equivalent to its index in dag_nodes_. The first entries are
  // reserved for leaf subsplits. The last entries are reserved for rootsplits. The DAG
  // root node has the highest node id.
  BitsetSizeMap subsplit_to_id_;
  // TODO: PROPOSED ADDITION
  // - Map of all Clades from DAG Nodes:
  //    - [ Node Clade (Bitset) => Id Vector of Nodes whose Subsplits contain Clade ]
  std::map<Bitset, SizeVector> clade_to_ids_;
  // - Map of all DAG Nodes:
  //    - [ Node Subsplit (Bitset) => [begin, end) Range of Child Ids ]
  // This indexer is an expanded version of parent_to_child_range_ in sbn_instance:
  // It includes single element range for leaf subsplits.
  BitsetSizePairMap parent_to_child_range_;
  // The number of taxa in the DAG. This is equivalent to the size of the clades in each
  // subsplit. Also equivalent to the number of leaf nodes in the DAG.
  size_t taxon_count_;
  // Number of internal edges in the DAG (excludes all edges that go to a root or
  // leaf).
  size_t edge_count_without_leaf_subsplits_;
  // Total number of tree topologies spanned by the DAG.
  double topology_count_;
  // Storage for the number of topologies possible below for each node. Each index maps
  // to the topology count for the corresponding node_id.
  EigenVectorXd topology_count_below_;
};

#endif  // SRC_SUBSPLIT_DAG_HPP_

#endif  // __SUBSPLIT_DAG_H__
