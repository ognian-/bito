// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag.hpp"

#include "combinatorics.hpp"
#include "numerical_utils.hpp"
#include "sbn_probability.hpp"

// ** Constructor methods:

SubsplitDAG::SubsplitDAG()
    : taxon_count_(0), gpcsp_count_without_fake_subsplits_(0), topology_count_(0.) {}

// TODO: Add Taxon names to the constructor.
SubsplitDAG::SubsplitDAG(size_t taxon_count,
                         const Node::TopologyCounter &topology_counter,
                         const TagStringMap &tag_taxon_map)
    : taxon_count_(taxon_count) {
  Assert(topology_counter.size() > 0, "Empty topology counter given to SubsplitDAG.");
  Assert(topology_counter.begin()->first->LeafCount() == taxon_count,
         "Taxon count mismatch in SubsplitDAG constructor.");
  auto [gpcsp_indexer, index_to_child, rootsplits] =
      ProcessTopologyCounter(topology_counter);
  BuildNodes(index_to_child, rootsplits);
  BuildEdges(index_to_child);
  BuildDAGEdgesFromGPCSPIndexer(gpcsp_indexer);
  AddFakeSubsplitsToDAGEdgesAndParentToRange();
  CountTopologies();

  dag_taxons_ = std::map<std::string, size_t>();
  // Insert all taxons from tree_collections's map to SubsplitDAG map.
  for (const auto &[tag, name] : tag_taxon_map) {
    // The "tag" key of the tree_collection's taxon_map is 2 bitpacked ints: [id,
    // topology count]. We only care about the id.
    size_t id = static_cast<size_t>(UnpackFirstInt(tag));
    dag_taxons_.insert(std::make_pair(name, id));
  }
}

SubsplitDAG::SubsplitDAG(const RootedTreeCollection &tree_collection)
    : SubsplitDAG(tree_collection.TaxonCount(), tree_collection.TopologyCounter(),
                  tree_collection.TagTaxonMap()) {}

// ** Comparator methods:

int SubsplitDAG::Compare(const SubsplitDAG &lhs, const SubsplitDAG &rhs) {
  // (1) Compare Taxon Sizes.
  if (lhs.TaxonCount() != rhs.TaxonCount()) {
    return -100;
  }
  // Create translation map (lhs->rhs) for bitset clades.
  auto taxon_map = SubsplitDAG::BuildTaxonTranslationMap(lhs, rhs);
  // (2) Compare Subsplit Nodes.
  auto lhs_nodes = lhs.GetSortedVectorOfNodeBitsets();
  auto rhs_nodes = rhs.GetSortedVectorOfNodeBitsets();
  // Translate to account for different Taxon mappings and sort output.
  for (size_t i = 0; i < lhs_nodes.size(); i++) {
    lhs_nodes[i] =
        SubsplitDAG::BitsetTranslateViaTaxonTranslationMap(lhs_nodes[i], taxon_map);
    // Verify subsplit clades are in 'sorted' order.
    lhs_nodes[i] = lhs_nodes[i].SubsplitSort();
  }
  std::sort(lhs_nodes.begin(), lhs_nodes.end());
  if (lhs_nodes != rhs_nodes) {
    return -200;
  }
  // (3) Compare PCSP Edges.
  auto lhs_edges = lhs.GetSortedVectorOfEdgeBitsets();
  auto rhs_edges = rhs.GetSortedVectorOfEdgeBitsets();
  // Translate to account for different Taxon mappings and sort output.
  for (size_t i = 0; i < lhs_edges.size(); i++) {
    lhs_edges[i] =
        SubsplitDAG::BitsetTranslateViaTaxonTranslationMap(lhs_edges[i], taxon_map);
    lhs_edges[i] = lhs_edges[i].PCSPSort();
  }
  std::sort(lhs_edges.begin(), lhs_edges.end());
  if (lhs_edges != rhs_edges) {
    return -300;
  }
  return 0;
}

bool operator==(const SubsplitDAG &lhs, const SubsplitDAG &rhs) {
  return (SubsplitDAG::Compare(lhs, rhs) == 0);
}

bool operator!=(const SubsplitDAG &lhs, const SubsplitDAG &rhs) {
  return (SubsplitDAG::Compare(lhs, rhs) != 0);
}

// ** Count methods:

void SubsplitDAG::CountTopologies() {
  topology_count_below_ = EigenVectorXd::Ones(NodeCount());

  for (const auto &node_id : RootwardPassTraversal(true)) {
    const auto &node = GetDAGNode(node_id);
    for (const bool rotated : {true, false}) {
      // When there are no leafward nodes in the `rotated` direction, we set the number
      // of topologies for the rotation of the node to be 1.
      double per_rotated_count = node->GetLeafward(rotated).empty() ? 1. : 0.;
      // Sum options across the possible children.
      for (const auto &child_id : node->GetLeafward(rotated)) {
        per_rotated_count += topology_count_below_[child_id];
      }
      // Take the product across the number of options for the left and right branches
      // of the tree.
      topology_count_below_[node_id] *= per_rotated_count;
    }
  }
  topology_count_ = topology_count_below_[DAGRootNodeId()];
}

size_t SubsplitDAG::TaxonCount() const { return dag_taxons_.size(); }

size_t SubsplitDAG::NodeCount() const { return dag_nodes_.size(); }

size_t SubsplitDAG::NodeCountWithoutDAGRoot() const { return NodeCount() - 1; }

double SubsplitDAG::TopologyCount() const { return topology_count_; }

size_t SubsplitDAG::RootsplitCount() const { return RootsplitIds().size(); }

size_t SubsplitDAG::GPCSPCount() const { return gpcsp_count_without_fake_subsplits_; }

size_t SubsplitDAG::GPCSPCountWithFakeSubsplits() const { return dag_edges_.size(); }

// ** Output methods:

void SubsplitDAG::Print() const {
  for (const auto &dag_node : dag_nodes_) {
    std::cout << dag_node->ToString() << std::endl;
  }
}

void SubsplitDAG::PrintNodes() const {
  for (const auto &dag_node : dag_nodes_) {
    std::cout << dag_node->Id() << ": " << dag_node->GetBitset().SubsplitToString()
              << std::endl;
  }
}

void SubsplitDAG::PrintGPCSPIndexer() const {
  for (const auto &[pcsp, idx] : BuildGPCSPIndexer()) {
    std::cout << idx << ": " << pcsp.PCSPToString() << std::endl;
  }
}

void SubsplitDAG::PrintDAGEdges() const {
  for (const auto &[parent_child_id, gpcsp_idx] : dag_edges_) {
    const auto &[parent_id, child_id] = parent_child_id;
    std::cout << "[" << parent_id << ", " << child_id << "], " << gpcsp_idx
              << std::endl;
  }
}

void SubsplitDAG::PrintParentToRange() const {
  for (const auto &[subsplit, pcsp_range] : parent_to_child_range_) {
    std::cout << subsplit.SubsplitToString() << ": [" << pcsp_range.first << ", "
              << pcsp_range.second << "]" << std::endl;
  }
}

void SubsplitDAG::ToDot(const std::string file_path, bool show_index_labels) const {
  std::ofstream out_file(file_path);
  out_file << ToDot(show_index_labels);
  out_file.close();
}

std::string SubsplitDAG::ToDot(bool show_index_labels) const {
  std::stringstream string_stream;
  string_stream << "digraph g {\n";
  string_stream << "node [shape=record];\n";
  string_stream << "edge [colorscheme=dark23];\n";
  DepthFirstWithAction(
      {DAGRootNodeId()},
      SubsplitDAGTraversalAction(
          // BeforeNode
          [this, &string_stream, &show_index_labels](size_t node_id) {
            const auto &node = GetDAGNode(node_id);
            if (node->IsDAGRootNode()) {
              string_stream << node_id << " [label=\"<f0>&rho;\"]\n";
              return;
            }
            auto bs = node->GetBitset();
            string_stream << node_id << " [label=\"<f0>"
                          << bs.SubsplitGetClade(0).ToVectorOfSetBitsAsString()
                          << "|<f1>";
            if (show_index_labels) {
              string_stream << node_id;
            }
            string_stream << "|<f2>"
                          << bs.SubsplitGetClade(1).ToVectorOfSetBitsAsString()
                          << "\"]\n";
          },
          // AfterNode
          [](size_t node_id) {},
          // BeforeNodeClade
          [](size_t node_id, bool rotated) {},
          // VisitEdge
          [this, &string_stream, &show_index_labels](size_t node_id, size_t child_id,
                                                     bool rotated) {
            if (GetDAGNode(child_id)->IsLeaf()) {
              string_stream << child_id << " [label=\"<f1>" << child_id << "\"]\n";
            }
            string_stream << "\"" << node_id << "\":";
            string_stream << (rotated ? "f0" : "f2");
            string_stream << "->\"";
            string_stream << child_id << "\":f1";
            if (show_index_labels) {
              string_stream << " [label=\"" << GetGPCSPEdgeIdx(node_id, child_id);
              if (rotated) {
                string_stream << "\", color=1, fontcolor=1";
              } else {
                string_stream << "\", color=3, fontcolor=3";
              }
              if (GetDAGNode(node_id)->IsDAGRootNode()) {
                string_stream << ",style=dashed]";
              } else {
                string_stream << "]";
              }
            } else {
              if (GetDAGNode(node_id)->IsDAGRootNode()) {
                string_stream << "[style=dashed]";
              }
            }
            string_stream << "\n";
          }));
  string_stream << "}";
  return string_stream.str();
}

// ** Build output indexes/vectors methods:

BitsetSizeMap SubsplitDAG::BuildGPCSPIndexer() const {
  auto gpcsp_indexer = BitsetSizeMap();
  ReversePostorderIndexTraversal([this, &gpcsp_indexer](size_t parent_id, bool rotated,
                                                        size_t child_id,
                                                        size_t gpcsp_idx) {
    const auto parent_subsplit = GetDAGNode(parent_id)->GetBitset(rotated);
    const auto child_subsplit = GetDAGNode(child_id)->GetBitset();
    SafeInsert(gpcsp_indexer, Bitset::PCSP(parent_subsplit, child_subsplit), gpcsp_idx);
  });
  return gpcsp_indexer;
}

size_t SubsplitDAG::GetTaxonId(const std::string &name) const {
  Assert(
      ContainsTaxon(name),
      "SubsplitDAG::GetTaxonId(): Taxon with given name is not contained in the DAG.");
  return dag_taxons_.find(name)->second;
}

SubsplitDAGNode *SubsplitDAG::GetDAGNode(const size_t node_id) const {
  Assert(ContainsNode(node_id),
         "Node with the given node_id does not exist in SubsplitDAG::GetDAGNode().");
  return dag_nodes_.at(node_id).get();
}

size_t SubsplitDAG::GetDAGNodeId(const Bitset &subsplit) const {
  Assert(ContainsNode(subsplit),
         "Node with the given subsplit does not exist in SubsplitDAG::GetDAGNodeId().");
  return subsplit_to_id_.at(subsplit);
}

size_t SubsplitDAG::DAGRootNodeId() const { return NodeCount() - 1; }

const SizeVector &SubsplitDAG::RootsplitIds() const {
  return GetDAGNode(DAGRootNodeId())->GetLeafwardRotated();
}

size_t SubsplitDAG::GetGPCSPEdgeIdx(const Bitset &parent_subsplit,
                                    const Bitset &child_subsplit) const {
  return GetGPCSPEdgeIdx(GetDAGNodeId(parent_subsplit), GetDAGNodeId(child_subsplit));
}

size_t SubsplitDAG::GetGPCSPEdgeIdx(size_t parent_id, size_t child_id) const {
  return dag_edges_.at({parent_id, child_id});
}

SizePair SubsplitDAG::GetEdgeRange(const Bitset &subsplit, const bool rotated) const {
  Assert(ContainsNode(subsplit),
         "Node with the given subsplit does not exist in SubsplitDAG::GetEdgeRange.");
  return parent_to_child_range_.at(SubsplitToSortedOrder(subsplit, rotated));
}

std::vector<std::string> SubsplitDAG::GetSortedVectorOfTaxonNames() const {
  std::vector<std::string> taxons;
  for (const auto &[name, id] : dag_taxons_) {
    taxons.push_back(name);
  }
  std::sort(taxons.begin(), taxons.end());
  return taxons;
}

std::vector<Bitset> SubsplitDAG::GetSortedVectorOfNodeBitsets() const {
  std::vector<Bitset> nodes;
  for (size_t i = 0; i < NodeCount(); i++) {
    Bitset node_bitset = GetDAGNode(i)->GetBitset();
    nodes.push_back(node_bitset);
  }
  std::sort(nodes.begin(), nodes.end());
  return nodes;
}

std::vector<Bitset> SubsplitDAG::GetSortedVectorOfEdgeBitsets() const {
  std::vector<Bitset> edges;
  for (const auto &key_value_pair : dag_edges_) {
    auto id_pair = key_value_pair.first;
    auto parent_bitset = GetDAGNode(id_pair.first)->GetBitset();
    auto child_bitset = GetDAGNode(id_pair.second)->GetBitset();
    Bitset edge_bitset = Bitset::PCSP(parent_bitset, child_bitset);
    edges.push_back(edge_bitset);
  }
  std::sort(edges.begin(), edges.end());
  return edges;
}

std::map<std::string, size_t> &SubsplitDAG::GetTaxonMap() { return dag_taxons_; }

EigenVectorXd SubsplitDAG::BuildUniformOnTopologicalSupportPrior() const {
  EigenVectorXd q = EigenVectorXd::Ones(GPCSPCountWithFakeSubsplits());

  for (const auto &node_id : RootwardPassTraversal(true)) {
    const auto &node = GetDAGNode(node_id);
    for (const bool rotated : {false, true}) {
      if (!node->GetLeafward(rotated).empty()) {
        double per_rotated_count = 0.;
        for (const auto &child_id : node->GetLeafward(rotated)) {
          per_rotated_count += topology_count_below_[child_id];
        }
        for (const auto &child_id : node->GetLeafward(rotated)) {
          size_t gpcsp_idx = GetGPCSPEdgeIdx(node->Id(), child_id);
          q(gpcsp_idx) = topology_count_below_(child_id) / per_rotated_count;
        }
      }
    }
  }
  return q;
}

Node::NodePtrVec SubsplitDAG::GenerateAllTopologies() const {
  std::vector<Node::NodePtrVec> topology_below(NodeCount());

  auto GetSubtopologies = [&topology_below](const SubsplitDAGNode *node) {
    Node::NodePtrVec rotated_subtopologies, sorted_subtopologies;
    for (const bool rotated : {false, true}) {
      for (const auto &child_id : node->GetLeafward(rotated)) {
        for (const auto &subtopology : topology_below.at(child_id)) {
          rotated ? rotated_subtopologies.push_back(subtopology)
                  : sorted_subtopologies.push_back(subtopology);
        }
      }
    }
    return std::make_pair(rotated_subtopologies, sorted_subtopologies);
  };

  auto MergeTopologies = [this](size_t node_id, Node::NodePtrVec &rotated_subtopologies,
                                Node::NodePtrVec &sorted_subtopologies) {
    Node::NodePtrVec topologies;
    for (const auto &rotated_subtopology : rotated_subtopologies) {
      for (const auto &sorted_subtopology : sorted_subtopologies) {
        Node::NodePtr new_topology =
            Node::Join(sorted_subtopology, rotated_subtopology, node_id);
        topologies.push_back(new_topology);
      }
    }
    if (node_id == DAGRootNodeId()) {
      // DAG root node has no `sorted_subtopologies`, so loop above yields empty
      // `topologies` vector.
      return rotated_subtopologies;
    }
    return topologies;
  };

  for (const auto &node_id : RootwardPassTraversal(true)) {
    const auto &node = GetDAGNode(node_id);
    if (node->IsLeaf()) {
      topology_below.at(node_id).push_back(Node::Leaf(node_id));
    } else {
      auto [rotated_topologies, sorted_topologies] = GetSubtopologies(node);
      topology_below[node_id] =
          MergeTopologies(node_id, rotated_topologies, sorted_topologies);
    }
  }

  const auto &topologies = topology_below.at(DAGRootNodeId());
  Assert(topologies.size() == TopologyCount(),
         "The realized number of topologies does not match the expected count.");

  // We return a deep copy of every Polished topology to avoid loops in the pointer
  // structure. Such loops can create problems when we Polish the topologies one at a
  // time: polishing a second topology can change the numbering of a} previous topology.
  // This is checked for in the "GPInstance: GenerateCompleteRootedTreeCollection" test.
  Node::NodePtrVec final_topologies;
  final_topologies.reserve(topologies.size());
  for (auto &topology : topologies) {
    topology->Polish();
    final_topologies.push_back(topology->DeepCopy());
  }

  return final_topologies;
}

EigenVectorXd SubsplitDAG::BuildUniformOnAllTopologiesPrior() const {
  EigenVectorXd result = EigenVectorXd::Zero(GPCSPCountWithFakeSubsplits());
  for (const auto &[parent_child_id, gpcsp_idx] : dag_edges_) {
    const auto &[parent_id, child_id] = parent_child_id;
    size_t child0_taxon_count =
        GetDAGNode(child_id)->GetBitset().SubsplitGetCladeByBinaryOrder(0).Count();
    size_t child1_taxon_count =
        GetDAGNode(child_id)->GetBitset().SubsplitGetCladeByBinaryOrder(1).Count();
    //
    result(gpcsp_idx) = Combinatorics::LogChildSubsplitCountRatio(child0_taxon_count,
                                                                  child1_taxon_count);
  }
  NumericalUtils::Exponentiate(result);
  return result;
}

// ** DAG Lambda Iterator methods:

void SubsplitDAG::IterateOverRealNodes(const NodeLambda &f) const {
  Assert(taxon_count_ < NodeCount(), "No real DAG nodes!");
  for (auto it = dag_nodes_.cbegin() + taxon_count_; it < dag_nodes_.cend() - 1; it++) {
    f((*it).get());
  }
}

void SubsplitDAG::IterateOverLeafwardEdges(const SubsplitDAGNode *node, bool rotated,
                                           const NodeLambda &f) const {
  for (const size_t child_id : node->GetLeafward(rotated)) {
    f(GetDAGNode(child_id));
  }
}

void SubsplitDAG::IterateOverLeafwardEdges(const SubsplitDAGNode *node,
                                           const EdgeDestinationLambda &f) const {
  for (bool rotated : {false, true}) {
    for (const size_t child_id : node->GetLeafward(rotated)) {
      f(rotated, GetDAGNode(child_id));
    }
  }
}

void SubsplitDAG::IterateOverLeafwardEdgesAndChildren(
    const SubsplitDAGNode *node, const EdgeAndNodeLambda &f) const {
  IterateOverLeafwardEdges(
      node, [this, &node, &f](bool rotated, const SubsplitDAGNode *child) {
        f(GetGPCSPEdgeIdx(node->Id(), child->Id()), rotated, child->Id());
      });
}

void SubsplitDAG::IterateOverRootwardEdges(const SubsplitDAGNode *node,
                                           const EdgeDestinationLambda &f) const {
  if (!node->IsRootsplit()) {
    for (bool rotated : {false, true}) {
      for (const size_t parent_id : node->GetRootward(rotated)) {
        f(rotated, GetDAGNode(parent_id));
      }
    }
  }
}

void SubsplitDAG::IterateOverRootwardEdgesAndParents(const SubsplitDAGNode *node,
                                                     const EdgeAndNodeLambda &f) const {
  IterateOverRootwardEdges(
      node, [this, &node, &f](bool rotated, const SubsplitDAGNode *parent) {
        f(GetGPCSPEdgeIdx(parent->Id(), node->Id()), rotated, parent->Id());
      });
}

void SubsplitDAG::IterateOverParentAndChildAndLeafwardEdges(
    const SubsplitDAGNode *node, const ParentRotationChildEdgeLambda &f) const {
  IterateOverLeafwardEdges(
      node, [this, &node, &f](bool rotated, const SubsplitDAGNode *child) {
        f(node->Id(), rotated, child->Id(), GetGPCSPEdgeIdx(node->Id(), child->Id()));
      });
}

RootedIndexerRepresentation SubsplitDAG::IndexerRepresentationOf(
    const BitsetSizeMap &indexer, const Node::NodePtr &topology,
    size_t default_index) const {
  return RootedSBNMaps::IndexerRepresentationOf(indexer, topology, default_index);
}

EigenVectorXd SubsplitDAG::UnconditionalNodeProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters) const {
  EigenVectorXd node_probabilities(NodeCount());
  node_probabilities.setZero();
  node_probabilities[DAGRootNodeId()] = 1.;

  ReversePostorderIndexTraversal([this, &node_probabilities,
                                  &normalized_sbn_parameters](
                                     const size_t parent_id, const bool,
                                     const size_t child_id, const size_t gpcsp_idx) {
    const double child_probability_given_parent = normalized_sbn_parameters[gpcsp_idx];
    Assert(child_probability_given_parent >= 0. && child_probability_given_parent <= 1.,
           "UnconditionalNodeProbabilities: got an out-of-range probability. Are these "
           "normalized and in linear space?");
    const double parent_probability = node_probabilities[parent_id];
    node_probabilities[child_id] += parent_probability * child_probability_given_parent;
  });

  return node_probabilities;
}

BitsetDoubleMap SubsplitDAG::UnconditionalSubsplitProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters) const {
  auto node_probabilities = UnconditionalNodeProbabilities(normalized_sbn_parameters);
  BitsetDoubleMap subsplit_probability_map;
  for (size_t node_id = 0; node_id < node_probabilities.size(); node_id++) {
    const auto &subsplit_bitset = GetDAGNode(node_id)->GetBitset();
    if (node_id != DAGRootNodeId() && !subsplit_bitset.SubsplitIsLeaf()) {
      SafeInsert(subsplit_probability_map, subsplit_bitset,
                 node_probabilities[node_id]);
    }
  }
  return subsplit_probability_map;
}

EigenVectorXd SubsplitDAG::InvertedGPCSPProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters,
    EigenConstVectorXdRef node_probabilities) const {
  EigenVectorXd inverted_probabilities =
      EigenVectorXd(normalized_sbn_parameters.size());
  inverted_probabilities.setOnes();
  ReversePostorderIndexTraversal(
      [this, &node_probabilities, &normalized_sbn_parameters, &inverted_probabilities](
          const size_t parent_id, const bool, const size_t child_id,
          const size_t gpcsp_idx) {
        // The traversal doesn't set the rootsplit probabilities, but those are always 1
        // (there is only one "parent" of a rootsplit).
        if (parent_id != DAGRootNodeId()) {
          // For a PCSP t -> s:
          inverted_probabilities[gpcsp_idx] =         // P(t|s)
              node_probabilities[parent_id] *         // P(t)
              normalized_sbn_parameters[gpcsp_idx] /  // P(s|t)
              node_probabilities[child_id];           // P(s)
        }
      });
  return inverted_probabilities;
}

std::vector<Bitset> SubsplitDAG::GetChildSubsplits(const SizeBitsetMap &index_to_child,
                                                   const Bitset &parent_subsplit,
                                                   bool include_fake_subsplits) {
  std::vector<Bitset> children_subsplits;
  // Add all child subsplit bitsets to vector.
  if (parent_to_child_range_.count(parent_subsplit) > 0) {
    const auto [start, stop] = parent_to_child_range_.at(parent_subsplit);
    for (auto idx = start; idx < stop; idx++) {
      children_subsplits.push_back(index_to_child.at(idx));
    }
  } else if (include_fake_subsplits) {
    // This method is designed to be called before calling
    // AddFakeSubsplitsToDAGEdgesAndParentToRange. In that case, if the second clade
    // of the subsplit is just a single taxon, the subsplit will not map to any value in
    // parent_to_child_range_.
    //
    // But we still need to create and connect to fake subsplits in the DAG. So, here we
    // make a fake child subsplit.
    children_subsplits.push_back(Bitset::FakeChildSubsplit(parent_subsplit));
  }

  return children_subsplits;
}

// ** Modify DAG methods:

std::tuple<BitsetSizeMap, SizeBitsetMap, BitsetVector>
SubsplitDAG::ProcessTopologyCounter(const Node::TopologyCounter &topology_counter) {
  BitsetSizeMap gpcsp_indexer;
  SizeBitsetMap index_to_child;
  BitsetVector rootsplits;
  std::tie(rootsplits, gpcsp_indexer, index_to_child, parent_to_child_range_,
           gpcsp_count_without_fake_subsplits_) =
      SBNMaps::BuildIndexerBundle(RootedSBNMaps::RootsplitCounterOf(topology_counter),
                                  RootedSBNMaps::PCSPCounterOf(topology_counter));
  return {gpcsp_indexer, index_to_child, rootsplits};
}

void SubsplitDAG::CreateAndInsertNode(const Bitset &subsplit) {
  size_t node_id = NodeCount();
  dag_nodes_.push_back(std::make_unique<SubsplitDAGNode>(node_id, subsplit));
  SafeInsert(subsplit_to_id_, subsplit, node_id);
}

void SubsplitDAG::CreateAndInsertEdge(const size_t parent_id, const size_t child_id,
                                      bool rotated) {
  Assert(ContainsNode(parent_id), "Node with the given parent_id does not exist.");
  Assert(ContainsNode(child_id), "Node with the given child_id does not exist.");
  // Insert edge between parent and child.
  ConnectGivenNodes(parent_id, child_id, rotated);
  SafeInsert(dag_edges_, {parent_id, child_id}, GPCSPCountWithFakeSubsplits());
}

void SubsplitDAG::DeleteNode(const Bitset &node_subsplit) {
  Assert(ContainsNode(node_subsplit),
         "SubsplitDAG::DeleteNode(): Given edge does not exist");
  size_t node_id = GetDAGNodeId(node_subsplit);
  // dag_nodes_.erase(node_id);
}

void SubsplitDAG::DeleteEdge(const size_t parent_id, const size_t child_id,
                             bool rotated) {
  auto edge = std::make_pair(parent_id, child_id);
  Assert(ContainsEdge(parent_id, child_id),
         "SubsplitDAG::DeleteEdge(): Given edge does not exist.");
  // auto edge_idx = dag_edges_.at(edge);
  // dag_edges_.erase(edge);
}

void SubsplitDAG::ConnectGivenNodes(const size_t parent_id, const size_t child_id,
                                    bool rotated) {
  const auto parent_node = GetDAGNode(parent_id);
  const auto child_node = GetDAGNode(child_id);
  if (rotated) {
    parent_node->AddLeafwardRotated(child_node->Id());
    child_node->AddRootwardRotated(parent_node->Id());
  } else {
    parent_node->AddLeafwardSorted(child_node->Id());
    child_node->AddRootwardSorted(parent_node->Id());
  }
}

void SubsplitDAG::ConnectNodes(const SizeBitsetMap &index_to_child, size_t node_id,
                               bool rotated) {
  // Get bitset of parent node according to its rotation.
  const Bitset subsplit = GetDAGNode(node_id)->GetBitset(rotated);
  // Build vector of child node's subsplits.
  const auto children = GetChildSubsplits(index_to_child, subsplit, true);
  // Connect parent node to all child nodes.
  for (const auto &child_subsplit : children) {
    ConnectGivenNodes(node_id, GetDAGNodeId(child_subsplit), rotated);
  }
}

void SubsplitDAG::BuildNodesDepthFirst(const SizeBitsetMap &index_to_child,
                                       const Bitset &subsplit,
                                       std::unordered_set<Bitset> &visited_subsplits) {
  visited_subsplits.insert(subsplit);
  for (bool rotated : {false, true}) {
    for (const auto &child_subsplit : GetChildSubsplits(
             index_to_child, SubsplitToSortedOrder(subsplit, rotated), false)) {
      if (visited_subsplits.count(child_subsplit) == 0) {
        BuildNodesDepthFirst(index_to_child, child_subsplit, visited_subsplits);
      }
    }
  }
  CreateAndInsertNode(subsplit);
}

void SubsplitDAG::BuildNodes(const SizeBitsetMap &index_to_child,
                             const BitsetVector &rootsplits) {
  std::unordered_set<Bitset> visited_subsplits;
  // We will create fake subsplits and insert to dag_nodes_.
  // Fake subsplits (i.e. leaf nodes) will take IDs in [0, taxon_count_).
  for (size_t taxon_idx = 0; taxon_idx < taxon_count_; taxon_idx++) {
    CreateAndInsertNode(
        Bitset::FakeSubsplit(Bitset::Singleton(taxon_count_, taxon_idx)));
  }
  // Then we add the remaining nodes.
  // The root splits will take on the higher IDs compared to the non-rootsplits.
  for (const auto &rootsplit : rootsplits) {
    BuildNodesDepthFirst(index_to_child, rootsplit, visited_subsplits);
  }
  // Finally, we add the DAG root node.
  CreateAndInsertNode(Bitset::DAGRootSubsplitOfTaxonCount(taxon_count_));
}

void SubsplitDAG::BuildEdges(const SizeBitsetMap &index_to_child) {
  // Connect every node except for the DAG root node.
  for (size_t node_id = taxon_count_; node_id < DAGRootNodeId(); node_id++) {
    ConnectNodes(index_to_child, node_id, true);
    ConnectNodes(index_to_child, node_id, false);
  }
  // Connect the DAG root node.
  ConnectNodes(index_to_child, DAGRootNodeId(), true);
}

void SubsplitDAG::BuildDAGEdgesFromGPCSPIndexer(BitsetSizeMap &gpcsp_indexer) {
  for (const auto &[gpcsp, index] : gpcsp_indexer) {
    Assert(gpcsp.size() == 3 * taxon_count_,
           "All edges should be bitsets with size 3 times taxon_count_.");
    const auto parent_id = GetDAGNodeId(gpcsp.PCSPGetParentSubsplit());
    const auto child_id = GetDAGNodeId(gpcsp.PCSPGetChildSubsplit());
    SafeInsert(dag_edges_, {parent_id, child_id}, index);
  }
}

void SubsplitDAG::AddFakeSubsplitsToDAGEdgesAndParentToRange() {
  for (size_t node_id = 0; node_id < taxon_count_; node_id++) {
    const auto current_bitset(GetDAGNode(node_id)->GetBitset());
    IterateOverRootwardEdges(GetDAGNode(node_id), [this, current_bitset](
                                                      const bool rotated,
                                                      const SubsplitDAGNode *node) {
      SafeInsert(parent_to_child_range_, node->GetBitset(rotated),
                 {GPCSPCountWithFakeSubsplits(), GPCSPCountWithFakeSubsplits() + 1});
      SafeInsert(dag_edges_, {node->Id(), GetDAGNodeId(current_bitset)},
                 GPCSPCountWithFakeSubsplits());
    });
  }
}

void SubsplitDAG::RootwardDepthFirst(size_t node_id, SizeVector &visit_order,
                                     std::unordered_set<size_t> &visited_nodes) const {
  // Add to set of all visited nodes.
  SafeInsert(visited_nodes, node_id);
  // Recurse on sorted children.
  const auto &node = GetDAGNode(node_id);
  for (size_t parent_id : node->GetRootwardSorted()) {
    if (visited_nodes.count(parent_id) == 0) {
      RootwardDepthFirst(parent_id, visit_order, visited_nodes);
    }
  }
  // Recurse on rotated children.
  for (size_t parent_id : node->GetRootwardRotated()) {
    if (visited_nodes.count(parent_id) == 0) {
      RootwardDepthFirst(parent_id, visit_order, visited_nodes);
    }
  }
  // Append to vector post-order (after all children have been visited).
  visit_order.push_back(node_id);
}

void SubsplitDAG::LeafwardDepthFirst(size_t node_id, SizeVector &visit_order,
                                     std::unordered_set<size_t> &visited_nodes) const {
  // Add to set of all visited nodes.
  SafeInsert(visited_nodes, node_id);
  // Recurse on sorted children.
  for (size_t child_id : GetDAGNode(node_id)->GetLeafwardSorted()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, visit_order, visited_nodes);
    }
  }
  // Recurse on rotated children.
  for (size_t child_id : GetDAGNode(node_id)->GetLeafwardRotated()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, visit_order, visited_nodes);
    }
  }
  // Append to vector post-order (after all children have been visited).
  visit_order.push_back(node_id);
}

SizeVector SubsplitDAG::LeafwardPassTraversal(bool include_dag_root_node) const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  if (!include_dag_root_node) {
    SafeInsert(visited_nodes, DAGRootNodeId());
  }
  for (size_t leaf_id = 0; leaf_id < taxon_count_; leaf_id++) {
    RootwardDepthFirst(leaf_id, visit_order, visited_nodes);
  }
  return visit_order;
}

SizeVector SubsplitDAG::RootwardPassTraversal(bool include_dag_root_node) const {
  SizeVector visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (const auto &rootsplit_id : RootsplitIds()) {
    LeafwardDepthFirst(rootsplit_id, visit_order, visited_nodes);
  }
  if (include_dag_root_node) {
    visit_order.push_back(DAGRootNodeId());
  }
  return visit_order;
}

SizeVector SubsplitDAG::ReversePostorderTraversal() const {
  auto visit_order = RootwardPassTraversal(true);
  std::reverse(visit_order.begin(), visit_order.end());
  return visit_order;
}

void SubsplitDAG::ReversePostorderIndexTraversal(
    ParentRotationChildEdgeLambda f) const {
  for (const auto node_id : ReversePostorderTraversal()) {
    IterateOverLeafwardEdgesAndChildren(
        GetDAGNode(node_id), [&f, &node_id](const size_t gpcsp_idx, const bool rotated,
                                            const size_t child_id) {
          f(node_id, rotated, child_id, gpcsp_idx);
        });
  }
}

// ** Miscellaneous methods:

Bitset SubsplitDAG::SubsplitToSortedOrder(const Bitset &subsplit, bool rotated) const {
  return rotated ? subsplit.SubsplitRotate() : subsplit;
}

SizeVector SubsplitDAG::BuildTaxonTranslationMap(const SubsplitDAG &dag_a,
                                                 const SubsplitDAG &dag_b) {
  auto names_a = dag_a.GetSortedVectorOfTaxonNames();
  auto names_b = dag_b.GetSortedVectorOfTaxonNames();
  Assert(names_a == names_b,
         "SubsplitDAG::BuildTaxonTranslationMap(): SubsplitDAGs do not cover the same "
         "taxon set.");
  Assert(names_a.size() == dag_a.TaxonCount(),
         "SubsplitDAG::BuildTaxonTranslationMap(): Number of taxon names does not "
         "match the number of taxons in the DAG.");

  SizeVector taxon_map(names_a.size());
  for (const auto &name : names_a) {
    taxon_map[dag_a.GetTaxonId(name)] = dag_b.GetTaxonId(name);
  }
  return taxon_map;
}

int SubsplitDAG::BitsetCompareViaTaxonTranslationMap(const Bitset &bitset_a,
                                                     const Bitset &bitset_b,
                                                     const SizeVector &taxon_map) {
  Bitset bitset_a_translated_to_b =
      BitsetTranslateViaTaxonTranslationMap(bitset_a, taxon_map);
  return Bitset::Compare(bitset_a_translated_to_b, bitset_b);
}

Bitset SubsplitDAG::BitsetTranslateViaTaxonTranslationMap(
    const Bitset &bitset, const SizeVector &taxon_map, const bool forward_translate) {
  Bitset translated_bitset(bitset.size());
  size_t clade_size = taxon_map.size();
  size_t clade_count = bitset.size() / taxon_map.size();
  // Remap each clade individually.
  for (size_t i = 0; i < clade_count; i++) {
    // Remap each bit in clade according to a new position according to the taxon map.
    for (size_t j = 0; j < clade_size; j++) {
      size_t input_offset = (i * clade_size) + j;
      size_t output_offset = (i * clade_size) + taxon_map[j];
      if (forward_translate) {
        translated_bitset.set(input_offset, bitset[output_offset]);
      } else {
        translated_bitset.set(output_offset, bitset[input_offset]);
      }
    }
  }
  return translated_bitset;
}

bool SubsplitDAG::ContainsTaxon(const std::string &name) const {
  return dag_taxons_.find(name) != dag_taxons_.end();
}

bool SubsplitDAG::ContainsNode(const Bitset &subsplit) const {
  return subsplit_to_id_.find(subsplit) != subsplit_to_id_.end();
}

bool SubsplitDAG::ContainsNode(const size_t node_id) const {
  return node_id < NodeCount();
}

bool SubsplitDAG::ContainsEdge(const size_t parent_id, const size_t child_id) const {
  return dag_edges_.find({parent_id, child_id}) != dag_edges_.end();
}

std::pair<SizeVector, SizeVector> SubsplitDAG::BuildParentIdVector(
    const Bitset &subsplit) const {
  SizeVector rotated_parents, sorted_parents;
  // Linear search for all parents.
  for (const auto &[potential_parent_subsplit, node_id] : subsplit_to_id_) {
    if (subsplit.SubsplitIsRotatedChildOf(potential_parent_subsplit)) {
      rotated_parents.push_back(node_id);
    } else if (subsplit.SubsplitIsSortedChildOf(potential_parent_subsplit)) {
      sorted_parents.push_back(node_id);
    }
  }
  return {rotated_parents, sorted_parents};
}

std::pair<SizeVector, SizeVector> SubsplitDAG::BuildChildIdVector(
    const Bitset &subsplit) const {
  SizeVector rotated_children, sorted_children;
  // Linear search for all children.
  for (const auto &[potential_child_subsplit, node_id] : subsplit_to_id_) {
    if (potential_child_subsplit.SubsplitIsRotatedChildOf(subsplit)) {
      rotated_children.push_back(node_id);
    } else if (potential_child_subsplit.SubsplitIsSortedChildOf(subsplit)) {
      sorted_children.push_back(node_id);
    }
  }
  return {rotated_children, sorted_children};
}

void SubsplitDAG::ConnectChildToAllChildren(const Bitset &child_subsplit,
                                            SizeVector &added_edge_idxs) {
  // Get all children of new subsplit.
  const auto [rotated_children_of_child, sorted_children_of_child] =
      BuildChildIdVector(child_subsplit);

  // Iterate over sorted and rotated children of subsplit.
  for (const auto &[children_of_child, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{rotated_children_of_child, true},
                                                {sorted_children_of_child, false}}) {
    // Add child nodes in range to parent->child_range map.
    SafeInsert(parent_to_child_range_, SubsplitToSortedOrder(child_subsplit, rotated),
               {GPCSPCountWithFakeSubsplits(),
                GPCSPCountWithFakeSubsplits() + children_of_child.size()});
    // Add all child_ids
    for (const size_t child_of_child_id : children_of_child) {
      added_edge_idxs.push_back(GPCSPCountWithFakeSubsplits());
      CreateAndInsertEdge(GetDAGNodeId(child_subsplit), child_of_child_id, rotated);
    }
  }
}

void SubsplitDAG::ConnectParentToAllChildrenExcept(const Bitset &parent_subsplit,
                                                   const Bitset &child_subsplit,
                                                   SizeVector &added_edge_idxs) {
  const auto [rotated_children_of_parent, sorted_children_of_parent] =
      BuildChildIdVector(parent_subsplit);
  for (const auto &[children_of_parent, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{rotated_children_of_parent, true},
                                                {sorted_children_of_parent, false}}) {
    SafeInsert(parent_to_child_range_, SubsplitToSortedOrder(parent_subsplit, rotated),
               {GPCSPCountWithFakeSubsplits(),
                GPCSPCountWithFakeSubsplits() + children_of_parent.size()});
    for (const size_t child_of_parent_id : children_of_parent) {
      //
      if (child_of_parent_id != GetDAGNodeId(child_subsplit)) {
        added_edge_idxs.push_back(GPCSPCountWithFakeSubsplits());
        CreateAndInsertEdge(GetDAGNodeId(parent_subsplit), child_of_parent_id, rotated);
      }
    }
  }
}

void SubsplitDAG::ConnectChildToAllParentsExcept(const Bitset &parent_subsplit,
                                                 const Bitset &child_subsplit,
                                                 SizeVector &added_edge_idxs) {
  const auto [rotated_parents_of_child, sorted_parents_of_child] =
      BuildParentIdVector(child_subsplit);
  for (const auto &[parents_of_child, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{rotated_parents_of_child, true},
                                                {sorted_parents_of_child, false}}) {
    for (const size_t parent_of_child_id : parents_of_child) {
      // if parent_of_child is the same
      if (parent_of_child_id != GetDAGNodeId(parent_subsplit)) {
        added_edge_idxs.push_back(GPCSPCountWithFakeSubsplits());
        CreateAndInsertEdge(parent_of_child_id, GetDAGNodeId(child_subsplit), rotated);
      }
    }
  }
}

void SubsplitDAG::ConnectParentToAllParents(const Bitset &parent_subsplit,
                                            SizeVector &added_edge_idxs) {
  const auto [rotated_parents_of_parent, sorted_parents_of_parent] =
      BuildParentIdVector(parent_subsplit);
  for (const auto &[parents_of_parent, rotated] :
       std::vector<std::pair<SizeVector, bool>>{{rotated_parents_of_parent, true},
                                                {sorted_parents_of_parent, false}}) {
    for (const size_t parent_of_parent_id : parents_of_parent) {
      added_edge_idxs.push_back(GPCSPCountWithFakeSubsplits());
      CreateAndInsertEdge(parent_of_parent_id, GetDAGNodeId(parent_subsplit), rotated);
    }
  }
}

SubsplitDAG::ModificationResult SubsplitDAG::AddNodePair(const Bitset &parent_subsplit,
                                                         const Bitset &child_subsplit) {
  // Check that node pair will create a valid SubsplitDAG.
  Assert(
      IsValidAddNodePair(parent_subsplit, child_subsplit),
      "The given pair of nodes is incompatible with DAG in SubsplitDAG::AddNodePair.");
  // Initialize output vectors.
  SizeVector added_node_ids, added_edge_idxs, node_reindexer, edge_reindexer;
  // Check if either parent or child don't already exist in the DAG.
  const bool parent_is_new = !ContainsNode(parent_subsplit);
  const bool child_is_new = !ContainsNode(child_subsplit);
  // Soft assert: This allows for parent-child pair to exist in the DAG, but no work is
  // done. If both the parent and child already exist in DAG, return added_node_ids and
  // added_edge_idxs as empty, and node_reindexer and edge_reindexer as identity
  // reindexers.
  if (!parent_is_new && !child_is_new) {
    // Return default reindexers if both nodes already exist.
    node_reindexer = Reindexer::IdentityReindexer(NodeCount());
    edge_reindexer = Reindexer::IdentityReindexer(GPCSPCountWithFakeSubsplits());
    return {added_node_ids, added_edge_idxs, node_reindexer, edge_reindexer};
  }

  // Note: `prev_node_count` acts as a place marker. We know what the DAG root node id
  // is
  // (`prev_node_count - 1`).
  const size_t prev_node_count = NodeCount();

  // Add parent/child nodes and connect them to their children
  // If child node is new, add node and connect it to all its children.
  if (child_is_new) {
    CreateAndInsertNode(child_subsplit);
    added_node_ids.push_back(GetDAGNodeId(child_subsplit));
    // Don't reindex these edges.
    ConnectChildToAllChildren(child_subsplit, added_edge_idxs);
  }
  // If parent node is new, add node it to all its children (except )
  if (parent_is_new) {
    CreateAndInsertNode(parent_subsplit);
    added_node_ids.push_back(GetDAGNodeId(parent_subsplit));
    // Don't reindex these edges.
    ConnectParentToAllChildrenExcept(parent_subsplit, child_subsplit, added_edge_idxs);
  }
  // Note: `prev_edge_count` is a marker conveying where we need to start
  // reindexing edge idxs.
  // Edges are only reindexed if the parent node already existed in the DAG
  // (so as to ensure that edges descending from the same node clade have
  // contiguous idxs).
  size_t prev_edge_count = GPCSPCountWithFakeSubsplits();
  // Connect the given parent node to the given child node.
  added_edge_idxs.push_back(GPCSPCountWithFakeSubsplits());
  CreateAndInsertEdge(GetDAGNodeId(parent_subsplit), GetDAGNodeId(child_subsplit),
                      child_subsplit.SubsplitIsRotatedChildOf(parent_subsplit));
  // Don't reindex the edge between the given parent and child if the parent is new.
  if (parent_is_new) {
    prev_edge_count = GPCSPCountWithFakeSubsplits();
  }
  if (child_is_new) {
    // Reindex these edges.
    ConnectChildToAllParentsExcept(parent_subsplit, child_subsplit, added_edge_idxs);
  }
  if (parent_is_new) {
    // Reindex these edges.
    ConnectParentToAllParents(parent_subsplit, added_edge_idxs);
  }
  // Create reindexers.
  node_reindexer = BuildNodeReindexer(prev_node_count);
  edge_reindexer = BuildEdgeReindexer(prev_edge_count);
  // Update the ids in added_node_ids and added_edge_idxs according to the reindexers.
  Reindexer::RemapIdVector(added_node_ids, node_reindexer);
  Reindexer::RemapIdVector(added_edge_idxs, edge_reindexer);
  // Update fields in the Subsplit DAG according to the reindexers.
  RemapNodeIds(node_reindexer);
  RemapEdgeIdxs(edge_reindexer);
  // Recount topologies to update `topology_count_below_`.
  CountTopologies();

  return {added_node_ids, added_edge_idxs, node_reindexer, edge_reindexer};
}

SubsplitDAG::ModificationResult SubsplitDAG::RemoveNodePair(
    const Bitset &parent_subsplit, const Bitset &child_subsplit) {
  // Check that removal of node pair will maintain a valid DAG.
  Assert(IsValidRemoveNodePair(parent_subsplit, child_subsplit),
         "The given pair of nodes is incompatible with DAG in "
         "SubsplitDAG::RemoveNodePair.");
  // Output vectors.
  SizeVector removed_node_ids, removed_edge_idxs, node_reindexer, edge_reindexer;

  auto parent_id = GetDAGNodeId(parent_subsplit);
  auto child_id = GetDAGNodeId(child_subsplit);
  removed_node_ids = {parent_id, child_id};
  // Disconnect nodes: remove all edges connected to parent node and child node.
  for (const auto node_id : removed_node_ids) {
    const auto is_parent = (node_id == parent_id);
    const auto node = GetDAGNode(node_id);
    // Remove parents and children of given node.
    for (const bool &leafward : {true, false}) {
      // Remove rotated and sorted adjacent nodes of given node.
      for (const bool &rotated : {true, false}) {
        const auto adjacent_ids = node->GetLeafwardOrRootward(leafward, rotated);
        if (is_parent) {
          for (const auto adjacent_id : adjacent_ids) {
            DeleteEdge(node_id, adjacent_id, rotated);
            auto edge_idx = GetGPCSPEdgeIdx(node_id, adjacent_id);
            removed_edge_idxs.push_back(edge_idx);
          }
        } else {
          for (const auto adjacent_id : adjacent_ids) {
            DeleteEdge(adjacent_id, node_id, rotated);
            auto edge_idx = GetGPCSPEdgeIdx(adjacent_id, node_id);
            removed_edge_idxs.push_back(edge_idx);
          }
        }
      }
    }
  }
  // Remove pair of nodes.
  DeleteNode(parent_subsplit);
  DeleteNode(child_subsplit);
}

// ** Validation Test methods:

bool SubsplitDAG::IsConsistent() const {
  Failwith("SubsplitDAG::IsConsistent() is not yet implemented.");
  return false;
}

bool SubsplitDAG::IsValid() const {
  IterateOverRealNodes([this](const SubsplitDAGNode *node) {
    if (node->IsValid() == false) {
      return false;
    }
  });
  return true;
}

bool SubsplitDAG::IsValidAddNodePair(const Bitset &parent_subsplit,
                                     const Bitset &child_subsplit) const {
  auto GetParentNodeCounts = [this](Bitset subsplit) {
    const auto [rotated_parents, sorted_parents] = BuildParentIdVector(subsplit);
    SizePair parents = {rotated_parents.size(), sorted_parents.size()};
    return parents;
  };
  auto GetChildNodeCounts = [this](Bitset subsplit) {
    const auto [rotated_children, sorted_children] = BuildChildIdVector(subsplit);
    SizePair children = {rotated_children.size(), sorted_children.size()};
    return children;
  };

  // Get all adjacent nodes, not including parent and child.
  auto [rotated_parents_of_parent, sorted_parents_of_parent] =
      GetParentNodeCounts(parent_subsplit);
  auto [rotated_children_of_parent, sorted_children_of_parent] =
      GetChildNodeCounts(parent_subsplit);
  auto [rotated_children_of_child, sorted_children_of_child] =
      GetChildNodeCounts(child_subsplit);
  // Add child to parent's adjacent nodes.
  const bool is_child_rotated =
      child_subsplit.SubsplitIsRotatedChildOf(parent_subsplit);
  if (is_child_rotated) {
    rotated_children_of_parent++;
  } else {
    sorted_children_of_child++;
  }
  const bool is_parent_root = parent_subsplit.SubsplitIsRoot();
  const bool is_child_leaf = child_subsplit.SubsplitIsLeaf();

  // (1) Added nodes are parent/child pair.
  if (Bitset::SubsplitIsParentChildPair(parent_subsplit, child_subsplit) == false) {
    return false;
  }
  // (2) Nodes do not add or remove taxa.
  if ((parent_subsplit.size() != 2 * taxon_count_) ||
      (child_subsplit.size() != 2 * taxon_count_)) {
    std::cout << "BITSETS: " << parent_subsplit.SubsplitToString() << ", "
              << child_subsplit.SubsplitToString() << std::endl;
    std::cout << "NOT VALID => taxon_count" << std::endl;
    return false;
  }
  // (3) The parent node has at least one parent, and at least one rotated and sorted
  // child.
  // (* including the added child node)
  const bool parent_has_parent =
      (rotated_parents_of_parent > 0) || (sorted_parents_of_parent > 0);
  const bool parent_has_children =
      (rotated_children_of_parent > 0) && (sorted_children_of_parent > 0);
  if ((parent_has_parent && parent_has_children) == false) {
    std::cout << "BITSETS: " << parent_subsplit.SubsplitToString() << ", "
              << child_subsplit.SubsplitToString() << std::endl;
    std::cout << "NOT VALID => parent_node" << std::endl;
    std::cout << "PARENT_COUNTS: " << rotated_parents_of_parent << " "
              << sorted_parents_of_parent << " " << rotated_children_of_parent << " "
              << sorted_children_of_parent << std::endl;
    return false;
  }
  // (4) The child node has at least one parent, and at least one rotated and sorted
  // child.
  // (* we know child node has a parent node, so only need to check children)
  const bool child_has_children =
      (rotated_children_of_child > 0) && (sorted_children_of_child > 0);
  if (child_has_children == false) {
    std::cout << "BITSETS: " << parent_subsplit.SubsplitToString() << ", "
              << child_subsplit.SubsplitToString() << std::endl;
    std::cout << "NOT VALID => child_node" << std::endl;
    std::cout << "CHILD_COUNTS: " << rotated_children_of_child << " "
              << sorted_children_of_child << std::endl;
    return false;
  }
  return true;
}

bool SubsplitDAG::IsValidRemoveNodePair(const Bitset &parent_subsplit,
                                        const Bitset &child_subsplit) const {
  // (1) Check that parent and child are both in the DAG.
  if ((ContainsNode(parent_subsplit) && ContainsNode(child_subsplit)) == false) {
    return false;
  }
  // (2) Check that they are a valid parent/child pair.
  if (child_subsplit.SubsplitIsRotatedChildOf(parent_subsplit) ||
      child_subsplit.SubsplitIsSortedChildOf(parent_subsplit)) {
  }
  // (3) Check that parent is not the root node and child is not a leaf node (at least
  // for now, considering removing leaves illegal).
  if (parent_subsplit.SubsplitIsRoot() || child_subsplit.SubsplitIsLeaf()) {
    return false;
  }
  auto parent_id = GetDAGNodeId(parent_subsplit);
  auto child_id = GetDAGNodeId(child_subsplit);
  auto parent_node = GetDAGNode(parent_id);
  auto child_node = GetDAGNode(child_id);
  // (4) Check that all parents and children of removed nodes have alternate parents and
  // children (more than one).
  for (const auto &node : {parent_node, child_node}) {
    const bool is_parent = (node == parent_node);
    const size_t node_id = (is_parent ? parent_id : child_id);
    for (const auto &leafward : {true, false}) {
      for (const auto &rotated : {true, false}) {
        // All nodes adjacent to the nodes to be removed.
        const auto adjacent_ids = node->GetLeafwardOrRootward(leafward, rotated);
        for (const auto &adjacent_id : adjacent_ids) {
          // Check if adjacent nodes have alternate paths (greater than two edges).
          // (4a) If looking leafward, then the adjacent node is the child.
          // Therefore, we are looking at the adjacent node's rootward edges.
          // In the case of rootward edges, we look at the sum of sorted and rotated
          // edges. (4b) If looking rootward, then the adjacent node is the parent.
          // Therefore, we are looking at the adjacent node's leafward edges.
          // In that case of leafward edges, we look at the sorted and rotated edges
          // separately.
          const auto adjacent_node = GetDAGNode(adjacent_id);
          auto adjacent_node_count =
              adjacent_node->GetLeafwardOrRootward(!leafward, rotated).size();
          if (leafward) {
            adjacent_node_count +=
                adjacent_node->GetLeafwardOrRootward(!leafward, !rotated).size();
          }
          if (adjacent_node_count < 2) {
            // If there are no alternate edges, then fails, except in the case of the
            // edge connecting the parent node to child node.
            if (is_parent) {
              if (adjacent_id != child_id) {
                return false;
              }
            } else {
              if (adjacent_id != parent_id) {
                return false;
              }
            }
          }
        }
      }
    }
  }
  return true;
}

// ** Reindexer methods:

// NOTE: To be performed *after* DAG modification.
SizeVector SubsplitDAG::BuildNodeReindexer(const size_t prev_node_count) {
  SizeVector node_reindexer = Reindexer::IdentityReindexer(NodeCount());
  // Begin reindex values at taxon count to account for ...leaves?
  size_t running_traversal_idx = taxon_count_;
  size_t dag_root_node_id = prev_node_count - 1;
  // Build node_reindexer by using post-order traversal (topological sort) of entire DAG
  // to assign new ids, where the index is the "before" node_id (stored in the node
  // object), and the value is the "after" node_id.
  DepthFirstWithAction(
      {dag_root_node_id},
      SubsplitDAGTraversalAction(
          // BeforeNode
          [this](size_t node_id) {},
          // AfterNode
          [this, &node_reindexer, &running_traversal_idx](size_t node_id) {
            node_reindexer.at(node_id) = running_traversal_idx;
            running_traversal_idx++;
          },
          // BeforeNodeClade
          [this](size_t node_id, bool rotated) {},
          // VisitEdge
          [this](size_t node_id, size_t child_id, bool rotated) {}));
  return node_reindexer;
}

SizeVector SubsplitDAG::BuildEdgeReindexer(const size_t prev_edge_count) {
  SizeVector edge_reindexer =
      Reindexer::IdentityReindexer(GPCSPCountWithFakeSubsplits());
  // Only edges from an existing parent node to a new child node need to be reindexed.
  // See SubsplitDAG::AddNodePair().
  for (size_t edge_idx = prev_edge_count; edge_idx < GPCSPCountWithFakeSubsplits();
       edge_idx++) {
    // Find edge with given idx.
    const auto &element =
        std::find_if(dag_edges_.begin(), dag_edges_.end(),
                     [edge_idx](auto pair) { return pair.second == edge_idx; });
    Assert(element != dag_edges_.end(),
           "An edge with given edge_idx did not exist in "
           "SubsplitDAG::BuildEdgeReindexer.");
    //
    const auto &[node_pair, idx] = *element;
    const auto &[parent_id, child_id] = node_pair;
    const Bitset parent_subsplit = GetDAGNode(parent_id)->GetBitset();
    const Bitset child_subsplit = GetDAGNode(child_id)->GetBitset();
    const auto idx_range = GetEdgeRange(
        parent_subsplit, child_subsplit.SubsplitIsRotatedChildOf(parent_subsplit));
    // New edge is added to the end of the range.
    const size_t new_idx = edge_reindexer.at(idx_range.second);
    Reindexer::ReassignAndShift(edge_reindexer, edge_idx, new_idx);
  }
  return edge_reindexer;
}

void SubsplitDAG::RemapNodeIds(const SizeVector &node_reindexer) {
  // Update `dag_nodes_`.
  auto dag_nodes_copy = Reindexer::Reindex(dag_nodes_, node_reindexer);
  dag_nodes_.swap(dag_nodes_copy);
  // Update each node's id and leafward/rootward ids.
  for (int node_id = 0; node_id < NodeCount(); node_id++) {
    GetDAGNode(node_id)->RemapNodeIds(node_reindexer);
  }
  // Update `subsplit_to_id_`.
  for (const auto &[subsplit, node_id] : subsplit_to_id_) {
    subsplit_to_id_.at(subsplit) = node_reindexer.at(node_id);
  }
  // Update `dag_edges_`.
  std::map<SizePair, size_t> dag_edges_copy;
  for (const auto &[node_pair, edge_idx] : dag_edges_) {
    SafeInsert(
        dag_edges_copy,
        {node_reindexer.at(node_pair.first), node_reindexer.at(node_pair.second)},
        edge_idx);
  }
  dag_edges_.swap(dag_edges_copy);
}

void SubsplitDAG::RemapEdgeIdxs(const SizeVector &edge_reindexer) {
  // Update `dag_edges_`.
  for (const auto &[node_pair, edge_idx] : dag_edges_) {
    dag_edges_.at(node_pair) = edge_reindexer.at(edge_idx);
  }
  // Update `parent_to_child_range_`.
  for (const auto &[subsplit, idx_range] : parent_to_child_range_) {
    parent_to_child_range_.at(subsplit) = {edge_reindexer.at(idx_range.first),
                                           edge_reindexer.at(idx_range.second)};
  }
}
