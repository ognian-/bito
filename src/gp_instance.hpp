// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// This is the comprehensive unit for performing Generalized Pruning.
//

#ifndef SRC_GP_INSTANCE_HPP_
#define SRC_GP_INSTANCE_HPP_

#include "gp_dag.hpp"
#include "gp_engine.hpp"
#include "nni_evaluation_engine.hpp"
#include "rooted_tree_collection.hpp"
#include "site_pattern.hpp"

class GPInstance {
 public:
  // Create instance with filepath for virtual memory map.
  explicit GPInstance(const std::string &mmap_file_path)
      : mmap_file_path_(mmap_file_path) {
    if (mmap_file_path.empty()) {
      Failwith("GPInstance needs a legal path as a constructor argument.");
    }
  };

  // Print following GPInstance stats to console:
  //  - Number of trees, Number of taxa/leaves.
  //  - Number of alignment sequences.
  //  - Number of SubsplitDAG nodes, Number of possible tree topologies in DAG.
  //  - Number of continuous parameters (this is the same as branch lengths/edges).
  //  - Amount of virtual memory used by Engine (or whether engine has been made).
  void PrintStatus();

  // Load fasta file and overwrites/adds all sequences to alignment_.
  void ReadFastaFile(const std::string &fname);
  // Load newick file and overwrites/adds all trees to tree_collection_.
  void ReadNewickFile(const std::string &fname);
  // Load compressed newick file and overwrites/adds all trees to tree_collection_.
  void ReadNewickFileGZ(const std::string &fname);
  // Load nexus file and overwrites/adds all trees to tree_collection_.
  void ReadNexusFile(const std::string &fname);
  // Load compressed nexus file and overwrites/adds all trees to tree_collection_.
  void ReadNexusFileGZ(const std::string &fname);

  // Initialize DAG and GPEngine. Assumes uniform distribution on topological support.
  // NOTE: Requires that sequences and trees have already been loaded.
  void MakeEngine(double rescaling_threshold = GPEngine::default_rescaling_threshold_);
  GPEngine *GetEngine() const;
  bool HasEngine() const;

  GPDAG &GetDAG();
  void PrintDAG();
  void PrintEdgeIndexer();
  void ProcessOperations(const GPOperationVector &operations);
  void HotStartBranchLengths();
  void EstimateSBNParameters();
  void EstimateBranchLengths(double tol, size_t max_iter, bool quiet = false);
  void PopulatePLVs();
  void ComputeLikelihoods();
  void ComputeMarginalLikelihood();
  void CalculateHybridMarginals();
  RootedTreeCollection GenerateCompleteRootedTreeCollection();

  // #348: A lot of code duplication here with things in SBNInstance.
  StringVector PrettyIndexer() const;
  EigenConstVectorXdRef GetSBNParameters();
  StringDoubleVector PrettyIndexedSBNParameters();
  StringDoubleVector PrettyIndexedBranchLengths();
  StringDoubleVector PrettyIndexedPerGPCSPLogLikelihoods();
  StringDoubleVector PrettyIndexedPerGPCSPComponentsOfFullLogMarginal();

  void SBNParametersToCSV(const std::string &file_path);
  void SBNPriorToCSV(const std::string &file_path);
  void BranchLengthsToCSV(const std::string &file_path);

  // Generate a version of the topologies in the current tree collection that use the
  // current GP branch lengths.
  RootedTreeCollection CurrentlyLoadedTreesWithGPBranchLengths();

  // Subset the currently loaded topologies to those that have a given PCSP, and equip
  // them with current GP branch lengths.
  RootedTreeCollection CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths(
      const std::string &pcsp_string);

  // Run CurrentlyLoadedTreesWithGPBranchLengths and export to a Newick file.
  void ExportTrees(const std::string &out_path);
  // Run CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths and export to a Newick
  // file.
  void ExportTreesWithAPCSP(const std::string &pcsp_string,
                            const std::string &newick_path);
  // Run CurrentlyLoadedTreesWithGPBranchLengths and export to a Newick file.
  // Export all trees in the span of the subsplit DAG (with GP branch lengths) to a
  // Newick file.
  void ExportAllGeneratedTrees(const std::string &out_path);
  // Generate all trees spanned by the DAG and load them into the instance.
  void LoadAllGeneratedTrees();

  // Export the subsplit DAG as a DOT file.
  void SubsplitDAGToDot(const std::string &out_path, bool show_index_labels = true);

  // TODO: Work in Progress
  // ** NNI Evaluation Engine

  // Initialize NNI Evaluation Engine.
  void MakeNNIEngine();
  // Get NNI Evaluation Engine.
  NNIEvaluationEngine &GetNNIEngine();
  // Get the number of adjacent NNI's for current DAG.
  size_t GetNNICount();

  // Add NNI pair.
  SubsplitDAG::ModificationResult AddNodePair(const Bitset &parent_bitset,
                                              const Bitset &child_bitset);
  SubsplitDAG::ModificationResult AddGraftNodePair(const Bitset &parent_bitset,
                                                   const Bitset &child_bitset);
  // Update engine.
  void UpdateEngineAfterModifyingDAG(const SizeVector &node_reindexer,
                                     const SizeVector &edge_reindexer);
  void UpdateEngineAfterGraftingDAG(const SizeVector &node_reindexer,
                                    const SizeVector &edge_reindexer);
  // Compute Marginal Likelihood for NNI.
  void ComputePerNNIPerPCSPLikelihood(const Bitset &parent_bitset,
                                      const Bitset &child_bitset);
  void ComputePerNNIPerPCSPLikelihood(const SubsplitDAGNode &parent_node,
                                      const SubsplitDAGNode &child_node);

  // Initialize graft dag.
  void MakeGraftDAG();

  // ** Alignment

  // Get taxon names.
  StringVector GetTaxonNames() { return tree_collection_.TaxonNames(); }

 private:
  void ClearTreeCollectionAssociatedState();
  // Verify that sequences and trees are nonempty.
  void CheckSequencesAndTreesLoaded() const;

  size_t GetEdgeIndexForLeafNode(const Bitset &parent_subsplit,
                                 const Node *leaf_node) const;
  RootedTreeCollection TreesWithGPBranchLengthsOfTopologies(
      Node::NodePtrVec &&topologies) const;
  StringDoubleVector PrettyIndexedVector(EigenConstVectorXdRef v);

 private:
  std::string mmap_file_path_;
  Alignment alignment_;
  std::unique_ptr<GPEngine> engine_;
  RootedTreeCollection tree_collection_;
  GPDAG dag_;
  static constexpr size_t plv_count_per_node_ = 6;

  // TODO: NNI Evaluation Engine
  std::optional<NNIEvaluationEngine> nni_engine_ = std::nullopt;
  std::optional<SubsplitDAGGraft> graft_dag_ = std::nullopt;
  std::string graft_mmap_file_path_;
};

#endif  // SRC_GP_INSTANCE_HPP_
