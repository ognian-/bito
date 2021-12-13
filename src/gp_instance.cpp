// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "gp_instance.hpp"

#include <chrono>
#include <iomanip>
#include <string>

#include "csv.hpp"
#include "driver.hpp"
#include "gp_operation.hpp"
#include "numerical_utils.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_probability.hpp"

using namespace GPOperations;  // NOLINT

void GPInstance::PrintStatus() {
  const auto tree_count = tree_collection_.TreeCount();
  const auto taxon_count = tree_collection_.TaxonCount();
  if (tree_count > 0) {
    std::cout << tree_count << " trees loaded on " << taxon_count << " leaves.\n";
  } else {
    std::cout << "No trees loaded.\n";
  }
  std::cout << alignment_.Data().size() << " sequences loaded.\n";
  std::cout << dag_.NodeCount() << " DAG nodes representing " << dag_.TopologyCount()
            << " trees.\n";
  std::cout << dag_.EdgeCountWithLeafSubsplits() << " continuous parameters.\n";
  if (HasEngine()) {
    std::cout << "Engine available using " << GetEngine()->PLVByteCount() / 1e9
              << "G virtual memory.\n";
  } else {
    std::cout << "Engine has not been made.\n";
  }
}

void GPInstance::ReadFastaFile(const std::string &fname) {
  alignment_ = Alignment::ReadFasta(fname);
}

void GPInstance::ReadNewickFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile(fname));
}

void GPInstance::ReadNewickFileGZ(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFileGZ(fname));
}

void GPInstance::ReadNexusFile(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFile(fname));
}

void GPInstance::ReadNexusFileGZ(const std::string &fname) {
  Driver driver;
  tree_collection_ =
      RootedTreeCollection::OfTreeCollection(driver.ParseNexusFileGZ(fname));
}

void GPInstance::CheckSequencesAndTreesLoaded() const {
  if (alignment_.SequenceCount() == 0) {
    Failwith(
        "Load an alignment into your GPInstance on which you wish to "
        "calculate phylogenetic likelihoods.");
  }
  if (tree_collection_.TreeCount() == 0) {
    Failwith(
        "Load some trees into your GPInstance on which you wish to "
        "calculate phylogenetic likelihoods.");
  }
}

void GPInstance::MakeEngine(double rescaling_threshold) {
  CheckSequencesAndTreesLoaded();
  // TODO: Explain site patterns.
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());
  // Initialize the DAG (uses TidySubsplitDAG constructor).
  dag_ = GPDAG(tree_collection_);
  // Initialize the GPEngine.
  // Uniform support assigns each possible topology expressed by the DAG to have equal
  // likelihood.
  auto sbn_prior = dag_.BuildUniformOnTopologicalSupportPrior();
  // TODO:
  auto unconditional_node_probabilities =
      dag_.UnconditionalNodeProbabilities(sbn_prior);
  // TODO:
  auto inverted_sbn_prior =
      dag_.InvertedGPCSPProbabilities(sbn_prior, unconditional_node_probabilities);
  //
  engine_ = std::make_unique<GPEngine>(
      std::move(site_pattern), plv_count_per_node_ * (dag_.NodeCountWithoutDAGRoot()),
      dag_.EdgeCountWithLeafSubsplits(), mmap_file_path_, rescaling_threshold,
      std::move(sbn_prior), std::move(unconditional_node_probabilities),
      std::move(inverted_sbn_prior));
}

GPEngine *GPInstance::GetEngine() const {
  if (engine_ != nullptr) {
    return engine_.get();
  }
  // else
  Failwith(
      "Engine not available. Call MakeEngine to make an engine for phylogenetic "
      "likelihood computation.");
}

bool GPInstance::HasEngine() const { return engine_ != nullptr; }

GPDAG &GPInstance::GetDAG() { return dag_; }

void GPInstance::PrintDAG() { dag_.Print(); }

void GPInstance::PrintEdgeIndexer() {
  std::cout << "Vector of taxon names: " << tree_collection_.TaxonNames() << std::endl;
  dag_.PrintEdgeIndexer();
}

void GPInstance::ProcessOperations(const GPOperationVector &operations) {
  GetEngine()->ProcessOperations(operations);
}

void GPInstance::ClearTreeCollectionAssociatedState() { dag_ = GPDAG(); }

void GPInstance::HotStartBranchLengths() {
  if (HasEngine()) {
    GetEngine()->HotStartBranchLengths(tree_collection_, dag_.BuildEdgeIndexer());
  } else {
    Failwith(
        "Please load and process some trees before calling HotStartBranchLengths.");
  }
}

void GPInstance::PopulatePLVs() { ProcessOperations(dag_.PopulatePLVs()); }

void GPInstance::ComputeLikelihoods() { ProcessOperations(dag_.ComputeLikelihoods()); }

void GPInstance::ComputeMarginalLikelihood() {
  ProcessOperations(dag_.MarginalLikelihood());
}

void GPInstance::EstimateBranchLengths(double tol, size_t max_iter, bool quiet) {
  std::stringstream dev_null;
  auto &our_ostream = quiet ? dev_null : std::cout;
  auto now = std::chrono::high_resolution_clock::now;
  auto t_start = now();
  our_ostream << "Begin branch optimization\n";
  GPOperationVector branch_optimization_operations = dag_.BranchLengthOptimization();
  GPOperationVector marginal_lik_operations = dag_.MarginalLikelihood();
  GPOperationVector populate_plv_operations = dag_.PopulatePLVs();

  our_ostream << "Populating PLVs\n";
  PopulatePLVs();
  std::chrono::duration<double> warmup_duration = now() - t_start;
  t_start = now();
  our_ostream << "Computing initial likelihood\n";
  ProcessOperations(marginal_lik_operations);
  double current_marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();
  std::chrono::duration<double> initial_likelihood_duration = now() - t_start;
  t_start = now();

  for (size_t i = 0; i < max_iter; i++) {
    our_ostream << "Iteration: " << (i + 1) << std::endl;
    ProcessOperations(branch_optimization_operations);
    // #321 Replace with a cleaned up traversal.
    ProcessOperations(populate_plv_operations);
    ProcessOperations(marginal_lik_operations);
    double marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();
    our_ostream << "Current marginal log likelihood: ";
    our_ostream << std::setprecision(9) << current_marginal_log_lik << std::endl;
    our_ostream << "New marginal log likelihood: ";
    our_ostream << std::setprecision(9) << marginal_log_lik << std::endl;
    if (marginal_log_lik < current_marginal_log_lik) {
      our_ostream << "Marginal log likelihood decreased.\n";
    }
    if (abs(current_marginal_log_lik - marginal_log_lik) < tol) {
      our_ostream << "Converged.\n";
      break;
    }
    current_marginal_log_lik = marginal_log_lik;
  }
  std::chrono::duration<double> optimization_duration = now() - t_start;
  our_ostream << "\n# Timing Report\n";
  our_ostream << "warmup: " << warmup_duration.count() << "s\n";
  our_ostream << "initial likelihood: " << initial_likelihood_duration.count() << "s\n";
  our_ostream << "optimization: " << optimization_duration.count() << "s or "
              << optimization_duration.count() / 60 << "m\n";
}

void GPInstance::EstimateSBNParameters() {
  std::cout << "Begin SBN parameter optimization\n";
  PopulatePLVs();
  ComputeLikelihoods();
  ProcessOperations(dag_.OptimizeSBNParameters());
}

void GPInstance::CalculateHybridMarginals() {
  std::cout << "Calculating hybrid marginals\n";
  PopulatePLVs();
  dag_.ReversePostorderIndexTraversal([this](const size_t parent_id, const bool rotated,
                                             const size_t child_id, const size_t) {
    this->GetEngine()->ProcessQuartetHybridRequest(
        dag_.QuartetHybridRequestOf(parent_id, rotated, child_id));
  });
}

size_t GPInstance::GetEdgeIndexForLeafNode(const Bitset &parent_subsplit,
                                           const Node *leaf_node) const {
  Assert(leaf_node->IsLeaf(), "Only leaf node is permitted.");
  return dag_.GetEdgeIdx(parent_subsplit, Bitset::LeafSubsplit(leaf_node->Leaves()));
}

RootedTreeCollection GPInstance::TreesWithGPBranchLengthsOfTopologies(
    Node::NodePtrVec &&topologies) const {
  const EigenVectorXd edge_indexed_branch_lengths = engine_->GetBranchLengths();
  RootedTree::RootedTreeVector tree_vector;

  for (auto &root_node : topologies) {
    size_t node_count = 2 * root_node->LeafCount() - 1;
    std::vector<double> branch_lengths(node_count);

    root_node->RootedPCSPPreorder([this, &branch_lengths, &edge_indexed_branch_lengths](
                                      const Node *sister, const Node *focal,
                                      const Node *child0, const Node *child1) {
      Bitset parent_subsplit = Bitset::Subsplit(sister->Leaves(), focal->Leaves());
      Bitset child_subsplit = Bitset::Subsplit(child0->Leaves(), child1->Leaves());
      size_t edge_idx = dag_.GetEdgeIdx(parent_subsplit, child_subsplit);
      branch_lengths[focal->Id()] = edge_indexed_branch_lengths[edge_idx];

      if (sister->IsLeaf()) {
        edge_idx = GetEdgeIndexForLeafNode(parent_subsplit, sister);
        branch_lengths[sister->Id()] = edge_indexed_branch_lengths[edge_idx];
      }
      if (child0->IsLeaf()) {
        edge_idx = GetEdgeIndexForLeafNode(child_subsplit, child0);
        branch_lengths[child0->Id()] = edge_indexed_branch_lengths[edge_idx];
      }
      if (child1->IsLeaf()) {
        edge_idx = GetEdgeIndexForLeafNode(child_subsplit, child1);
        branch_lengths[child1->Id()] = edge_indexed_branch_lengths[edge_idx];
      }
    });

    tree_vector.emplace_back(root_node, std::move(branch_lengths));
  }

  return RootedTreeCollection(tree_vector, tree_collection_.TagTaxonMap());
}

RootedTreeCollection GPInstance::GenerateCompleteRootedTreeCollection() {
  return TreesWithGPBranchLengthsOfTopologies(dag_.GenerateAllTopologies());
}

StringVector GPInstance::PrettyIndexer() const {
  StringVector pretty_representation(dag_.BuildEdgeIndexer().size());
  for (const auto &[pcsp, idx] : dag_.BuildEdgeIndexer()) {
    pretty_representation[idx] = pcsp.EdgeToString();
  }
  return pretty_representation;
}

StringDoubleVector GPInstance::PrettyIndexedVector(EigenConstVectorXdRef v) {
  StringDoubleVector result;
  result.reserve(v.size());
  const auto pretty_indexer = PrettyIndexer();
  Assert(v.size() <= pretty_indexer.size(), "v is too long in PrettyIndexedVector");
  for (size_t i = 0; i < v.size(); i++) {
    result.push_back({pretty_indexer.at(i), v(i)});
  }
  return result;
}

EigenConstVectorXdRef GPInstance::GetSBNParameters() {
  return engine_->GetSBNParameters();
}

StringDoubleVector GPInstance::PrettyIndexedSBNParameters() {
  return PrettyIndexedVector(GetSBNParameters());
}

StringDoubleVector GPInstance::PrettyIndexedBranchLengths() {
  return PrettyIndexedVector(GetEngine()->GetBranchLengths());
}

StringDoubleVector GPInstance::PrettyIndexedPerGPCSPLogLikelihoods() {
  return PrettyIndexedVector(GetEngine()->GetPerGPCSPLogLikelihoods());
}

StringDoubleVector GPInstance::PrettyIndexedPerGPCSPComponentsOfFullLogMarginal() {
  return PrettyIndexedVector(GetEngine()->GetPerGPCSPComponentsOfFullLogMarginal());
}

void GPInstance::SBNParametersToCSV(const std::string &file_path) {
  CSV::StringDoubleVectorToCSV(PrettyIndexedSBNParameters(), file_path);
}

void GPInstance::SBNPriorToCSV(const std::string &file_path) {
  CSV::StringDoubleVectorToCSV(
      PrettyIndexedVector(dag_.BuildUniformOnTopologicalSupportPrior()), file_path);
}

void GPInstance::BranchLengthsToCSV(const std::string &file_path) {
  CSV::StringDoubleVectorToCSV(PrettyIndexedBranchLengths(), file_path);
}

RootedTreeCollection GPInstance::CurrentlyLoadedTreesWithGPBranchLengths() {
  Node::NodePtrVec topologies;
  for (const auto &tree : tree_collection_.Trees()) {
    topologies.push_back(tree.Topology()->DeepCopy());
  }
  return TreesWithGPBranchLengthsOfTopologies(std::move(topologies));
}

RootedTreeCollection GPInstance::CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths(
    const std::string &pcsp_string) {
  const BitsetSizeMap &indexer = dag_.BuildEdgeIndexer();
  Bitset pcsp(pcsp_string);
  auto search = indexer.find(pcsp);
  if (search == indexer.end()) {
    Failwith("Don't have " + pcsp_string + " as a PCSP in the instance!");
  }
  auto pcsp_index = search->second;

  Node::NodePtrVec topologies;
  for (const auto &tree : tree_collection_.Trees()) {
    auto indexer_representation = dag_.IndexerRepresentationOf(
        indexer, tree.Topology(), std::numeric_limits<size_t>::max());
    if (std::find(indexer_representation.begin(), indexer_representation.end(),
                  pcsp_index) != indexer_representation.end()) {
      topologies.push_back(tree.Topology()->DeepCopy());
    }
  }
  return TreesWithGPBranchLengthsOfTopologies(std::move(topologies));
}

void GPInstance::ExportTrees(const std::string &out_path) {
  auto trees = CurrentlyLoadedTreesWithGPBranchLengths();
  trees.ToNewickFile(out_path);
}

void GPInstance::ExportTreesWithAPCSP(const std::string &pcsp_string,
                                      const std::string &out_path) {
  auto trees = CurrentlyLoadedTreesWithAPCSPStringAndGPBranchLengths(pcsp_string);
  trees.ToNewickFile(out_path);
}

void GPInstance::ExportAllGeneratedTrees(const std::string &out_path) {
  auto trees = GenerateCompleteRootedTreeCollection();
  trees.ToNewickFile(out_path);
}

void GPInstance::LoadAllGeneratedTrees() {
  tree_collection_ = GenerateCompleteRootedTreeCollection();
}

void GPInstance::SubsplitDAGToDot(const std::string &out_path, bool show_index_labels) {
  std::ofstream out_stream(out_path);
  out_stream << dag_.ToDot(show_index_labels) << std::endl;
  if (out_stream.bad()) {
    Failwith("Failure writing to " + out_path);
  }
  out_stream.close();
}

void GPInstance::AddNodePair(const Bitset &parent_bitset, const Bitset &child_bitset) {
  // Add node pair to SubsplitDAG.
  auto node_addition_result = dag_.AddNodePair(parent_bitset, child_bitset);
  // Fetch the reorganized locations of nodes and edges (relative to before AddNodePair)
  auto node_reindexer = node_addition_result.node_reindexer;
  auto edge_reindexer = node_addition_result.edge_reindexer;
  // Reindex data
  // Reindexer::Reindex(, node_reindexer);

  // If NNI Evaluation Engine has been initialized, then update NNIs.
  if (nni_engine_) {
    nni_engine_->UpdateSetOfNNIsAfterDAGAddNodePair(parent_bitset, child_bitset);
  }
}

void GPInstance::MakeNNIEngine() {
  nni_engine_ = std::make_unique<NNIEvaluationEngine>(dag_);
  nni_engine_->SyncSetOfNNIsWithDAG();
}

size_t GPInstance::GetNNICount() { return nni_engine_->GetAdjacentNNICount(); }
