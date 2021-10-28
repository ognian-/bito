// Copyright 2019-2021 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

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
  std::cout << dag_.GPCSPCountWithFakeSubsplits() << " continuous parameters.\n";
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
  SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());

  dag_ = GPDAG(tree_collection_);
  auto sbn_prior = dag_.BuildUniformOnTopologicalSupportPrior();
  auto unconditional_node_probabilities =
      dag_.UnconditionalNodeProbabilities(sbn_prior);
  auto inverted_sbn_prior =
      dag_.InvertedGPCSPProbabilities(sbn_prior, unconditional_node_probabilities);
  engine_ = std::make_unique<GPEngine>(
      std::move(site_pattern), plv_count_per_node_ * (dag_.NodeCountWithoutDAGRoot()),
      dag_.GPCSPCountWithFakeSubsplits(), mmap_file_path_, rescaling_threshold,
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

void GPInstance::PrintGPCSPIndexer() {
  std::cout << "Vector of taxon names: " << tree_collection_.TaxonNames() << std::endl;
  dag_.PrintGPCSPIndexer();
}

void GPInstance::ProcessOperations(const GPOperationVector &operations) {
  GetEngine()->ProcessOperations(operations);
}

void GPInstance::ClearTreeCollectionAssociatedState() { dag_ = GPDAG(); }

void GPInstance::HotStartBranchLengths() {
  if (HasEngine()) {
    GetEngine()->HotStartBranchLengths(tree_collection_, dag_.BuildGPCSPIndexer());
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
  GPOperationVector compute_lik_operations = dag_.ComputeLikelihoods();

  our_ostream << "Populating PLVs\n";
  PopulatePLVs();
  std::chrono::duration<double> warmup_duration = now() - t_start;
  t_start = now();
  our_ostream << "Computing initial likelihood\n";
  ProcessOperations(marginal_lik_operations);
  ProcessOperations(compute_lik_operations);

  size_t gpcsp_count = dag_.GPCSPCountWithFakeSubsplits();
  per_pcsp_branch_lengths_ = EigenMatrixXd(gpcsp_count, 1);
  per_pcsp_marg_lik_ = EigenMatrixXd(gpcsp_count, 1);

  per_pcsp_branch_lengths_.col(0) = GetEngine()->GetBranchLengths();
  per_pcsp_marg_lik_.col(0) = GetEngine()->GetPerGPCSPLogLikelihoods();

  double current_marginal_log_lik = GetEngine()->GetLogMarginalLikelihood();
  std::chrono::duration<double> initial_likelihood_duration = now() - t_start;

  t_start = now();

  for (size_t i = 0; i < max_iter; i++) {
    our_ostream << "Iteration: " << (i + 1) << std::endl;
    ProcessOperations(branch_optimization_operations);
    // #321 Replace with a cleaned up traversal.
    ProcessOperations(populate_plv_operations);
    ProcessOperations(marginal_lik_operations);
    ProcessOperations(compute_lik_operations);

    per_pcsp_branch_lengths_.conservativeResize(Eigen::NoChange, i + 2);
    per_pcsp_marg_lik_.conservativeResize(Eigen::NoChange, i + 2);

    per_pcsp_branch_lengths_.col(i + 1) = GetEngine()->GetBranchLengths();
    per_pcsp_marg_lik_.col(i + 1) = GetEngine()->GetPerGPCSPLogLikelihoods();

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

size_t GPInstance::GetGPCSPIndexForLeafNode(const Bitset &parent_subsplit,
                                            const Node *leaf_node) const {
  Assert(leaf_node->IsLeaf(), "Only leaf node is permitted.");
  return dag_.GetGPCSPIndex(parent_subsplit, Bitset::FakeSubsplit(leaf_node->Leaves()));
}

RootedTreeCollection GPInstance::TreesWithGPBranchLengthsOfTopologies(
    Node::NodePtrVec &&topologies) const {
  const EigenVectorXd gpcsp_indexed_branch_lengths = engine_->GetBranchLengths();
  RootedTree::RootedTreeVector tree_vector;

  for (auto &root_node : topologies) {
    size_t node_count = 2 * root_node->LeafCount() - 1;
    std::vector<double> branch_lengths(node_count);

    root_node->RootedPCSPPreorder(
        [this, &branch_lengths, &gpcsp_indexed_branch_lengths](
            const Node *sister, const Node *focal, const Node *child0,
            const Node *child1) {
          Bitset parent_subsplit =
              Bitset::SubsplitOfPair(sister->Leaves(), focal->Leaves());
          Bitset child_subsplit =
              Bitset::SubsplitOfPair(child0->Leaves(), child1->Leaves());
          size_t gpcsp_idx = dag_.GetGPCSPIndex(parent_subsplit, child_subsplit);
          branch_lengths[focal->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];

          if (sister->IsLeaf()) {
            gpcsp_idx = GetGPCSPIndexForLeafNode(parent_subsplit, sister);
            branch_lengths[sister->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];
          }
          if (child0->IsLeaf()) {
            gpcsp_idx = GetGPCSPIndexForLeafNode(child_subsplit, child0);
            branch_lengths[child0->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];
          }
          if (child1->IsLeaf()) {
            gpcsp_idx = GetGPCSPIndexForLeafNode(child_subsplit, child1);
            branch_lengths[child1->Id()] = gpcsp_indexed_branch_lengths[gpcsp_idx];
          }
        });

    tree_vector.emplace_back(root_node, std::move(branch_lengths));
  }

  return RootedTreeCollection(tree_vector, tree_collection_.TagTaxonMap());
}

RootedTreeCollection GPInstance::GenerateCompleteRootedTreeCollection() {
  return TreesWithGPBranchLengthsOfTopologies(dag_.GenerateAllTopologies());
}

void GPInstance::GetPerGPCSPLogLikelihoodSurfaces(int steps) {
  const EigenVectorXd optimized_branch_lengths = GetEngine()->GetBranchLengths();

  size_t gpcsp_count = optimized_branch_lengths.size();
  const EigenVectorXd scaling_vector = EigenVectorXd::LinSpaced(steps, 0.01, 10.0);
  per_pcsp_lik_surfaces_ = EigenMatrixXd(gpcsp_count * steps, 2);

  for (int gpcsp_idx = 0; gpcsp_idx < gpcsp_count; gpcsp_idx++) {
    EigenVectorXd gpcsp_new_branch_lengths =
        scaling_vector * optimized_branch_lengths[gpcsp_idx];
    EigenVectorXd new_branch_length_vector = optimized_branch_lengths;

    for (int i = 0; i < steps; i++) {
      new_branch_length_vector[gpcsp_idx] = gpcsp_new_branch_lengths[i];
      GetEngine()->SetBranchLengths(new_branch_length_vector);
      PopulatePLVs();
      ComputeLikelihoods();

      int matrix_position = gpcsp_count * i;
      per_pcsp_lik_surfaces_(matrix_position + gpcsp_idx, 0) =
          gpcsp_new_branch_lengths[i];
      per_pcsp_lik_surfaces_(matrix_position + gpcsp_idx, 1) =
          GetEngine()->GetPerGPCSPLogLikelihoods(gpcsp_idx, 1)(0, 0);
    }
  }
}

StringEigenVectorXdVector GPInstance::TrackValuesFromOptimization() {
  // We assume this is run after a proper branch length optimization
  const EigenVectorXd optimized_branch_lengths = GetEngine()->GetBranchLengths();
  const EigenVectorXd optimized_per_pcsp_llhs =
      GetEngine()->GetPerGPCSPLogLikelihoods();

  size_t gpcsp_count = optimized_branch_lengths.size();
  EigenVectorXi run_counts = EigenVectorXi::Zero(gpcsp_count);
  EigenMatrixXd tracked_optimization_values(1, 2);

  const auto pretty_indexer = PrettyIndexer();
  StringVector pretty_index_vector;
  GPOperationVector branch_optimization_operations = dag_.BranchLengthOptimization();

  for (int gpcsp_idx = 0; gpcsp_idx < gpcsp_count; gpcsp_idx++) {
    double optimized_llh = optimized_per_pcsp_llhs[gpcsp_idx];
    double current_branch_length = 0.1;

    while (true) {
      run_counts[gpcsp_idx]++;

      EigenVectorXd new_branch_length_vector = optimized_branch_lengths;
      new_branch_length_vector[gpcsp_idx] = current_branch_length;
      GetEngine()->SetBranchLengths(new_branch_length_vector);

      PopulatePLVs();
      ComputeLikelihoods();

      double current_llh = GetEngine()->GetPerGPCSPLogLikelihoods(gpcsp_idx, 1)(0, 0);

      tracked_optimization_values.row(tracked_optimization_values.rows() - 1)
          << current_branch_length,
          current_llh;
      tracked_optimization_values.conservativeResize(
          tracked_optimization_values.rows() + 1, Eigen::NoChange);

      if (fabs(current_llh - optimized_llh) < 1e-5 | run_counts[gpcsp_idx] > 4) {
        break;
      } else {
        ProcessOperations(branch_optimization_operations);
        current_branch_length = GetEngine()->GetBranchLengths()[gpcsp_idx];
      }
    }
    pretty_index_vector.insert(pretty_index_vector.end(), run_counts[gpcsp_idx],
                               pretty_indexer.at(gpcsp_idx));
  }

  tracked_optimization_values.conservativeResize(tracked_optimization_values.rows() - 1,
                                                 Eigen::NoChange);
  StringEigenVectorXdVector result;
  result.reserve(tracked_optimization_values.rows());
  for (size_t i = 0; i < tracked_optimization_values.rows(); i++) {
    result.push_back({pretty_index_vector.at(i), tracked_optimization_values.row(i)});
  }
  return result;
}

StringVector GPInstance::PrettyIndexer() const {
  StringVector pretty_representation(dag_.BuildGPCSPIndexer().size());
  for (const auto &[pcsp, idx] : dag_.BuildGPCSPIndexer()) {
    pretty_representation[idx] = pcsp.PCSPToString();
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

StringEigenVectorXdVector GPInstance::PrettyIndexedMatrix(EigenConstMatrixXdRef m) {
  StringEigenVectorXdVector result;
  result.reserve(m.rows());
  const auto pretty_indexer = PrettyIndexer();
  // Assert(m.rows() <= pretty_indexer.size(), "m is too long in
  // PrettyIndexedMatrix");
  for (size_t i = 0; i < m.rows(); i++) {
    int idx = i % pretty_indexer.size();
    result.push_back({pretty_indexer.at(idx), m.row(i)});
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

StringEigenVectorXdVector
GPInstance::PrettyIndexedPerGPCSPBranchLengthsFromOptimization() {
  return PrettyIndexedMatrix(per_pcsp_branch_lengths_);
}

StringEigenVectorXdVector
GPInstance::PrettyIndexedPerGPCSPLogLikelihoodsFromOptimization() {
  return PrettyIndexedMatrix(per_pcsp_marg_lik_);
}

StringEigenVectorXdVector GPInstance::PrettyIndexedPerGPCSPLogLikelihoodSurfaces() {
  return PrettyIndexedMatrix(per_pcsp_lik_surfaces_);
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

void GPInstance::PerGPCSPLogLikelihoodsToCSV(const std::string &file_path) {
  CSV::StringDoubleVectorToCSV(PrettyIndexedPerGPCSPLogLikelihoods(), file_path);
}

void GPInstance::PerPCSPIndexedMatrixToCSV(
    StringEigenVectorXdVector per_pcsp_indexed_matrix, const std::string &file_path) {
  std::ofstream out_stream(file_path);

  for (const auto &[s, eigen] : per_pcsp_indexed_matrix) {
    out_stream << s;
    for (const auto &value : eigen) {
      out_stream << "," << std::setprecision(9) << value;
    }
    out_stream << std::endl;
  }
  if (out_stream.bad()) {
    Failwith("Failure writing to " + file_path);
  }
  out_stream.close();
}

void GPInstance::PerGPCSPBranchLengthsFromOptimizationToCSV(
    const std::string &file_path) {
  return PerPCSPIndexedMatrixToCSV(PrettyIndexedPerGPCSPBranchLengthsFromOptimization(),
                                   file_path);
}

void GPInstance::PerGPCSPLogLikelihoodsFromOptimizationToCSV(
    const std::string &file_path) {
  return PerPCSPIndexedMatrixToCSV(
      PrettyIndexedPerGPCSPLogLikelihoodsFromOptimization(), file_path);
}

void GPInstance::PerGPCSPLogLikelihoodSurfacesToCSV(const std::string &file_path) {
  std::ofstream out_stream(file_path);
  StringEigenVectorXdVector vect = PrettyIndexedPerGPCSPLogLikelihoodSurfaces();

  for (const auto &[s, eigen] : vect) {
    out_stream << s;
    for (const auto &value : eigen) {
      out_stream << "," << std::setprecision(9) << value;
    }
    out_stream << std::endl;
  }
  if (out_stream.bad()) {
    Failwith("Failure writing to " + file_path);
  }
  out_stream.close();
}

void GPInstance::TrackValuesFromOptimizationToCSV(const std::string &file_path) {
  return PerPCSPIndexedMatrixToCSV(TrackValuesFromOptimization(), file_path);
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
  const BitsetSizeMap &indexer = dag_.BuildGPCSPIndexer();
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
