// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// A visitor for GPOperations. See
// https://arne-mertz.de/2018/05/modern-c-features-stdvariant-and-stdvisit/

#ifndef SRC_GP_ENGINE_HPP_
#define SRC_GP_ENGINE_HPP_

#include "eigen_sugar.hpp"
#include "gp_operation.hpp"
#include "mmapped_plv.hpp"
#include "numerical_utils.hpp"
#include "quartet_hybrid_request.hpp"
#include "reindexer.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"
#include "site_pattern.hpp"
#include "subsplit_dag_graft.hpp"
#include "substitution_model.hpp"

class GPEngine {
 public:
  GPEngine(SitePattern site_pattern, size_t plv_count, size_t gpcsp_count,
           const std::string& mmap_file_path, double rescaling_threshold,
           EigenVectorXd sbn_prior, 
           EigenVectorXd unconditional_node_probabilities,
           EigenVectorXd inverted_sbn_prior);

  // These operators mean that we can invoke this class on each of the operations.
  void operator()(const GPOperations::ZeroPLV& op);
  void operator()(const GPOperations::SetToStationaryDistribution& op);
  void operator()(const GPOperations::IncrementWithWeightedEvolvedPLV& op);
  void operator()(const GPOperations::ResetMarginalLikelihood& op);
  void operator()(const GPOperations::IncrementMarginalLikelihood& op);
  void operator()(const GPOperations::Multiply& op);
  void operator()(const GPOperations::Likelihood& op);
  void operator()(const GPOperations::OptimizeBranchLength& op);
  void operator()(const GPOperations::UpdateSBNProbabilities& op);
  void operator()(const GPOperations::PrepForMarginalization& op);

  // Apply all operations in vector in order from beginning to end.
  void ProcessOperations(GPOperationVector operations);
  //
  void SetTransitionMatrixToHaveBranchLength(double branch_length);
  void SetTransitionAndDerivativeMatricesToHaveBranchLength(double branch_length);
  void SetTransitionMatrixToHaveBranchLengthAndTranspose(double branch_length);
  const Eigen::Matrix4d& GetTransitionMatrix() { return transition_matrix_; };
  //
  void SetBranchLengths(EigenVectorXd branch_lengths);
  void SetBranchLengthsToConstant(double branch_length);
  void ResetLogMarginalLikelihood();
  double GetLogMarginalLikelihood() const;
  EigenVectorXd GetBranchLengths() const;
  // This function returns a vector indexed by GPCSP such that the i-th entry
  // stores the log of the across-sites product of
  // (the marginal likelihood conditioned on a given GPCSP) *
  //     (the unconditional probability of i's parent subsplit).
  // That is, it's sum_m r^m(t) P(t -> s) p^m(s).
  // See lem:PerPCSPMarginalLikelihood.
  // #288 rename?
  EigenVectorXd GetPerGPCSPLogLikelihoods() const;
  // This override of GetPerGPCSPLogLikelihoods computes the marginal log
  // likelihood for GPCSPs in the range [start, start + length).
  EigenVectorXd GetPerGPCSPLogLikelihoods(size_t start, size_t length) const;
  // This is the full marginal likelihood sum restricted to trees containing a PCSP.
  // When we sum the log of eq:PerGPCSPComponentsOfFullMarginal over the sites, we get
  // out a term that is the number of sites times the log of the prior conditional PCSP
  // probability.
  EigenVectorXd GetPerGPCSPComponentsOfFullLogMarginal() const;
  // #288 reconsider this name
  EigenConstMatrixXdRef GetLogLikelihoodMatrix() const;
  //
  EigenConstVectorXdRef GetHybridMarginals() const;
  //
  EigenConstVectorXdRef GetSBNParameters() const;

  // Calculate a vector of likelihoods, one for each summand of the hybrid marginal.
  EigenVectorXd CalculateQuartetHybridLikelihoods(const QuartetHybridRequest& request);
  // Calculate the actual hybrid marginal and store it in the corresponding entry of
  // hybrid_marginal_log_likelihoods_.
  void ProcessQuartetHybridRequest(const QuartetHybridRequest& request);

  // Calculate a vector of likelihoods, one for each summand of the hybrid marginal.
  EigenVectorXd CalculateQuartetHybridLikelihoodsWithGraft(
      const SubsplitDAGGraft& graft, const QuartetHybridRequest& request);
  // Calculate the actual hybrid marginal and store it in the corresponding entry of
  // hybrid_marginal_log_likelihoods_.
  void ProcessQuartetHybridRequestWithGraft(const SubsplitDAGGraft& graft,
                                            const QuartetHybridRequest& request);

  //
  void PrintPLV(size_t plv_idx);

  // Use branch lengths from loaded sample as a starting point for optimization.
  void HotStartBranchLengths(const RootedTreeCollection& tree_collection,
                             const BitsetSizeMap& indexer);

  DoublePair LogLikelihoodAndDerivative(const GPOperations::OptimizeBranchLength& op);

  double PLVByteCount() const { return mmapped_master_plv_.ByteCount(); };

 private:
  //
  void InitializePLVsWithSitePatterns();

  void RescalePLV(size_t plv_idx, int amount);
  void AssertPLVIsFinite(size_t plv_idx, const std::string& message) const;
  std::pair<double, double> PLVMinMax(size_t plv_idx) const;
  // If a PLV all entries smaller than rescaling_threshold_ then rescale it up and
  // increment the corresponding entry in rescaling_counts_.
  void RescalePLVIfNeeded(size_t plv_idx);
  double LogRescalingFor(size_t plv_idx);

  void BrentOptimization(const GPOperations::OptimizeBranchLength& op);
  void GradientAscentOptimization(const GPOperations::OptimizeBranchLength& op);

  inline void PrepareUnrescaledPerPatternLikelihoodDerivatives(size_t src1_idx,
                                                               size_t src2_idx) {
    per_pattern_likelihood_derivatives_ =
        (plvs_.at(src1_idx).transpose() * derivative_matrix_ * plvs_.at(src2_idx))
            .diagonal()
            .array();
  }

  inline void PrepareUnrescaledPerPatternLikelihoods(size_t src1_idx, size_t src2_idx) {
    per_pattern_likelihoods_ =
        (plvs_.at(src1_idx).transpose() * transition_matrix_ * plvs_.at(src2_idx))
            .diagonal()
            .array();
  }

  // This function is used to compute the marginal log likelihood over all trees that
  // have a given PCSP. We assume that transition_matrix_ is as desired, and src1_idx
  // and src2_idx are the two PLV indices on either side of the PCSP.
  inline void PreparePerPatternLogLikelihoodsForGPCSP(size_t src1_idx,
                                                      size_t src2_idx) {
    per_pattern_log_likelihoods_ =
        (plvs_.at(src1_idx).transpose() * transition_matrix_ * plvs_.at(src2_idx))
            .diagonal()
            .array()
            .log() +
        LogRescalingFor(src1_idx) + LogRescalingFor(src2_idx);
  }

 public:
  static constexpr double default_rescaling_threshold_ = 1e-40;
  // Initial branch length during first branch length opimization.
  static constexpr double default_branch_length_ = 0.1;

 private:
  // Absolute lower bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double min_log_branch_length_ = -13.9;
  // Absolute upper bound for possible branch lengths during optimization (in log
  // space).
  static constexpr double max_log_branch_length_ = 1.1;
  // Precision used for checking convergence of branch length optimization.
  int significant_digits_for_optimization_ = 6;
  //
  double relative_tolerance_for_optimization_ = 1e-2;
  // Step size used for gradient-based branch length optimization.
  double step_size_for_optimization_ = 5e-4;
  // Number of iterations allowed for branch length optimization.
  size_t max_iter_for_optimization_ = 1000;

  // ** Per-Node Data

  // Total number of PLVs across entire DAG. Proportional to the number of nodes in DAG.
  // plv_per_node * node_count_without_root.
  size_t plv_count_;
  // Master PLV: Large data block of virtual memory for Partial Likelihood Vectors.
  // Subdivided into sections for plvs_.
  std::unique_ptr<MmappedNucleotidePLV> mmapped_master_plv_ptr_ = nullptr;
  MmappedNucleotidePLV mmapped_master_plv_;
  // Partial Likelihood Vectors.
  // plvs_ store the following (see GPDAG::GetPLVIndexStatic):
  // [0, num_nodes): p(s).
  // [num_nodes, 2*num_nodes): phat(s).
  // [2*num_nodes, 3*num_nodes): phat(s_tilde).
  // [3*num_nodes, 4*num_nodes): rhat(s) = rhat(s_tilde).
  // [4*num_nodes, 5*num_nodes): r(s).
  // [5*num_nodes, 6*num_nodes): r(s_tilde).
  NucleotidePLVRefVector plvs_;
  // Rescaling count for each plv.
  EigenVectorXi rescaling_counts_;
  // For hybrid marginal calculations. #328
  // The PLV coming down from the root.
  EigenMatrixXd quartet_root_plv_;
  // The R-PLV pointing leafward from s.
  EigenMatrixXd quartet_r_s_plv_;
  // The Q-PLV pointing leafward from s.
  EigenMatrixXd quartet_q_s_plv_;
  // The R-PLV pointing leafward from t.
  EigenMatrixXd quartet_r_sorted_plv_;

  // ** Per-Edge Data

  // Total number of edges in DAG.
  size_t gpcsp_count_;
  // branch_lengths_, q_, etc. are indexed in the same way as sbn_parameters_ in
  // gp_instance.
  EigenVectorXd branch_lengths_;
  // TODO: Add data array to store the sbn_prior (in linear space?)
  // During initialization, stores the SBN prior.
  // After UpdateSBNProbabilities(), stores the SBN probabilities.
  // Stored in log space.
  EigenVectorXd q_;
  //
  EigenVectorXd unconditional_node_probabilities_;
  //
  EigenVectorXd inverted_sbn_prior_;
  // The number of rows is equal to the number of GPCSPs.
  // The number of columns is equal to the number of site patterns.
  // The rows are indexed in the same way as branch_lengths_ and q_.
  // Entry (i,j) stores the marginal log likelihood over all trees that include
  // a GPCSP corresponding to index i at site j.
  EigenMatrixXd log_likelihoods_;
  // The length of this vector is equal to the number of site patterns.
  // Entry j stores the marginal log likelihood over all trees at site pattern
  // j.
  EigenVectorXd log_marginal_likelihood_;
  // This vector is indexed by the GPCSPs and stores the hybrid marginals if they are
  // available.
  EigenVectorXd hybrid_marginal_log_likelihoods_;
  // Descriptor containing all taxons and sequence alignments.
  SitePattern site_pattern_;
  // Rescaling threshold factor to prevent under/overflow errors.
  const double rescaling_threshold_;
  // Rescaling threshold in log space.
  const double log_rescaling_threshold_;
  // Internal "temporaries" useful for likelihood and derivative calculation.
  EigenVectorXd per_pattern_log_likelihoods_;
  EigenVectorXd per_pattern_likelihoods_;
  EigenVectorXd per_pattern_likelihood_derivatives_;
  EigenVectorXd per_pattern_likelihood_derivative_ratios_;

  // ** Substitution Model

  // When we change from JC69Model, check that we are actually doing transpose in
  // leafward calculations.
  JC69Model substitution_model_;
  Eigen::Matrix4d eigenmatrix_ = substitution_model_.GetEigenvectors().reshaped(4, 4);
  Eigen::Matrix4d inverse_eigenmatrix_ =
      substitution_model_.GetInverseEigenvectors().reshaped(4, 4);
  Eigen::Vector4d eigenvalues_ = substitution_model_.GetEigenvalues();
  Eigen::Vector4d diagonal_vector_;
  Eigen::DiagonalMatrix<double, 4> diagonal_matrix_;
  Eigen::Matrix4d transition_matrix_;
  Eigen::Matrix4d derivative_matrix_;
  Eigen::Vector4d stationary_distribution_ = substitution_model_.GetFrequencies();
  EigenVectorXd site_pattern_weights_;

 public:
  // TODO:
  // ** Calculate Hybrid Likelihoods with Graft
  // A SubsplitDAGGraft is a proposed (graft) set of nodes and edges to be added to the
  // (host) SubsplitDAG but have not been formally indexed into the data structure.
  // These functions can operate on the DAG and the graft as if they were a singular
  // object.

  // TODO: Update after modify/grafting DAG.
  // Update stats and data to reflect current DAG state.
  void InitEngine(size_t plv_count, size_t gpcsp_count,
                  const std::string& mmap_file_path, double rescaling_threshold,
                  EigenVectorXd sbn_prior,
                  EigenVectorXd unconditional_node_probabilities,
                  EigenVectorXd inverted_sbn_prior);

  // Initialize engine's graftDAG data with current size.
  void InitEngineForDAG(const size_t plv_count, const size_t gpcsp_count,
                        const std::string& mmap_file_path);
  //
  void ResizeAfterModifyingDAG(const size_t old_plv_count, 
                               const size_t new_plv_count,
                               const size_t old_gpcsp_count,
                               const size_t new_gpcsp_count);
  // Resize data members to store SubsplitDAG after modification.
  void UpdateAfterModifyingDAG(const size_t old_plv_count, 
                               const size_t new_plv_count,
                               const size_t old_gpcsp_count,
                               const size_t new_gpcsp_count,
                               const std::string& mmap_file_path,
                               const SizeVector& node_reindexer,
                               const SizeVector& edge_reindexer);

  // Initialize engine's graftDAG data with current size.
  void InitEngineForGraftDAG(const size_t plv_count, const size_t gpcsp_count,
                             const std::string& graft_mmap_file_path);
  //
  void ResizeAfterGraftingDAG(const size_t old_plv_count,
                              const size_t new_plv_count,
                              const size_t old_gpcsp_count,
                              const size_t new_gpcsp_count);
  // Append new data members to store SubsplitDAG growth after modification.
  void UpdateAfterGraftingDAG(const size_t old_plv_count, const size_t new_plv_count,
                              const size_t old_gpcsp_count,
                              const size_t new_gpcsp_count,
                              const std::string& mmap_file_path,
                              const SizeVector& node_reindexer,
                              const SizeVector& edge_reindexer);

  //
  EigenVectorXd ComputeAllNNIsPerPCSPLikelihood();
  // Compute PerPCSP Likelihoods for given NNI pair.
  EigenVectorXd ComputePerNNIPerPCSPLikelihood(
      const size_t parent_node_id, const size_t child_node_id,
      SizeVector& parents_of_parent_nodes, SizeVector& children_of_parent_nodes,
      SizeVector& parents_of_child_nodes, SizeVector& children_of_child_nodes,
      SizeVector& parents_of_parent_edges, SizeVector& children_of_parent_edges,
      SizeVector& parents_of_child_edges, SizeVector& children_of_child_edges);

 private:
  // ** NNI Engine for GraftedDAG.
  // // NNI per-node data.
  // size_t graft_plv_count_;
  // MmappedNucleotidePLV graft_mmapped_master_plv_;
  // NucleotidePLVRefVector graft_plvs_;
  // EigenMatrixXd graft_quartet_root_plv_;
  // EigenMatrixXd graft_quartet_r_s_plv_;
  // EigenMatrixXd graft_quartet_q_s_plv_;
  // EigenMatrixXd graft_quartet_r_sorted_plv_;
  // EigenVectorXi graft_rescaling;
  // // NNI per-edge data.
  // size_t graft_gpcsp_count;
  // EigenVectorXd graft_branch_lengths_;
  // EigenVectorXd graft_log_likelihoods_;
  // EigenVectorXd graft_log_marginal_likelihoods_;
  // EigenVectorXd graft_hybrid_marginal_log_likelihoods_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("GPEngine") {
  EigenVectorXd empty_vector;
  SitePattern hello_site_pattern = SitePattern::HelloSitePattern();
  GPEngine engine(hello_site_pattern, 6 * 5, 5, "_ignore/mmapped_plv.data",
                  GPEngine::default_rescaling_threshold_, empty_vector, empty_vector,
                  empty_vector);
  engine.SetTransitionMatrixToHaveBranchLength(0.75);
  // Computed directly:
  // https://en.wikipedia.org/wiki/Models_of_DNA_evolution#JC69_model_%28Jukes_and_Cantor_1969%29
  CHECK(fabs(0.52590958087 - engine.GetTransitionMatrix()(0, 0)) < 1e-10);
  CHECK(fabs(0.1580301397 - engine.GetTransitionMatrix()(0, 1)) < 1e-10);
}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_GP_ENGINE_HPP_
