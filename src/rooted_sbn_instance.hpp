// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include "csv.hpp"
#include "generic_sbn_instance.hpp"
#include "phylo_flags.hpp"
#include "rooted_gradient_transforms.hpp"
#include "rooted_sbn_support.hpp"

using PreRootedSBNInstance = GenericSBNInstance<RootedTreeCollection, RootedSBNSupport,
                                                RootedIndexerRepresentation>;
template class GenericSBNInstance<RootedTreeCollection, RootedSBNSupport,
                                  RootedIndexerRepresentation>;

class RootedSBNInstance : public PreRootedSBNInstance {
 public:
  using PreRootedSBNInstance::PreRootedSBNInstance;

  // ** SBN-related items

  // Turn an IndexerRepresentation into a string representation of the underlying
  // bitsets. This is really just so that we can make a test of indexer
  // representations.
  StringSet StringIndexerRepresentationOf(
      const RootedIndexerRepresentation& indexer_representation) const;
  StringSet StringIndexerRepresentationOf(const Node::NodePtr& topology,
                                          size_t out_of_sample_index) const;

  // Make a map from each subsplit to its overall probability when we sample a tree from
  // the SBN.
  BitsetDoubleMap UnconditionalSubsplitProbabilities() const;
  void UnconditionalSubsplitProbabilitiesToCSV(const std::string& csv_path) const;

  // ** Phylogenetic likelihood

  std::vector<double> LogLikelihoods();
  std::vector<double> LogLikelihoods(std::optional<PhyloFlags> flags);
  std::vector<double> LogLikelihoods(StringVector& flags, bool use_defaults = true);
  std::vector<double> LogLikelihoods(StringDoubleVector& flags,
                                     bool use_defaults = true);
  std::vector<double> UnrootedLogLikelihoods();
  // For each loaded tree, return the phylogenetic gradient.
  std::vector<PhyloGradient> PhyloGradients();
  std::vector<PhyloGradient> PhyloGradients(std::optional<PhyloFlags> flags);
  std::vector<PhyloGradient> PhyloGradients(StringVector& flags,
                                            bool use_defaults = true);
  std::vector<PhyloGradient> PhyloGradients(StringDoubleVector& flags,
                                            bool use_defaults = true);
  // For each loaded tree, return the log determinant jacobian of the height transform.
  std::vector<double> LogDetJacobianHeightTransform();

  // ** I/O

  void ReadNewickFile(const std::string& fname);
  void ReadNexusFile(const std::string& fname);

  void SetDatesToBeConstant(bool initialize_time_trees_using_branch_lengths);
  void ParseDatesFromTaxonNames(bool initialize_time_trees_using_branch_lengths);
  void ParseDatesFromCSV(const std::string& csv_path,
                         bool initialize_time_trees_using_branch_lengths);
};

#ifdef DOCTEST_LIBRARY_INCLUDED

#include "doctest_constants.hpp"

// Centered finite difference approximation of the derivative wrt rate.
std::vector<double> DerivativeStrictClock(RootedSBNInstance& inst) {
  double eps = 0.00000001;
  std::vector<double> rates;
  std::vector<double> gradients;

  for (auto& tree : inst.tree_collection_.trees_) {
    rates.push_back(tree.rates_[0]);
    tree.rates_.assign(tree.rates_.size(), rates.back() - eps);
  }
  auto lm = inst.LogLikelihoods();

  int i = 0;
  for (auto& tree : inst.tree_collection_.trees_) {
    tree.rates_.assign(tree.rates_.size(), rates[i++] + eps);
  }
  auto lp = inst.LogLikelihoods();

  for (size_t index = 0; index < lm.size(); index++) {
    gradients.push_back((lp[index] - lm[index]) / (2. * eps));
  }
  return gradients;
}

// Centered finite difference approximation of the derivative wrt to each rate.
std::vector<std::vector<double>> DerivativeRelaxedClock(RootedSBNInstance& inst) {
  double eps = 0.00000001;
  std::vector<std::vector<double>> gradients;
  std::vector<double> lp;
  std::vector<double> lm;
  size_t edge_count = inst.TaxonCount() * 2 - 2;

  for (size_t index = 0; index < edge_count; index++) {
    std::vector<double> gradient;
    std::vector<double> rates;
    for (size_t i = 0; i < inst.tree_collection_.TreeCount(); i++) {
      double value = inst.tree_collection_.trees_[i].rates_[index];
      rates.push_back(value);
      inst.tree_collection_.trees_[i].rates_[index] = rates.back() - eps;
    }
    lm = inst.LogLikelihoods();

    for (size_t i = 0; i < inst.tree_collection_.TreeCount(); i++) {
      inst.tree_collection_.trees_[i].rates_[index] = rates[i] + eps;
    }
    lp = inst.LogLikelihoods();

    for (size_t i = 0; i < inst.tree_collection_.TreeCount(); i++) {
      inst.tree_collection_.trees_[i].rates_[index] = rates[i];
      gradient.push_back((lp[i] - lm[i]) / (2. * eps));
    }

    gradients.push_back(gradient);
  }
  return gradients;
}

RootedSBNInstance MakeFiveTaxonRootedInstance() {
  RootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon_rooted.nwk");
  inst.ProcessLoadedTrees();
  return inst;
}

TEST_CASE("RootedSBNInstance: subsplit support and TrainSimpleAverage") {
  auto inst = MakeFiveTaxonRootedInstance();
  auto pretty_indexer = inst.PrettyIndexer();
  StringSet pretty_indexer_set{pretty_indexer.begin(), pretty_indexer.end()};
  // The indexer_ is to index the sbn_parameters_. Note that neither of these
  // data structures attempt to catalog the complete collection of rootsplits or
  // PCSPs, but just those that are present in the the input trees.
  //
  // The indexer_ and sbn_parameters_ are laid out as follows (I'll just call it
  // the "index" in what follows). Say there are rootsplit_count rootsplits in
  // the support.
  // The first rootsplit_count entries of the index are assigned to the
  // rootsplits (again, those rootsplits that are present for some rooting of
  // the unrooted input trees). The rest of the entries of the index are laid out as
  // blocks of parameters for PCSPs that share the same parent. Take a look at the
  // description of PCSP bitsets (and the unit tests) in bitset.hpp to understand the
  // notation used here.
  //
  // In contrast to the unrooted case, we can write out the pretty indexer here and
  // verify it by hand. There is the block structure in which the two children of
  // 10000|01111 are grouped together.
  StringSet correct_pretty_indexer_set{
      "00000|11111|00111",  // ((x0,x1),(x2,(x3,x4)))
      "00000|11111|01111",  // (x0,(((x1,x3),x2),x4)) and ((x1,((x2,x4),x3)),x0)
      "00000|11111|00010",  // (x3,((x0,(x4,x1)),x2))
      "00100|01010|00010",  // ((x1,x3),x2)
      "00111|11000|01000",  // ((x0,x1),(x2,(x3,x4)))
      "00100|00011|00001",  // (x2,(x3,x4))
      "11000|00111|00011",  // ((x0,x1),(x2,(x3,x4)))
      "00100|11001|01001",  // ((x0,(x4,x1)),x2)
      "10000|01001|00001",  // (x0,(x4,x1))
      "01000|00111|00010",  // (x1,((x2,x4),x3))
      "10000|01111|00001",  // (x0,(((x1,x3),x2),x4))
      "10000|01111|00111",  // ((x1,((x2,x4),x3)),x0)
      "00010|00101|00001",  // ((x2,x4),x3)
      "00001|01110|00100",  // (((x1,x3),x2),x4)
      "00010|11101|00100"   // (x3,((x0,(x4,x1)),x2))
  };
  CHECK_EQ(pretty_indexer_set, correct_pretty_indexer_set);

  // Test of rooted IndexerRepresentationOf.
  // Topology is ((0,1),(2,(3,4)));, or with internal nodes ((0,1)5,(2,(3,4)6)7)8;
  auto indexer_test_rooted_topology = Node::OfParentIdVector({5, 5, 7, 6, 6, 8, 7, 8});
  auto correct_rooted_indexer_representation =
      StringSet({"00000|11111|00111", "11000|00111|00011", "00100|00011|00001",
                 "00111|11000|01000"});
  CHECK_EQ(inst.StringIndexerRepresentationOf(indexer_test_rooted_topology,
                                              out_of_sample_index),
           correct_rooted_indexer_representation);

  inst.TrainSimpleAverage();
  StringVector correct_taxon_names({"x0", "x1", "x2", "x3", "x4"});
  CHECK_EQ(inst.SBNSupport().TaxonNames(), correct_taxon_names);
  StringDoubleVector correct_parameters({{"00000|11111|00111", 0.25},
                                         {"00000|11111|01111", 0.5},
                                         {"00000|11111|00010", 0.25},
                                         {"00100|01010|00010", 1},
                                         {"00111|11000|01000", 1},
                                         {"00100|00011|00001", 1},
                                         {"11000|00111|00011", 1},
                                         {"00100|11001|01001", 1},
                                         {"10000|01001|00001", 1},
                                         {"01000|00111|00010", 1},
                                         {"10000|01111|00001", 0.5},
                                         {"10000|01111|00111", 0.5},
                                         {"00010|00101|00001", 1},
                                         {"00001|01110|00100", 1},
                                         {"00010|11101|00100", 1}});
  std::sort(correct_parameters.begin(), correct_parameters.end());
  auto parameters = inst.PrettyIndexedSBNParameters();
  std::sort(parameters.begin(), parameters.end());
  CHECK_EQ(correct_parameters.size(), parameters.size());
  for (size_t i = 0; i < correct_parameters.size(); i++) {
    CHECK_EQ(correct_parameters[i].first, parameters[i].first);
    CHECK_LT(fabs(correct_parameters[i].second - parameters[i].second), 1e-8);
  }
}

TEST_CASE("RootedSBNInstance: UnconditionalSubsplitProbabilities") {
  RootedSBNInstance inst("rooted instance");
  inst.ReadNewickFile("data/five_taxon_rooted_more.nwk");
  inst.ProcessLoadedTrees();
  inst.TrainSimpleAverage();
  // See diagram at https://github.com/phylovi/bito/issues/349#issuecomment-898022916
  // Numbering in comments is...                               node: subsplit.
  StringDoubleMap correct_parameters({{"1100000111", 0.5},  // 10: 01|234
                                      {"1000001111", 0.3},  // 15: 0|1234
                                      {"1110100010", 0.2},  // 19: 0124|3
                                      {"1100100100", 0.2},  // 18: 014|2
                                      {"0100000111", 0.1},  // 14: 1|234
                                      {"0111000001", 0.2},  // 13: 123|4
                                      {"0101000100", 0.2},  // 12: 13|2
                                      {"1000001001", 0.2},  // 17: 0|14
                                      {"0010000011", 0.4},  //  8: 2|34
                                      {"0011000001", 0.2},  //  6: 23|4
                                      {"1000001000", 0.5},  //  9: 0|1
                                      {"0100000010", 0.2},  // 11: 1|3
                                      {"0100000001", 0.2},  // 16: 1|4
                                      {"0010000010", 0.2},  //  5: 2|3
                                      {"0001000001", 0.4}}  //  7: 3|4
  );

  auto subsplit_probabilities = inst.UnconditionalSubsplitProbabilities();
  CHECK_EQ(correct_parameters.size(), subsplit_probabilities.size());
  for (const auto& [subsplit, probability] : subsplit_probabilities) {
    CHECK_LT(fabs(correct_parameters.at(subsplit.ToString()) - probability), 1e-8);
  }
}

// Instance SA-trained on a sample of 20-taxon trees.
RootedSBNInstance MakeRootedSimpleAverageInstance() {
  RootedSBNInstance inst("rooted instance");
  inst.ReadNewickFile("data/rooted_simple_average.nwk");
  inst.ProcessLoadedTrees();
  inst.TrainSimpleAverage();
  return inst;
}

TEST_CASE("RootedSBNInstance: TrainSimpleAverage on 20 taxa") {
  auto inst = MakeRootedSimpleAverageInstance();
  auto results = inst.PrettyIndexedSBNParameters();
  // Values confirmed with
  // https://github.com/mdkarcher/vbsupertree/commit/b7f87f711e8a1044b7c059b5a92e94c117d8cee1
  auto correct_map =
      CSV::StringDoubleMapOfCSV("data/rooted_simple_average_results.csv");
  for (const auto& [found_string, found_probability] : results) {
    CHECK(fabs(found_probability - correct_map.at(found_string)) < 1e-6);
  }
}

RootedSBNInstance MakeFluInstance(bool initialize_time_trees) {
  RootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/fluA.tree");
  inst.ParseDatesFromTaxonNames(initialize_time_trees);
  inst.ReadFastaFile("data/fluA.fa");
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  inst.PrepareForPhyloLikelihood(simple_specification, 1);
  return inst;
}

TEST_CASE("RootedSBNInstance: gradients") {
  auto inst = MakeFluInstance(true);
  for (auto& tree : inst.tree_collection_.trees_) {
    tree.rates_.assign(tree.rates_.size(), 0.001);
  }

  auto likelihood = inst.LogLikelihoods();
  double physher_ll = -4777.616349;
  double physher_jacobian = -9.25135166;
  double physher_ll_jacobian = physher_ll + physher_jacobian;
  CHECK_LT(fabs(likelihood[0] - physher_ll_jacobian), 0.0001);

  auto gradients = inst.PhyloGradients();
  std::vector<double> physher_gradients = {
      -0.593654, 6.441290,   11.202945, 5.173924,  -0.904631, 2.731402,   3.157131,
      7.082914,  10.305417,  13.988206, 20.709336, 48.897993, 99.164949,  130.205747,
      17.314019, 21.033290,  -1.336335, 12.259822, 22.887291, 27.176564,  47.487426,
      3.637276,  12.955169,  15.315953, 83.254605, -3.806996, 105.385095, 4.874023,
      22.754466, 6.036534,   25.651478, 29.535185, 29.598789, 1.817247,   10.598685,
      76.259248, 56.481423,  10.679778, 6.587179,  3.330556,  -4.622247,  33.417304,
      63.415767, 188.809515, 23.540875, 17.421076, 1.222568,  22.372012,  34.239511,
      3.486115,  4.098873,   13.200954, 19.726890, 96.808738, 4.240029,   7.414585,
      48.871694, 3.488516,   82.969065, 9.009334,  8.032474,  3.981016,   6.543650,
      53.702423, 37.835952,  2.840831,  7.517186,  19.936861};
  for (size_t i = 0; i < physher_gradients.size(); i++) {
    CHECK_LT(fabs(gradients[0].gradient_[PhyloGradient::ratios_root_height_key_][i] -
                  physher_gradients[i]),
             0.0001);
  }
  CHECK_LT(fabs(gradients[0].log_likelihood_ - physher_ll), 0.0001);
}

TEST_CASE("RootedSBNInstance: clock gradients") {
  auto inst = MakeFluInstance(true);
  for (auto& tree : inst.tree_collection_.trees_) {
    tree.rates_.assign(tree.rates_.size(), 0.001);
  }

  auto likelihood = inst.LogLikelihoods();
  double physher_ll = -4777.616349;
  double physher_jacobian = -9.25135166;
  double physher_ll_jacobian = physher_ll + physher_jacobian;
  CHECK_LT(fabs(likelihood[0] - physher_ll_jacobian), 0.0001);

  // Gradient with a strict clock.
  auto gradients_strict = inst.PhyloGradients();
  std::vector<double> gradients_strict_approx = DerivativeStrictClock(inst);
  CHECK_LT(fabs(gradients_strict[0].gradient_[PhyloGradient::clock_model_key_][0] -
                gradients_strict_approx[0]),
           0.001);
  CHECK_LT(fabs(gradients_strict[0].log_likelihood_ - physher_ll), 0.001);

  // Gradient with a "relaxed" clock.
  auto& tree = inst.tree_collection_.trees_[0];
  // Make a clock with some rate variation.
  for (size_t i = 0; i < tree.rates_.size(); i++) {
    tree.rates_[i] *= i % 3 + 1.0;
  }
  tree.rate_count_ = tree.rates_.size();

  auto gradients_relaxed = inst.PhyloGradients();
  auto gradients_relaxed_approx = DerivativeRelaxedClock(inst);

  for (size_t j = 0; j < gradients_relaxed_approx.size(); j++) {
    CHECK_LT(fabs(gradients_relaxed[0].gradient_[PhyloGradient::clock_model_key_][j] -
                  gradients_relaxed_approx[j][0]),
             0.001);
  }
}

TEST_CASE("RootedSBNInstance: GTR gradients") {
  auto inst = MakeFluInstance(true);
  PhyloModelSpecification gtr_specification{"GTR", "constant", "strict"};
  inst.PrepareForPhyloLikelihood(gtr_specification, 1);
  for (auto& tree : inst.tree_collection_.trees_) {
    tree.rates_.assign(tree.rates_.size(), 0.001);
  }
  auto param_block_map = inst.GetPhyloModelParamBlockMap();
  EigenVectorXdRef frequencies = param_block_map.at(GTRModel::frequencies_key_);
  EigenVectorXdRef rates = param_block_map.at(GTRModel::rates_key_);
  frequencies << 0.1, 0.2, 0.3, 0.4;
  rates << 0.05, 0.1, 0.15, 0.20, 0.25, 0.25;

  auto likelihood = inst.LogLikelihoods();
  double phylotorch_ll = -5221.438941335706;
  double physher_jacobian = -9.25135166;
  double expected_ll_jacobian = phylotorch_ll + physher_jacobian;
  CHECK_LT(fabs(likelihood[0] - expected_ll_jacobian), 0.001);

  auto gradients = inst.PhyloGradients();
  std::vector<double> phylotorch_gradients = {49.06451538, 151.83105912, 26.40235659,
                                              -8.25135661, 75.29759338,  352.56545247,
                                              90.07046995, 30.12301652};
  for (size_t i = 0; i < phylotorch_gradients.size(); i++) {
    CHECK_LT(fabs(gradients[0].gradient_[PhyloGradient::substitution_model_key_][i] -
                  phylotorch_gradients[i]),
             0.001);
  }
  CHECK_LT(fabs(gradients[0].log_likelihood_ - phylotorch_ll), 0.001);
}

TEST_CASE("RootedSBNInstance: HKY gradients") {
  auto inst = MakeFluInstance(true);
  PhyloModelSpecification specification{"HKY", "constant", "strict"};
  inst.PrepareForPhyloLikelihood(specification, 1);
  for (auto& tree : inst.tree_collection_.trees_) {
    tree.rates_.assign(tree.rates_.size(), 0.001);
  }
  auto param_block_map = inst.GetPhyloModelParamBlockMap();
  EigenVectorXdRef frequencies =
      param_block_map.at(SubstitutionModel::frequencies_key_);
  EigenVectorXdRef rates = param_block_map.at(SubstitutionModel::rates_key_);
  frequencies << 0.1, 0.2, 0.3, 0.4;
  rates << 3.0;

  auto likelihood = inst.LogLikelihoods();
  double phylotorch_ll = -4931.770106816288;
  double physher_jacobian = -9.25135166;
  double expected_ll_jacobian = phylotorch_ll + physher_jacobian;
  CHECK_LT(fabs(expected_ll_jacobian - likelihood[0]), 0.001);

  auto gradients = inst.PhyloGradients();
  std::vector<double> phylotorch_gradients = {18.218397759598506, 309.56536079428355,
                                              47.15713892857574, 42.98132033283943};
  for (size_t i = 0; i < phylotorch_gradients.size(); i++) {
    CHECK_LT(
        fabs(gradients[0].gradient_["substitution_model"][i] - phylotorch_gradients[i]),
        0.001);
  }
  CHECK_LT(fabs(phylotorch_ll - gradients[0].log_likelihood_), 0.0001);
}

TEST_CASE("RootedSBNInstance: Weibull gradients") {
  auto inst = MakeFluInstance(true);
  PhyloModelSpecification weibull_specification{"JC69", "weibull+4", "strict"};
  inst.PrepareForPhyloLikelihood(weibull_specification, 1);
  for (auto& tree : inst.tree_collection_.trees_) {
    tree.rates_.assign(tree.rates_.size(), 0.001);
  }
  auto param_block_map = inst.GetPhyloModelParamBlockMap();
  param_block_map.at(WeibullSiteModel::shape_key_).setConstant(0.1);

  auto likelihood = inst.LogLikelihoods();
  double physher_ll = -4618.2062529058;
  double physher_jacobian = -9.25135166;
  double physher_ll_jacobian = physher_ll + physher_jacobian;
  CHECK_LT(fabs(likelihood[0] - physher_ll_jacobian), 0.0001);

  // Gradient wrt Weibull site model.
  auto gradients = inst.PhyloGradients();
  double physher_gradient = -5.231329;
  CHECK_LT(fabs(gradients[0].gradient_["site_model"][0] - physher_gradient), 0.001);
  CHECK_LT(fabs(gradients[0].log_likelihood_ - physher_ll), 0.001);
}

TEST_CASE("RootedSBNInstance: parsing dates") {
  RootedSBNInstance inst("charlie");
  inst.ReadNexusFile("data/test_beast_tree_parsing.nexus");
  inst.ParseDatesFromTaxonNames(true);
  std::vector<double> dates;
  for (const auto& [tag, date] : inst.tree_collection_.GetTagDateMap()) {
    std::ignore = tag;
    dates.push_back(date);
  }
  std::sort(dates.begin(), dates.end());
  CHECK_EQ(dates[0], 0);
  CHECK_EQ(dates.back(), 80.0);

  RootedSBNInstance alt_inst("betty");
  alt_inst.ReadNexusFile("data/test_beast_tree_parsing.nexus");
  alt_inst.tree_collection_.ParseDatesFromCSV("data/test_beast_tree_parsing.csv", true);
  CHECK_EQ(inst.tree_collection_.GetTagDateMap(),
           alt_inst.tree_collection_.GetTagDateMap());
}

TEST_CASE("RootedSBNInstance: uninitialized time trees raise an exception") {
  auto inst = MakeFluInstance(false);
  CHECK_THROWS(inst.PhyloGradients());
}

TEST_CASE("RootedSBNInstance: reading SBN parameters from a CSV") {
  auto inst = MakeFiveTaxonRootedInstance();
  inst.ReadSBNParametersFromCSV("data/test_modifying_sbn_parameters.csv");
  auto pretty_indexer = inst.PrettyIndexer();
  auto gpcsp_it =
      std::find(pretty_indexer.begin(), pretty_indexer.end(), "10000|01111|00001");
  CHECK(gpcsp_it != pretty_indexer.end());
  auto gpcsp_idx = std::distance(pretty_indexer.begin(), gpcsp_it);
  CHECK_LT(fabs(inst.sbn_parameters_[gpcsp_idx] - log(0.15)), 1e-8);
  inst.SetSBNParameters({}, false);
  CHECK_EQ(inst.sbn_parameters_[gpcsp_idx], DOUBLE_MINIMUM);
  CHECK_THROWS(inst.SetSBNParameters({{"10000|01111|00001", -5.}}, false));
}

TEST_CASE("RootedSBNInstance: SBN parameter round trip") {
  std::string csv_test_file_path = "_ignore/for_sbn_parameter_round_trip.csv";
  auto inst = MakeRootedSimpleAverageInstance();
  auto original_normalized_sbn_parameters = inst.NormalizedSBNParameters();
  inst.SBNParametersToCSV(csv_test_file_path);
  inst.ReadSBNParametersFromCSV(csv_test_file_path);
  auto reloaded_normalized_sbn_parameters = inst.NormalizedSBNParameters();
  CheckVectorXdEquality(original_normalized_sbn_parameters,
                        reloaded_normalized_sbn_parameters, 1e-6);
}

TEST_CASE("RootedSBNInstance: gradient requests") {
  // GP Instance defaults for gradients.
  auto CreateNewInstance = []() {
    auto inst = MakeFluInstance(true);
    PhyloModelSpecification gtr_specification{"GTR", "constant", "strict"};
    inst.PrepareForPhyloLikelihood(gtr_specification, 1);
    for (auto& tree : inst.tree_collection_.trees_) {
      tree.rates_.assign(tree.rates_.size(), 0.001);
    }
    auto param_block_map = inst.GetPhyloModelParamBlockMap();
    EigenVectorXdRef frequencies = param_block_map.at(GTRModel::frequencies_key_);
    EigenVectorXdRef rates = param_block_map.at(GTRModel::rates_key_);
    frequencies << 0.1, 0.2, 0.3, 0.4;
    rates << 0.05, 0.1, 0.15, 0.20, 0.25, 0.25;
    return inst;
  };
  // "Golden" instance for determining correctness.
  auto gold_inst = CreateNewInstance();
  auto gold_likelihoods = gold_inst.LogLikelihoods();
  auto gold_gradients = gold_inst.PhyloGradients();
  size_t num_trees = gold_inst.tree_collection_.trees_.size();

  // Test All "include" options for gradient flags.
  StringVector gradient_flags = {
      PhyloFlags::gradient_options_.clock_model_rates_.flag_,
      PhyloFlags::gradient_options_.ratios_root_height_.flag_,
      PhyloFlags::gradient_options_.site_model_parameters_.flag_,
      PhyloFlags::gradient_options_.substitution_model_rates_.flag_,
      PhyloFlags::gradient_options_.substitution_model_frequencies_.flag_,
  };

  // Iterate through all "include" gradient flag combinations.
  auto IterateOverAllCombinations =
      [&](StringVector& all_flags,
          std::function<void(StringVector&, StringVector&, StringVector&)> func) {
        size_t num_flags = all_flags.size();
        size_t num_combinations = pow(2, num_flags);
        for (size_t i = 0; i < num_combinations; i++) {
          auto used_flags = StringVector();
          auto unused_flags = StringVector();
          // Split between groups of used and unused flags.
          for (size_t j = 1, k = 0; j < num_combinations; j <<= 1, k += 1) {
            if ((j & i)) {
              used_flags.push_back(all_flags[k]);
            } else {
              unused_flags.push_back(all_flags[k]);
            }
          }
          func(used_flags, unused_flags, all_flags);
        }
      };

  // Test that expected flagged keys are populated with correct data,
  // and that unflagged keys are not stored in map.
  auto ComparePhyloGradients = [&](StringVector& used_flags, StringVector& unused_flags,
                                   StringVector& all_flags) {
    // Create instance and run phylogradients with used_flags.
    auto inst = CreateNewInstance();
    auto gradients = inst.PhyloGradients(used_flags, false);
    // Check that proper fields are populated correctly.
    for (size_t i = 0; i < num_trees; i++) {
      auto& gold_grad_map = gold_gradients[i].gradient_;
      auto& grad_map = gradients[i].gradient_;
      for (auto& flag : used_flags) {
        for (size_t i = 0; i < gold_grad_map[flag].size(); i++) {
          Assert(
              std::abs(gold_grad_map[flag][i] - grad_map[flag][i]) < 0.01,
              std::string(
                  "gold_grad_map and grad_map do not have the same value for flag: " +
                  flag));
        }
      }
      for (auto& flag : unused_flags) {
        Assert(grad_map.find(flag) == grad_map.end(),
               std::string("grad_map has a key that should not exist: " + flag));
      }
    }
  };
  IterateOverAllCombinations(gradient_flags, ComparePhyloGradients);

  // Test likelihood "exclude" options.
  auto LikelihoodExcludeLogDeterminant = [&]() {
    auto inst = CreateNewInstance();
    StringVector flags = {
        PhyloFlags::loglikelihood_options_.run_defaults_.flag_,
        PhyloFlags::loglikelihood_options_.exclude_log_det_jacobian_likelihood_.flag_};
    double likelihood_exclude_log_det = inst.LogLikelihoods(flags)[0];
    double log_det_jacobian_likelihood =
        LogDetJacobianHeightTransform(inst.tree_collection_.trees_[0]);
    double gold_likelihood = gold_likelihoods[0];
    Assert(gold_likelihood == likelihood_exclude_log_det + log_det_jacobian_likelihood,
           "LogLikelihood not equal to excluded LogLikelihood plus "
           "LogDetJacobianHeightTransform.");
  };
  LikelihoodExcludeLogDeterminant();

  // Test gradient "exclude" options.
  auto GradientExcludeLogDeterminant = [&]() {
    auto inst = CreateNewInstance();
    StringVector flags = {
        PhyloFlags::gradient_options_.run_defaults_.flag_,
        PhyloFlags::gradient_options_.exclude_log_det_jacobian_gradient_.flag_};
    GradientMap grad_map = inst.PhyloGradients(flags, true)[0].gradient_;
    DoubleVector grad_exclude_log_det =
        grad_map[PhyloGradient::ratios_root_height_key_];
    DoubleVector log_det_jacobian_grad =
        GradientLogDeterminantJacobian(inst.tree_collection_.trees_[0]);
    GradientMap gold_grad_map = gold_gradients[0].gradient_;
    DoubleVector gold_grad = gold_grad_map[PhyloGradient::ratios_root_height_key_];
    for (size_t i = 0; i < gold_grad.size(); i++) {
      CHECK_MESSAGE(
          gold_grad[i] == grad_exclude_log_det[i] + log_det_jacobian_grad[i],
          "Gradient not equal to excluded Gradient plus GradientLogDetJacobian.");
    }
  };
  GradientExcludeLogDeterminant();

  // Test gradient "set" options.
  auto GradientSetDelta = [&]() {
    auto inst = CreateNewInstance();
    StringDoubleVector flags = {
        {PhyloFlags::gradient_options_.run_defaults_.flag_, 1.0},
        {PhyloFlags::gradient_options_.set_gradient_delta_.flag_, 1.0e1}};
    GradientMap grad_map = inst.PhyloGradients(flags, true)[0].gradient_;
    DoubleVector subst_grad =
        grad_map[PhyloFlags::gradient_mapkeys_.substitution_model_.flag_];
    // delta = 1.0e-6 (default)
    DoubleVector gold_subst_grad_1e6 = {49.0649, 151.831, 26.4022, -8.25114,
                                        75.2975, 352.565, 90.0701, 30.1228};
    // delta = 1.0e1
    DoubleVector gold_subst_grad_1e1 = {-73.2611, 25.4074,  -33.2865, -54.0479,
                                        47.9938,  -2696.06, -84.2954, 6.0563};
    for (size_t i = 0; i < gold_subst_grad_1e1.size(); i++) {
      CHECK_MESSAGE(
          std::abs(subst_grad[i] - gold_subst_grad_1e1[i]) < 0.01,
          "Delta value set by flag did not result in correct gradient values.");
    }
  };
  GradientSetDelta();

  // Test options that have children.
  auto GradientCompareParentChildren = [&]() {
    auto parent_inst = CreateNewInstance();
    StringVector parent_flags = {
        PhyloFlags::gradient_options_.substitution_model_.flag_};
    GradientMap parent_map =
        parent_inst.PhyloGradients(parent_flags, true)[0].gradient_;

    auto child_inst = CreateNewInstance();
    StringVector child_flags = {
        PhyloFlags::gradient_options_.substitution_model_rates_.flag_,
        PhyloFlags::gradient_options_.substitution_model_frequencies_.flag_};
    GradientMap child_map = child_inst.PhyloGradients(child_flags, true)[0].gradient_;

    CHECK_MESSAGE(parent_map == child_map, "Parent and Child-based Maps not equal.");
  };
  GradientCompareParentChildren();
}

#endif  // DOCTEST_LIBRARY_INCLUDED
