#include "unrooted_sbn_instance.hpp"
#include "rooted_sbn_instance.hpp"
#include "subsplit_dag.hpp"
#include "topology_sampler.hpp"

// This is just a place to muck around, and check out performance.

auto now = std::chrono::high_resolution_clock::now;

// To valgrind (you can pip install gprof2dot):
// valgrind --tool=callgrind ./_build/noodle
// gprof2dot -f callgrind callgrind.out.16763 | dot -Tpng -o ~/output.png

constexpr size_t out_of_sample_index = 99999999;

void testDAGSampling() {
  RootedSBNInstance sbninst("charlie");
  sbninst.ReadNewickFile("data/five_taxon_rooted.nwk");
  sbninst.ProcessLoadedTrees();
  sbninst.TrainSimpleAverage();

  Driver driver;
  auto tree_collection =
      RootedTreeCollection::OfTreeCollection(driver.ParseNewickFile("data/five_taxon_rooted.nwk"));
  SubsplitDAG dag(tree_collection);
  EigenVectorXd normalized_sbn_parameters = dag.BuildUniformOnTopologicalSupportPrior();

  // Count the frequencies of rooted trees in a file.
  size_t rooted_tree_count_from_file = 0;
  RootedIndexerRepresentationSizeDict counter_from_file(0);
  for (const auto &indexer_representation : sbninst.MakeIndexerRepresentations()) {
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(counter_from_file,
                                                                indexer_representation);
    rooted_tree_count_from_file += indexer_representation.size();
  }
  // Count the frequencies of trees when we sample after training with
  // SimpleAverage.
  size_t sampled_tree_count = 1'000'000;
  RootedIndexerRepresentationSizeDict counter_from_sampling(0);
  ProgressBar progress_bar(sampled_tree_count / 1000);
  TopologySampler sampler;
  SubsplitDAGSamplerInput input{dag, normalized_sbn_parameters};
  for (size_t sample_idx = 0; sample_idx < sampled_tree_count; ++sample_idx) {
    const auto rooted_topology = sampler.SampleTopology(input, true);
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(
        counter_from_sampling,
        RootedSBNMaps::IndexerRepresentationOf(dag.BuildGPCSPIndexer(),
                                               rooted_topology, out_of_sample_index));
    if (sample_idx % 1000 == 0) {
      ++progress_bar;
      progress_bar.display();
    }
  }
  // These should be equal in the limit when we're training with SA.
  for (const auto &[key, _] : counter_from_file) {
    double observed =
        static_cast<double>(counter_from_sampling.at(key)) / sampled_tree_count;
    double expected =
        static_cast<double>(counter_from_file.at(key)) / rooted_tree_count_from_file;
    std::cout << "Observed: " << observed << "  Expected: " << expected << "\n";
    Assert(fabs(observed - expected) <= 5e-3, "Failed equality");
  }
  progress_bar.done();
}

[[maybe_unused]]
void testSBNInstanceSampling() {
  UnrootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon_unrooted.nwk");
  inst.ProcessLoadedTrees();
  inst.TrainSimpleAverage();
  // Count the frequencies of rooted trees in a file.
  size_t rooted_tree_count_from_file = 0;
  RootedIndexerRepresentationSizeDict counter_from_file(0);
  for (const auto &indexer_representation : inst.MakeIndexerRepresentations()) {
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(counter_from_file,
                                                                indexer_representation);
    rooted_tree_count_from_file += indexer_representation.size();
  }
  // Count the frequencies of trees when we sample after training with
  // SimpleAverage.
  size_t sampled_tree_count = 1'000'000;
  RootedIndexerRepresentationSizeDict counter_from_sampling(0);
  ProgressBar progress_bar(sampled_tree_count / 1000);
  TopologySampler sampler;
  UnrootedSBNInstanceSamplerInput input{inst};
  for (size_t sample_idx = 0; sample_idx < sampled_tree_count; ++sample_idx) {
    const auto rooted_topology = sampler.SampleTopology(input, true);
    RootedSBNMaps::IncrementRootedIndexerRepresentationSizeDict(
        counter_from_sampling,
        RootedSBNMaps::IndexerRepresentationOf(inst.SBNSupport().Indexer(),
                                               rooted_topology, out_of_sample_index));
    if (sample_idx % 1000 == 0) {
      ++progress_bar;
      progress_bar.display();
    }
  }
  // These should be equal in the limit when we're training with SA.
  for (const auto &[key, _] : counter_from_file) {
    std::ignore = _;
    double observed =
        static_cast<double>(counter_from_sampling.at(key)) / sampled_tree_count;
    double expected =
        static_cast<double>(counter_from_file.at(key)) / rooted_tree_count_from_file;
    Assert(fabs(observed - expected) < 5e-3, "Sampling failed");
  }
  progress_bar.done();
}

int main() {

  testDAGSampling();
  //testSBNInstanceSampling();
  return 0;

  uint32_t leaf_count = 10000;

  Node::NodePtr topology = Node::Ladder(leaf_count);

  std::vector<size_t> ids;
  ids.reserve(1 + 2 * leaf_count);

  auto t_start = now();
  for (int i = 0; i < 100; i++) {
    ids.clear();
    topology->Preorder([&ids](const Node* node) { ids.push_back(node->Id()); });
  }

  std::chrono::duration<double> duration = now() - t_start;
  std::cout << "time: " << duration.count() << " seconds\n";
}
