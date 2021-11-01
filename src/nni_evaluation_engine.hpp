// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// NNI Evaluation Engine
//

#ifndef SRC_NNI_EVALUATION_ENGINE_HPP_
#define SRC_NNI_EVALUATION_ENGINE_HPP_

#include "bitset.hpp"
#include "subsplit_dag.hpp"
#include "subsplit_dag_nni.hpp"
#include "sugar.hpp"

class NNIEvaluationEngine {
 public:
  NNIEvaluationEngine(SubsplitDAG &dag, SetOfNNIs &set_of_nnis)
      : dag_(dag), set_of_nnis_(set_of_nnis){};

  // Maintainence Methods: These maintain SetOfNNIs to stay consistent with the state of
  // associated DAG.
  //
  // Freshly synchonizes SetOfNNIs to match the current state of its DAG. Wipes old NNI
  // data and finds all all parent/child pairs adjacent to DAG by iterating over all
  // internal edges in DAG.
  void SyncSetOfNNIsWithDAG(SetOfNNIs &set_of_nnis, const SubsplitDAG &dag);
  // Updates NNI Set after given parent/child node pair have been added to the DAG.
  // Removes pair from NNI Set and adds adjacent pairs coming from newly created edges.
  void UpdateSetOfNNIsAfterDAGAddNodePair(SetOfNNIs &set_of_nnis,
                                          const SubsplitDAG &dag,
                                          const Bitset &parent_bitset,
                                          const Bitset &child_bitset);
  // Maintainence Helper Methods:
  //
  // Adds all NNIs from all (node_id, other_id) pairs, where other_id's are elements of
  // the node_id_vector.
  void AddAllNNIsFromNodeVectorToSetOfNNIs(SetOfNNIs &set_of_nnis,
                                           const SubsplitDAG &dag,
                                           const size_t &node_id,
                                           const SizeVector &adjacent_node_ids,
                                           const bool is_edge_rotated,
                                           const bool is_edge_leafward);
  // Based on given input NNIOperation, produces the two possible output NNIOperations
  // and adds those results to the NNI Set (if results are not a member of the DAG).
  void SafeAddOutputNNIsToSetOfNNIs(SetOfNNIs &set_of_nnis, const SubsplitDAG &dag,
                                    const Bitset &parent_bitset,
                                    const Bitset &child_bitset, const bool rotated);

 private:
  std::reference_wrapper<SubsplitDAG> dag_;
  std::reference_wrapper<SetOfNNIs> set_of_nnis_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_NNI_EVALUATION_ENGINE_HPP_
