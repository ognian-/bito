// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// NNI Evaluation Engine
//

#ifndef SRC_NNI_EVALUATION_ENGINE_HPP_
#define SRC_NNI_EVALUATION_ENGINE_HPP_

#include "bitset.hpp"
#include "gp_dag.hpp"
#include "subsplit_dag.hpp"
#include "subsplit_dag_nni.hpp"
#include "sugar.hpp"

class NNIEvaluationEngine {
 public:
  // Constructors:
  // Use externally sourced DAG, internally sourced NNI.
  NNIEvaluationEngine(GPDAG &dag_src);

  // ** Getter/Setters:

  //
  GPDAG *GetReferenceGPDAG() { return dag_; };
  //
  SetOfNNIs &GetAdjacentNNIs() { return adjacent_nnis_; };

  // ** Runner Methods:
  // These start the engine, which procedurally ranks and adds (and
  // maybe removes) NNIs to the DAG, until some termination criteria has been satisfied.

  // Primary Runner Selector.
  void Runner();

  // ** Evaluation/Ranking Methods:
  // These evaluate all NNIs from the SetOfNNIs by a criterion in order to place
  // them in a relative ordering.

  // Primary Full Evaluation Selector.
  void Evaluate();
  // Naive Implementation of Evaluation. Adds NNI to DAG, computes likelihood, then
  // removes NNI from DAG.
  void NaiveEvaluate();
  //
  // Individual NNI Evaluation.
  double EvaluateNNI(NNIOperation &proposed_nni);
  // Naive Implementation of NNI Evaluation.
  // Simply adds node pair to DAG, computes NNI score, then removes NNI.
  double NaiveEvaluateNNI(NNIOperation &proposed_nni);

  // ** Maintainence Methods:
  // These maintain SetOfNNIs to stay consistent with the state of associated DAG.

  // Freshly synchonizes SetOfNNIs to match the current state of its DAG. Wipes old NNI
  // data and finds all all parent/child pairs adjacent to DAG by iterating over all
  // internal edges in the DAG. (For each internal edges, two NNI's are possible.)
  void SyncSetOfNNIsWithDAG();
  // Updates NNI Set after given parent/child node pair have been added to the DAG.
  // Removes pair from NNI Set and adds adjacent pairs coming from newly created edges.
  void UpdateSetOfNNIsAfterDAGAddNodePair(const Bitset &parent_bitset,
                                          const Bitset &child_bitset);
  // Adds all NNIs from all (node_id, other_id) pairs, where other_id's are elements of
  // the adjacent_node_ids vector. is_edge_leafward tells whether node_id is the child
  // or parent. is_edge_rotated determines which side of parent the child descends from.
  void AddAllNNIsFromNodeVectorToSetOfNNIs(const size_t &node_id,
                                           const SizeVector &adjacent_node_ids,
                                           const bool is_edge_rotated,
                                           const bool is_edge_leafward);
  // Based on given input NNIOperation, produces the two possible output NNIOperations
  // and adds those results to the NNI Set (if results are not a member of the DAG).
  void SafeAddOutputNNIsToSetOfNNIs(const Bitset &parent_bitset,
                                    const Bitset &child_bitset, const bool rotated);

 private:
  GPDAG dag_src_;
  GPDAG *dag_;

  SetOfNNIs adjacent_nnis_;
  GPDAG nni_partial_dag_;
  RankedSetOfNNIs ranked_nnis_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_NNI_EVALUATION_ENGINE_HPP_
