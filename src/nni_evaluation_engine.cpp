// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//

#include "nni_evaluation_engine.hpp"

#include "subsplit_dag.hpp"
#include "subsplit_dag_nni.hpp"

// ** NNIEvaluationEngine Methods

NNIEvaluationEngine::NNIEvaluationEngine(GPDAG &dag_src) {
  dag_ = &dag_src;

  adjacent_nnis_ = SetOfNNIs();
  nni_partial_dag_ = GPDAG();
  ranked_nnis_ = RankedSetOfNNIs();
}

// ** Runner Methods

// TODO:
void NNIEvaluationEngine::Runner() {
  bool termination_criteria = true;

  // (0a) Initialize NNI Set based on state of DAG.
  SyncSetOfNNIsWithDAG();
  // (0b) Evaluate likelihood of all possible NNIs.
  // for (auto &nni : adjacent_nnis_) {
  //   ranked_nnis_.insert(nni, EvaluateNNI(nni));
  // }

  while (termination_criteria) {
    // If there are no remaining NNIs, then NNI Evaluation is over.
    if (ranked_nnis_.Empty()) {
      break;
    }
    // (1) Add the best NNI from adjacent NNIs to the DAG.
    // NNIOperation best_nni = ranked_nnis_.GetMaxNNI();
    // dag_->AddNodePair(best_nni.parent_, best_nni.child_);
    // (2) Update set of all NNIs to reflect added NNI.
    // UpdateSetOfNNIsAfterDAGAddNodePair(best_nni.parent_, best_nni.child_);
    // (3) Evaluate all newly added NNIs to the set.
    // for (auto &nni : adjacent_nnis_) {
    //   EvaluateNNI(nni);
    // }
    termination_criteria = false;
  }
}

// ** Evaluate Methods

void Evaluate() {}

void NaiveEvaluate() {}

// ** Individual NNI Evaluation Methods

double NNIEvaluationEngine::EvaluateNNI(NNIOperation &proposed_nni) {
  return NaiveEvaluateNNI(proposed_nni);
}

double NNIEvaluationEngine::NaiveEvaluateNNI(NNIOperation &proposed_nni) {
  dag_->AddNodePair(proposed_nni.parent_, proposed_nni.child_);
  // dag_->QuartetHybridRequestOf(proposed_nni.parent_, proposed_nni.child_);
  // dag_->RemoveNodePair(proposed_nni.parent_, proposed_nni.child_);

  return 0.0f;
}

double NNIEvaluationEngine::GraftEvaluateNNI(NNIOperation &proposed_nni) {}

// ** Maintainence Methods

void NNIEvaluationEngine::SyncSetOfNNIsWithDAG() {
  adjacent_nnis_.Clear();
  // Only real node pairs are viable NNIs.
  dag_->IterateOverRealNodes([this](const SubsplitDAGNode *node) {
    dag_->IterateOverParentAndChildAndLeafwardEdges(
        node, [this](const size_t parent_id, const bool is_rotated,
                     const size_t child_id, const size_t edge_idx) {
          // Only internal node pairs are viable NNIs.
          Bitset parent_bitset = dag_->GetDAGNode(parent_id)->GetBitset();
          Bitset child_bitset = dag_->GetDAGNode(child_id)->GetBitset();
          if (!(parent_bitset.SubsplitIsRoot() || child_bitset.SubsplitIsLeaf())) {
            SafeAddOutputNNIsToSetOfNNIs(parent_bitset, child_bitset, is_rotated);
          }
        });
  });
}

void NNIEvaluationEngine::UpdateSetOfNNIsAfterDAGAddNodePair(
    const Bitset &parent_bitset, const Bitset &child_bitset) {
  size_t parent_id = dag_->GetNodeId(parent_bitset);
  size_t child_id = dag_->GetNodeId(child_bitset);
  // Every new edge added is a potential new NNI.
  // Iterate over the parent and child node of the new pair.
  for (const size_t &node_id : {parent_id, child_id}) {
    // Get nodes adjacent to current node from both sorted and rotated edges.
    for (const bool is_edge_leafward : {true, false}) {
      // Get nodes adjacent to current node from both leafward and rootward directions.
      for (const bool is_edge_rotated : {true, false}) {
        SizeVector adjacent_node_ids = dag_->GetDAGNode(node_id)->GetLeafwardOrRootward(
            is_edge_leafward, is_edge_rotated);
        AddAllNNIsFromNodeVectorToSetOfNNIs(node_id, adjacent_node_ids, is_edge_rotated,
                                            is_edge_leafward);
      }
    }
  }
  // Remove the pair that was just added to the DAG from NNI Set.
  NNIOperation new_nni = NNIOperation(parent_bitset, child_bitset);
  adjacent_nnis_.Erase(new_nni);
}

void NNIEvaluationEngine::AddAllNNIsFromNodeVectorToSetOfNNIs(
    const size_t &node_id, const SizeVector &adjacent_node_ids,
    const bool is_edge_rotated, const bool is_edge_leafward) {
  Bitset node_bitset = dag_->GetDAGNode(node_id)->GetBitset();
  // Determine whether node_id corresponds to parent or child of the pair.
  // Add every edge's NNI to NNI Set.
  // If edges are leafward, node_id is the parent to all vector nodes.
  // If edges are rootward, node_id is the child to all vector nodes.
  if (is_edge_leafward) {
    const Bitset &parent_bitset = node_bitset;
    for (const auto &adjacent_node_id : adjacent_node_ids) {
      const Bitset child_bitset = dag_->GetDAGNode(adjacent_node_id)->GetBitset();
      SafeAddOutputNNIsToSetOfNNIs(parent_bitset, child_bitset, is_edge_rotated);
    }
  } else {
    const Bitset &child_bitset = node_bitset;
    for (const auto &adjacent_node_id : adjacent_node_ids) {
      const Bitset parent_bitset = dag_->GetDAGNode(adjacent_node_id)->GetBitset();
      SafeAddOutputNNIsToSetOfNNIs(parent_bitset, child_bitset, is_edge_rotated);
    }
  }
}

void NNIEvaluationEngine::SafeAddOutputNNIsToSetOfNNIs(const Bitset &parent_bitset,
                                                       const Bitset &child_bitset,
                                                       const bool is_edge_rotated) {
  // Soft assert that parent is not the root and child is not a leaf.
  if (parent_bitset.SubsplitIsRoot() || child_bitset.SubsplitIsLeaf()) {
    return;
  }
  // Input pair is in the DAG, so remove it from the Set if it exists.
  adjacent_nnis_.Erase(parent_bitset, child_bitset);
  // Add NNI for sorted clade swap and rotated clade swap.
  for (bool is_swap_with_sorted_child : {true, false}) {
    bool is_in_dag = false;
    const auto new_nni = NNIOperation::NNIOperationFromNeighboringSubsplits(
        parent_bitset, child_bitset, is_swap_with_sorted_child, !is_edge_rotated);
    // If DAG already contains output parent and child nodes, and an edge between them,
    // then don't add it to the adjacent_nnis.
    if (dag_->ContainsNode(new_nni.parent_) && dag_->ContainsNode(new_nni.child_)) {
      const size_t parent_id = dag_->GetNodeId(new_nni.parent_);
      const size_t child_id = dag_->GetNodeId(new_nni.child_);
      is_in_dag = dag_->ContainsEdge(parent_id, child_id);
    }
    if (is_in_dag == false) {
      adjacent_nnis_.Insert(new_nni);
    }
  }
}
