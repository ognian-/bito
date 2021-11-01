// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag_nni.hpp"

#include "bitset.hpp"
#include "subsplit_dag.hpp"

// ** NNIOperation Methods

int NNIOperation::Compare(const NNIOperation &nni_a, const NNIOperation &nni_b) {
  auto compare_parent = Bitset::Compare(nni_a.parent_, nni_b.parent_);
  if (compare_parent != 0) {
    return compare_parent;
  }
  auto compare_child = Bitset::Compare(nni_a.child_, nni_b.child_);
  return compare_child;
};

int NNIOperation::Compare(const NNIOperation &nni_b) const {
  const NNIOperation &nni_a = *this;
  return Compare(nni_a, nni_b);
}

bool operator<(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) < 0;
}
bool operator>(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) > 0;
}
bool operator==(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) == 0;
}
bool operator!=(const NNIOperation &lhs, const NNIOperation &rhs) {
  return NNIOperation::Compare(lhs, rhs) != 0;
}

NNIOperation NNIOperation::NNIOperationFromNeighboringSubsplits(
    const Bitset parent_in, const Bitset child_in,
    const bool swap_which_child_clade_with_sister, const bool which_child_of_parent) {
  // Input: Parent(X,YZ) -> Child(Y,Z).
  Bitset X = parent_in.SubsplitGetClade(!which_child_of_parent);
  // "Y" clade can be chosen arbitrarily from (Y,Z), so "Y" is chosen based on which
  // we want to swap with "X".
  Bitset Y = child_in.SubsplitGetClade(swap_which_child_clade_with_sister);
  Bitset Z = child_in.SubsplitGetClade(!swap_which_child_clade_with_sister);
  // Output: Parent(Y,XZ) -> Child(X,Z).
  Bitset parent_out = Bitset::Subsplit(Y, X | Z);
  Bitset child_out = Bitset::Subsplit(X, Z);
  return NNIOperation(parent_out, child_out);
}

NNIOperation NNIOperation::NNIOperationFromNeighboringSubsplits(
    const Bitset parent_in, const Bitset child_in,
    const bool swap_which_child_clade_with_sister) {
  bool which_clade_of_parent = Bitset::SubsplitIsWhichChildOf(parent_in, child_in);
  return NNIOperationFromNeighboringSubsplits(
      parent_in, child_in, swap_which_child_clade_with_sister, which_clade_of_parent);
}

// ** SetOfNNIs Methods

bool operator==(const SetOfNNIs &lhs, const SetOfNNIs &rhs) {
  return lhs.set_ == rhs.set_;
}

bool operator!=(const SetOfNNIs &lhs, const SetOfNNIs &rhs) {
  return lhs.set_ != rhs.set_;
}

void SetOfNNIs::Insert(NNIOperation nni_op) { set_.insert(nni_op); };
void SetOfNNIs::Insert(Bitset parent, Bitset child) {
  Insert(NNIOperation(parent, child));
}

void SetOfNNIs::Erase(NNIOperation nni_op) { set_.erase(nni_op); };
void SetOfNNIs::Erase(Bitset parent, Bitset child) {
  Erase(NNIOperation(parent, child));
}

void SetOfNNIs::Clear() { set_.clear(); }

size_t SetOfNNIs::GetSize() const { return set_.size(); }
