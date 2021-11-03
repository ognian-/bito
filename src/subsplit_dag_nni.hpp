// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Subsplit DAG NNI (Nearest Neighbor Interchange):
//
// An `NNIOperation` contains an output parent/child Subsplit pair which is the
// result of an NNI operation on an input parent/child pair. An NNI operation can be
// seen as a swapping of the branches in the SubsplitDAG, or alternatively as a a
// reordering of the set of clades in an input parent/child pair: the parent's sister
// clade, the child's left clade, and the child's right clade. For any given
// parent/child pair, there are two possible NNIs: swapping the sister clade with the
// left child clade, or swapping the sister clade with the right child clade.
//
// The `SetOfNNIs` is a set of `NNIOperations` used to account for all "adjacent" NNIs
// to a SubsplitDAG.  That is, all output parent/child pairs which can be generated from
// a single NNI operation on an input parent/child pair taken from the set of all the
// parent/child pairs currently in the SubsplitDAG, where the output parent/child pair
// is not also already in the SubpslitDAG.

#ifndef SRC_SUBSPLIT_DAG_NNI_HPP_
#define SRC_SUBSPLIT_DAG_NNI_HPP_

#include <numeric>
#include "bitset.hpp"
#include "sugar.hpp"

// Nearest Neighbor Interchange Operations
// NNIOperation stores output parent/child pair which are the product of an NNI.
class NNIOperation {
 public:
  NNIOperation(Bitset parent, Bitset child) : parent_(parent), child_(child){};

  // Comparator:
  // NNIOperations are ordered according to the std::bitset ordering of their parent
  // subsplit, then the std::bitset order their child subsplit.
  static int Compare(const NNIOperation &nni_a, const NNIOperation &nni_b);
  int Compare(const NNIOperation &nni_b) const;

  friend bool operator<(const NNIOperation &lhs, const NNIOperation &rhs);
  friend bool operator>(const NNIOperation &lhs, const NNIOperation &rhs);
  friend bool operator==(const NNIOperation &lhs, const NNIOperation &rhs);
  friend bool operator!=(const NNIOperation &lhs, const NNIOperation &rhs);

  // Special Constructors:
  // Produces the output NNIOperation from an input subsplit pair that results from an
  // NNI Swap according to `swap_which_child_clade_with_sister`.
  static NNIOperation NNIOperationFromNeighboringSubsplits(
      const Bitset parent_in, const Bitset child_in,
      const bool swap_which_child_clade_with_sister, const bool which_child_of_parent);
  // If it is not known whether the child is sorted/rotated, it can be inferred by
  // this overload.
  static NNIOperation NNIOperationFromNeighboringSubsplits(
      const Bitset parent_in, const Bitset child_in,
      const bool swap_which_child_clade_with_sister);

  Bitset parent_;
  Bitset child_;
};

// SetOfNNIs: 
// Contain all NNI output parent/child pairs which are "adjacent" to the
// current SubsplitDAG. That is, pairs which are the result of an NNI on an input bitset
// pair that are currently present in the SubsplitDAG.
class SetOfNNIs : protected std::set<NNIOperation> {
 public:
  friend bool operator==(const SetOfNNIs &lhs, const SetOfNNIs &rhs);
  friend bool operator!=(const SetOfNNIs &lhs, const SetOfNNIs &rhs);

  void Insert(NNIOperation nni_op) { insert(nni_op); };
  void Insert(Bitset parent, Bitset child) { Insert(NNIOperation(parent, child)); };
  void Erase(NNIOperation nni_op) { erase(nni_op); };
  void Erase(Bitset parent, Bitset child) { Erase(NNIOperation(parent, child)); };

  void Clear() { return clear(); };
  size_t GetSize() const { return size(); };

  std::set<NNIOperation>::iterator Begin() const { return begin(); };
  std::set<NNIOperation>::iterator End() const { return end(); };
};

// RankedSetOfNNIs: 
// Bi-directional Map between an NNI and a associated score for ranking the quality of NNI.
class RankedSetOfNNIs {
 public:
  void Insert(double, NNIOperation);
  void Insert(std::pair<double, NNIOperation>);
  void Remove(double, NNIOperation);
  void Remove(std::pair<double, NNIOperation>);
  bool Empty() const { return score_to_nni_.empty(); };

  double GetMaxScore() const { return score_to_nni_.rbegin()->first; };
  NNIOperation GetMaxNNI() const { return score_to_nni_.rbegin()->second; };

  std::set<std::pair<double, NNIOperation>> score_to_nni_;
  std::set<std::pair<double, NNIOperation>> nni_to_score_;
};

// Sorts a positional index array [0,1,2,3,...] with respect to input data array.
template<typename T>
std::vector<size_t> ArgSort(const std::vector<T> &data) {
  std::vector<size_t> sorted_index(data.size());
  std::iota(data.begin(), data.end(), 0);
  std::sort(data.begin(), data.end(),
    // Sort indices of sorted_index according to their index in input_array.
    [&sorted_index, &data](int left, int right) {
      return data[left] < data[right];
  });
  return sorted_index;
};

// // Maintains an sorted index vector with respect to a reference data array.
// template<typename T, std::function<int(T,T)>
// class ArgsortVector {
//  public:

//  private:
//   std::vector<size_t> argsort_;
//   std::vector<T> *data_;
// };

#ifdef DOCTEST_LIBRARY_INCLUDED

// See tree diagram at:
// https://user-images.githubusercontent.com/31897211/136849710-de0dcbe3-dc2b-42b7-b3de-dd9b1a60aaf4.gif
TEST_CASE("NNIOperation") {
  // Clades for NNI.
  Bitset X("100");
  Bitset Y("010");
  Bitset Z("001");
  // Initial Child and Parent.
  Bitset parent_in = Bitset::Subsplit(X, Y | Z);
  Bitset child_in = Bitset::Subsplit(Y, Z);
  // Correct Solutions.
  Bitset correct_parent_xy = Bitset::Subsplit(Y, X | Z);
  Bitset correct_child_xy = Bitset::Subsplit(X, Z);
  NNIOperation correct_nni_xy = NNIOperation(correct_parent_xy, correct_child_xy);
  Bitset correct_parent_xz = Bitset::Subsplit(Z, Y | X);
  Bitset correct_child_xz = Bitset::Subsplit(Y, X);
  NNIOperation correct_nni_xz = NNIOperation(correct_parent_xz, correct_child_xz);

  // Swap X and Y
  auto nni_xy =
      NNIOperation::NNIOperationFromNeighboringSubsplits(parent_in, child_in, 0);
  CHECK_EQ(correct_nni_xy, nni_xy);
  // Swap X and Z
  auto nni_xz =
      NNIOperation::NNIOperationFromNeighboringSubsplits(parent_in, child_in, 1);
  CHECK_EQ(correct_nni_xz, nni_xz);

  // Relationship is known (child_in is the rotated clade of parent_in)
  auto nni_xy_2 =
      NNIOperation::NNIOperationFromNeighboringSubsplits(parent_in, child_in, 0, 1);
  CHECK_EQ(correct_nni_xy, nni_xy_2);
  CHECK_THROWS(
      NNIOperation::NNIOperationFromNeighboringSubsplits(parent_in, child_in, 0, 0));
};

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_SUBSPLIT_DAG_NNI_HPP_
