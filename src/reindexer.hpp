// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Reindexers are argsorted index arrays that describe reordering of data due to
// modifications of data vectors, such as when adding NNI pairs to the SubsplitDAG.
// These reindexers can then be used to reorder associated data arrays whose indices
// correspond the ordering of SubsplitDAG data arrays, such as with the node or edge
// arrays.
//
// A reindexer is a SizeVector that can be thought of as a one-to-one function that maps
// from an old indexing scheme to a new indexing scheme. In other words, if old index
// `i` maps to new index `j`, then reindexer[`i`] = `j`. For example, if old_vector =
// [A, B, C] and reindexer = [1, 2, 0], then new_vector = [C, A, B]. Note that
// old_vector and reindexer must have the same size.

#ifndef SRC_REINDEX_HPP_
#define SRC_REINDEX_HPP_

#include <numeric>

#include "eigen_sugar.hpp"
#include "sugar.hpp"
#include "bitset.hpp"

using ReindexVector = SizeVector;

namespace Reindexer {
// These operations are for working with reindexers.

// For each position in a identity reindexer, reindexer[`i`] = `i`.
// E.g. for size = 5, reindexer = [0, 1, 2, 3, 4].
inline SizeVector IdentityReindexer(const size_t size) {
  SizeVector reindexer(size);
  std::iota(reindexer.begin(), reindexer.end(), 0);
  return reindexer;
}

// Check if a SizeVector is a valid reindexer (contains every index exactly once).
inline bool IsValidReindexer(const SizeVector &reindexer) {
  std::vector<bool> already_used(reindexer.size(), false);
  for (size_t idx = 0; idx < reindexer.size(); idx++) {
    if (reindexer[idx] >= reindexer.size() || already_used[reindexer[idx]]) {
      return false;
    }
    already_used[reindexer[idx]] = true;
  }
  return true;
}

// Reindexes the given vector, according to the reindexer.
// NOTE: Is it worth having the option to perform in-place?
template <typename VectorType>
inline VectorType Reindex(VectorType &old_vector, const SizeVector &reindexer) {
  Assert(IsValidReindexer(reindexer), "Reindexer must be valid in Reindexer::Reindex.");
  Assert(old_vector.size() == reindexer.size(),
         "The vector and reindexer must have the same size in Reindexer::Reindex.");
  VectorType new_vector(reindexer.size());
  for (size_t idx = 0; idx < old_vector.size(); idx++) {
    new_vector[reindexer[idx]] = std::move(old_vector[idx]);
  }
  return new_vector;
}

// Reindexes the given vector concatenated with some additional values.
template <typename VectorType>
inline VectorType Reindex(VectorType &old_vector, const SizeVector &reindexer,
                          VectorType &additional_values) {
  Assert(IsValidReindexer(reindexer), "Reindexer must be valid in Reindexer::Reindex.");
  Assert(old_vector.size() + additional_values.size() == reindexer.size(),
         "Size of the vector and additional values must add up to the reindexer size "
         "in Reindexer::Reindex.");
  VectorType new_vector(reindexer.size());
  for (size_t idx = 0; idx < old_vector.size(); idx++) {
    new_vector[reindexer[idx]] = std::move(old_vector[idx]);
  }
  for (size_t idx = 0; idx < additional_values.size(); idx++) {
    new_vector[reindexer[old_vector.size() + idx]] = std::move(additional_values[idx]);
  }
  return new_vector;
}

// TODO: Work in Progress
// Reindexes the given vector of given datatype according to the reindexer, done in place with bool checking.
template <typename VectorType, typename DataType>
inline void ReindexInPlace(VectorType &data_vector, const SizeVector &reindexer, std::optional<DataType> fill_value = std::nullopt) {
  size_t old_size = data_vector.size();
  size_t new_size = reindexer.size();
  Bitset is_reindexed(new_size);
  data_vector.resize(new_size);
  if (fill_value.has_value()) {
    for (size_t i = old_size; i < new_size; i++) {
      data_vector[i] = *fill_value;
    }
  }
  for (size_t i = 0; i < new_size; i++) {
    size_t old_index;
    size_t new_index = i;
    DataType old_swap_value = data_vector[new_index];
    DataType new_swap_value;
    while (is_reindexed[new_index] != true) {
      old_index = new_index;
      new_index = reindexer[old_index];
      // Push current value from old index to new index, and save overwritten value.
      new_swap_value = data_vector[new_index];
      data_vector[new_index] = old_swap_value;
      old_swap_value = new_swap_value;
      is_reindexed[old_index] = true;
    } 
  }
}

// Gets the inverse of a given reindexer.
inline SizeVector InvertReindexer(const SizeVector &reindexer) {
  Assert(IsValidReindexer(reindexer),
         "Reindexer must be valid in Reindexer::InvertReindexer.");
  SizeVector inverted_reindexer(reindexer.size());
  for (size_t idx = 0; idx < reindexer.size(); idx++) {
    inverted_reindexer[reindexer[idx]] = idx;
  }
  return inverted_reindexer;
}

// Remaps each of the ids in the vector according to the reindexer.
inline void RemapIdVector(SizeVector &vector, const SizeVector &reindexer) {
  Assert(IsValidReindexer(reindexer),
         "Reindexer must be valid in Reindexer::RemapIdVector.");
  for (size_t id : vector) {
    Assert(id < reindexer.size(),
           "The vector cannot contain an id out of bounds of the reindexer in "
           "Reindexer::RemapIdVector.");
  }
  std::transform(vector.begin(), vector.end(), vector.begin(),
                 [reindexer](size_t id) { return reindexer[id]; });
}

// In a given reindexer, take the old_id in the reindexer and reassign it to the new_id
// and shift over the ids strictly between old_id and new_id to ensure that the
// reindexer remains valid. For example if old_id = 1 and new_id = 4, this method would
// shift 1 -> 4, 4 -> 3, 3 -> 2, and 2 -> 1.
inline void ReassignAndShift(SizeVector &reindexer, const size_t old_id,
                             const size_t new_id) {
  Assert(old_id < reindexer.size() && new_id < reindexer.size(),
         "The given ids must be within the bounds of the reindexer in "
         "Reindexer::ReassignAndShift.");
  Assert(IsValidReindexer(reindexer),
         "Reindexer must be valid in Reindexer::ReassignAndShift.");
  if (old_id == new_id) {
    return;
  }
  // Find position with value old_id.
  const size_t old_id_position =
      std::find(reindexer.begin(), reindexer.end(), old_id) - reindexer.begin();
  // Shift.
  if (old_id > new_id) {
    for (size_t &id : reindexer) {
      if (id < old_id && id >= new_id) {
        id++;
      }
    }
  } else {
    for (size_t &id : reindexer) {
      if (id > old_id && id <= new_id) {
        id--;
      }
    }
  }
  // Reassign old_id to new_id.
  reindexer[old_id_position] = new_id;
}

}  // namespace Reindexer

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("IdentityReindexer") {
  // Check that IdentityReindexer returns correctly.
  SizeVector correct_default{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  CHECK_EQ(correct_default, Reindexer::IdentityReindexer(10));
}

TEST_CASE("IsValidReindexer") {
  // Index appears more than once.
  CHECK(!Reindexer::IsValidReindexer({1, 3, 0, 0}));
  // Missing an index.
  CHECK(!Reindexer::IsValidReindexer({1, 3, 4, 2}));
  // Valid reindexer.
  CHECK(Reindexer::IsValidReindexer({1, 3, 0, 2}));
}

TEST_CASE("Reindex") {
  // Check that Reindex throws if the given vector and reindexer have different sizes.
  SizeVector old_size_vector{7, 8, 9};
  SizeVector reindexer{2, 0, 3, 1};
  CHECK_THROWS(Reindexer::Reindex(old_size_vector, reindexer));
  // Check that Reindex returns correctly.
  reindexer = {2, 0, 1};
  SizeVector new_size_vector = Reindexer::Reindex(old_size_vector, reindexer);
  SizeVector correct_new_size_vector{8, 9, 7};
  CHECK_EQ(new_size_vector, correct_new_size_vector);
  // Check that Reindex also works with EigenVectorXd and additional values.
  EigenVectorXd old_eigen_vector(3);
  old_eigen_vector << 7, 8, 9;
  EigenVectorXd additional_values(2);
  additional_values << 10, 11;
  reindexer = {2, 4, 0, 3, 1};
  EigenVectorXd new_eigen_vector =
      Reindexer::Reindex(old_eigen_vector, reindexer, additional_values);
  EigenVectorXd correct_new_eigen_vector(5);
  correct_new_eigen_vector << 9, 11, 7, 10, 8;
  CHECK_EQ(new_eigen_vector, correct_new_eigen_vector);
}

TEST_CASE("InvertReindexer") {
  // Check that inverting a vector twice results in the original vector.
  SizeVector reindexer{1, 3, 0, 2};
  SizeVector correct_inverted_reindexer{2, 0, 3, 1};
  SizeVector inverted_reindexer = Reindexer::InvertReindexer(reindexer);
  CHECK_EQ(inverted_reindexer, correct_inverted_reindexer);
  SizeVector correct_reindexer{1, 3, 0, 2};
  reindexer = Reindexer::InvertReindexer(inverted_reindexer);
  CHECK_EQ(reindexer, correct_reindexer);
}

TEST_CASE("RemapIdVector") {
  // Check that Reindex throws if the given vector has an index out of bounds of the
  // reindexer.
  SizeVector size_vector{3, 5};
  SizeVector reindexer{2, 0, 3, 1};
  CHECK_THROWS(Reindexer::RemapIdVector(size_vector, reindexer));
  // Check that RemapIdVector returns correctly.
  size_vector = {3, 5};
  reindexer = {2, 0, 3, 1, 6, 4, 5};
  Reindexer::RemapIdVector(size_vector, reindexer);
  SizeVector correct_size_vector{1, 4};
  CHECK_EQ(size_vector, correct_size_vector);
}

TEST_CASE("ReassignAndShift") {
  // Check that ReassignAndShift returns correctly when old_id > new_id.
  SizeVector reindexer{0, 1, 2, 3, 4, 5, 6};
  Reindexer::ReassignAndShift(reindexer, 4, 1);
  SizeVector correct_reindexer{0, 2, 3, 4, 1, 5, 6};
  CHECK_EQ(reindexer, correct_reindexer);
  Reindexer::ReassignAndShift(reindexer, 5, 2);
  correct_reindexer = {0, 3, 4, 5, 1, 2, 6};
  CHECK_EQ(reindexer, correct_reindexer);
  Reindexer::ReassignAndShift(reindexer, 1, 3);
  correct_reindexer = {0, 2, 4, 5, 3, 1, 6};
  CHECK_EQ(reindexer, correct_reindexer);
  // Check that ReassignAndShift returns correctly when old_id = new_id.
  reindexer = {1, 0, 4, 6, 5, 3, 2};
  Reindexer::ReassignAndShift(reindexer, 4, 4);
  correct_reindexer = {1, 0, 4, 6, 5, 3, 2};
  CHECK_EQ(reindexer, correct_reindexer);
  // Check that ReassignAndShift returns correctly when old_id < new_id.
  reindexer = {6, 0, 4, 1, 5, 3, 2};
  Reindexer::ReassignAndShift(reindexer, 1, 5);
  correct_reindexer = {6, 0, 3, 5, 4, 2, 1};
  CHECK_EQ(reindexer, correct_reindexer);
}

#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_REINDEX_HPP_
