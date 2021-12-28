// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// TODO: Implement SortedVector, SortedUniqueVector, BidirectionalMap, ArgsortVector,
// BisortVector This defines various specialized containers:
// - SortedVector, SortedUniqueVector:
//    These two data structures maintain sorted vectors. There one difference is that
//    SortedVector allows for duplicate data, SortedUniqueVector does not. SortedVectors
//    are an alternative to std::set. std::set is useful when you have a large dataset
//    that you are frequently inserting into.  But vectors have superior cache
//    performance when during access than the Red/Black Trees of std::set.

#include "sugar.hpp"

#ifndef SRC_CONTAINERS_HPP_
#define SRC_CONTAINERS_HPP_

// TODO: Implementation
// ** Comparator for generic vectors
// This generally compares any two vectors containing a comparable datatype.
// Uses a lexicographic-like sorting scheme: looks at each element from beginning to end
// in the vector.  The first element where its ith element is less than the other ith
// element is the lessor vector.  If they are equal to the end of one vectors, the
// shorter vector is lessor. template<class T, class Compare> int Compare(const
// std::vector<T> lhs, const std::vector<T> rhs);

// ** Sorted Vector
// Preserves sort for all operations except where explicitly indicated otherwise.
// O(n) insertion, O(n) deletion, O(log n) find.
template <class T, class Compare>
class SortedVector {
 public:
  // Constructor.
  SortedVector() : data_(){};
  SortedVector(std::vector &data) : data_(data){};

  // ** Sorted methods:
  // These methods preserve the sort in the vector.

  // Perfom a fresh sort on all values in the vector.
  void Sort();
  // Insert value into vector. Finds insertion point of value and makes gap by right
  // shift.
  void Insert(T value);
  // Insert multiple values into vector.
  void Insert(std::vector<T> values);
  // Delete value from vector. Finds value to delete and closes gap with left shift.
  void Delete(const T value);
  // Lookup by value.
  size_t Find(const T value) const;
  // Check if vector contains value.
  bool Contains(const T value) const;
  // Lookup by value if exists, otherwise return std::nullopt.
  std::optional<size_t> FindIf(T value) const;

  // ** Nonsorting methods:
  // These methods do NOT preserve sort in the vector.

  // Insert value into vector. Does not maintain sort.
  void InsertWithoutSorting(T value);
  // Delete value from vector. Does not maintain sort.
  void DeleteWithoutSorting(T value);

 protected:
  // Insert value into given index. Does not maintain sort.
  void InsertByIndex(size_t index, T value);
  // Delete data at given index. Does not maintain sort.
  void DeleteByIndex(size_t index);

 protected:
  // Sorted Data Array.
  std::vector data_;
  // Some functions allow for the postponing of the sort.  This verifies vector is
  // currently sorted.
  bool is_sorted_;
};

// ** Sorted Unique Vector
// Vector does not allow for equal values. Inserting duplicate is a allowed, acts as a
// null operation.
template <typename T>
class SortedUniqueVector : public SortedVector {
 public:
  // NOTE: If vector already contains a value, then it not inserted.
  // Insert value into vector.
  Insert(T value);
  // Insert value into vector without maintaining sort.
  void InsertNoSort(T value);
  // Remove ith element from vector.
  void Remove(size_t i);
  // Get ith element from vector.
  T Get(size_t i)
};

// TODO: Implement class that maintains a argsort vector/set on a data vector?
// Sorts a positional index array [0,1,2,...] with respect to input data array.
template <typename T>
std::vector<size_t> ArgSort(const std::vector<T> &data_vector,
                            const std::function<int(T, T)> compare) {
  std::vector<size_t> argsort_vector(data_vector.size());
  std::iota(data_vector.begin(), data_vector.end(), 0);
  std::sort(data_vector.begin(), data_vector.end(),
            // Sort indices of sorted_index according to their index in input_array.
            [&argsort_vector, &data_vector](int left, int right) {
              return Compare(data_vector[left], data_vector[right]);
            });
  return argsort_vector;
};

template <typename T>
class ArgSortVector {};

// This maintains a vector of pairs that are sorted by both pairs.
template <typename T1, typename T2>
class TwoSortVector {
 private:
  // Sorts arg vector by data vector, using the specified data field.
  void ArgsortByField(std::vector<std::pair<T1, T2>> data,
                      std::vector<size_t> argvector, size_t sort_by_field) {}

  std::vector<std::pair<T1, T2>> data;
  std::vector<size_t> sort_by_first;
  std::vector<size_t> sort_by_second;
};

#endif  // SRC_CONTAINERS_HPP_
