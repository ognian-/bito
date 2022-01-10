// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// TODO: Implement SortedVector, SortedUniqueVector, BidirectionalMap, ArgsortVector,
// BisortVector This defines various specialized containers:
// - SortedVector, SortedUniqueVector:
//    These two data structures maintain sorted vectors. There one difference is that
//    SortedVector allows for duplicate elements in vector, SortedUniqueVector does not. SortedVectors
//    are an alternative to std::set. std::set is useful when you have a large dataset
//    that you are frequently inserting into.  But vectors have superior cache
//    performance during access than the Red/Black Trees of std::set.

// TODO: Work in Progress.

#include "sugar.hpp"

#pragma once

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
template <class T>
class SortedVector {
 public:
  // ** Constructors:

  // Empty vector.
  SortedVector() : is_sorted_(true), data_() {};
  // Build sorted vector from normal vector.
  SortedVector(const std::vector<T> &data, const bool is_vector_sorted = false)
      : is_sorted_(false), data_(data) {
    Sort();
  };

  // ** Access:

  T at(size_t idx) { return data_.at(idx); }

  // ** Sorted methods:
  // These methods preserve the sort in the vector.
  // Insert in O(n), Delete in O(n).

  // Perfom a fresh sort on all values in the vector.
  void Sort() {
    std::sort(data_.begin(), data_.end());
    is_sorted_ = true;
  };
  // Checks whether the vector is currently in a sorted state.
  bool IsSorted() const { return is_sorted_; }

  // Insert value into sorted position in the vector.
  // Finds insertion point of value and makes gap by right shift.
  void Insert(const T value) { 
    data_.push_back(std::move(value));
    // Find last instance of value or first value greater than value.
    size_t pos = BinarySearchLast(value, true);
    // Right shift from position.
    for (size_t i = pos; i < data_.size(); i++) {
      data_[i+1] = std::move(data_[i]);
    }
    data_[pos] = std::move(value);
  };

  // Insert multiple values into vector. Puts off sorting until all new values inserted.  
  void Insert(const std::vector<T> values){
    for (size_t i = 0; i < values.size(); i++) {
      InsertWithoutSorting(values[i]);
    }
    Sort();
  };

  // Delete value from vector. Finds value to delete and closes gap with left shift.
  void Delete(const T value) {
    std::optional<size_t> idx = FindIf(value);
    Assert(idx.has_value(),
           "SortedVector::DeleteByValue(): Value does not exist in vector.");
    return Delete(*idx);
  };

  void Delete(const size_t idx){
    for (size_t i = idx + 1; i < data_.size(); i++) {
      data_[i-1] = std::move(data_[i]);
    }
    data_.pop_back();
  };

  // ** Nonsorting methods:
  // These methods do NOT preserve sort in the vector.

  // Insert value into vector. Does not maintain sort.
  void InsertWithoutSorting(T value) {
    data_.push_back(value);
    is_sorted_ = false;
  };

  // Delete value from vector. Does not maintain sort.
  void DeleteWithoutSorting(T value) {
    std::optional<size_t> idx = FindIf(value);
    Assert(idx.has_value(),
           "SortedVector::DeleteByValue(): Value does not exist in vector.");
    return DeleteWithoutSorting(*idx);
  };

  void DeleteWithoutSorting(size_t idx) {
    data_[idx] = data_[data_.size() - 1];
    data_.pop_back();
    is_sorted_ = false;
  };

  // ** Search:
  // Lookup by value. If sorted, performs binary search.
  // If unsorted, performs linear scan.

  // Returns position of value
  size_t Find(const T value) const {
    std::optional<size_t> pos = FindIf(value);
    Assert(pos.has_value(), "SortedVector::Find(): Value does not exist in vector.");
  };

  // Check if vector contains value.
  bool Contains(const T value) const { return (FindIf(value) != std::nullopt); };

  // Lookup by value if exists, otherwise return std::nullopt.
  std::optional<size_t> FindIf(T value) const {
    if (is_sorted_) {
      return BinarySearchFirst(value);
    }
    return LinearSearch(value);
  };

  // Get data vector
  std::vector<T> &GetData() {
    return data_;
  };

 protected:
  // Checks if vector is currently in a sorted state.
  bool IsValidSort() const {
    for (size_t i = 1; i < data_.size(); i++) {
      if (CompareFn(data_[i - 1], data_[i]) > 0) {
        return false;
      }
    }
    return true;
  }

  // Perform binary search to find position of first occuraance of the greatest element
  // that is less than or equal to target element.
  std::optional<size_t> BinarySearch(T value) const {
    size_t pos = data_.size() / 2;
    for (size_t i = data_.size() / 4; i > 1; i /= 2) {
      int cmp_val = CompareFn(value, data_[pos]);
      if (cmp_val > 0) {
        pos -= i;
      } else if (cmp_val < 0) {
        pos += i;
      }
      else {
        return pos;
      }
    }
    return std::nullopt;
  }

  // Perform binary search to find position of first occuraance of the greatest element
  // that is less than or equal to target element.
  std::optional<size_t> BinarySearchFirst(const T value, const bool accept_nonequal = true) const {
    size_t pos = data_.size() / 2;
    for (size_t i = data_.size() / 4; i > 1; i /= 2) {
      int cmp_val = CompareFn(value, data_[pos]);
      if (cmp_val >= 0) {
        pos -= i;
      } else if (cmp_val < 0) {
        pos += i;
      }
    }
    if (accept_nonequal || (CompareFn(value, data_[pos]) == 0)) {
      return pos;
    }
    return std::nullopt;
  }

  // Perform binary search to find position of last occurance of the least element that
  // is greater than or equal to target element.
  std::optional<size_t> BinarySearchLast(const T value, const bool accept_nonequal = true) const {
    size_t pos = data_.size() / 2;
    for (size_t i = data_.size() / 4; i > 1; i /= 2) {
      int cmp_val = CompareFn(value, data_[pos]);
      if (cmp_val > 0) {
        pos -= i;
      } else if (cmp_val <= 0) {
        pos += i;
      }
    }
    if (accept_nonequal || (CompareFn(value, data_[pos]) == 0)) {
      return pos;
    }
    return std::nullopt;
  }

  // Perform linear search to find value.
  std::optional<size_t> LinearSearch(T value) const {
    for (size_t i = 1; i < data_.size(); i++) {
      if (CompareFn(data_[i], value) == 0) {
        return i;
      }
    }
    return std::nullopt;
  }

 protected:
  // Sorted Data Array.
  std::vector<T> data_;
  // Some functions allow for the postponing of the sort.  This verifies vector is
  // currently sorted.
  bool is_sorted_;
};

// // ** Sorted Unique Vector
// // Vector does not allow for equal values. Inserting duplicate is a allowed, acts as
// a
// // null operation.
// template <typename T>
// class SortedUniqueVector : public SortedVector {
//  public:
//   // NOTE: If vector already contains a value, then it not inserted.
//   // Insert value into vector.
//   Insert(T value);
//   // Insert value into vector without maintaining sort.
//   void InsertNoSort(T value);
//   // Remove ith element from vector.
// };

// // TODO: Implement class that maintains a argsort vector/set on a data vector?
// // Sorts a positional index array [0,1,2,...] with respect to input data array.
// template <typename T>
// std::vector<size_t> ArgSort(const std::vector<T> &data_vector,
//                             const std::function<int(T, T)> compare) {
//   std::vector<size_t> argsort_vector(data_vector.size());
//   std::iota(data_vector.begin(), data_vector.end(), 0);
//   std::sort(data_vector.begin(), data_vector.end(),
//             // Sort indices of sorted_index according to their index in input_array.
//             [&argsort_vector, &data_vector](int left, int right) {
//               return Compare(data_vector[left], data_vector[right]);
//             });
//   return argsort_vector;
// };

// template <typename T>
// class ArgSortVector {};

// // This maintains a vector of pairs that are sorted by both pairs.
// template <typename T1, typename T2>
// class TwoSortVector {
//  private:
//   // Sorts arg vector by data vector, using the specified data field.
//   void ArgsortByField(std::vector<std::pair<T1, T2>> data,
//                       std::vector<size_t> argvector, size_t sort_by_field) {}

//   std::vector<std::pair<T1, T2>> data;
//   std::vector<size_t> sort_by_first;
//   std::vector<size_t> sort_by_second;
// };

#ifdef DOCTEST_LIBRARY_INCLUDED

TEST_CASE("Sorted Vectors") { 
  std::cout << "== SORTED_VECTORS" << std::endl;
  // Test with int vector.
  std::vector<int> int_vec = { -2, -2, 0, 0, 0, 1, 2, 3, 5, 7, 10, 11 };
  std::vector<int> int_vec_shuf = { -5, 5, 3, 2, 0, 4, 5, -1 };
  // std::cout << "Vector: " << int_vec_shuf << std::endl;
  // SortedVector int_sortvec = SortedVector(int_vec);
  // std::cout << "SortedVector: " << int_sortvec.GetData() << std::endl;

  // Test with unique pointers.

}

#endif  // DOCTEST_LIBRARY_INCLUDED
