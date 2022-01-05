// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// Note that the default constructor for std::vector<bool> is filled with false:
// https://stackoverflow.com/a/22984114/467327

#include "bitset.hpp"

#include <utility>

#include "sugar.hpp"

Bitset::Bitset(std::vector<bool> value) : value_(std::move(value)) {}

Bitset::Bitset(const size_t n, const bool initial_value) : value_(n, initial_value) {}

Bitset::Bitset(const std::string str) : Bitset(str.length()) {
  for (size_t i = 0; i < value_.size(); i++) {
    if (str[i] == '0') {
      value_[i] = false;
    } else if (str[i] == '1') {
      value_[i] = true;
    } else {
      Failwith("String constructor for Bitset must use only 0s or 1s; found '" +
               std::string(1, str[i]) + "'.");
    }
  }
}

Bitset::Bitset(const SizeVector bits_on, const size_t n) : Bitset(n, false) {
  for (auto i : bits_on) {
    Assert(i < n, "Bitset SizeVector constructor has values out of range.");
    value_[i] = true;
  }
}

// ** std::bitset Interface Methods

bool Bitset::operator[](size_t i) const { return value_[i]; }

size_t Bitset::size() const { return value_.size(); }

void Bitset::set(size_t i, bool value) {
  Assert(i < value_.size(), "i out of range in Bitset::set.");
  value_[i] = value;
}

void Bitset::reset(size_t i) {
  Assert(i < value_.size(), "i out of range in Bitset::reset.");
  value_[i] = false;
}

void Bitset::flip() { value_.flip(); }

int Bitset::Compare(const Bitset& bitset_a, const Bitset& bitset_b) {
  Assert(bitset_a.size() == bitset_b.size(),
         "Bitsets must be same size for Bitset::Compare.");
  for (size_t i = 0; i < bitset_a.size(); i++) {
    if (bitset_a[i] != bitset_b[i]) {
      return static_cast<int>(bitset_a[i]) - static_cast<int>(bitset_b[i]);
    }
  }
  return 0;
}

int Bitset::Compare(const Bitset& that_bitset) const {
  const Bitset& this_bitset = *this;
  return Compare(this_bitset, that_bitset);
}

bool Bitset::operator==(const Bitset& other) const { return value_ == other.value_; }
bool Bitset::operator!=(const Bitset& other) const { return value_ != other.value_; }
bool Bitset::operator<(const Bitset& other) const { return value_ < other.value_; }
bool Bitset::operator<=(const Bitset& other) const { return value_ <= other.value_; }
bool Bitset::operator>(const Bitset& other) const { return value_ > other.value_; }
bool Bitset::operator>=(const Bitset& other) const { return value_ >= other.value_; }

Bitset Bitset::operator&(const Bitset& other) const {
  Assert(value_.size() == other.size(), "Size mismatch in Bitset::operator&.");
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] && other.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

Bitset Bitset::operator|(const Bitset& other) const {
  Assert(value_.size() == other.size(), "Size mismatch in Bitset::operator|.");
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] || other.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

Bitset Bitset::operator^(const Bitset& other) const {
  Assert(value_.size() == other.size(), "Size mismatch in Bitset::operator^.");
  Bitset r(value_.size());
  for (size_t i = 0; i < value_.size(); i++) {
    if (value_[i] != other.value_[i]) {
      r.set(i);
    }
  }
  return r;
}

Bitset Bitset::operator~() const {
  Bitset r(value_);
  r.value_.flip();
  return r;
}

Bitset Bitset::operator+(const Bitset& other) const {
  Bitset sum(value_.size() + other.size());
  sum.CopyFrom(*this, 0, false);
  sum.CopyFrom(other, value_.size(), false);
  return sum;
}

void Bitset::operator&=(const Bitset& other) {
  Assert(value_.size() == other.size(), "Size mismatch in Bitset::operator&=.");
  for (size_t i = 0; i < value_.size(); i++) {
    value_[i] = value_[i] && other[i];
  }
}

void Bitset::operator|=(const Bitset& other) {
  Assert(value_.size() == other.size(), "Size mismatch in Bitset::operator|=.");
  for (size_t i = 0; i < value_.size(); i++) {
    value_[i] = value_[i] || other[i];
  }
}

std::ostream& operator<<(std::ostream& os, const Bitset& bitset) {
  os << bitset.SubsplitToString();
  return os;
}

// ** Bitset Methods

void Bitset::Zero() { std::fill(value_.begin(), value_.end(), false); }

size_t Bitset::Hash() const { return std::hash<std::vector<bool>>{}(value_); }

std::string Bitset::ToString() const {
  std::string str;
  for (const auto& bit : value_) {
    str += (bit ? '1' : '0');
  }
  return str;
}

std::vector<size_t> Bitset::ToVectorOfSetBits() const {
  std::vector<size_t> vec;
  for (size_t i = 0; i < size(); i++) {
    if (value_[i]) {
      vec.push_back(i);
    }
  }
  return vec;
}

bool Bitset::All() const {
  for (const auto& bit : value_) {
    if (!bit) {
      return false;
    }
  }
  return true;
}

bool Bitset::Any() const {
  for (const auto& bit : value_) {
    if (bit) {
      return true;
    }
  }
  return false;
}

bool Bitset::None() const { return !Any(); }

bool Bitset::IsSingleton() const { return SingletonOption().has_value(); }

bool Bitset::IsDisjoint(const Bitset& other) const {
  Assert(size() == other.size(), "Size mismatch in Bitset::IsDisjoint.");
  for (size_t i = 0; i < size(); i++) {
    if (value_[i] && other.value_[i]) {
      return false;
    }
  }
  return true;
}

void Bitset::Minorize() {
  Assert(!(value_.empty()), "Can't Bitset::Minorize an empty bitset.");
  if (value_[0]) {
    value_.flip();
  }
}

void Bitset::CopyFrom(const Bitset& other, size_t begin, bool flip) {
  Assert(begin + other.size() <= size(), "Can't fit copy in Bitset::CopyFrom.");
  if (flip) {
    for (size_t i = 0; i < other.size(); i++) {
      value_[i + begin] = !other[i];
    }
  } else {
    for (size_t i = 0; i < other.size(); i++) {
      value_[i + begin] = other[i];
    }
  }
}

std::optional<uint32_t> Bitset::SingletonOption() const {
  bool found_already = false;
  uint32_t found_index;
  for (uint32_t i = 0; i < size(); i++) {
    if (value_[i]) {
      if (found_already) {
        // We previously found an index, so this isn't a singleton.
        return std::nullopt;
      }  // else
      found_already = true;
      found_index = i;
    }
  }
  if (found_already) {
    return found_index;
  }  // else
  return std::nullopt;
}

size_t Bitset::Count() const { return std::count(value_.begin(), value_.end(), true); }

std::string Bitset::ToVectorOfSetBitsAsString() const {
  std::string str;
  for (size_t i = 0; i < size(); i++) {
    if (value_[i]) {
      str += std::to_string(i);
      str += ",";
    }
  }
  if (!str.empty()) {
    str.pop_back();
  }
  return str;
}

// ** SBN-related functions

// ** Clade / MultiClade functions

int Bitset::CladeCompare(const Bitset& bitset_a, const Bitset& bitset_b) {
  // Comparing by lexigraphical taxon representation is the precise opposite of
  // comparing by their binary representation.  See header file for details.
  return (-1 * Bitset::Compare(bitset_a, bitset_b));
}

int Bitset::CladeCompare(const Bitset& that_bitset) const {
  const Bitset& this_bitset = *this;
  return CladeCompare(this_bitset, that_bitset);
}

size_t Bitset::MultiCladeGetCladeSize(const size_t clade_count) const {
  Assert(size() % clade_count == 0,
         "Bitset::MultiCladeGetCladeSize: Size isn't evenly divisible by clade_count.");
  return size() / clade_count;
}

Bitset Bitset::MultiCladeGetClade(const size_t i, const size_t clade_count) const {
  Assert(i < clade_count, "Bitset::MultiCladeGetClade: index is too large.");
  size_t clade_size = MultiCladeGetCladeSize(clade_count);
  std::vector<bool> new_value(
      value_.begin() + static_cast<std::vector<bool>::difference_type>(i * clade_size),
      value_.begin() +
          static_cast<std::vector<bool>::difference_type>((i + 1) * clade_size));
  return Bitset(std::move(new_value));
}

std::string Bitset::MultiCladeToString(const size_t clade_count) const {
  Assert(size() % clade_count == 0,
         "Size isn't a multiple of clade_count in Bitset::MultiCladeToString.");
  size_t clade_size = size() / clade_count;
  std::string str;
  for (size_t i = 0; i < value_.size(); ++i) {
    str += (value_[i] ? '1' : '0');
    if ((i + 1) % clade_size == 0 && i + 1 < value_.size()) {
      // The next item will start a new clade, so add a separator.
      str += '|';
    }
  }
  return str;
}

// ** Subsplit functions

Bitset Bitset::Subsplit(const Bitset& clade_0, const Bitset& clade_1) {
  // This asserts that clades are disjoint and equal-sized.
  Assert(clade_0.IsDisjoint(clade_1),
         "SubsplitOfPair: given bitsets are not a valid clade pair.");
  return SubsplitFromUnorderedClades(clade_0, clade_1);
}

Bitset Bitset::Subsplit(const std::string clade_0, const std::string clade_1) {
  return Bitset::Subsplit(Bitset(clade_0), Bitset(clade_1));
}

Bitset Bitset::Subsplit(const SizeVector clade_0, const SizeVector clade_1,
                        const size_t n) {
  return Bitset::Subsplit(Bitset(clade_0, n), Bitset(clade_1, n));
}

Bitset Bitset::SubsplitFromUnorderedClades(const Bitset& clade_0,
                                           const Bitset& clade_1) {
  Assert(clade_0.size() == clade_1.size(),
         "Bitset::SubsplitOrderClades requires Bitsets be the same size.");
  return CladeCompare(clade_0, clade_1) < 0 ? clade_0 + clade_1 : clade_1 + clade_0;
}

int Bitset::SubsplitCompare(const Bitset& subsplit_a, const Bitset& subsplit_b) {
  Assert(subsplit_a.size() == subsplit_b.size(),
         "Bitset::SubsplitCompare requires Bitsets be the same size.");
  // (1) Compare the number of taxa of the Subsplits.
  auto count_a = subsplit_a.Count();
  auto count_b = subsplit_b.Count();
  if (count_a != count_b) {
    return count_a - count_b;
  }
  // (2) Compare their respective union Bitsets.
  auto union_a = subsplit_a.SubsplitCladeUnion();
  auto union_b = subsplit_b.SubsplitCladeUnion();
  auto compare_union = Bitset::Compare(union_a, union_b);
  if (compare_union != 0) {
    return compare_union;
  }
  // (3) Compare the subsplit Bitsets.
  auto compare_subsplit = Bitset::Compare(subsplit_a, subsplit_b);
  return compare_subsplit;
}

int Bitset::SubsplitCompare(const Bitset& subsplit_b) const {
  const Bitset& subsplit_a = *this;
  return SubsplitCompare(subsplit_a, subsplit_b);
}

Bitset Bitset::SubsplitRotate() const {
  Assert(size() % 2 == 0, "Bitset::SubsplitRotate requires an even-size bitset.");
  Bitset clade_0 = SubsplitGetClade(0);
  Bitset clade_1 = SubsplitGetClade(1);
  return clade_1 + clade_0;
}

Bitset Bitset::SubsplitSort() const {
  Assert(size() % 2 == 0, "Bitset::SubsplitRotate requires an even-size bitset.");
  Bitset clade_0 = SubsplitGetClade(0);
  Bitset clade_1 = SubsplitGetClade(1);
  return SubsplitFromUnorderedClades(clade_0, clade_1);
}

std::string Bitset::SubsplitToString() const { return MultiCladeToString(2); }

std::string Bitset::SubsplitToVectorOfSetBitsAsString() const {
  std::string str;
  str += SubsplitGetClade(0).ToVectorOfSetBitsAsString();
  str += "|";
  str += SubsplitGetClade(1).ToVectorOfSetBitsAsString();
  return str;
}

size_t Bitset::SubsplitGetCladeSize() const {
  return MultiCladeGetCladeSize(SubsplitCladeCount);
}

Bitset Bitset::SubsplitGetClade(const size_t which_clade) const {
  return MultiCladeGetClade(which_clade, SubsplitCladeCount);
}

Bitset Bitset::SubsplitGetCladeByBinaryOrder(const size_t i) const {
  Assert(i < 2, "SubsplitGetClade index too large.");
  // The clades appear in taxon ordering (opposite of binary ordering); see "Clade
  // Methods" in header file for details.
  Bitset clade_0 = SubsplitGetClade(1);
  Bitset clade_1 = SubsplitGetClade(0);
  Assert(clade_1 > clade_0,
         "Bitset::SubsplitGetClade: Subsplit clade_1 should be larger than clade_0 (in "
         "binary rep).");
  return i == 0 ? clade_0 : clade_1;
}

bool Bitset::SubsplitIsLeaf() const {
  // A subsplit is a leaf if left clade has a
  bool is_left_clade_nonempty = SubsplitGetClade(0).IsSingleton();
  bool is_right_clade_empty = SubsplitGetClade(1).None();
  return is_left_clade_nonempty && is_right_clade_empty;
}

bool Bitset::SubsplitIsRoot() const {
  // A subsplit is a root if the left clade contains all taxons and the right clade
  // contains no taxons. If subsplit is valid, then we can assume the right clade is
  // empty.
  bool is_left_clade_full = SubsplitGetClade(0).All();
  return is_left_clade_full;
}

bool Bitset::SubsplitIsRootsplit() const {
  // A subsplit is a rootsplit if the union of the clades contain all clades.
  // But is also not the root, meaning both clades are nonempty.
  bool is_union_of_clades_full = SubsplitCladeUnion().All();
  bool is_left_clade_nonempty = !SubsplitGetClade(0).None();
  bool is_right_clade_nonempty = !SubsplitGetClade(1).None();
  return is_union_of_clades_full && is_left_clade_nonempty && is_right_clade_nonempty;
}

bool Bitset::SubsplitIsLeftChildOf(const Bitset& parent) const {
  return (size() == parent.size()) &&
         (SubsplitCladeUnion() == parent.SubsplitGetClade(0));
}

bool Bitset::SubsplitIsRightChildOf(const Bitset& parent) const {
  return (size() == parent.size()) &&
         (SubsplitCladeUnion() == parent.SubsplitGetClade(1));
}

Bitset Bitset::SubsplitCladeUnion() const {
  Assert(size() % SubsplitCladeCount == 0,
         "Size isn't 0 mod 2 in Bitset::SubsplitCladeUnion.");
  return SubsplitGetClade(0) | SubsplitGetClade(1);
}

bool Bitset::SubsplitIsWhichChildOf(const Bitset& parent, const Bitset& child) {
  Assert(parent.size() == child.size(),
         "Bitset::SubsplitIsWhichChildOf() bitsets are different sizes.");
  Bitset child_union = child.SubsplitCladeUnion();
  for (bool clade : {0, 1}) {
    if (child_union == parent.SubsplitGetClade(clade)) {
      return clade;
    }
  }
  // If it reaches the end, then it is not a parent.
  Failwith(
      "Bitset::SubsplitIsWhichChildOf(): given parent is not a parent of given child.");
}

bool Bitset::SubsplitIsParentChildPair(const Bitset& parent, const Bitset& child) {
  return child.SubsplitIsLeftChildOf(parent) || child.SubsplitIsRightChildOf(parent);
}

bool Bitset::SubsplitIsAdjacent(const Bitset& subsplit_a, const Bitset& subsplit_b) {
  return SubsplitIsParentChildPair(subsplit_a, subsplit_b) ||
         SubsplitIsParentChildPair(subsplit_b, subsplit_a);
}

bool Bitset::SubsplitIsValid() const {
  return SubsplitGetClade(0).IsDisjoint(SubsplitGetClade(1));
}

// ** Edge functions

Bitset Bitset::Edge(const Bitset& parent_subsplit, const Bitset& child_subsplit) {
  // Assert that:
  // - child_subsplit is either a sorted or rotated child of parent_subsplit.
  // - child_subsplit forms a valid subsplit.
  bool is_parent_valid = parent_subsplit.SubsplitIsValid();
  bool is_child_valid = child_subsplit.SubsplitIsValid();
  bool is_pair_valid = child_subsplit.SubsplitIsLeftChildOf(parent_subsplit) ||
                       child_subsplit.SubsplitIsRightChildOf(parent_subsplit);
  Assert(is_parent_valid && is_child_valid && is_pair_valid,
         "Edge(): given bitsets are not a valid parent/child pair.");

  if (child_subsplit.SubsplitIsLeftChildOf(parent_subsplit)) {
    return parent_subsplit.SubsplitRotate() +
           child_subsplit.SubsplitGetCladeByBinaryOrder(0);
  } else {
    return parent_subsplit + child_subsplit.SubsplitGetCladeByBinaryOrder(0);
  }
}

Bitset Bitset::Edge(const Bitset& sister_clade, const Bitset& focal_clade,
                    const Bitset& sorted_child_clade) {
  Assert(sister_clade.size() == focal_clade.size() &&
             focal_clade.size() == sorted_child_clade.size(),
         "Edge(): all clades must be of equal size.");
  Bitset pcsp = sister_clade + focal_clade + sorted_child_clade;
  Assert(pcsp.EdgeIsValid(), "Edge(): given clades form an invalid Edge.");
  return pcsp;
}

Bitset Bitset::Edge(const std::string sister_clade, const std::string focal_clade,
                    const std::string sorted_child_clade) {
  return Edge(Bitset(sister_clade), Bitset(focal_clade), Bitset(sorted_child_clade));
}

size_t Bitset::EdgeGetCladeSize() const { return MultiCladeGetCladeSize(3); }

Bitset Bitset::EdgeGetClade(const size_t i) const { return MultiCladeGetClade(i, 3); }

Bitset Bitset::EdgeGetParentSubsplit() const {
  Bitset sister = EdgeGetClade(0);
  Bitset focal = EdgeGetClade(1);
  return Bitset::Subsplit(sister, focal);
}

Bitset Bitset::EdgeGetChildSubsplit() const {
  Bitset focal = EdgeGetClade(1);
  Bitset child_0 = EdgeGetClade(2);
  Bitset child_1 = focal & ~child_0;
  return Bitset::Subsplit(child_0, child_1);
}

std::string Bitset::EdgeToString() const { return MultiCladeToString(EdgeCladeCount); }

bool Bitset::EdgeIsValid() const {
  if (size() % EdgeCladeCount != 0) {
    return false;
  }
  Bitset sister = EdgeGetClade(0);
  Bitset focal = EdgeGetClade(1);
  Bitset child_0 = EdgeGetClade(2);
  // The parent clades should be disjoint.
  if (!sister.IsDisjoint(focal)) {
    return false;
  }
  // The clade should split the focal clade of the parent,
  // so the taxa of child_0 should be a subset of those of focal clade.
  if (!child_0.IsDisjoint(~focal)) {
    return false;
  }
  // Something has to be set in each clade.
  if (sister.None() || focal.None() || child_0.None()) {
    return false;
  }
  return true;
}

bool Bitset::EdgeIsLeaf() const {
  Assert(size() % EdgeCladeCount == 0, "Size isn't 0 mod 3 in Bitset::EdgeIsLeaf.");
  // If third clade of Edge is empty, that means that the associated clade's sorted
  // subsplit is empty, so it is leaf.
  return EdgeGetClade(2).None();
}

Bitset Bitset::EdgeSort() const {
  Bitset parent = EdgeGetParentSubsplit();
  Bitset child = EdgeGetChildSubsplit();
  return Bitset::Edge(parent, child);
}

bool Bitset::EdgeIsParentRootsplit() const {
  Assert(size() % EdgeCladeCount == 0,
         "Size isn't 0 mod 3 in Bitset::EdgeIsRootsplit.");
  return EdgeGetParentSubsplit().SubsplitIsRootsplit();
}

SizePair Bitset::EdgeGetChildSubsplitTaxonCounts() const {
  auto clade_size = EdgeGetCladeSize();
  auto total_clade_taxon_count =
      std::count(value_.begin() + clade_size,
                 value_.begin() + SubsplitCladeCount * clade_size, true);
  auto clade0_taxon_count =
      std::count(value_.begin() + SubsplitCladeCount * clade_size, value_.end(), true);
  Assert(clade0_taxon_count < total_clade_taxon_count,
         "EdgeGetChildSubsplitTaxonCounts: not a proper Edge bitset.");
  return {static_cast<size_t>(clade0_taxon_count),
          static_cast<size_t>(total_clade_taxon_count - clade0_taxon_count)};
}

Bitset Bitset::Singleton(size_t n, size_t which_on) {
  Assert(which_on < n, "which_on too big in Bitset::Singleton.");
  Bitset singleton(n);
  singleton.set(which_on);
  return singleton;
}

Bitset Bitset::LeafSubsplit(const Bitset& nonzero_clade) {
  // Leaf pairs a nonempty left clade and an empty right clade.
  Bitset leaf_subsplit =
      Bitset::Subsplit(nonzero_clade, Bitset(nonzero_clade.size(), false));
  return leaf_subsplit;
}

void AssertSubsplitIsLeafAdjacent(const Bitset& subsplit) {
  // For subsplit to be adjacent to
  bool is_left_clade_nonempty = subsplit.SubsplitGetClade(0).Any();
  bool is_right_clade_singleton = subsplit.SubsplitGetClade(1).IsSingleton();
  Assert(is_left_clade_nonempty && is_right_clade_singleton,
         "Assertion SisterAndLeafSubsplit failed: we want the left-hand clade of the "
         "subsplit be "
         "non-empty and the right-hand clade be a singleton.");
}

Bitset Bitset::LeafChildSubsplit(const Bitset& parent_subsplit) {
  AssertSubsplitIsLeafAdjacent(parent_subsplit);
  // Put the right-hand clade of the subsplit as the nonzero contents of the leaf
  // subsplit.
  return LeafSubsplit(parent_subsplit.SubsplitGetClade(1));
}

Bitset Bitset::LeafEdge(const Bitset& parent_subsplit) {
  AssertSubsplitIsLeafAdjacent(parent_subsplit);
  const auto taxon_count = parent_subsplit.SubsplitGetCladeSize();
  Bitset leaf(3 * taxon_count);
  // Put the nonzero contents on the left of the leaf subsplit.
  leaf.CopyFrom(parent_subsplit, 0, false);
  return leaf;
}

Bitset Bitset::DAGRootSubsplitOfTaxonCount(const size_t taxon_count) {
  Bitset zeros(taxon_count);
  return ~zeros + zeros;
}

Bitset Bitset::RootsplitOfHalf(const Bitset& subsplit_half) {
  Bitset half = subsplit_half;
  half.Minorize();
  return ~half + half;
}

Bitset Bitset::EdgeOfRootsplit(const Bitset& rootsplit) {
  Assert(rootsplit.SubsplitIsRootsplit(),
         "Given subsplit is not rootsplit in Bitset::EdgeOfRootsplit.");
  return Edge(DAGRootSubsplitOfTaxonCount(rootsplit.size() / 2), rootsplit);
}

Bitset Remap(const Bitset& bitset, const SizeOptionVector& idx_table) {
  Bitset result(idx_table.size(), false);
  for (size_t i = 0; i < idx_table.size(); ++i) {
    if (idx_table[i].has_value() && bitset[idx_table[i].value()]) {
      result.set(i, true);
    }
  }
  return result;
}
