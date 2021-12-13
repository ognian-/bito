// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "sbn_support.hpp"

StringVector SBNSupport::PrettyIndexer() const {
  StringVector pretty_representation(indexer_.size());
  for (const auto& [key, idx] : indexer_) {
    pretty_representation[idx] = key.EdgeToString();
  }
  return pretty_representation;
}

void SBNSupport::PrettyPrintIndexer() const {
  auto pretty_representation = PrettyIndexer();
  for (size_t i = 0; i < pretty_representation.size(); i++) {
    std::cout << i << "\t" << pretty_representation[i] << std::endl;
  }
}

// Return indexer_ and parent_to_range_ converted into string-keyed maps.
std::tuple<StringSizeMap, StringSizePairMap> SBNSupport::GetIndexers() const {
  auto str_indexer = StringifyMap(indexer_);
  auto str_parent_to_range = StringifyMap(parent_to_range_);
  std::string rootsplit("DAG Root Node");
  SafeInsert(str_parent_to_range, rootsplit, {0, rootsplits_.size()});
  return {str_indexer, str_parent_to_range};
}

// Get the indexer, but reversed and with bitsets appropriately converted to
// strings.
StringVector SBNSupport::StringReversedIndexer() const {
  StringVector reversed_indexer(indexer_.size());
  for (const auto& [key, idx] : indexer_) {
    reversed_indexer[idx] = key.EdgeToString();
  }
  return reversed_indexer;
}

void SBNSupport::ProbabilityNormalizeSBNParametersInLog(
    EigenVectorXdRef sbn_parameters) const {
  SBNProbability::ProbabilityNormalizeParamsInLog(sbn_parameters, rootsplits_.size(),
                                                  parent_to_range_);
}
