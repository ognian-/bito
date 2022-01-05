// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <map>
#include <vector>

using GradientMap = std::map<std::string, std::vector<double>>;

struct PhyloGradient {
  PhyloGradient() = default;
  PhyloGradient(double log_likelihood, GradientMap& gradient)
      : log_likelihood_(log_likelihood), gradient_(gradient){};

  double log_likelihood_;

  GradientMap gradient_;

  // GradientMap keys
  inline const static std::string site_model_key_ = "site_model";
  inline const static std::string clock_model_key_ = "clock_model";
  inline const static std::string substitution_model_key_ = "substitution_model";
  inline const static std::string substitution_model_rates_key_ =
      SubstitutionModel::rates_key_;
  inline const static std::string substitution_model_frequencies_key_ =
      SubstitutionModel::frequencies_key_;
  inline const static std::string branch_lengths_key_ = "branch_lengths";
  inline const static std::string ratios_root_height_key_ = "ratios_root_height";
};
