// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "phylo_flags.hpp"

#include "phylo_model.hpp"
#include "substitution_model.hpp"
#include "sugar.hpp"
#include "tree_gradient.hpp"

// ** PhyloFlagName

PhyloFlagName::PhyloFlagName()
    : name_(),
      flag_(),
      is_set_when_running_defaults_(false),
      flag_type_(FlagType::None),
      data_type_(DataType::None),
      child_flags_(){};

PhyloFlagName::PhyloFlagName(std::string name, std::string flag,
                             bool is_set_when_running_defaults, FlagType flag_type,
                             DataType data_type, StringVector child_flags)
    : name_(name),
      flag_(flag),
      is_set_when_running_defaults_(is_set_when_running_defaults),
      flag_type_(flag_type),
      data_type_(data_type),
      child_flags_(child_flags){};

PhyloFlagName::PhyloFlagName(std::string name, std::string flag,
                             bool is_set_when_running_defaults, FlagType flag_type,
                             DataType data_type)
    : name_(name),
      flag_(flag),
      is_set_when_running_defaults_(is_set_when_running_defaults),
      flag_type_(flag_type),
      data_type_(data_type),
      child_flags_(){};

PhyloFlagName::PhyloFlagName(std::string name, std::string flag,
                             bool is_set_when_running_defaults, FlagType flag_type)
    : name_(name),
      flag_(flag),
      is_set_when_running_defaults_(is_set_when_running_defaults),
      flag_type_(flag_type),
      data_type_(DataType::None),
      child_flags_(){};

void PhyloFlagName::AddChild(const PhyloFlagName &child) { AddChild(child.flag_); }

void PhyloFlagName::AddChild(const std::string child_flag) {
  child_flags_.push_back(child_flag);
}

void PhyloFlagName::AddChildren(const StringVector child_flags) {
  for (const auto child_flag : child_flags) {
    child_flags_.push_back(child_flag);
  }
}

StringVector PhyloFlagName::GetChildren() const { return child_flags_; }

std::string PhyloFlagName::ToString() const {
  std::stringstream str;
  str << "{ " << name_ << ": " << flag_ << " }";
  return str.str();
};

int PhyloFlagName::Compare(const PhyloFlagName &flag_a, const PhyloFlagName &flag_b) {
  return StringCompare(flag_a.flag_, flag_b.flag_);
}

int PhyloFlagName::StringCompare(const std::string &str_a, const std::string &str_b) {
  for (size_t i = 0; i < std::min(str_a.size(), str_b.size()); i++) {
    if (str_a[i] - str_b[i] != 0) {
      return str_a[i] - str_b[i];
    }
  }
  return str_a.size() - str_b.size();
}

bool PhyloFlagName::operator==(const PhyloFlagName &other) {
  return Compare(*this, other) == 0;
}

bool operator==(const PhyloFlagName &lhs, const PhyloFlagName &rhs) {
  return PhyloFlagName::Compare(lhs, rhs) == 0;
}

bool PhyloFlagName::operator<(const PhyloFlagName &other) {
  return Compare(*this, other) < 0;
}

bool operator<(const PhyloFlagName &lhs, const PhyloFlagName &rhs) {
  return PhyloFlagName::Compare(lhs, rhs) < 0;
}

bool PhyloFlagName::operator==(const std::string &other_name) {
  return flag_ == other_name;
}

bool PhyloFlagName::operator<(const std::string &other_name) {
  return flag_ < other_name;
}

bool PhyloFlagName::FlagTypeDefaultBehavior() const {
  switch (flag_type_) {
    case FlagType::Include:
      return false;
    case FlagType::Exclude:
      return false;
    case FlagType::SetValue:
      return false;
    case FlagType::MapKey:
      return false;
    case FlagType::RunAll:
      return false;
    default:
      Failwith("Invalid flag_type_ set.");
  }
}

// ** PhyloFlagOptions

std::string PhyloFlagOptions::ToString() {
  std::stringstream str;
  for (const auto option : all_options_) {
    str << option.name_ << " | " << option.flag_ << " | " << option.child_flags_
        << std::endl;
  }
  return str.str();
}

PhyloFlagGradientOptions PhyloFlagGradientOptions::CreatePhyloFlagGradientOptions() {
  PhyloFlagGradientOptions opt;
  // Include options
  opt.clock_model_rates_ = {"CLOCK_MODEL_RATES", StrictClockModel::rate_key_, true,
                            PhyloFlagName::FlagType::Include};
  opt.ratios_root_height_ = {"RATIOS_ROOT_HEIGHT",
                             PhyloGradient::ratios_root_height_key_, true,
                             PhyloFlagName::FlagType::Include};
  opt.site_model_parameters_ = {"SITE_MODEL_PARAMETERS", WeibullSiteModel::shape_key_,
                                true, PhyloFlagName::FlagType::Include};
  opt.substitution_model_rates_ = {"SUBSTITUTION_MODEL_RATES",
                                   SubstitutionModel::rates_key_, true,
                                   PhyloFlagName::FlagType::Include};
  opt.substitution_model_frequencies_ = {"SUBSTITUTION_MODEL_FREQUENCIES",
                                         SubstitutionModel::frequencies_key_, true,
                                         PhyloFlagName::FlagType::Include};
  opt.substitution_model_ = {"SUBSTITUTION_MODEL",
                             PhyloGradient::substitution_model_key_, true,
                             PhyloFlagName::FlagType::Include};
  // Exclude options
  opt.exclude_log_det_jacobian_gradient_ = {"EXCLUDE_LOG_DET_JACOBIAN_GRADIENT",
                                            "exclude_log_det_jacobian_gradient", false,
                                            PhyloFlagName::FlagType::Exclude};
  // SetValue options
  opt.set_gradient_delta_ = {"SET_GRADIENT_DELTA", "set_gradient_delta", false,
                             PhyloFlagName::FlagType::SetValue,
                             PhyloFlagName::DataType::Double};
  // Add child flags.
  opt.substitution_model_.AddChild(opt.substitution_model_rates_);
  opt.substitution_model_.AddChild(opt.substitution_model_frequencies_);
  // Set of All options
  opt.all_options_ = {opt.clock_model_rates_,
                      opt.ratios_root_height_,
                      opt.site_model_parameters_,
                      opt.substitution_model_,
                      opt.substitution_model_rates_,
                      opt.substitution_model_frequencies_,
                      opt.exclude_log_det_jacobian_gradient_,
                      opt.set_gradient_delta_,
                      opt.run_defaults_};

  return opt;
};

PhyloFlagLoglikelihoodOptions
PhyloFlagLoglikelihoodOptions::CreatePhyloFlagLoglikelihoodOptions() {
  PhyloFlagLoglikelihoodOptions opt;
  // Exclude options
  opt.exclude_log_det_jacobian_likelihood_ = {"EXCLUDE_LOG_DET_JACOBIAN_LIKELIHOOD",
                                              "exclude_log_det_jacobian_likelihood",
                                              false, PhyloFlagName::FlagType::Exclude};
  // All LogLikelihood options
  opt.all_options_ = {opt.exclude_log_det_jacobian_likelihood_, opt.run_defaults_};

  return opt;
};

PhyloFlagGradientMapKeys PhyloFlagGradientMapKeys::CreatePhyloFlagGradientMapKeys() {
  PhyloFlagGradientMapKeys opt;
  opt.site_model_ = {"SITE_MODEL", PhyloGradient::site_model_key_, true,
                     PhyloFlagName::FlagType::MapKey};
  opt.ratios_root_height_ = {"RATIOS_ROOT_HEIGHT",
                             PhyloGradient::ratios_root_height_key_, true,
                             PhyloFlagName::FlagType::MapKey};
  opt.substitution_model_ = {"SUBSTITUTION_MODEL",
                             PhyloGradient::substitution_model_key_, true,
                             PhyloFlagName::FlagType::MapKey};
  opt.substitution_model_rates_ = {"SUBSTITUTION_MODEL_RATES",
                                   SubstitutionModel::rates_key_, true,
                                   PhyloFlagName::FlagType::MapKey};
  opt.substitution_model_frequencies_ = {"SUBSTITUTION_MODEL_FREQUENCIES",
                                         SubstitutionModel::frequencies_key_, true,
                                         PhyloFlagName::FlagType::MapKey};
  opt.site_model_parameters_ = {"SITE_MODEL_PARAMETERS", WeibullSiteModel::shape_key_,
                                true, PhyloFlagName::FlagType::MapKey};
  opt.clock_model_rates_ = {"CLOCK_MODEL_RATES", StrictClockModel::rate_key_, true,
                            PhyloFlagName::FlagType::MapKey};
  opt.entire_substitution_model_ = {"ENTIRE_SUBSTITUTION_MODEL",
                                    PhyloModel::entire_substitution_key_, true,
                                    PhyloFlagName::FlagType::MapKey};
  opt.entire_site_model_ = {"ENTIRE_SITE_MODEL", PhyloModel::entire_site_key_, true,
                            PhyloFlagName::FlagType::MapKey};
  opt.entire_clock_model_ = {"ENTIRE_CLOCK_MODEL", PhyloModel::entire_clock_key_, true,
                             PhyloFlagName::FlagType::MapKey};
  // Add child flags.
  opt.substitution_model_.AddChild(opt.substitution_model_rates_);
  opt.substitution_model_.AddChild(opt.substitution_model_frequencies_);
  // All Gradient map keys.
  opt.all_keys_ = {opt.site_model_,
                   opt.ratios_root_height_,
                   opt.substitution_model_,
                   opt.substitution_model_rates_,
                   opt.substitution_model_frequencies_,
                   opt.site_model_parameters_,
                   opt.clock_model_rates_,
                   opt.entire_substitution_model_,
                   opt.entire_site_model_,
                   opt.entire_clock_model_};

  return opt;
};

// ** PhyloFlags

PhyloFlags::PhyloFlags(const bool is_run_defaults)
    : map_(), is_run_defaults_(is_run_defaults){};

PhyloFlags::PhyloFlags(const StringVector &key_vec, const bool is_run_defaults)
    : map_(), is_run_defaults_(is_run_defaults) {
  for (auto &key : key_vec) {
    SetFlag(key);
  }
};

PhyloFlags::PhyloFlags(const StringDoubleVector &key_vec, const bool is_run_defaults)
    : map_(), is_run_defaults_(is_run_defaults) {
  for (auto &key : key_vec) {
    SetFlag(key.first, true, key.second);
  }
}

std::optional<PhyloFlagName> PhyloFlags::FindFlagByName(
    const std::set<PhyloFlagName> &options, const std::string &flag_name,
    std::optional<PhyloFlagName> flag) const {
  if (flag.has_value()) {
    return flag;
  }
  for (const auto &flag : options) {
    if (flag.flag_ == flag_name) {
      return flag;
    }
  }
  return std::nullopt;
}

void PhyloFlags::SetFlag(const std::string &flag_name, const bool set,
                         const double value) {
  // Find flag.
  std::optional<PhyloFlagName> flag = std::nullopt;
  flag = FindFlagByName(gradient_options_.all_options_, flag_name, flag);
  flag = FindFlagByName(loglikelihood_options_.all_options_, flag_name, flag);
  Assert(flag.has_value(),
         "Attempted to set a option flag by name that does not exist: \"" + flag_name +
             "\"");
  SetFlag(*flag, set, value);
};

void PhyloFlags::SetFlag(const PhyloFlagName &flag, const bool set,
                         const double value) {
  // Add given flag.
  AddFlagToMap(flag, set, value);
  // Add all child flags of given flag.
  for (const auto child_flag : flag.child_flags_) {
    SetFlag(child_flag, value);
  }
  // If flag being set is the special run_defaults_ flag.
  if (PhyloFlagOptions::run_defaults_.name_ == flag.name_) {
    SetRunDefaultsFlag(true);
  }
};

void PhyloFlags::SetFlagWithValue(const std::string &flag_name, const double value) {
  SetFlag(flag_name, true, value);
};

void PhyloFlags::SetFlagWithValue(const PhyloFlagName &flag, const double value) {
  SetFlag(flag, true, value);
};

void PhyloFlags::AddFlagToMap(const PhyloFlagName &flag, const bool set,
                              const double value) {
  map_.insert(std::make_pair(flag.flag_, std::make_pair(set, value)));
}

bool PhyloFlags::IsFlagSet(const PhyloFlagName &flag) const {
  return IsFlagSet(flag.flag_);
}

bool PhyloFlags::IsFlagSet(const std::string &flag_name) const {
  for (auto map_flag = map_.begin(); map_flag != map_.end(); map_flag++) {
    if (map_flag->first == flag_name) {
      return map_flag->second.first;
    }
  }
  return false;
};

void PhyloFlags::SetRunDefaultsFlag(bool is_set) { is_run_defaults_ = is_set; }

bool PhyloFlags::IsRunDefaultsFlagSet() const { return is_run_defaults_; }

bool PhyloFlags::IsFlagNotSet(const PhyloFlagName &flag) const {
  return IsFlagNotSet(flag.flag_);
}

bool PhyloFlags::IsFlagNotSet(const std::string &flag_name) const {
  return (!IsFlagSet(flag_name));
};

bool PhyloFlags::IsFlagInMap(const std::string &flag_name) const {
  return (map_.find(flag_name) != map_.end());
};

std::optional<double> PhyloFlags::GetFlagValue(const PhyloFlagName &flag) const {
  return GetFlagValue(flag.flag_);
}

std::optional<double> PhyloFlags::GetFlagValue(const std::string &flag_name) const {
  if (IsFlagInMap(flag_name)) {
    return map_.at(flag_name).second;
  }
  return std::nullopt;
};

// Returns the value of the flag if set, otherwise returns default value.
double PhyloFlags::GetFlagValueIfSet(const std::string &flag_name,
                                     const double default_value) const {
  auto opt_value = GetFlagValue(flag_name);
  if (opt_value.has_value()) {
    return *opt_value;
  }
  return default_value;
};
double PhyloFlags::GetFlagValueIfSet(const PhyloFlagName &flag,
                                     const double default_value) const {
  return GetFlagValueIfSet(flag.flag_, default_value);
};
double PhyloFlags::GetFlagValueIfSet(const std::optional<PhyloFlags> phylo_flags,
                                     const PhyloFlagName &flag, double default_value) {
  if (phylo_flags.has_value()) {
    return phylo_flags->GetFlagValueIfSet(flag, default_value);
  }
  return default_value;
};

// Pass an vector of keys, and returns an vector of booleans whether their flag is
// set.
BoolVector PhyloFlags::ToBoolFlags(const StringVector &vec_of_keys) const {
  BoolVector bool_flags = BoolVector(vec_of_keys.size());
  for (size_t i = 0; i < vec_of_keys.size(); i++) {
    bool_flags[i] = IsFlagSet(vec_of_keys[i]);
  }
  return bool_flags;
};

std::string PhyloFlags::ToString() const {
  std::ostringstream rep;
  rep << "{ ";
  for (const auto &[key, value] : map_) {
    rep << "(" << key << ": " << value.first << ") ";
  }
  rep << "}";
  return rep.str();
};

bool PhyloFlags::IsFlagChecked(const PhyloFlagName &flag) const {
  // (1) Priority is given to explicitly set options.
  if (IsFlagSet(flag)) {
    return true;
  }
  // (2) If is_run_default_ option is set, use given individual flag's defined default
  // behavior.
  if (is_run_defaults_) {
    return flag.is_set_when_running_defaults_;
  }
  // (3) Otherwise, use flag type-based's default behavior.
  return flag.FlagTypeDefaultBehavior();
};

bool PhyloFlags::IsFlagNotChecked(const PhyloFlagName &flag) const {
  return !IsFlagChecked(flag);
};

bool PhyloFlags::IsFlagChecked(const std::optional<PhyloFlags> phylo_flags,
                               const PhyloFlagName &flag) {
  // (1) If user has not passed any flags, then fall back to default behavior.
  if (!phylo_flags.has_value()) {
    return flag.is_set_when_running_defaults_;
  };
  // (2) If user passed flags, then check if option is set.
  bool is_checked = phylo_flags->IsFlagChecked(flag);
  return is_checked;
};

bool PhyloFlags::IsFlagNotChecked(const std::optional<PhyloFlags> phylo_flags,
                                  const PhyloFlagName &flag) {
  return !PhyloFlags::IsFlagChecked(phylo_flags, flag);
};
