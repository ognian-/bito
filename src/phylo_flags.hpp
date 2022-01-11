// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// PhyloFlags are used for adding optional arguments to functions that are specified by
// the user. Additionally, it contains the keys used for accessing members of the
// GradientMap returned by PhyloGradients.
//

#pragma once

#include "phylo_model.hpp"
#include "substitution_model.hpp"
#include "sugar.hpp"
#include "tree_gradient.hpp"

// ** Single Option
// This is the base option flag type.  Also specifies default behaviour for flags.
class PhyloFlagName {
 public:
  enum class FlagType { None, Include, Exclude, SetValue, MapKey, RunAll };
  enum class DataType { None, Bool, Double };

  PhyloFlagName();
  PhyloFlagName(std::string name, std::string flag, bool is_set_when_running_defaults,
                FlagType flag_type, DataType data_type, StringVector child_flags);
  PhyloFlagName(std::string name, std::string flag, bool is_set_when_running_defaults,
                FlagType flag_type, DataType data_type);
  PhyloFlagName(std::string name, std::string flag, bool is_set_when_running_defaults,
                FlagType flag_type);

  // This is a descriptive name of the flag option that will be visible to the user in
  // bito pybind interface.
  std::string name_;
  // This is the uniquely identifiable string that is used for setting/adding flag
  // options.  For convenience, this should match the output mapkey if that mapkey
  // directly corresponds to a function option for the computation of underlying data.
  // (e.g. in FatBeagle::Gradient, we have an option to compute
  // `substitution_model_rates`. If flag is set, then the return map will contain a
  // `substitution_model_rates` key).
  std::string flag_;
  // Determines whether to consider this option set or unset if `is_run_defaults` flag
  // is set, but this option has not explicitly been set. If option is not set and
  // `is_run_defaults` is not set, then all flags are considered false, so includes are
  // not run and excludes are run.
  bool is_set_when_running_defaults_;
  // This gives the type of flag. There are:
  // - Include: these options are run only when specified.
  // - Exclude: these options are NOT run only when specified.
  // - SetValue: these options override some input value when specified.
  // - Mapkey: these options correspond to an key for an output map.
  FlagType flag_type_;
  // This gives the underlying datatype of the flag.
  // Datatype is boolean if anything other than a SetValue.  Currently only bool and
  // double are supported.
  DataType data_type_;
  // These allow for subflags, corresponding to subroutines of given superflag routine.
  // (e.g. in FatBeagle::Gradient, `substitution_model` flag has two subflags,
  // `substitution_model_rates` and `substitution_model_frequencies`. If both subflags
  // are set, we should consider the superflag set as well.)
  StringVector child_flags_;

  // Child Ids: Ids for the subflags.
  void AddChild(const PhyloFlagName &child);
  void AddChild(const std::string child_flag);
  void AddChildren(const StringVector child_flags);
  StringVector GetChildren() const;
  // Output to String.
  virtual std::string ToString() const;
  // Comparators
  static int Compare(const PhyloFlagName &flag_a, const PhyloFlagName &flag_b);
  static int StringCompare(const std::string &str_a, const std::string &str_b);
  // General compare.
  bool operator==(const PhyloFlagName &other);
  friend bool operator==(const PhyloFlagName &lhs, const PhyloFlagName &rhs);
  bool operator<(const PhyloFlagName &other);
  friend bool operator<(const PhyloFlagName &lhs, const PhyloFlagName &rhs);
  // Compare against String
  bool operator==(const std::string &other_name);
  bool operator<(const std::string &other_name);
  // This establishes how flag type default value if not using `is_run_defaults_`.
  bool FlagTypeDefaultBehavior() const;
};

// ** Base options for all functions.
class PhyloFlagOptions {
 public:
  // This is a special option that is contained in all option sets.
  static inline PhyloFlagName run_defaults_ = {"RUN_DEFAULTS", "run_defaults", false,
                                               PhyloFlagName::FlagType::RunAll};
  // List of all possible options user can set.
  std::set<PhyloFlagName> all_options_ = {run_defaults_};
  // Flag to String.
  std::string ToString();
};

// ** Function options for FatBeagle::Gradient().
class PhyloFlagGradientOptions : public PhyloFlagOptions {
 public:
  static PhyloFlagGradientOptions CreatePhyloFlagGradientOptions();

 public:
  // Include options.
  PhyloFlagName clock_model_rates_;
  PhyloFlagName ratios_root_height_;
  PhyloFlagName site_model_parameters_;
  PhyloFlagName substitution_model_rates_;
  PhyloFlagName substitution_model_frequencies_;
  PhyloFlagName substitution_model_;
  // Exclude options.
  PhyloFlagName exclude_log_det_jacobian_gradient_;
  // SetValue options.
  PhyloFlagName set_gradient_delta_;
  // Set of All Gradient options.
  std::set<PhyloFlagName> all_options_;
};

// ** Function options for FatBeagle::LogLikelihood().
class PhyloFlagLoglikelihoodOptions : public PhyloFlagOptions {
 public:
  static PhyloFlagLoglikelihoodOptions CreatePhyloFlagLoglikelihoodOptions();

 public:
  // Exclude options.
  PhyloFlagName exclude_log_det_jacobian_likelihood_;
  // Set of All LogLikelihood options.
  std::set<PhyloFlagName> all_options_;
};

// ** Map Keys for FatBeagle::Gradient().
class PhyloFlagGradientMapKeys : public PhyloFlagOptions {
 public:
  static PhyloFlagGradientMapKeys CreatePhyloFlagGradientMapKeys();

 public:
  // Gradient Mapkeys.
  PhyloFlagName site_model_;
  PhyloFlagName ratios_root_height_;
  PhyloFlagName substitution_model_;
  PhyloFlagName substitution_model_rates_;
  PhyloFlagName substitution_model_frequencies_;
  PhyloFlagName site_model_parameters_;
  PhyloFlagName clock_model_rates_;
  PhyloFlagName entire_substitution_model_;
  PhyloFlagName entire_site_model_;
  PhyloFlagName entire_clock_model_;
  // Vector of All Gradient map keys.
  std::set<PhyloFlagName> all_keys_;
};

// ** All User Options
class PhyloFlags {
 public:
  PhyloFlags(bool is_run_defaults = true);
  PhyloFlags(const StringVector &key_vec, bool is_run_defaults = true);
  PhyloFlags(const StringDoubleVector &key_vec, bool is_run_defaults = true);

  // ** Find Flag
  // Find flag by name.
  std::optional<PhyloFlagName> FindFlagByName(const std::set<PhyloFlagName> &options,
                                              const std::string &flag_name,
                                              std::optional<PhyloFlagName> flag) const;

  // ** FlagCheck
  // These functions are the ones to be called by utilizing functions (FatBeagle, etc).
  // FlagCheck functions determine whether the associated flag will be evaluated as true
  // or false. (1) First, returns the flag's value if it has been explicitly set by a
  // SetFlag function. (2) If not, then it checks whether the RunDefault flag has been
  // set, in which case the flag's default behavior is returned. (3) If neither have
  // been set, then it returns false.
  bool IsFlagChecked(const PhyloFlagName &flag) const;
  bool IsFlagNotChecked(const PhyloFlagName &flag) const;
  // Checks if a flag if user may or may not have passed any options.
  // If options have not been passed, uses flag's default behavior.
  static bool IsFlagChecked(const std::optional<PhyloFlags> phylo_flags,
                            const PhyloFlagName &flag);
  static bool IsFlagNotChecked(const std::optional<PhyloFlags> phylo_flags,
                               const PhyloFlagName &flag);

  // Returns the value of the flag.
  std::optional<double> GetFlagValue(const std::string &flag_name) const;
  std::optional<double> GetFlagValue(const PhyloFlagName &flag) const;
  // Returns the flag's value if set, otherwise returns default value.
  double GetFlagValueIfSet(const std::string &flag_name, double default_value) const;
  double GetFlagValueIfSet(const PhyloFlagName &flag, double default_value) const;
  static double GetFlagValueIfSet(const std::optional<PhyloFlags> phylo_flags,
                                  const PhyloFlagName &flag, double default_value);

  // ** Misc
  // Convert flag name vector to boolean vector.
  BoolVector ToBoolFlags(const StringVector &vec_of_keys) const;
  // Interprets flags as a string.
  std::string ToString() const;

  // Default Option
  static inline PhyloFlagName run_defaults_ = PhyloFlagOptions::run_defaults_;
  // All Options
  static inline PhyloFlagGradientOptions gradient_options_ =
      PhyloFlagGradientOptions::CreatePhyloFlagGradientOptions();
  static inline PhyloFlagLoglikelihoodOptions loglikelihood_options_ =
      PhyloFlagLoglikelihoodOptions::CreatePhyloFlagLoglikelihoodOptions();
  // All MapKeys
  static inline PhyloFlagGradientMapKeys gradient_mapkeys_ =
      PhyloFlagGradientMapKeys::CreatePhyloFlagGradientMapKeys();
  // Current Options
  PhyloFlagOptions *options_ = nullptr;

 private:
  // ** FlagSet
  // FlagSet functions add or return a boolean and associated value to/from the map.

  // Set flag option.
  void SetFlag(const std::string &flag_name, const bool set = true,
               const double value = 1.0f);
  void SetFlag(const PhyloFlagName &flag, const bool set = true,
               const double value = 1.0f);
  // Set flag option with specified value.
  void SetFlagWithValue(const std::string &flag_name, const double value);
  void SetFlagWithValue(const PhyloFlagName &flag, const double value);
  // Check if flag option has been explicitly set.
  bool IsFlagSet(const std::string &flag_name) const;
  bool IsFlagSet(const PhyloFlagName &flag) const;
  bool IsFlagNotSet(const std::string &flag_name) const;
  bool IsFlagNotSet(const PhyloFlagName &flag) const;
  // Add flag to map.
  void AddFlagToMap(const PhyloFlagName &flag, const bool set = true,
                    const double value = 1.0f);
  // Check if flag has been added to the map.
  bool IsFlagInMap(const std::string &flag) const;

  // ** RunDefaults
  // RunDefaults flag triggers each flag's default behavior, if flag has no been
  // explicitly set.

  // Special case for the is_run_defaults key
  void SetRunDefaultsFlag(bool is_set);
  // Special flag for indicating all options to be considered true.
  bool IsRunDefaultsFlagSet() const;
  // Stores all option flags that have been manually modified, with a bool whether the
  // flag has been set, and an associated data value.
  std::map<std::string, std::pair<bool, std::optional<double>>> map_;
  // This is a special flag that determines behavior if option is not explicitly set.
  // If is_run_defaults_ is false, all flags are treated as if unset.
  // Otherwise, all flags are treated as their default.
  bool is_run_defaults_;
};

