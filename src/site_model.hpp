// Copyright 2019-2021 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#pragma once

#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "block_model.hpp"

class SiteModel : public BlockModel {
 public:
  SiteModel(const BlockSpecification::ParamCounts& param_counts)
      : BlockModel(param_counts) {}
  virtual ~SiteModel() = default;

  virtual size_t GetCategoryCount() const = 0;
  virtual const EigenVectorXd& GetCategoryRates() const = 0;
  virtual const EigenVectorXd& GetCategoryProportions() const = 0;
  virtual const EigenVectorXd& GetRateGradient() const = 0;

  static std::unique_ptr<SiteModel> OfSpecification(const std::string& specification);
};

class ConstantSiteModel : public SiteModel {
 public:
  ConstantSiteModel()
      : SiteModel({}), zero_(EigenVectorXd::Zero(1)), one_(EigenVectorXd::Ones(1)) {}

  size_t GetCategoryCount() const override { return 1; }

  const EigenVectorXd& GetCategoryRates() const override { return one_; }

  const EigenVectorXd& GetCategoryProportions() const override { return one_; }

  const EigenVectorXd& GetRateGradient() const override { return zero_; };

  void SetParameters(const EigenVectorXdRef param_vector) override{};

 private:
  EigenVectorXd zero_;
  EigenVectorXd one_;
};

class WeibullSiteModel : public SiteModel {
 public:
  explicit WeibullSiteModel(size_t category_count, double shape)
      : SiteModel({{shape_key_, 1}}),
        category_count_(category_count),
        shape_(shape),
        rate_derivatives_(category_count) {
    category_rates_.resize(category_count);
    category_proportions_.resize(category_count);
    for (size_t i = 0; i < category_count; i++) {
      category_proportions_[i] = 1.0 / category_count;
    }
    UpdateRates();
  }

  size_t GetCategoryCount() const override;
  const EigenVectorXd& GetCategoryRates() const override;
  const EigenVectorXd& GetCategoryProportions() const override;
  const EigenVectorXd& GetRateGradient() const override;

  void SetParameters(const EigenVectorXdRef param_vector) override;

  inline const static std::string shape_key_ = "Weibull_shape";

 private:
  void UpdateRates();

  size_t category_count_;
  double shape_;  // shape of the Weibull distribution
  EigenVectorXd rate_derivatives_;
  EigenVectorXd category_rates_;
  EigenVectorXd category_proportions_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
#include <algorithm>
TEST_CASE("SiteModel") {
  // Test 1: First we test using the "built in" default values.
  auto weibull_model = std::make_unique<WeibullSiteModel>(4, 1.0);
  const EigenVectorXd rates = weibull_model->GetCategoryRates();
  EigenVectorXd rates_r(4);
  rates_r << 0.1457844, 0.5131316, 1.0708310, 2.2702530;
  CheckVectorXdEquality(rates, rates_r, 0.0001);

  // Test 2: Now set param_vector using SetParameters.
  weibull_model = std::make_unique<WeibullSiteModel>(4, 1.0);
  EigenVectorXd param_vector(1);
  param_vector << 0.1;
  weibull_model->SetParameters(param_vector);
  rates_r << 4.766392e-12, 1.391131e-06, 2.179165e-03, 3.997819e+00;
  const EigenVectorXd rates2 = weibull_model->GetCategoryRates();
  CheckVectorXdEquality(rates2, rates_r, 0.0001);

  // Test 3: Check proportions.
  const EigenVectorXd proportions = weibull_model->GetCategoryProportions();
  CheckVectorXdEquality(0.25, proportions, 0.0001);

  // Test 4: Check sum rates[i]*proportions[i]==1.
  CHECK_LT(fabs(rates.dot(proportions) - 1.), 0.0001);
  CHECK_LT(fabs(rates2.dot(proportions) - 1.), 0.0001);
}
#endif  // DOCTEST_LIBRARY_INCLUDED
