#ifndef GBOP2_TYPES_H
#define GBOP2_TYPES_H

#include <vector>
#include <RcppArmadillo.h>

struct Data {
  std::vector<int> num_event;
  std::vector<double> prob;
};

struct Result_singleE {
  Data nonstopping_df_null;
  Data nonstopping_df_alter;  // You can leave this unused in some files
  double nonstop_prob_null;
  double nonstop_prob_alter;
  double expected_size;
};

struct Result_dualE {
  Data effective_df;
  double effective_prob;
  double cummu_effective;
  Data nonstop_df;
  double expected_size;
  double stopping_prob;
  double nonstopping_prob;
  Rcpp::NumericVector prob_size_futile;
  Rcpp::NumericVector prob_size_effective;
  Rcpp::NumericVector prob_size_nonstop;
};

#endif
