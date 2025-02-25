#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>
#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>

using namespace Rcpp;


//#include <RcppEigen.h>
//using namespace Eigen;
using namespace std;

// set seed
// [[Rcpp::export]]
void set_seed_dual(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}



// Define a struct to hold data
struct Data {
  std::vector<int> num_event;
  std::vector<double> prob;
};

// Define a struct to hold result
struct Result {
  Data effective_df;
  double effective_prob;
  double cummu_effective;
  Data nonstop_df;
  double expected_size;
  double stopping_prob;
  double nonstopping_prob;
  NumericVector prob_size_futile;
  NumericVector prob_size_effective;
  NumericVector prob_size_nonstop;
};

// Function to calculate exact error recursively
Result exacterror_recursive(Rcpp::NumericVector nobs, Rcpp::NumericVector ncut1, Rcpp::NumericVector ncut2, double test_prob, int ntotal) {
  int n_c = nobs.size();
  double cummu_effective = 0;
  double ptsa = 0;
  double stopping_prob = 0;
  double nonstopping_prob = 0;
  std::vector<double> prob_sizes_futile;
  std::vector<double> prob_sizes_effective;
  std::vector<double> prob_sizes_nonstop;

  if (n_c == 1) {
    Data effective_df;
    std::vector<double> effective_prob_vec;
    for (int i = ncut2[0]; i < nobs[0] + 1; i++) {
      effective_df.num_event.push_back(i);
      effective_prob_vec.push_back(R::dbinom(i, nobs[0], test_prob, false));
    }
    effective_df.prob = effective_prob_vec;

    Data nonstop_df;
    std::vector<double> nonstop_prob_vec;
    for (int i = ncut1[0] + 1; i < ncut2[0]; i++) {
      nonstop_df.num_event.push_back(i);
      nonstop_prob_vec.push_back(R::dbinom(i, nobs[0], test_prob, false));
    }
    nonstop_df.prob = nonstop_prob_vec;

    Data futile_df;
    std::vector<double> futile_prob_vec;
    for (int i = 0; i < ncut1[0] + 1; i++) {
      futile_df.num_event.push_back(i);
      futile_prob_vec.push_back(R::dbinom(i, nobs[0], test_prob, false));
    }
    futile_df.prob = futile_prob_vec;

    double futile_prob = std::accumulate(futile_df.prob.begin(), futile_df.prob.end(), 0.0);
    double effective_prob = std::accumulate(effective_df.prob.begin(), effective_df.prob.end(), 0.0);
    double nonstop_prob = std::accumulate(nonstop_df.prob.begin(), nonstop_df.prob.end(), 0.0);

    if (n_c == ntotal) {
      ptsa = ptsa + (futile_prob + effective_prob + nonstop_prob) * nobs[0];
    } else {
      ptsa = ptsa + futile_prob * nobs[0] + effective_prob * nobs[0];
      stopping_prob = stopping_prob + futile_prob + effective_prob;
      nonstopping_prob = nonstopping_prob + nonstop_prob;
    }

    // Rcout << "Round" << n_c << ptsa << std::endl;
    // Rcout << "Effprob, fut prob" << futile_prob << effective_prob << std::endl;
    // Rcout << "Stopping prob: " << stopping_prob << std::endl;
    // Rcout << "1-Nonstopping prob " << 1-nonstopping_prob << std::endl;
    prob_sizes_futile.insert(prob_sizes_futile.end(), futile_prob);
    prob_sizes_effective.insert(prob_sizes_effective.end(), effective_prob);
    prob_sizes_nonstop.insert(prob_sizes_nonstop.end(), nonstop_prob);

    Result res_list;
    res_list.effective_df = effective_df;
    res_list.effective_prob = effective_prob;
    res_list.cummu_effective = cummu_effective + effective_prob;
    res_list.nonstop_df = nonstop_df;
    res_list.expected_size = ptsa;
    res_list.stopping_prob = stopping_prob;
    res_list.nonstopping_prob = nonstopping_prob;
    res_list.prob_size_futile = wrap(prob_sizes_futile);
    res_list.prob_size_effective = wrap(prob_sizes_effective);
    res_list.prob_size_nonstop = wrap(prob_sizes_nonstop);

    return res_list;
  } else {
    Rcpp::NumericVector nobs_new(nobs.begin(), nobs.end() - 1);
    // nobs_new.erase(nobs_new.end() - 1);
    Result res_before = exacterror_recursive(nobs_new, ncut1, ncut2, test_prob, ntotal);
    Data effective_df = res_before.effective_df;
    Data nonstop_df = res_before.nonstop_df;
    cummu_effective = res_before.cummu_effective;
    ptsa = res_before.expected_size;
    stopping_prob = res_before.stopping_prob;
    nonstopping_prob = res_before.nonstopping_prob;
    NumericVector prob_sizes_futile = res_before.prob_size_futile;
    NumericVector prob_sizes_effective = res_before.prob_size_effective;
    NumericVector prob_sizes_nonstop = res_before.prob_size_nonstop;

    Data effective_df_new;
    for (int i = 0; i < nonstop_df.num_event.size(); i++) {
      std::vector<int> newevent;
      for (int j = ncut2[n_c - 1]; j <= nobs[n_c - 1]; j++) {
        int temp_newevent = j - nonstop_df.num_event[i];
        if (temp_newevent <= nobs[n_c - 1] - nobs[n_c - 2] && temp_newevent >= 0) {
          newevent.push_back(temp_newevent);
        }
      }
      // for (int j = ncut2[n_c - 1]; j <= nobs[n_c - 1]; j++) {
      //   int adjusted_x = j - nonstop_df.num_event[i];
      //   if (adjusted_x >= 0) {
      //     newevent.push_back(adjusted_x);
      //   } else {
      //     newevent.push_back(0);
      //   }
      // }

      std::sort(newevent.begin(), newevent.end());
      auto last = std::unique(newevent.begin(), newevent.end());
      newevent.erase(last, newevent.end());

      double prob = nonstop_df.prob[i];
      Data newdf;
      for (int k = 0; k < newevent.size(); k++) {
        int new_event_value = nonstop_df.num_event[i] + newevent[k];
        double dbinom_result = R::dbinom(newevent[k], nobs[n_c - 1] - nobs[n_c - 2], test_prob, false);

        newdf.num_event.push_back(new_event_value);
        newdf.prob.push_back(prob * dbinom_result);
      }

      effective_df_new.num_event.insert(effective_df_new.num_event.end(), newdf.num_event.begin(), newdf.num_event.end());
      effective_df_new.prob.insert(effective_df_new.prob.end(), newdf.prob.begin(), newdf.prob.end());
    }

    Data nonstop_df_new;
    for (int i = 0; i < nonstop_df.num_event.size(); i++) {
      std::vector<int> newevent;
      for (int j = ncut1[n_c-1] + 1; j < ncut2[n_c-1]; j++) {
        int temp_newevent = j - nonstop_df.num_event[i];
        if (temp_newevent <= nobs[n_c - 1] - nobs[n_c - 2] && temp_newevent >= 0) {
          newevent.push_back(temp_newevent);
        }
      }
      // for (int j = ncut1[n_c-1] + 1; j < ncut2[n_c-1]; j++) {
      //   int adjusted_x = j - nonstop_df.num_event[i];
      //   if (adjusted_x >= 0) {
      //     newevent.push_back(adjusted_x);
      //   } else {
      //     newevent.push_back(0);
      //   }
      // }

      std::sort(newevent.begin(), newevent.end());
      auto last = std::unique(newevent.begin(), newevent.end());
      newevent.erase(last, newevent.end());

      double prob = nonstop_df.prob[i];
      Data newdf;
      for (int k = 0; k < newevent.size(); k++) {
        int new_event_value = nonstop_df.num_event[i] + newevent[k];
        double dbinom_result = R::dbinom(newevent[k], nobs[n_c - 1] - nobs[n_c - 2], test_prob, false);

        newdf.num_event.push_back(new_event_value);
        newdf.prob.push_back(prob * dbinom_result);
      }

      nonstop_df_new.num_event.insert(nonstop_df_new.num_event.end(), newdf.num_event.begin(), newdf.num_event.end());
      nonstop_df_new.prob.insert(nonstop_df_new.prob.end(), newdf.prob.begin(), newdf.prob.end());
    }

    // aggregate non_df_new:
    // Create a hashmap to store num_event as keys and their corresponding index in nonstop_df_new_sum as values
    Data nonstop_df_new_sum;
    std::unordered_map<int, int> num_event_index_map_null;

    for (int i = 0; i < nonstop_df_new.num_event.size(); i++) {
      int current_num_event = nonstop_df_new.num_event[i];
      auto it = num_event_index_map_null.find(current_num_event);

      if (it != num_event_index_map_null.end()) {
        // If num_event is already in nonstop_df_new_sum, update the probability
        nonstop_df_new_sum.prob[it->second] += nonstop_df_new.prob[i];
      } else {
        // If num_event is not found, add it to nonstop_df_new_sum and update the hashmap
        nonstop_df_new_sum.num_event.push_back(current_num_event);
        nonstop_df_new_sum.prob.push_back(nonstop_df_new.prob[i]);
        num_event_index_map_null[current_num_event] = nonstop_df_new_sum.num_event.size() - 1;
      }
    }

    // futile probability:
    Data futile_df_new;
    for (int i = 0; i < nonstop_df.num_event.size(); i++) {
      std::vector<int> newevent;
      for (int j = 0; j <= ncut1[n_c-1]; j++) {
        int temp_newevent = j - nonstop_df.num_event[i];
        if (temp_newevent <= nobs[n_c - 1] - nobs[n_c - 2] && temp_newevent >= 0) {
          newevent.push_back(temp_newevent);
        }
      }
      // for (int j = 0; j <= ncut1[n_c-1]; j++) {
      //   int adjusted_x = j - nonstop_df.num_event[i];
      //   if (adjusted_x >= 0) {
      //     newevent.push_back(adjusted_x);
      //   } else {
      //     newevent.push_back(0);
      //   }
      // }

      std::sort(newevent.begin(), newevent.end());
      auto last = std::unique(newevent.begin(), newevent.end());
      newevent.erase(last, newevent.end());

      double prob = nonstop_df.prob[i];
      Data newdf;
      for (int k = 0; k < newevent.size(); k++) {
        int new_event_value = nonstop_df.num_event[i] + newevent[k];
        double dbinom_result = R::dbinom(newevent[k], nobs[n_c - 1] - nobs[n_c - 2], test_prob, false);

        newdf.num_event.push_back(new_event_value);
        newdf.prob.push_back(prob * dbinom_result);
      }

      futile_df_new.num_event.insert(futile_df_new.num_event.end(), newdf.num_event.begin(), newdf.num_event.end());
      futile_df_new.prob.insert(futile_df_new.prob.end(), newdf.prob.begin(), newdf.prob.end());
    }

    double futile_prob = std::accumulate(futile_df_new.prob.begin(), futile_df_new.prob.end(), 0.0);
    double effective_prob = std::accumulate(effective_df_new.prob.begin(), effective_df_new.prob.end(), 0.0);
    double nonstop_prob = std::accumulate(nonstop_df_new_sum.prob.begin(), nonstop_df_new_sum.prob.end(), 0.0);

    // double test_prob2 = std::accumulate(nonstop_df.prob.begin(), nonstop_df.prob.end(), 0.0);

    if (n_c == ntotal) {
      ptsa = ptsa + (futile_prob + effective_prob + nonstop_prob) * nobs[n_c - 1];
      // Rcout << "Final Round Prob. " << futile_prob + effective_prob + nonstop_prob << std::endl;
      // Rcout << "Final Round Prob. (plus) " << futile_prob + effective_prob + nonstop_prob << std::endl;
      // Rcout << "Final Round Prob. (nonstop) " << test_prob2 << std::endl;

    } else {
      ptsa = ptsa + (effective_prob + futile_prob) * nobs[n_c - 1];
      stopping_prob = stopping_prob + futile_prob + effective_prob;
      nonstopping_prob = nonstopping_prob + nonstop_prob;
    }

    // Rcout << "Round " << n_c << ptsa << std::endl;
    // Rcout << "Stopping prob: " << stopping_prob << std::endl;
    // Rcout << "1-Nonstopping prob: " << 1-nonstopping_prob << std::endl;

    cummu_effective = cummu_effective + effective_prob;
    prob_sizes_futile.insert(prob_sizes_futile.end(), futile_prob);
    prob_sizes_effective.insert(prob_sizes_effective.end(), effective_prob);
    prob_sizes_nonstop.insert(prob_sizes_nonstop.end(), nonstop_prob);

    Result res_list;
    res_list.effective_df = effective_df_new;
    res_list.effective_prob = effective_prob;
    res_list.cummu_effective = cummu_effective;
    res_list.nonstop_df = nonstop_df_new_sum;
    res_list.expected_size = ptsa;
    res_list.stopping_prob = stopping_prob;
    res_list.nonstopping_prob = nonstopping_prob;
    res_list.prob_size_futile = wrap(prob_sizes_futile);
    res_list.prob_size_effective = wrap(prob_sizes_effective);
    res_list.prob_size_nonstop = wrap(prob_sizes_nonstop);

    return res_list;
  }
}


// [[Rcpp::export]]
List exact_error_recursive_Rcpp(Rcpp::NumericVector nobs, Rcpp::NumericVector ncut1, Rcpp::NumericVector ncut2, double pnull, double palter, int ntotal) {
  // Call the C++ function
  Result result_null = exacterror_recursive(nobs, ncut1, ncut2, pnull, ntotal);
  Result result_alter = exacterror_recursive(nobs, ncut1, ncut2, palter, ntotal);

  // Convert Result to Rcpp List
  List result_list = List::create(
    Named("t1err") = result_null.cummu_effective,
    Named("power") = result_alter.cummu_effective,
    Named("pts_H0") = result_null.expected_size,
    Named("pts_Ha") = result_alter.expected_size,
    Named("prob_futile_null") = result_null.prob_size_futile,
    Named("prob_effective_null") = result_null.prob_size_effective,
    Named("prob_nonstop_null") = result_null.prob_size_nonstop,
    Named("prob_futile_alter") = result_alter.prob_size_futile,
    Named("prob_effective_alter") = result_alter.prob_size_effective,
    Named("prob_nonstop_alter") = result_alter.prob_size_nonstop
  );

  return result_list;
}


//Binary Efficacy
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List GetocBiRcpp_dual(int seed, double nsim, NumericMatrix contrast, NumericVector nobs, double b,
                      double b_grad1, double b_grad2, double pow1, double pow2, double pow3,
                      NumericVector dprior, double ptrue, double phi, double delta0, double delta1, Function fff)
{

  int n1 =nobs.size();
  NumericVector stopbound1(2);
  // std::map<std::string, int> stopbound1;
  NumericVector stopbound_fut(n1);
  NumericVector stopbound_grad(n1);
  double nmax=max(nobs);

  List resultslist(5);

  set_seed_dual(seed);

  for (int i=0; i<n1; i++)
  {
    double n=nobs[i];
    double lambda2 = b_grad1*pow(n/nmax,pow3); //cutoff1,cutoff2,cutoff2_lambda1,cutoff2_lambda2
    stopbound1 = fff(b*pow(n/nmax,pow1),pow(n/nmax,pow2),lambda2,b_grad2,dprior,n,phi,delta0,delta1,contrast);

    if (stopbound1[0] < 0 || stopbound1[1] < 0 ) { // check if there is -Inf bound
      resultslist[0]= 0; // edrop
      resultslist[1]= 0; // t1err
      resultslist[2]= 0; // power
      resultslist[3]= 1; // pts
      resultslist[4]= 1;
      return resultslist;
    }

    stopbound_fut[i]=stopbound1[0];
    stopbound_grad[i]=stopbound1[1];
  }

  for (int i=0; i<n1; i++){
    nobs[i] = std::ceil(nobs[i]);
  }
  // Rcpp::Rcout << "nobs " << nobs << std::endl;
  // Rcpp::Rcout << "Fut Boundary " << stopbound_fut << std::endl;
  // Rcpp::Rcout << "Grad Boundary " << stopbound_grad << std::endl;

  List temp(4);
  temp = exact_error_recursive_Rcpp(nobs, stopbound_fut, stopbound_grad, phi, ptrue, n1);

  resultslist[0]= 0; // edrop
  resultslist[1]= temp[0]; // t1err
  resultslist[2]= temp[1]; // power
  resultslist[3]= temp[2]; // pts
  resultslist[4]= temp[3]; // pts_Ha

  return resultslist;
}

