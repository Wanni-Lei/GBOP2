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



//
// // set seed
// // [[Rcpp::export]]
// void set_seed(unsigned int seed) {
//   Rcpp::Environment base_env("package:base");
//   Rcpp::Function set_seed_r = base_env["set.seed"];
//   set_seed_r(seed);
// }
//
// double dbinom_product( const NumericVector & vXs,  NumericVector & nobs, double prob);


// // [[Rcpp::export]]
// List Exacterror(NumericVector nobs, NumericVector ncut, double pnull, double palter)
// {
//   double dropn(0);
//   double dropa(0);
//   double prea(0);
//   double ptsa(0);
//
//   if (nobs.size() == 1){
//     dropn = R::pbinom(ncut[0],nobs[0],pnull,true,false);
//     dropa = R::pbinom(ncut[0],nobs[0],palter,true,false);
//     ptsa = nobs[0];
//   }else if(nobs.size() == 2){
//     for(int i=0;i<=nobs[0];i++){
//       if(i<=ncut[0]){
//           dropn = dropn + R::dbinom(i,nobs[0],pnull,false);
//           dropa = dropa + R::dbinom(i,nobs[0],palter,false);
//           prea = prea + R::dbinom(i,nobs[0],palter,false);
//           ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*nobs[0];
//       }else{
//         for(int j=0;j<=nobs[1]-nobs[0];j++){
//           if(i+j<=ncut[1]){
//             dropn = dropn + R::dbinom(i,nobs[0],pnull,false)*R::dbinom(j,nobs[1]-nobs[0],pnull,false);
//             dropa = dropa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false);
//             ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*nobs[1];
//           }else{
//             ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*nobs[1];
//           }
//         }
//       }
//     }
//   }else if(nobs.size() == 3){
//     for(int i=0;i<=nobs[0];i++){
//       if(i<=ncut[0]){
//         dropn = dropn + R::dbinom(i,nobs[0],pnull,false);
//         dropa = dropa + R::dbinom(i,nobs[0],palter,false);
//         prea = prea + R::dbinom(i,nobs[0],palter,false);
//         ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*nobs[0];
//       }else{
//         for(int j=0;j<=nobs[1]-nobs[0];j++){
//           if(i+j<=ncut[1]){
//             dropn = dropn + R::dbinom(i,nobs[0],pnull,false)*R::dbinom(j,nobs[1]-nobs[0],pnull,false);
//             dropa = dropa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false);
//             prea = prea + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false);
//             ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*nobs[1];
//           }else{
//             for(int k=0;k<=nobs[2]-nobs[1];k++){
//               if(i+j+k<=ncut[2]){
//                 dropn = dropn + R::dbinom(i,nobs[0],pnull,false)*R::dbinom(j,nobs[1]-nobs[0],pnull,false)*R::dbinom(k,nobs[2]-nobs[1],pnull,false);
//                 dropa = dropa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false);
//                 ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false)*nobs[2];
//               }else{
//                 ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false)*nobs[2];
//               }
//             }
//           }
//         }
//       }
//     }
//   }else if(nobs.size() == 4){
//     for(int i=0;i<=nobs[0];i++){
//       if(i<=ncut[0]){
//         dropn = dropn + R::dbinom(i,nobs[0],pnull,false);
//         dropa = dropa + R::dbinom(i,nobs[0],palter,false);
//         prea = prea + R::dbinom(i,nobs[0],palter,false);
//         ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*nobs[0];
//       }else{
//         for(int j=0;j<=nobs[1]-nobs[0];j++){
//           if(i+j<=ncut[1]){
//               dropn = dropn + R::dbinom(i,nobs[0],pnull,false)*R::dbinom(j,nobs[1]-nobs[0],pnull,false);
//               dropa = dropa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false);
//               prea = prea + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false);
//               ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*nobs[1];
//           }else{
//             for(int k=0;k<=nobs[2]-nobs[1];k++){
//               if(i+j+k<=ncut[2]){
//                 dropn = dropn + R::dbinom(i,nobs[0],pnull,false)*R::dbinom(j,nobs[1]-nobs[0],pnull,false)*R::dbinom(k,nobs[2]-nobs[1],pnull,false);
//                 dropa = dropa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false);
//                 prea = prea + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false);
//                 ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false)*nobs[2];
//               }else{
//                 for(int m=0;m<=nobs[3]-nobs[2];m++){
//                   if(i+j+k+m<=ncut[3]) {
//                     dropn = dropn + R::dbinom(i,nobs[0],pnull,false)*R::dbinom(j,nobs[1]-nobs[0],pnull,false)*R::dbinom(k,nobs[2]-nobs[1],pnull,false)*R::dbinom(m,nobs[3]-nobs[2],pnull,false);
//                     dropa = dropa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false)*R::dbinom(m,nobs[3]-nobs[2],palter,false);
//                     ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false)*R::dbinom(m,nobs[3]-nobs[2],palter,false)*nobs[3];
//                   }else{
//                     ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false)*R::dbinom(m,nobs[3]-nobs[2],palter,false)*nobs[3];
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }else if(nobs.size() == 5){
//     for(int i=0;i<=nobs[0];i++){
//       if(i<=ncut[0]){
//         dropn = dropn + R::dbinom(i,nobs[0],pnull,false);
//         dropa = dropa + R::dbinom(i,nobs[0],palter,false);
//         prea  = prea  + R::dbinom(i,nobs[0],palter,false);
//         ptsa  = ptsa  + R::dbinom(i,nobs[0],palter,false)*nobs[0];
//       }else{
//         for(int j=0;j<=nobs[1]-nobs[0];j++){
//           if(i+j<=ncut[1]){
//             //dropn = dropn + R::dbinom(i,nobs[0],pnull, false) * R::dbinom(j,nobs[1]-nobs[0],pnull, false);
//             dropn = dropn + dbinom_product( NumericVector::create(i,j), nobs, pnull );
//             dropa = dropa + dbinom_product( NumericVector::create(i,j), nobs, palter);
//             prea  = prea  + dbinom_product( NumericVector::create(i,j), nobs, palter);
//             ptsa  = ptsa  + dbinom_product( NumericVector::create(i,j), nobs, palter) * nobs[1];
//           }else{
//             for(int k=0;k<=nobs[2]-nobs[1];k++){
//               if(i+j+k<=ncut[2]){
//                  dropn = dropn + dbinom_product( NumericVector::create(i,j,k), nobs, pnull );
//                  dropa = dropa + dbinom_product( NumericVector::create(i,j,k), nobs, palter);
//                  prea  = prea  + dbinom_product( NumericVector::create(i,j,k), nobs, palter);
//                  ptsa  = ptsa  + dbinom_product( NumericVector::create(i,j,k), nobs, palter) * nobs[2];
//               }else{
//                 for(int m=0;m<=nobs[3]-nobs[2];m++){
//                   if(i+j+k+m<=ncut[3]){
//                      dropn = dropn + dbinom_product( NumericVector::create(i,j,k,m), nobs, pnull );
//                      dropa = dropa + dbinom_product( NumericVector::create(i,j,k,m), nobs, palter);
//                      prea  = prea  + dbinom_product( NumericVector::create(i,j,k,m), nobs, palter);
//                      ptsa  = ptsa  + dbinom_product( NumericVector::create(i,j,k,m), nobs, palter) * nobs[3];
//                   }else{
//                     for(int n=0;n<=nobs[4]-nobs[3];n++){
//                       if(i+j+k+m+n<=ncut[4]){
//                         dropn = dropn + R::dbinom(i,nobs[0],pnull ,false)*R::dbinom(j,nobs[1]-nobs[0],pnull ,false)*R::dbinom(k,nobs[2]-nobs[1],pnull ,false)*R::dbinom(m,nobs[3]-nobs[2],pnull ,false)*R::dbinom(n,nobs[4]-nobs[3],pnull ,false);
//                         dropa = dropa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false)*R::dbinom(m,nobs[3]-nobs[2],palter,false)*R::dbinom(n,nobs[4]-nobs[3],palter,false);
//                         ptsa  = ptsa  + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false)*R::dbinom(m,nobs[3]-nobs[2],palter,false)*R::dbinom(n,nobs[4]-nobs[3],palter,false)*nobs[4];
//                       }else{
//                         ptsa = ptsa + R::dbinom(i,nobs[0],palter,false)*R::dbinom(j,nobs[1]-nobs[0],palter,false)*R::dbinom(k,nobs[2]-nobs[1],palter,false)*R::dbinom(m,nobs[3]-nobs[2],palter,false)*R::dbinom(n,nobs[4]-nobs[3],palter,false)*nobs[4];
//                       }
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
//   List result = List::create(Named("edrop")=prea, Named("pts")=ptsa, Named("t1err")=1-dropn, Named("power")=1-dropa);
//   return result;
// }


// double dbinom_product(  const NumericVector & vXs,  NumericVector & nobs, double prob)
// {
//     int iMax = vXs.size();   // nobs.size() must be at least iMax also
//    double dResult = R::dbinom(vXs[0], nobs[0], prob, false);
//    for (int i = 1; i < iMax; ++i)
//    {
//        dResult *= R::dbinom(vXs[i], nobs[i] - nobs[i-1], prob, false);
//    }
//    return dResult;
// }










// [[Rcpp::export]]
DataFrame create_grid(int n) {
  int total_rows = (n + 1) * (n + 1) * (n + 1);
  NumericMatrix grid(total_rows, 3);

  int row = 0;
  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j <= n; ++j) {
      for (int k = 0; k <= n; ++k) {
        grid(row, 0) = i;
        grid(row, 1) = j;
        grid(row, 2) = k;
        row++;
      }
    }
  }

  // Convert NumericMatrix to DataFrame and assign column names
  DataFrame grid_df = DataFrame::create(
    Named("x") = grid(_, 0),
    Named("y") = grid(_, 1),
    Named("z") = grid(_, 2)
  );

  return grid_df;

}

// Helper function to calculate factorial
double factorial(int n) {
  if (n == 0) return 1;
  double result = 1;
  for (int i = 2; i <= n; ++i) {
    result *= i;
  }
  return result;
}

// Helper function to calculate multinomial coefficient
double multinomial_coeff(NumericVector x) {
  double n = sum(x);
  double coeff = factorial(n);
  for (int i = 0; i < x.size(); ++i) {
    coeff /= factorial(x[i]);
  }
  return coeff;
}

// Helper function to calculate the multinomial probability
double dmultinom(NumericVector x, NumericVector p) {
  double prob = multinomial_coeff(x);
  for (int i = 0; i < x.size(); ++i) {
    prob *= pow(p[i], x[i]);
  }
  return prob;
}

// [[Rcpp::export]]
DataFrame filter_and_calculate(DataFrame grid, int n) {
  // Extract the columns as IntegerVectors
  IntegerVector x = grid["x"];
  IntegerVector y = grid["y"];
  IntegerVector z = grid["z"];

  // Create vectors to hold filtered results
  std::vector<int> x_filtered, y_filtered, z_filtered, a_filtered;

  // Loop over the rows and apply the condition
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] + y[i] + z[i] <= n) {
      x_filtered.push_back(x[i]);
      y_filtered.push_back(y[i]);
      z_filtered.push_back(z[i]);
      a_filtered.push_back(n - (x[i] + y[i] + z[i]));
    }
  }

  // Convert filtered results to IntegerVectors
  IntegerVector x_res = wrap(x_filtered);
  IntegerVector y_res = wrap(y_filtered);
  IntegerVector z_res = wrap(z_filtered);
  IntegerVector a_res = wrap(a_filtered);

  // Return the filtered DataFrame with the new column 'a'
  return DataFrame::create(Named("x") = x_res,
                           Named("y") = y_res,
                           Named("z") = z_res,
                           Named("a") = a_res);
}

// Example usage

// [[Rcpp::export]]
DataFrame calculate_and_add_probabilities(DataFrame nobs_grid, NumericVector p) {
  int n = nobs_grid.nrows();
  NumericVector prob(n);

  // Extract columns as vectors
  NumericVector nobs1 = nobs_grid["x"];
  NumericVector nobs2 = nobs_grid["y"];
  NumericVector nobs3 = nobs_grid["z"];
  NumericVector nobs4 = nobs_grid["a"];

  for (int i = 0; i < n; ++i) {
    NumericVector nobs = NumericVector::create(nobs4[i], nobs2[i], nobs1[i], nobs3[i]);
    prob[i] = dmultinom(nobs, p);
  }

  // Add the probability vector as a new column to the DataFrame
  nobs_grid["prob"] = prob;

  return nobs_grid;
}





// [[Rcpp::export]]
double den_cpp(int x, int y, int n, NumericVector p) {
  // p is a vector (p00, p01, p10, p11)
  // p11 is patients having response and toxicity
  // p01 = Pt - p11
  // p10 = Pr - p11
  double Pr = p[2] + p[3];  // total prob of response

  double Pt = p[1] + p[3];  // total prob of toxicity

  double Pt_1 = p[3] / Pr;  // condition toxicity from response patients p(0|1)

  double Pt_0 = p[1] / (1 - Pr);  // condition toxicity from non-response patients p(0|0)

  // This is done by iterating over all possible ways to distribute the y toxicities
  // among the x response patients and the remaining n-x non-response patients.
  double dx = R::dbinom(x, n, Pr, false);
  double sumb = 0.0;

  for (int j = 0; j <= std::min(x, y); j++) {
    sumb += R::dbinom(j, x, Pt_1, false) * R::dbinom((y - j), (n - x), Pt_0, false);
  }

  double den = dx * sumb;
  return den;
}




// [[Rcpp::export]]
List create_tt(int s, int npt, NumericVector t_bound) {
  // Create the hashmap
  std::unordered_map<int, std::vector<std::pair<int, int>>> hashmap;

  // Generate the combinations and populate the hashmap
  for(int j = 0; j < t_bound[s-1]; ++j) {
    for(int i = 0; i <= npt; ++i) {
      int sum = i + j;
      if (sum < t_bound[s]){
        hashmap[sum].push_back(std::make_pair(j, i));
      }
    }
  }

  // Prepare the result to return to R
  List result;

  // Convert the hashmap to a list where each element is a list of pairs
  for (auto& elem : hashmap) {
    int key = elem.first;
    std::vector<std::pair<int, int>>& values = elem.second;

    // Convert each pair to a two-column matrix (or DataFrame)
    IntegerMatrix mat(values.size(), 2);
    for (size_t i = 0; i < values.size(); ++i) {
      mat(i, 0) = values[i].first;
      mat(i, 1) = values[i].second;
    }

    // Add the matrix to the result list, with the sum as the name
    result[std::to_string(key)] = mat;
    // result[key] = mat;
  }

  return result;
}


// [[Rcpp::export]]
List create_rr(int s, int npt, NumericVector r_bound, NumericVector r_interm) {
  // Create the hashmap
  std::unordered_map<int, std::vector<std::pair<int, int>>> hashmap;

  // Generate the combinations and populate the hashmap
  for(int j = r_bound[s-1] + 1; j <= r_interm[s-1]; ++j) {
    for(int i = 0; i <= npt; ++i) {
      int sum = i + j;
      hashmap[sum].push_back(std::make_pair(j, i));
    }
  }

  // Prepare the result to return to R
  List result;

  // Convert the hashmap to a list where each element is a list of pairs
  for (auto& elem : hashmap) {
    int key = elem.first;
    std::vector<std::pair<int, int>>& values = elem.second;

    // Convert each pair to a two-column matrix (or DataFrame)
    IntegerMatrix mat(values.size(), 2);
    for (size_t i = 0; i < values.size(); ++i) {
      mat(i, 0) = values[i].first;
      mat(i, 1) = values[i].second;
    }

    // Add the matrix to the result list, with the sum as the name
    result[std::to_string(key)] = mat;
  }

  return result;
}


// [[Rcpp::export]]
List efftox_recursive_optim(NumericVector interm, IntegerVector npt, NumericVector p,
                      NumericVector bound_eff, NumericVector bound_tox) {

  int n_c = interm.size();
  double ptsa = 0;

  if (n_c == 1) {

    int n = npt[n_c - 1];
    int n_obs = interm[n_c - 1];
    // Rcpp::Rcout << "n_stage1: " << n << std::endl;

    // get nonstop/go event counts:
    int b_eff = bound_eff[n_c - 1];
    int b_tox = bound_tox[n_c - 1];

    NumericMatrix grid(interm[n_c - 1] - b_eff, b_tox); // eff counts x tox counts

    for (int i = b_eff + 1; i <= n; ++i){ // marginal efficacy
      for (int j = 0; j <= b_tox - 1; ++j){ // marginal toxicity
        grid(i - (b_eff + 1), j) = den_cpp(i, j, n, p);
      }
    }

    double nonstop_prob = sum(grid);
    // Calculate ptsa
    ptsa = n;
    // Rcpp::Rcout << "pass stage 1 "  << std::endl;
    return List::create(Named("nonstop_grid") = grid,
                        Named("nonstop_prob") = nonstop_prob,
                        Named("ptsa") = ptsa);


  } else {
    Rcpp::NumericVector interm_new(interm.begin(), interm.end() - 1);
    // Rcpp::Rcout << "interm: " << interm_new << std::endl;
    List res_prev = efftox_recursive_optim(interm_new, npt, p, bound_eff, bound_tox);

    NumericMatrix nonstop_grid_prev = res_prev["nonstop_grid"];

    double ptsa = res_prev["ptsa"];

    int n = npt[n_c - 1];
    int n_obs = interm[n_c - 1];
    int n_obs_prev = interm[n_c - 2];

    // get nonstop/go event counts:
    int b_eff = bound_eff[n_c - 1];
    int b_eff_prev = bound_eff[n_c - 2];
    int b_tox = bound_tox[n_c - 1];

    List hash_rr = create_rr(n_c - 1, n, bound_eff, interm);
    List hash_tt = create_tt(n_c - 1, n, bound_tox);
    // Rcpp::Rcout << "pass hashmap, stage" << "n_c"  << std::endl;

    NumericMatrix grid(n_obs - b_eff, b_tox); // eff counts x tox counts
    IntegerMatrix rc, tc;
    int a1, a2, b1, b2;

    // Rcpp::Rcout << "k_com size: " << (n_obs - b_eff) * b_tox << std::endl;
    // Rcpp::Rcout << "nonstop_grid_prev nrow: " << nonstop_grid_prev.nrow() << std::endl;
    // Rcpp::Rcout << "nonstop_grid_prev ncol: " << nonstop_grid_prev.ncol() << std::endl;

    for (int r = b_eff + 1; r <= n_obs; r++ ){ // efficacy count for go decision
      for (int t = 0; t <= b_tox - 1; t ++){ // toxic count for go decision

        // Check if the key exists in hash_tt
        std::string r_str = std::to_string(r);
        std::string t_str = std::to_string(t);
        if (!hash_rr.containsElementNamed(r_str.c_str()) || !hash_tt.containsElementNamed(t_str.c_str())) {
          // Rcpp::Rcout << "Key not found in hash_tt: t = " << t << " (key = '" << t_str << "')" << std::endl;
          continue; // Skip this iteration if the key is missing
        }

        rc = as<IntegerMatrix>(hash_rr[std::to_string(r)]);
        tc = as<IntegerMatrix>(hash_tt[std::to_string(t)]);

        // Rcpp::Rcout << "r, t: " << NumericVector::create(r,t) << std::endl;
        // Rcpp::Rcout << "rc: " << rc << std::endl;
        // Rcpp::Rcout << "tc: " << tc << std::endl;
        double temp = 0;
        for (int i = 0; i < rc.nrow(); i++){
          for (int j = 0; j < tc.nrow(); j++){
            a1 = rc(i,0) - (b_eff_prev + 1);
            a2 = rc(i,1);
            b1 = tc(j,0);
            b2 = tc(j,1);
            // Rcpp::Rcout << "a,b: " << NumericVector::create(a1,b1,a2,b2) << std::endl;

            temp += nonstop_grid_prev(a1, b1) * den_cpp(a2, b2, n, p);
            // Rcpp::Rcout << "nonstop_prev_prob " << nonstop_grid_prev(a1, b1) << std::endl;
            // Rcpp::Rcout << "density of this term " << den_cpp(a2, b2, n, p) << std::endl;
          }
        }


        grid(r - (b_eff + 1), t) = temp;
      }
    }

    // if (n_c == 4){
    //   Rcpp::Rcout << "nonstop_prev_prob " << grid << std::endl;
    // }

    // Calculate the sum of probabilities in nonstop_grid
    double nonstop_prob = sum(grid);
    // Calculate ptsa
    double sum_prob_prev = sum(nonstop_grid_prev);
    ptsa = ptsa + n * sum_prob_prev;

    // Rcpp::Rcout << "pass stage" << n_c << std::endl;
    // Rcpp::Rcout << "nobs_grid:" << nonstop_prob << std::endl;
    // Rcpp::Rcout << "ptsa:" << ptsa << std::endl;

    return List::create(Named("nonstop_grid") = grid,
                        Named("nonstop_prob") = nonstop_prob,
                        Named("ptsa") = ptsa);
  }
}

// /*** R
// n <- 3
// grid <- create_grid(n)
// nobs_grid <- filter_and_calculate(grid, n)
// nobs_grid
// # calculate_and_add_probabilities(nobs_grid, p)
// */
//



// // Define a struct to hold the data
// struct Data {
//   std::vector<int> num_event;
//   std::vector<double> prob;
// };
//
// // Define a struct to hold the result
// struct Result {
//   Data nonstopping_df_null;
//   Data nonstopping_df_alter;
//   double nonstop_prob_null;
//   double nonstop_prob_alter;
//   double expected_size;
// };
//
// // Function to calculate exact error recursively
// Result exact_error_recursive2(Rcpp::NumericVector nobs, Rcpp::NumericVector ncut, double pnull, double palter, int ntotal) {
//   int n_c = nobs.size();
//   double ptsa = 0;
//
//   if (n_c == 1) {
//     // pnull
//     Data nonstopping_df_null;
//     std::vector<double> probabilities_null;
//     for (int i = ncut[0] + 1; i <= nobs[0]; i++) {
//       probabilities_null.push_back(R::dbinom(i, nobs[0], pnull, false));
//     }
//     nonstopping_df_null.prob = probabilities_null;
//     std::vector<int> mySet_null;
//     for (int i = ncut[0] + 1; i <= nobs[0]; i++) {
//       mySet_null.push_back(i);
//     }
//     nonstopping_df_null.num_event = mySet_null;
//
//
//     // palter
//     Data nonstopping_df_alter;
//     std::vector<double> probabilities_alter;
//     for (int i = ncut[0] + 1; i <= nobs[0]; i++) {
//       probabilities_alter.push_back(R::dbinom(i, nobs[0], palter, false));
//     }
//     nonstopping_df_alter.prob = probabilities_alter;
//     std::vector<int> mySet_alter;
//     for (int i = ncut[0] + 1; i <= nobs[0]; i++) {
//       mySet_alter.push_back(i);
//     }
//     nonstopping_df_alter.num_event = mySet_alter;
//
//     // expected sample size under the H0:
//     Data stopping_df_null;
//     std::vector<double> probabilities_null2;
//     for (int i = 0; i <= ncut[0]; i++) {
//       probabilities_null2.push_back(R::dbinom(i, nobs[0], pnull, false));
//     }
//     stopping_df_null.prob = probabilities_null2;
//     std::vector<int> mySet_null2;
//     for (int i = 0; i <= ncut[0]; i++) {
//       mySet_null2.push_back(i);
//     }
//     stopping_df_null.num_event = mySet_null2;
//
//     double nonstop_prob_null = std::accumulate(nonstopping_df_null.prob.begin(), nonstopping_df_null.prob.end(), 0.0);
//     double nonstop_prob_alter = std::accumulate(nonstopping_df_alter.prob.begin(), nonstopping_df_alter.prob.end(), 0.0);
//     double stop_prob_null = std::accumulate(stopping_df_null.prob.begin(), stopping_df_null.prob.end(), 0.0);
//     ptsa = ptsa + stop_prob_null * nobs[0];
//     if (n_c == ntotal) {
//       ptsa = ptsa + nonstop_prob_null * nobs[0];
//     }
//
//     Result res_list;
//     res_list.nonstopping_df_null = nonstopping_df_null;
//     res_list.nonstopping_df_alter = nonstopping_df_alter;
//     res_list.nonstop_prob_null = nonstop_prob_null;
//     res_list.nonstop_prob_alter = nonstop_prob_alter;
//     res_list.expected_size = ptsa;
//
//     return res_list;
//   } else {
//     Rcpp::NumericVector nobs_new(nobs.begin(), nobs.end() - 1);
//     Result res_before = exact_error_recursive2(nobs_new, ncut, pnull, palter, ntotal);
//     Data nonstop_null = res_before.nonstopping_df_null;
//     Data nonstop_alter = res_before.nonstopping_df_alter;
//     ptsa = res_before.expected_size;
//
//     // nonstopping prob
//     // pnull
//     Data nonstop_df_new_null;
//     // std::map<int, double> sum_map_null;
//     for (int i = 0; i < nonstop_null.num_event.size(); i++) {
//       std::vector<int> newevent;
//       for (int x = ncut[n_c-1] + 1; x <= nobs[n_c-1]; x++) {
//         int adjusted_x = x - nonstop_null.num_event[i];
//         if (adjusted_x >= 0) {
//           newevent.push_back(adjusted_x); //
//         } else {
//           newevent.push_back(0);
//         }
//       }
//
//       // find unique values in newevent:
//       std::sort(newevent.begin(), newevent.end());
//       auto last = std::unique(newevent.begin(), newevent.end());
//       newevent.erase(last, newevent.end());
//
//       double prob = nonstop_null.prob[i];
//       int nobs_diff = nobs[n_c - 1] - nobs[n_c - 2];
//
//       Data newdf;
//       for (int j = 0; j < newevent.size(); j++) {
//         int new_event_value = nonstop_null.num_event[i] + newevent[j];
//         double dbinom_result = R::dbinom(newevent[j], nobs_diff, pnull, false);
//
//         newdf.num_event.push_back(new_event_value);
//         newdf.prob.push_back(prob * dbinom_result);
//
//         if (n_c == ntotal) {
//           ptsa = ptsa + nobs[n_c - 1] * newdf.prob[j];
//         }
//       }
//       nonstop_df_new_null.num_event.insert(nonstop_df_new_null.num_event.end(), newdf.num_event.begin(), newdf.num_event.end());
//       nonstop_df_new_null.prob.insert(nonstop_df_new_null.prob.end(), newdf.prob.begin(), newdf.prob.end());
//     }
//
//     Data nonstop_df_new_null_sum;
//
//     // Create a hashmap to store num_event as keys and their corresponding index in nonstop_df_new_sum as values
//     std::unordered_map<int, int> num_event_index_map_null;
//
//     for (int i = 0; i < nonstop_df_new_null.num_event.size(); i++) {
//       int current_num_event = nonstop_df_new_null.num_event[i];
//       auto it = num_event_index_map_null.find(current_num_event);
//
//       if (it != num_event_index_map_null.end()) {
//         // If num_event is already in nonstop_df_new_sum, update the probability
//         nonstop_df_new_null_sum.prob[it->second] += nonstop_df_new_null.prob[i];
//       } else {
//         // If num_event is not found, add it to nonstop_df_new_sum and update the hashmap
//         nonstop_df_new_null_sum.num_event.push_back(current_num_event);
//         nonstop_df_new_null_sum.prob.push_back(nonstop_df_new_null.prob[i]);
//         num_event_index_map_null[current_num_event] = nonstop_df_new_null_sum.num_event.size() - 1;
//       }
//     }
//
//     // nonstopping * stopping:
//     // Data nonstop_df_new_stop;
//     for (int i = 0; i < nonstop_null.num_event.size(); i++) {
//       std::vector<int> newevent;
//       for (int x = ncut[n_c - 2] + 1; x <= ncut[n_c-1]; x++) {
//         int adjusted_x = x - nonstop_null.num_event[i];
//         if (adjusted_x < 0) {
//           newevent.push_back(-1);
//         } else {
//           newevent.push_back(adjusted_x);
//         }
//       }
//
//       std::sort(newevent.begin(), newevent.end());
//       auto last = std::unique(newevent.begin(), newevent.end());
//       newevent.erase(last, newevent.end());
//
//       double prob = nonstop_null.prob[i];
//       int nobs_diff = nobs[n_c - 1] - nobs[n_c - 2];
//
//       Data newdf;
//       for (int j = 0; j < newevent.size(); j++) {
//         int new_event_value = nonstop_null.num_event[i] + newevent[j];
//         double dbinom_result = R::dbinom(newevent[j], nobs_diff, pnull, false);
//
//         newdf.num_event.push_back(new_event_value);
//         newdf.prob.push_back(prob * dbinom_result);
//
//         // Aggregate sum of probabilities
//         ptsa = ptsa + nobs[n_c - 1] * newdf.prob[j];
//       }
//     }
//
//     // palter
//     Data nonstop_df_new_alter;
//     // std::map<int, double> sum_map_alter;
//     for (int i = 0; i < nonstop_alter.num_event.size(); i++) {
//       std::vector<int> newevent;
//       for (int x = ncut[n_c - 1] + 1; x <= nobs[n_c - 1]; x++) {
//         int adjusted_x = x - nonstop_alter.num_event[i];
//         if (adjusted_x >= 0) {
//           newevent.push_back(adjusted_x);
//         } else {
//           newevent.push_back(0);
//         }
//       }
//
//       std::sort(newevent.begin(), newevent.end());
//       auto last = std::unique(newevent.begin(), newevent.end());
//       newevent.erase(last, newevent.end());
//
//       double prob = nonstop_alter.prob[i];
//       int nobs_diff = nobs[n_c - 1] - nobs[n_c - 2];
//
//       Data newdf;
//       for (int j = 0; j < newevent.size(); j++) {
//         int new_event_value = nonstop_alter.num_event[i] + newevent[j];
//         double dbinom_result = R::dbinom(newevent[j], nobs_diff, palter, false);
//
//         newdf.num_event.push_back(new_event_value);
//         newdf.prob.push_back(prob * dbinom_result);
//
//         // if (n_c == ntotal) {
//         //   ptsa = ptsa + nobs[n_c - 1] * newdf.prob[j];
//         // }
//       }
//       nonstop_df_new_alter.num_event.insert(nonstop_df_new_alter.num_event.end(), newdf.num_event.begin(), newdf.num_event.end());
//       nonstop_df_new_alter.prob.insert(nonstop_df_new_alter.prob.end(), newdf.prob.begin(), newdf.prob.end());
//     }
//
//     // // Convert map to vectors
//     // for (const auto& pair : sum_map_alter) {
//     //   nonstop_df_new_alter.num_event.push_back(pair.first);
//     //   nonstop_df_new_alter.prob.push_back(pair.second);
//     // }
//     Data nonstop_df_new_alter_sum;
//
//     // Create a hashmap to store num_event as keys and their corresponding index in nonstop_df_new_sum as values
//     std::unordered_map<int, int> num_event_index_map_alter;
//
//     for (int i = 0; i < nonstop_df_new_alter.num_event.size(); i++) {
//       int current_num_event = nonstop_df_new_alter.num_event[i];
//       auto it = num_event_index_map_alter.find(current_num_event);
//
//       if (it != num_event_index_map_alter.end()) {
//         // If num_event is already in nonstop_df_new_sum, update the probability
//         nonstop_df_new_alter_sum.prob[it->second] += nonstop_df_new_alter.prob[i];
//       } else {
//         // If num_event is not found, add it to nonstop_df_new_sum and update the hashmap
//         nonstop_df_new_alter_sum.num_event.push_back(current_num_event);
//         nonstop_df_new_alter_sum.prob.push_back(nonstop_df_new_alter.prob[i]);
//         num_event_index_map_alter[current_num_event] = nonstop_df_new_alter_sum.num_event.size() - 1;
//       }
//     }
//
//
//
//     // nonstop_df_new_stop.num_event.insert(nonstop_df_new_stop.num_event.end(), newdf.num_event.begin(), newdf.num_event.end());
//     // nonstop_df_new_stop.prob.insert(nonstop_df_new_stop.prob.end(), newdf.prob.begin(), newdf.prob.end());
//
//
//
//     double nonstop_prob_null = std::accumulate(nonstop_df_new_null_sum.prob.begin(), nonstop_df_new_null_sum.prob.end(), 0.0);
//     double nonstop_prob_alter = std::accumulate(nonstop_df_new_alter_sum.prob.begin(), nonstop_df_new_alter_sum.prob.end(), 0.0);
//
//     // ptsa = ptsa + nonstop_prob_alter * nobs[n_c];
//
//     Result res_list;
//     res_list.nonstopping_df_null = nonstop_df_new_null_sum;
//     res_list.nonstopping_df_alter = nonstop_df_new_alter_sum;
//     res_list.nonstop_prob_null = nonstop_prob_null;
//     res_list.nonstop_prob_alter = nonstop_prob_alter;
//     res_list.expected_size = ptsa;
//
//     return res_list;
//   }
// }
//
// // [[Rcpp::export]]
// List exact_error_recursive2_Rcpp(Rcpp::NumericVector nobs, Rcpp::NumericVector ncut, double pnull, double palter, int ntotal) {
//   // Call the C++ function
//   Result result = exact_error_recursive2(nobs, ncut, pnull, palter, ntotal);
//
//   // Convert Result to Rcpp List
//   List result_list = List::create(
//     // Named("nonstopping_df_null_num_event") = wrap(result.nonstopping_df_null.num_event),
//     // Named("nonstopping_df_null_prob") = wrap(result.nonstopping_df_null.prob),
//     // Named("nonstopping_df_alter_num_event") = wrap(result.nonstopping_df_alter.num_event),
//     // Named("nonstopping_df_alter_prob") = wrap(result.nonstopping_df_alter.prob),
//     Named("t1err") = result.nonstop_prob_null,
//     Named("power") = result.nonstop_prob_alter,
//     Named("pts") = result.expected_size
//   );
//
//   return result_list;
// }



