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
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

double dbinom_product( const NumericVector & vXs,  NumericVector & nobs, double prob);
double my_dbinom(double x, double size, double prob);
double my_pbinom(double q, double size, double prob);


// [[Rcpp::export]]
List Exacterror(NumericVector nobs, NumericVector ncut, double pnull, double palter)
{
  double dropn(0);
  double dropa(0);
  double prea(0);
  double ptsa(0);
  
  if (nobs.size() == 1){
    dropn = my_pbinom(ncut[0],nobs[0],pnull);
    dropa = my_pbinom(ncut[0],nobs[0],palter);
    ptsa = nobs[0];
  }else if(nobs.size() == 2){
    for(int i=0;i<=nobs[0];i++){
      if(i<=ncut[0]){
        dropn = dropn + my_dbinom(i,nobs[0],pnull);
        dropa = dropa + my_dbinom(i,nobs[0],palter);
        prea = prea + my_dbinom(i,nobs[0],palter);
        ptsa = ptsa + my_dbinom(i,nobs[0],palter)*nobs[0];
      }else{
        for(int j=0;j<=nobs[1]-nobs[0];j++){
          if(i+j<=ncut[1]){
            dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull);
            dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter);
            ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*nobs[1];
          }else{
            ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*nobs[1];
          }
        }
      }
    }
  }else if(nobs.size() == 3){
    for(int i=0;i<=nobs[0];i++){
      if(i<=ncut[0]){
        dropn = dropn + my_dbinom(i,nobs[0],pnull);
        dropa = dropa + my_dbinom(i,nobs[0],palter);
        prea = prea + my_dbinom(i,nobs[0],palter);
        ptsa = ptsa + my_dbinom(i,nobs[0],palter)*nobs[0];
      }else{
        for(int j=0;j<=nobs[1]-nobs[0];j++){
          if(i+j<=ncut[1]){
            dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull);
            dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter);
            prea = prea + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter);
            ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*nobs[1];
          }else{
            for(int k=0;k<=nobs[2]-nobs[1];k++){
              if(i+j+k<=ncut[2]){
                dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull)*my_dbinom(k,nobs[2]-nobs[1],pnull);
                dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter);
                ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*nobs[2];
              }else{
                ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*nobs[2];
              }
            }
          }
        }
      }
    }
  }else if(nobs.size() == 4){
    for(int i=0;i<=nobs[0];i++){
      if(i<=ncut[0]){
        dropn = dropn + my_dbinom(i,nobs[0],pnull);
        dropa = dropa + my_dbinom(i,nobs[0],palter);
        prea = prea + my_dbinom(i,nobs[0],palter);
        ptsa = ptsa + my_dbinom(i,nobs[0],palter)*nobs[0];
      }else{
        for(int j=0;j<=nobs[1]-nobs[0];j++){
          if(i+j<=ncut[1]){
            dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull);
            dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter);
            prea = prea + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter);
            ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*nobs[1];
          }else{
            for(int k=0;k<=nobs[2]-nobs[1];k++){
              if(i+j+k<=ncut[2]){
                dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull)*my_dbinom(k,nobs[2]-nobs[1],pnull);
                dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter);
                prea = prea + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter);
                ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*nobs[2];
              }else{
                for(int m=0;m<=nobs[3]-nobs[2];m++){
                  if(i+j+k+m<=ncut[3]) {
                    dropn = dropn + my_dbinom(i,nobs[0],pnull)*my_dbinom(j,nobs[1]-nobs[0],pnull)*my_dbinom(k,nobs[2]-nobs[1],pnull)*my_dbinom(m,nobs[3]-nobs[2],pnull);
                    dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter);
                    ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter)*nobs[3];
                  }else{
                    ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter)*nobs[3];
                  }
                }
              }
            }
          }
        }
      }
    }
  }else if(nobs.size() == 5){
    for(int i=0;i<=nobs[0];i++){
      if(i<=ncut[0]){
        dropn = dropn + my_dbinom(i,nobs[0],pnull);
        dropa = dropa + my_dbinom(i,nobs[0],palter);
        prea  = prea  + my_dbinom(i,nobs[0],palter);
        ptsa  = ptsa  + my_dbinom(i,nobs[0],palter)*nobs[0];
      }else{
        for(int j=0;j<=nobs[1]-nobs[0];j++){
          if(i+j<=ncut[1]){
            //dropn = dropn + my_dbinom(i,nobs[0],pnull, false) * my_dbinom(j,nobs[1]-nobs[0],pnull, false);
            dropn = dropn + dbinom_product( NumericVector::create(i,j), nobs, pnull );
            dropa = dropa + dbinom_product( NumericVector::create(i,j), nobs, palter);
            prea  = prea  + dbinom_product( NumericVector::create(i,j), nobs, palter);
            ptsa  = ptsa  + dbinom_product( NumericVector::create(i,j), nobs, palter) * nobs[1];
          }else{
            for(int k=0;k<=nobs[2]-nobs[1];k++){
              if(i+j+k<=ncut[2]){
                dropn = dropn + dbinom_product( NumericVector::create(i,j,k), nobs, pnull );
                dropa = dropa + dbinom_product( NumericVector::create(i,j,k), nobs, palter);
                prea  = prea  + dbinom_product( NumericVector::create(i,j,k), nobs, palter);
                ptsa  = ptsa  + dbinom_product( NumericVector::create(i,j,k), nobs, palter) * nobs[2];
              }else{
                for(int m=0;m<=nobs[3]-nobs[2];m++){
                  if(i+j+k+m<=ncut[3]){
                    dropn = dropn + dbinom_product( NumericVector::create(i,j,k,m), nobs, pnull );
                    dropa = dropa + dbinom_product( NumericVector::create(i,j,k,m), nobs, palter);
                    prea  = prea  + dbinom_product( NumericVector::create(i,j,k,m), nobs, palter);
                    ptsa  = ptsa  + dbinom_product( NumericVector::create(i,j,k,m), nobs, palter) * nobs[3];
                  }else{
                    for(int n=0;n<=nobs[4]-nobs[3];n++){
                      if(i+j+k+m+n<=ncut[4]){
                        dropn = dropn + my_dbinom(i,nobs[0],pnull )*my_dbinom(j,nobs[1]-nobs[0],pnull )*my_dbinom(k,nobs[2]-nobs[1],pnull )*my_dbinom(m,nobs[3]-nobs[2],pnull )*my_dbinom(n,nobs[4]-nobs[3],pnull );
                        dropa = dropa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter)*my_dbinom(n,nobs[4]-nobs[3],palter);
                        ptsa  = ptsa  + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter)*my_dbinom(n,nobs[4]-nobs[3],palter)*nobs[4];
                      }else{
                        ptsa = ptsa + my_dbinom(i,nobs[0],palter)*my_dbinom(j,nobs[1]-nobs[0],palter)*my_dbinom(k,nobs[2]-nobs[1],palter)*my_dbinom(m,nobs[3]-nobs[2],palter)*my_dbinom(n,nobs[4]-nobs[3],palter)*nobs[4];
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  List result = List::create(Named("edrop")=prea, Named("pts")=ptsa, Named("t1err")=1-dropn, Named("power")=1-dropa);
  return result;
}

// Define a struct to hold the data
struct Data {
  std::vector<int> num_event;
  std::vector<double> prob;
};

// Define a struct to hold the result
struct Result {
  Data nonstopping_df_null;
  Data nonstopping_df_alter;
  double nonstop_prob_null;
  double nonstop_prob_alter;
  double expected_size;
};

// Function to calculate exact error recursively
Result exact_error_recursive2(Rcpp::NumericVector nobs, Rcpp::NumericVector ncut, double pnull, int ntotal) {
  int n_c = nobs.size();
  double ptsa = 0;
  
  if (n_c == 1) {
    // pnull
    Data nonstopping_df_null;
    std::vector<double> probabilities_null;
    for (int i = ncut[0] + 1; i <= nobs[0]; i++) {
      probabilities_null.push_back(R::dbinom(i, nobs[0], pnull, false));
    }
    nonstopping_df_null.prob = probabilities_null;
    std::vector<int> mySet_null;
    for (int i = ncut[0] + 1; i <= nobs[0]; i++) {
      mySet_null.push_back(i);
    }
    nonstopping_df_null.num_event = mySet_null;
    
    // expected sample size under the H0:
    Data stopping_df_null;
    std::vector<double> probabilities_null2;
    for (int i = 0; i <= ncut[0]; i++) {
      probabilities_null2.push_back(R::dbinom(i, nobs[0], pnull, false));
    }
    stopping_df_null.prob = probabilities_null2;
    std::vector<int> mySet_null2;
    for (int i = 0; i <= ncut[0]; i++) {
      mySet_null2.push_back(i);
    }
    stopping_df_null.num_event = mySet_null2;
    
    double nonstop_prob_null = std::accumulate(nonstopping_df_null.prob.begin(), nonstopping_df_null.prob.end(), 0.0);
    double stop_prob_null = std::accumulate(stopping_df_null.prob.begin(), stopping_df_null.prob.end(), 0.0);
    ptsa = ptsa + stop_prob_null * nobs[0];
    if (n_c == ntotal) {
      ptsa = ptsa + nonstop_prob_null * nobs[0];
    }
    
    Result res_list;
    res_list.nonstopping_df_null = nonstopping_df_null;
    res_list.nonstop_prob_null = nonstop_prob_null;
    res_list.expected_size = ptsa;
    
    return res_list;
  } else {
    Rcpp::NumericVector nobs_new(nobs.begin(), nobs.end() - 1);
    Result res_before = exact_error_recursive2(nobs_new, ncut, pnull, ntotal);
    Data nonstop_null = res_before.nonstopping_df_null;
    ptsa = res_before.expected_size;
    
    // nonstopping prob
    // pnull
    Data nonstop_df_new_null;
    // std::map<int, double> sum_map_null;
    for (int i = 0; i < nonstop_null.num_event.size(); i++) {
      std::vector<int> newevent;
      for (int x = ncut[n_c-1] + 1; x <= nobs[n_c-1]; x++) {
        int adjusted_x = x - nonstop_null.num_event[i];
        if (adjusted_x >= 0) {
          newevent.push_back(adjusted_x); //
        } else {
          newevent.push_back(0);
        }
      }
      
      // find unique values in newevent:
      std::sort(newevent.begin(), newevent.end());
      auto last = std::unique(newevent.begin(), newevent.end());
      newevent.erase(last, newevent.end());
      
      double prob = nonstop_null.prob[i];
      int nobs_diff = nobs[n_c - 1] - nobs[n_c - 2];
      
      Data newdf;
      for (int j = 0; j < newevent.size(); j++) {
        int new_event_value = nonstop_null.num_event[i] + newevent[j];
        double dbinom_result = R::dbinom(newevent[j], nobs_diff, pnull, false);
        
        newdf.num_event.push_back(new_event_value);
        newdf.prob.push_back(prob * dbinom_result);
        
        if (n_c == ntotal) {
          ptsa = ptsa + nobs[n_c - 1] * newdf.prob[j];
        }
      }
      nonstop_df_new_null.num_event.insert(nonstop_df_new_null.num_event.end(), newdf.num_event.begin(), newdf.num_event.end());
      nonstop_df_new_null.prob.insert(nonstop_df_new_null.prob.end(), newdf.prob.begin(), newdf.prob.end());
    }
    
    Data nonstop_df_new_null_sum;
    
    // Create a hashmap to store num_event as keys and their corresponding index in nonstop_df_new_sum as values
    std::unordered_map<int, int> num_event_index_map_null;
    
    for (int i = 0; i < nonstop_df_new_null.num_event.size(); i++) {
      int current_num_event = nonstop_df_new_null.num_event[i];
      auto it = num_event_index_map_null.find(current_num_event);
      
      if (it != num_event_index_map_null.end()) {
        // If num_event is already in nonstop_df_new_sum, update the probability
        nonstop_df_new_null_sum.prob[it->second] += nonstop_df_new_null.prob[i];
      } else {
        // If num_event is not found, add it to nonstop_df_new_sum and update the hashmap
        nonstop_df_new_null_sum.num_event.push_back(current_num_event);
        nonstop_df_new_null_sum.prob.push_back(nonstop_df_new_null.prob[i]);
        num_event_index_map_null[current_num_event] = nonstop_df_new_null_sum.num_event.size() - 1;
      }
    }
    
    // nonstopping * stopping:
    // Data nonstop_df_new_stop;
    for (int i = 0; i < nonstop_null.num_event.size(); i++) {
      std::vector<int> newevent;
      for (int x = ncut[n_c - 2] + 1; x <= ncut[n_c-1]; x++) {
        int adjusted_x = x - nonstop_null.num_event[i];
        if (adjusted_x < 0) {
          newevent.push_back(-1);
        } else {
          newevent.push_back(adjusted_x);
        }
      }
      
      std::sort(newevent.begin(), newevent.end());
      auto last = std::unique(newevent.begin(), newevent.end());
      newevent.erase(last, newevent.end());
      
      double prob = nonstop_null.prob[i];
      int nobs_diff = nobs[n_c - 1] - nobs[n_c - 2];
      
      Data newdf;
      for (int j = 0; j < newevent.size(); j++) {
        int new_event_value = nonstop_null.num_event[i] + newevent[j];
        double dbinom_result = R::dbinom(newevent[j], nobs_diff, pnull, false);
        
        newdf.num_event.push_back(new_event_value);
        newdf.prob.push_back(prob * dbinom_result);
        
        // Aggregate sum of probabilities
        ptsa = ptsa + nobs[n_c - 1] * newdf.prob[j];
      }
    }
    
    double nonstop_prob_null = std::accumulate(nonstop_df_new_null_sum.prob.begin(), nonstop_df_new_null_sum.prob.end(), 0.0);
    
    // ptsa = ptsa + nonstop_prob_alter * nobs[n_c];
    
    Result res_list;
    res_list.nonstopping_df_null = nonstop_df_new_null_sum;
    res_list.nonstop_prob_null = nonstop_prob_null;
    res_list.expected_size = ptsa;
    
    return res_list;
  }
}

// [[Rcpp::export]]
List exact_error_recursive2_Rcpp(Rcpp::NumericVector nobs, Rcpp::NumericVector ncut, double pnull, double palter, int ntotal) {
  // Call the C++ function
  Result result = exact_error_recursive2(nobs, ncut, pnull, ntotal);
  Result result_H1 = exact_error_recursive2(nobs, ncut, palter, ntotal);
  
  // Convert Result to Rcpp List
  List result_list = List::create(
    // Named("nonstopping_df_null_num_event") = wrap(result.nonstopping_df_null.num_event),
    // Named("nonstopping_df_null_prob") = wrap(result.nonstopping_df_null.prob),
    // Named("nonstopping_df_alter_num_event") = wrap(result.nonstopping_df_alter.num_event),
    // Named("nonstopping_df_alter_prob") = wrap(result.nonstopping_df_alter.prob),
    Named("t1err") = result.nonstop_prob_null,
    Named("power") = result_H1.nonstop_prob_null,
    Named("pts") = result.expected_size,
    Named("pts_H1") = result_H1.expected_size
  );
  
  return result_list;
}



// [[Rcpp::export]]
double my_dbinom(double x, double size, double prob) {
  if (x < 0 || x > size) {
    return 0;
  }
  return R::choose(size, x) * std::pow(prob, x) * std::pow(1 - prob, size - x);
}

// [[Rcpp::export]]
double my_pbinom(double q, double size, double prob) {
  if (q < 0) {
    return 0;
  }
  if (q >= size) {
    return 1;
  }
  double sum = 0;
  double step = 1; // Increment step size
  for (double i = 0; i <= q; i += step) {
    Rcpp::Rcout << i << std::endl;
    sum += my_dbinom(i, size, prob) * step; // Multiply by step size
  }
  return sum;
}

// [[Rcpp::export]]
double dbinom_product(  const NumericVector & vXs,  NumericVector & nobs, double prob)
{
  int iMax = vXs.size();   // nobs.size() must be at least iMax also
  double dResult = my_dbinom(vXs[0], nobs[0], prob);
  for (int i = 1; i < iMax; ++i)
  {
    dResult *= my_dbinom(vXs[i], nobs[i] - nobs[i-1], prob);
  }
  return dResult;
}



//Binary Efficacy
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
List GetocBiRcpp(int seed, double nsim, NumericMatrix contrast, NumericVector nobs, double b, double b2, double pow2,
                 NumericVector dprior, double ptrue, double phi, Function fff)
  {

  int n1 =nobs.size();
  NumericVector stopbound1;
  NumericVector stopbound(n1);
  double nmax=max(nobs);
  
  List resultslist(5);

  set_seed(seed);

  for (int i=0; i<n1; i++)
  {
    double n=nobs[i];
    if (i == n1-1){
      stopbound1 = fff(b2,dprior,n,phi,contrast);
    } else {
      stopbound1 = fff(b*pow(n/nmax,pow2),dprior,n,phi,contrast);
    }
    
    // Rcpp::Rcout << "nobs " << n << std::endl;
    // Rcpp::Rcout << "ncut " << stopbound1[0] << std::endl;
    if (stopbound1[0] < 0) { // check if there is -Inf bound
      resultslist[0]= 0; // edrop
      resultslist[1]= 0; // t1err
      resultslist[2]= 0; // power
      resultslist[3]= 1; // pts
      return resultslist;
    }
    
    stopbound[i]=stopbound1[0];
  }
  
  double tolerance = 1e-5;
  for (int i=0; i<n1; i++){
    // Check if it is the last element
    if (i == n1 - 1) {
      // Check if the last element is effectively an integer
      if (std::abs(nobs[i] - std::floor(nobs[i])) > tolerance) {
        nobs[i] = std::ceil(nobs[i]);
      }
    } else {
      nobs[i] = std::ceil(nobs[i]);
    }
  }
  // Rcpp::Rcout << "Boundary " << nobs[n1-1] << std::endl;
  // Rcpp::Rcout << "Boundary " << stopbound[0] << ", " << stopbound[1] << std::endl;
  
  
  List temp(5);
  // Rcpp::Rcout << "exact method " << std::endl;
  // temp = Exacterror(nobs, stopbound, phi, ptrue);
  // resultslist[0]= temp[0];
  // resultslist[1]= temp[2];
  // resultslist[2]= temp[3];
  // resultslist[3]= temp[1];
  temp = exact_error_recursive2_Rcpp(nobs, stopbound, phi, ptrue, n1);
  resultslist[0]= 0; // edrop
  resultslist[1]= temp[0]; // t1err
  resultslist[2]= temp[1]; // power
  resultslist[3]= temp[2]; // pts
  resultslist[4]= temp[3]; // sample size under H1
  // Rcpp::Rcout << "t1e " << static_cast<int>(temp[0]) << std::endl;
  // auto&& item = temp[0];
  // using ItemType = std::remove_reference_t<decltype(item)>;
  // std::cout << "Type of temp[0]: " << typeid(ItemType).name() << std::endl;
  

  return resultslist;
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
NumericMatrix GridSearchBiRcpp(int seed, NumericMatrix contrast, NumericVector nobs, NumericVector dprior,
                               NumericVector b1, NumericVector b2, NumericVector pow, double pn, double pa,
                               double cutstart, double nsim, double err1, Function ffff)
{
  List temp1;
  List temp2;
  int k1=cutstart-1;
  NumericVector ocmatbinary_v(6);
  NumericMatrix ocmatbinary(6,1);
  
  // if nobs.size()<5, the calculation is based on exact computation.
  // Rcpp::Rcout << "exact method " << std::endl;
  for (int j=0; j<b1.size(); j++)
  {
    for (int l=0; l<b2.size(); l++){
      for (int k = 0; k<pow.size(); k++)
      {
        
        k1=k;
        double b_mid1=b1[j];
        double b_mid2=b2[l];
        double pow_mid=pow[k];
        
        temp1=GetocBiRcpp(seed, nsim, contrast, nobs, b_mid1, b_mid2, pow_mid, dprior, pa, pn, ffff);
        double temp1_accept = temp1[1];
        //double temp1_pts = temp1[3]/nobs[nobs.size()-1];
        double temp1_pts = static_cast<double>(temp1[3]) / static_cast<double>(nobs[nobs.size() - 1]);
        
        
        
        ocmatbinary_v[0]=b1[j];
        ocmatbinary_v[1]=pow[k];
        ocmatbinary_v[2]=temp1[1]; // t1err
        ocmatbinary_v[3]=temp1[2]; // power
        ocmatbinary_v[4]=temp1[3]; // pts
        ocmatbinary_v[5]=b2[l];
        //ocmatbinary(_,0)=ocmatbinary_v;
        ocmatbinary = cbind(ocmatbinary, ocmatbinary_v);
        
        
      }
    }
    
  }
  
  
  return(ocmatbinary);
}





