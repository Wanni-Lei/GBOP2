#' PSOGO for optimal/minimax design with dual bounderies
#'
#' @param design choose from "optimal" or "minimax"
#' @param nlooks number of interim looks
#' @param b1n Null hypothesis response rate
#' @param b1a Alternative hypothesis response rate
#' @param err1 Type I error rate
#' @param nParallel number of pso ensemble
#' @param minPower power
#' @param weight weight of sample size under null
#' @param maxPatients maximum number of patients
#' @param Nmin_cohort1 minimum number of first cohort
#' @param Nmin_increase minimum number of increase in each cohort
#' @param pso_method "all" for using three distinct pso, otherwise indicate single pso method
#' @param seed seed for pso
#' @param nSwarm nSwarm for pso
#' @param maxIter maxIter for pso
#' @param nCore number of core
#'
#' @return A list on design parameters and operating characteristics
#' @export
#' @import globpso R6 Rcpp RcppArmadillo doParallel
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @importFrom foreach %dopar% foreach
#' @importFrom parallel makePSOCKcluster 
#' @importFrom doParallel registerDoParallel

GBOP2_minSS_dual <- function(
    design = "optimal",
    nlooks = 1,
    b1n = 0.2,
    b1a = 0.4,
    err1 = 0.05,
    nParallel = 3,
    minPower = 0.8,
    weight = 1,
    maxPatients = 50,
    Nmin_cohort1 = 10,
    Nmin_increase = 5,
    pso_method = "all",
    seed = 123,
    nSwarm = 64,
    maxIter = 200,
    nCore = 4){

  # Load necessary libraries
  # library(foreach)
  # library(doParallel)
  # library(globpso)
  # library(R6)
  # library(Rcpp)
  # library(RcppArmadillo)
  # library(dplyr)

  # Set up parallel computing
  cl <- makePSOCKcluster(nCore)
  registerDoParallel(cl)

  # Define the seed list
  set.seed(seed)
  seeds_list <- round(runif(1000) * 1e4)

  if (pso_method == "all") {
    # Perform parallel computation using foreach and %dopar%
    res <- foreach(i = 1:nParallel,
                   .packages = c("dplyr", "globpso", "R6", "Rcpp", "RcppArmadillo"),
                   .combine = rbind) %dopar%  {
                     
                     # # Load necessary Rcpp and R scripts
                     # source("boundcode_equalrand_jsm.R")
                     # Rcpp::sourceCpp(file = "Calculation_twoboundaries_jsm.cpp", cacheDir = "cache")
                     # source('PSO_design_dual.gbop2.R')
                     
                     # Extract the seed for the current iteration
                     current_seed <- seeds_list[i]
                     
                     # Call PSO_design_dual with different methods
                     r1 <- PSO_design_dual(
                       design = design,
                       method = "default",
                       maxPatients = maxPatients,
                       nlooks = nlooks,
                       Nmin_cohort1 = Nmin_cohort1,
                       Nmin_increase = Nmin_increase,
                       weight = weight,
                       b1n = b1n,
                       b1a = b1a,
                       err1 = err1,
                       minPower = minPower,
                       seed = current_seed,
                       nSwarm = nSwarm,
                       maxIter = maxIter
                     )
                     
                     r2 <- PSO_design_dual(
                       design = design,
                       method = "quantum",
                       maxPatients = maxPatients,
                       nlooks = nlooks,
                       Nmin_cohort1 = Nmin_cohort1,
                       Nmin_increase = Nmin_increase,
                       weight = weight,
                       b1n = b1n,
                       b1a = b1a,
                       err1 = err1,
                       minPower = minPower,
                       seed = current_seed,
                       nSwarm = nSwarm,
                       maxIter = maxIter
                     )
                     
                     r3 <- PSO_design_dual(
                       design = design,
                       method = "dexp",
                       maxPatients = maxPatients,
                       nlooks = nlooks,
                       Nmin_cohort1 = Nmin_cohort1,
                       Nmin_increase = Nmin_increase,
                       weight = weight,
                       b1n = b1n,
                       b1a = b1a,
                       err1 = err1,
                       minPower = minPower,
                       seed = current_seed,
                       nSwarm = nSwarm,
                       maxIter = maxIter
                     )
                     r1 <- unclass(r1)
                     r1 <- as.data.frame(r1)
                     r2 <- unclass(r2)
                     r2 <- as.data.frame(r2)
                     r3 <- unclass(r3)
                     r3 <- as.data.frame(r3)
                     
                     cohort_name <- c()
                     boudary_name <- c()
                     for(i in 1:(nlooks+1)){
                       cohort_name[i] <- paste0("cohort", i)
                     }
                     
                     
                     listname <- c("function", "design", "weight"                 
                                   ,"method", "parameters.lambda1", "parameters.lambda_grad1", 
                                   "parameters.lambda_grad2", "parameters.Gamma_1", "parameters.Gamma_2", "parameters.Gamma_3", "parameters.delta0",       "parameters.delta1", cohort_name,   "boundary.1",  "boundary.2", "Type.I.Error", 
                                   "Power", "Expected.Sample.Size", "Utility"  )
                     
                     colnames(r1) <- listname
                     colnames(r2) <- listname
                     colnames(r3) <- listname
                     r_ensemble <- rbind(r1, r2,r3)
                     
                     r_ensemble1 <- r_ensemble |>
                       filter(Utility == min(Utility))
                     
                     boundary1 <- t(as.vector(r_ensemble1$boundary.1))
                     colnames(boundary1) <-c("cohort1bd1", "cohort2bd1")
                     boundary2 <- t(as.vector(r_ensemble1$boundary.2))
                     colnames(boundary2) <-c("cohort1bd2", "cohort2bd2")
                     
                     r_ensemble2 <- r_ensemble1 |>
                       select(-c("boundary.1", "boundary.2")) |>
                       distinct()
                     
                     r_ensemble1_final <- cbind(r_ensemble2, boundary1, boundary2)
                     
                     r_ensemble1_final[[1]] <- "GBOP2_maxP_dual"
                     
                     results <- r_ensemble1_final
                     
                     return(results)
                   }
    
    res_final <- res |>
      distinct(Utility, .keep_all = TRUE) |>
      filter(Utility == min(Utility))
    
  } else { # Single method
                     r <- PSO_design_dual(   design = design,
                                              nlooks = nlooks,
                                              b1n = b1n,
                                              b1a = b1a,
                                              err1 = err1,
                                              minPower = minPower,
                                              weight = weight,
                                              maxPatients = maxPatients,
                                              Nmin_cohort1 = Nmin_cohort1,
                                              Nmin_increase = Nmin_increase,
                                              method = pso_method,
                                              seed = seed,
                                              nSwarm = nSwarm,
                                              maxIter = maxIter
                     )
                     res_final <- r
                   }

                

  # Stop the cluster
  stopCluster(cl)
  registerDoSEQ()

  
  if (pso_method == "all"){
  # Return the final result as a list
  res_final <- as.list(res_final)
  res_final[[1]] <- "GBOP2_minSS_dual"
  } else{
  res_final[[1]] <- "PSO_power_dual" 
  }
  
  class(res_final)<-"gbop2"
  return(res_final)
}











   