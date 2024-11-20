#' PSOGO for optimal/minimax design with single boundary
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

PSOGO_Design <- function(design = "optimal",
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
                  nSwarm = 64,
                  maxIter = 200,
                  nCore = 4) {


  # library(foreach)
  # library(doParallel)
  # library(globpso)
  # library(R6)
  # library(Rcpp)
  # library(RcppArmadillo)
  # library(dplyr)

  # Set up parallel computing
  cl <- makePSOCKcluster(nCore)  # Define cluster with specified number of cores
  registerDoParallel(cl)

  # Define the seed list
  set.seed(123)
  seeds_list <- round(runif(1000) * 1e4)

  # Perform parallel computation using foreach and %dopar%
  res <- foreach(i = 1:nParallel,
                 .packages = c("dplyr", "globpso", "R6", "Rcpp", "RcppArmadillo"),
                 .combine = rbind) %dopar%  {

                   # # Load necessary Rcpp source and custom functions
                   # Rcpp::sourceCpp(file = "Calculation_minimizeN_twolambda_update.cpp", cacheDir = "cache")
                   # source('PSO_design.gbop2.R')

                   # Extract the seed for the current iteration
                   current_seed <- seeds_list[i]

                   if (pso_method == "all") {

                     # Call PSO_design with different methods
                     r1 <- PSO_design(
                       design = design, nlooks = nlooks, b1n = b1n, b1a = b1a, err1 = err1,
                       minPower = minPower, weight = weight, maxPatients = maxPatients,
                       Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                       method = "default", seed = current_seed, nSwarm = nSwarm, maxIter = maxIter
                     )

                     r2 <- PSO_design(
                       design = design, nlooks = nlooks, b1n = b1n, b1a = b1a, err1 = err1,
                       minPower = minPower, weight = weight, maxPatients = maxPatients,
                       Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                       method = "quantum", seed = current_seed, nSwarm = nSwarm, maxIter = maxIter
                     )

                     r3 <- PSO_design(
                       design = design, nlooks = nlooks, b1n = b1n, b1a = b1a, err1 = err1,
                       minPower = minPower, weight = weight, maxPatients = maxPatients,
                       Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                       method = "dexp", seed = current_seed, nSwarm = nSwarm, maxIter = maxIter
                     )

                     # Combine the results into a list and select the best based on Utility
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
                       boudary_name[i] <- paste0("boundary", i)
                     }

                     listname <- c("function", "design", "weight", "method", "cputime", "parameter.lambda1", "parameter.lambda2","parameter.gamma",cohort_name,boudary_name, "Type.I.Error","Power","Expected.Sample.Size", "Utility")

                     colnames(r1) <- listname
                     colnames(r2) <- listname
                     colnames(r3) <- listname
                     r_ensemble <- rbind(r1, r2,r3)

                     r_ensemble1 <- r_ensemble %>%
                       distinct(Utility, .keep_all = TRUE)

                     index <- which(r_ensemble1$Utility == min(r_ensemble1$Utility))
                     results <- r_ensemble1[index, ]

                   } else{
                     r <- PSO_design(
                       design = design,
                       nlooks = nlooks,
                       b1n = b1n,  # Null hypothesis response rate
                       b1a = b1a,  # Alternative hypothesis response rate
                       err1 = err1,  # Type I error rate
                       minPower = minPower,  # Power
                       weight = weight,  # Weight of sample size under null
                       maxPatients = maxPatients,  # Maximum number of patients
                       Nmin_cohort1 = Nmin_cohort1,
                       Nmin_increase = Nmin_increase,
                       method = pso_method,  # PSO method
                       seed = current_seed,  # Set seed to calculate OC
                       nSwarm = nSwarm,
                       maxIter = maxIter)

                     results <- r

                   }



                   return(results)
                 }

  ## res returns N PSO ensemble
  ## res_final chooses the PSO ensemble with minimum utility
  if (pso_method == "all"){
    res_final <- res |>
      distinct(Utility, .keep_all = TRUE) |>
      filter(Utility == min(Utility))
  } else{
    res <- unclass(res)
    res_final <- as.data.frame(res) |>
      distinct(Utility, .keep_all = TRUE) |>
      filter(Utility == min(Utility))
  }

  # Stop the cluster
  stopCluster(cl)
  registerDoSEQ()

  res_final <- as.list( res_final)
  res_final[[1]] <- "PSOGO_Design"
  class(res_final)<-"gbop2"
  return(res_final)
}






