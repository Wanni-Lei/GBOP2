#' PSOGO for maximizing power design with single boundary
#'
#' @param nlooks number of interim looks
#' @param b1n Null hypothesis response rate
#' @param b1a Alternative hypothesis response rate
#' @param err1 Type I error rate
#' @param nParallel number of pso ensemble
#' @param minPower power
#' @param weight weight of sample size under null
#' @param totalPatients total number of patients
#' @param Nmin_cohort1 minimum number of first cohort
#' @param Nmin_increase minimum number of increase in each cohort
#' @param pso_method "all" for using three distinct pso, otherwise indicate single pso method
#' @param nSwarm nSwarm for pso
#' @param maxIter maxIter for pso
#' @param nCore number of core
#'
#' @return a list
#' @export
#' @import globpso R6 Rcpp RcppArmadillo doParallel magrittr
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @importFrom foreach %dopar% foreach registerDoSEQ
#' @importFrom parallel makePSOCKcluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom magrittr %>%

GBOP2_maxP_single <- function(
    nlooks = 1,
    b1n = 0.2,  # Null hypothesis response rate
    b1a = 0.4,  # Alternative hypothesis response rate
    err1 = 0.05,  # Type I error rate
    nParallel = 3,
    minPower = 0.8, ## power
    weight = 1, ## weight of sample size under null
    totalPatients = 50,  ## maximum number of patients
    Nmin_cohort1 = 10,
    Nmin_increase = 5,
    pso_method = "all", ## three different pso or three single pso
    seed = 1024,
    nSwarm = 64,
    maxIter = 200,
    nCore = 4){ ## how many cores to use
  ## option for which pso to use



# library(foreach)
# library(doParallel)
# library(globpso)
# library(R6)
# library(Rcpp)
# library(RcppArmadillo)
# library(dplyr)

  # Set up parallel computing
  cl <- makePSOCKcluster(nCore)  # Define cluster with 4 cores
  registerDoParallel(cl)

  # Define the seed list
  set.seed(123)
  seeds_list <- round(runif(1000) * 1e4)

  # Perform parallel computation using foreach
  if (pso_method == "all") {
  res <- foreach(i = 1:nParallel,
                 .packages = c("dplyr", "globpso", "R6", "Rcpp", "RcppArmadillo"),
                 .combine = rbind) %dopar%  {

                   # Load necessary libraries for each worker
                   # library(globpso)
                   # library(R6)
                   # library(Rcpp)
                   # library(RcppArmadillo)
                   # library(dplyr)

                   # Rcpp::sourceCpp(file = "Calculation_minimizeN_twolambda_update.cpp", cacheDir = "cache")
                   # source('PSO_power.gbop2.R')

                   # Extract the seed for the current iteration
                   current_seed <- seeds_list[i]

                     # Call PSO_power with different methods
                     r1 <- PSO_power(
                       nlooks = nlooks,
                       totalPatients = totalPatients,
                       Nmin_cohort1 = Nmin_cohort1,
                       Nmin_increase = Nmin_increase,
                       method = "default",
                       b1n = b1n,
                       b1a = b1a,
                       err1 = err1,
                       minPower = minPower,
                       seed = current_seed,
                       nSwarm = nSwarm,
                       maxIter = maxIter
                     )

                     r2 <- PSO_power(
                       nlooks = nlooks,
                       totalPatients = totalPatients,
                       Nmin_cohort1 = Nmin_cohort1,
                       Nmin_increase = Nmin_increase,
                       method = "quantum",
                       b1n = b1n,
                       b1a = b1a,
                       err1 = err1,
                       minPower = minPower,
                       seed = current_seed,
                       nSwarm = nSwarm,
                       maxIter = maxIter
                     )

                     r3 <- PSO_power(
                       nlooks = nlooks,
                       totalPatients = totalPatients,
                       Nmin_cohort1 = Nmin_cohort1,
                       Nmin_increase = Nmin_increase,
                       method = "dexp",
                       b1n = b1n,
                       b1a = b1a,
                       err1 = err1,
                       minPower = minPower,
                       seed = current_seed,
                       nSwarm = nSwarm,
                       maxIter = maxIter
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

                     listname <- c("function", "design", "method", "cputime", "lambda1", "lambda2",
                                   "gamma", cohort_name, boudary_name, "TypeI", "Power", "EN.P0",      "Utility" )

                     colnames(r1) <- listname
                     colnames(r2) <- listname
                     colnames(r3) <- listname

                     r_ensemble <- rbind(r1, r2,r3)

                     r_ensemble1 <- r_ensemble %>% distinct(Utility, .keep_all = TRUE)

                     # Filter the rows with maximum absolute Utility
                     index <- which(abs(r_ensemble1$Utility) == max(abs(r_ensemble1$Utility)))
                     results <- r_ensemble1[index, ]
                   
                     return(results)
                    }  
    res_final <- res |>
    distinct(Utility, .keep_all = TRUE) |>
    filter(Utility == min(Utility))
  
      } else {
                     # Single PSO method
                     r <- PSO_power(
                       nlooks = nlooks,
                       totalPatients = totalPatients,
                       Nmin_cohort1 = Nmin_cohort1,
                       Nmin_increase = Nmin_increase,
                       method = pso_method,
                       b1n = b1n,
                       b1a = b1a,
                       err1 = err1,
                       minPower = minPower,
                       seed = seed,
                       nSwarm = nSwarm,
                       maxIter = maxIter
                     )
                     
                   

          r <- unclass(r)
          res_final <- as.data.frame(r) |>
          distinct(Utility, .keep_all = TRUE) |>
          filter(Utility == min(Utility))
      }

  # Stop the cluster
  stopCluster(cl)
  registerDoSEQ()

  # Return the final result as a list
  res_final <- as.list(res_final)
  res_final[[1]] <- "GBOP2_maxP_single"
  class(res_final)<-"gbop2"
  return(res_final)
}






