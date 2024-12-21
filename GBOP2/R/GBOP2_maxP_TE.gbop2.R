#' PSOGO for maximizing power design with efficacy and toxicity boundaries
#'
#' @param design choose from "optimal" or "minimax"
#' @param pso_method method for single PSO, choose from "default", "quantum" or "dexp"
#' @param nlooks number of interim looks
#' @param skip_efficacy default is NULL, indicate skip efficacy as 1 and not skip as 0 in a vector
#' @param skip_toxicity default is NULL, indicate skip toxicity as 1 and not skip as 0 in a vector
#' @param totalPatients number of total patients
#' @param Nmin_cohort1 maximum number of patients
#' @param Nmin_increase minimum number of first cohort
#' @param e1n H0 for efficacy
#' @param e2n H0 for toxicity
#' @param e3n H0 for Eff and Tox
#' @param e1a H1 for efficacy
#' @param e2a H1 for toxicity
#' @param e3a H1 for Eff and Tox
#' @param err_eff Type I error rate: Efficacious but toxic
#' @param err_tox Type I error rate: Safe but futile
#' @param err_all Type I error rate: Futile and toxic
#' @param power_eff power: Efficacious but toxic
#' @param power_tox power: Safe but futile
#' @param power_all power: Futile and toxic
#' @param nSwarm nSwarm in PSO
#' @param maxIter maxIter in PSO
#' @param nParallel number of PSO ensemble
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

GBOP2_maxP_TE <- function(design = "optimal",
                          nlooks = 4,
                          skip_efficacy = NULL,
                          skip_toxicity = NULL,
                          totalPatients = 50,
                          Nmin_cohort1 = 10,
                          Nmin_increase = 5,
                          nParallel = 3,
                          e1n = 0.15,  # H0 for Eff
                          e2n = 0.16,  # H0 for Tox
                          e3n = 0.024, # H0 for Eff and Tox
                          e1a = 0.4,  # Ha for Eff
                          e2a = 0.08,  # Ha for Tox
                          e3a = 0.032, # Ha for Eff and Tox
                          err_eff = 1,  # Type I error rate: Efficacious but toxic
                          err_tox = 1 , # Type I error rate: Safe but futile
                          err_all = 0.1,  # Type I error rate: Futile and toxic
                          power_eff = 0.8,
                          power_tox = 0.8,
                          power_all = 0.8,
                          seed = 1024,
                          pso_method = "all",
                          nSwarm = 32,
                          maxIter = 100,
                          nCore = 4) {




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
  set.seed(123)
  seeds_list <- round(runif(1000) * 1e4)

  if (pso_method == "all") {
  # Perform parallel computation using foreach
  res <- foreach(i = 1:nParallel, .packages = c("dplyr", "globpso", "R6", "Rcpp", "RcppArmadillo"),
                 .combine = rbind) %dopar% {

                   # source("BOP2_functions_EffTox.R")
                   # source("BOP2_TE_function.R")
                   # source("boundcode.R")
                   # Rcpp::sourceCpp(file = "Calculation2_original.cpp")
                   # source('PSO_power_TE.gbop2.R')
                   current_seed <- seeds_list[i]

                   
                     # Run PSO with three methods
  r1 <- PSO_power_TE(design = design, method = "default",nlooks = nlooks,
                  skip_efficacy = NULL, skip_toxicity = NULL ,totalPatients = totalPatients,
                  Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                  e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,  e3a = e3a, # Ha for Eff and Tox
                  err_eff = err_eff,  # Type I error rate: Efficacious but toxic
                  err_tox = err_tox,  err_all = err_all,  # Type I error rate: Futile and toxic
                  power_eff = power_eff,
                  power_tox = power_tox, power_all = power_all, nSwarm = nSwarm,maxIter = maxIter
     )



                     r2 <-  PSO_power_TE(design = design, method = "quantum",nlooks = nlooks,
                                         skip_efficacy = NULL, skip_toxicity = NULL ,totalPatients = totalPatients,
                                         Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                                         e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,  e3a = e3a, # Ha for Eff and Tox
                                        err_eff = err_eff,  # Type I error rate: Efficacious but toxic
                                         err_tox = err_tox,  err_all = err_all,  # Type I error rate: Futile and toxic
                                         power_eff = power_eff,
                                         power_tox = power_tox, power_all = power_all, nSwarm = nSwarm,maxIter = maxIter
                     )

                     r3 <- PSO_power_TE(design = design, method = "dexp",nlooks = nlooks,
                                        skip_efficacy = NULL, skip_toxicity = skip_toxicity ,totalPatients = totalPatients,
                                        Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                                        e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,  e3a = e3a, # Ha for Eff and Tox
                                        err_eff = err_eff,  # Type I error rate: Efficacious but toxic
                                        err_tox = err_tox,  err_all = err_all,  # Type I error rate: Futile and toxic
                                        power_eff = power_eff,
                                        power_tox = power_tox, power_all = power_all, nSwarm = nSwarm,maxIter = maxIter
                     )

                     # Combine the results and select best
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
                     
                     for(i in 1:(nlooks+1)){
                       boudary_name[i] <- paste0("boundary_effi", i)
                       boudary_name[i+ nlooks+1] <- paste0("boundary_toxi", i)
                     }
                     
                     
                     listname <- c("function","design", "method", "lambdae1",
                                   "lambdae2", "lambdat1", "lambdat2", "gamma" , cohort_name,      
                                   boudary_name,"typeI_01", "typeI_10"       
                                   ,"typeI_00", "power", "Utility" )
                     colnames(r1) <- listname
                     colnames(r2) <- listname
                     colnames(r3) <- listname
                     r_ensemble <- rbind(r1, r2,r3)
                     
                     r_ensemble <- r_ensemble |>
                       filter(Utility == min(Utility))
                     
                     r_ensemble[[1]] <- "GBOP2_minSS_TE"
                     results <- r_ensemble
                     return(results)
                 } 
  
  
  res_final <- res |>
    distinct(Utility, .keep_all = TRUE) |>
    filter(Utility == min(Utility))
  
} else {
                     r <- PSO_power_TE(design = design,
                                       method = pso_method,
                                       nlooks = nlooks,
                                       skip_efficacy = NULL,
                                       skip_toxicity = NULL,
                                       totalPatients = totalPatients,
                                       Nmin_cohort1 = Nmin_cohort1,
                                       Nmin_increase = Nmin_increase,
                                       e1n = e1n,  # H0 for Eff
                                       e2n = e2n,  # H0 for Tox
                                       e3n = e3n, # H0 for Eff and Tox
                                       e1a = e1a,  # Ha for Eff
                                       e2a = e2a,  # Ha for Tox
                                       e3a = e3a, # Ha for Eff and Tox
                                       err_eff = err_eff,  # Type I error rate: Efficacious but toxic
                                       err_tox = err_tox , # Type I error rate: Safe but futile
                                       err_all = err_all,  # Type I error rate: Futile and toxic
                                       power_eff = power_eff,
                                       power_tox = power_tox,
                                       power_all = power_all,
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
  res_final[[1]] <- "GBOP2_maxP_TE"
} else{
  res_final[[1]] <- "PSO_power_TE" 
}

class(res_final)<-"gbop2"
return(res_final)
}



