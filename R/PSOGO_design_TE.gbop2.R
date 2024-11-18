#' PSOGO for optimal/minimax design with efficacy and toxicity boundaries
#'
#' @param design choose from "optimal" or "minimax"
#' @param method method for single PSO, choose from "default", "quantum" or "dexp"
#' @param nlooks number of interim looks
#' @param skip_efficacy default is NULL, indicate skip efficacy as 1 and not skip as 0 in a vector
#' @param skip_toxicity default is NULL, indicate skip toxicity as 1 and not skip as 0 in a vector
#' @param maxPatients maximum number of patients
#' @param Nmin_cohort1  minimum number of first cohort
#' @param Nmin_increase minimum number of increase in each cohort
#' @param nParallel number of pso ensemble
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
#' @param pso_method "all" for using three distinct pso, otherwise indicate single pso method
#' @param nSwarm nSwarm in PSO
#' @param maxIter maxIter in PSO
#' @param nCore number of core
#'
#' @return A list on design parameters and operating characteristics
#' @export
#' @import globpso R6 Rcpp RcppArmadillo doParallel
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @importFrom foreach %dopar% foreach



PSOGO_design_TE <- function(design = "optimal",
                            method = "quantum",
                            nlooks = 2,
                            skip_efficacy = NULL,
                            skip_toxicity = NULL,
                            maxPatients = 50,
                            Nmin_cohort1 = 10,
                            Nmin_increase = 5,
                            nParallel = 3,
                            e1n = 0.3,
                            e2n = 0.4,
                            e3n = 0.2,
                            e1a = 0.6,
                            e2a = 0.2,
                            e3a = 0.15,
                            err_eff = 0.1,
                            err_tox = 0.1,
                            err_all = 0.05,
                            power_eff = 0.8,
                            power_tox = 0.8,
                            power_all = 0.8,
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
  cl <- makePSOCKcluster(nCore)
  registerDoParallel(cl)

  # Define the seed list
  set.seed(123)
  seeds_list <- round(runif(1000) * 1e4)

  # Perform parallel computation using foreach
  res <- foreach(i = 1:nParallel, .packages = c("dplyr", "globpso", "R6", "Rcpp", "RcppArmadillo"),
                 .combine = rbind) %dopar% {

                   # source("BOP2_functions_EffTox.R")
                   # source("BOP2_TE_function.R")
                   # source("boundcode.R")
                   # Rcpp::sourceCpp(file = "Calculation2_original.cpp")
                   # source('PSODesign_TE.gbop2.R')
                   current_seed <- seeds_list[i]

                   if (pso_method == "all") {
                     # Run PSO with three methods
                     r1 <- PSODesign_TE(design = design, method = "default", nlooks = nlooks, skip_efficacy = skip_efficacy,
                                        skip_toxicity = skip_toxicity, maxPatients = maxPatients, Nmin_cohort1 = Nmin_cohort1,
                                        Nmin_increase = Nmin_increase, e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,
                                        e3a = e3a, err_eff = err_eff, err_tox = err_tox, err_all = err_all, power_eff = power_eff,
                                        power_tox = power_tox, power_all = power_all, nSwarm = nSwarm, maxIter = maxIter)

                     r2 <- PSODesign_TE(design = design, method = "quantum", nlooks = nlooks, skip_efficacy = skip_efficacy,
                                        skip_toxicity = skip_toxicity, maxPatients = maxPatients, Nmin_cohort1 = Nmin_cohort1,
                                        Nmin_increase = Nmin_increase, e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,
                                        e3a = e3a, err_eff = err_eff, err_tox = err_tox, err_all = err_all, power_eff = power_eff,
                                        power_tox = power_tox, power_all = power_all, nSwarm = nSwarm, maxIter = maxIter)

                     r3 <- PSODesign_TE(design = design, method = "dexp", nlooks = nlooks, skip_efficacy = skip_efficacy,
                                        skip_toxicity = skip_toxicity, maxPatients = maxPatients, Nmin_cohort1 = Nmin_cohort1,
                                        Nmin_increase = Nmin_increase, e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,
                                        e3a = e3a, err_eff = err_eff, err_tox = err_tox, err_all = err_all, power_eff = power_eff,
                                        power_tox = power_tox, power_all = power_all, nSwarm = nSwarm, maxIter = maxIter)

                     # Combine the results and select best
                     r_ensemble <- rbind(r1, r2, r3) %>%
                       as.data.frame() %>%
                       distinct(Utility, .keep_all = TRUE) %>%
                       filter(Utility == min(Utility))

                     r_ensemble[[1]] <- "PSOGO_design_TE"
                     results <- r_ensemble
                   } else {
                     r <- PSODesign_TE(design = design, method = pso_method, nlooks = nlooks, skip_efficacy = skip_efficacy,
                                       skip_toxicity = skip_toxicity, maxPatients = maxPatients, Nmin_cohort1 = Nmin_cohort1,
                                       Nmin_increase = Nmin_increase, e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,
                                       e3a = e3a, err_eff = err_eff, err_tox = err_tox, err_all = err_all, power_eff = power_eff,
                                       power_tox = power_tox, power_all = power_all, nSwarm = nSwarm, maxIter = maxIter)

                     r[[1]] <- "PSOGO_design_TE"
                     results <- r
                   }

                   return(results)
                 }

  # Post-process results to choose best PSO ensemble based on utility
  if (pso_method == "all") {
    res_final <- res %>%
      distinct(Utility, .keep_all = TRUE) %>%
      filter(Utility == min(Utility))
  } else {
    res_final <- as.data.frame(res) %>%
      distinct(Utility, .keep_all = TRUE) %>%
      filter(Utility == min(Utility))
  }

  # Stop the cluster
  stopCluster(cl)
  registerDoSEQ()

  res_final <- as.list(res_final)
  class(res_final)<-"gbop2"
  return(res_final)
}












