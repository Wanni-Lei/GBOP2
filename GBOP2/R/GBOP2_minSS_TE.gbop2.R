#' PSOGO: Optimal/Minimax design with efficacy and toxicity boundaries
#'
#' This function implements PSOGO to find an optimal or minimax design with efficacy and toxicity boundaries.
#'
#' @param design choose from "optimal", "minimax", or "unified"
#' @param unified.u specify when design = "unified", u in zero to one
#' @param nlooks number of interim looks
#' @param skip_efficacy default is NULL, indicate skip efficacy as 1 and not skip as 0 in a vector
#' @param skip_toxicity default is NULL, indicate skip toxicity as 1 and not skip as 0 in a vector
#' @param maxPatients maximum number of patients
#' @param Nmin_cohort1  minimum number of first cohort
#' @param Nmin_increase minimum number of increase in each cohort
#' @param nParallel number of pso ensemble
#' @param p01 H0 for efficacy
#' @param p02 H0 for toxicity
#' @param p03 H0 for Eff and Tox
#' @param p11 H1 for efficacy
#' @param p12 H1 for toxicity
#' @param p13 H1 for Eff and Tox
#' @param err_eff Type I error rate: Efficacious but toxic
#' @param err_tox Type I error rate: Safe but futile
#' @param err_all Type I error rate: Futile and toxic
#' @param power_eff power: Efficacious but toxic
#' @param power_tox power: Safe but futile
#' @param power_all power: Futile and toxic
#' @param pso_method "all" for using three distinct pso, otherwise indicate single pso method
#' @param nSwarm nSwarm in PSO
#' @param seed  Random seed for reproducibility
#' @param maxIter maxIter in PSO
#'
#' @return A list on design parameters and operating characteristics
#' @examples
#' \donttest{
#' # init_cluster(2)
#' #  GBOP2_minSS_TE(
#' #   design = "optimal", 
#' #    unified.u = 1, 
#' #    nlooks = 1, 
#' #    skip_efficacy = NULL, 
#' #    skip_toxicity = NULL, 
#' #    maxPatients = 25, 
#' #    Nmin_cohort1 = 10, 
#' #    Nmin_increase = 5, 
#' #    p01 = 0.3, 
#' #    p02 = 0.4, 
#' #    p03 = 0.2, 
#' #    p11 = 0.6, 
#' #    p12 = 0.2, 
#' #    p13 = 0.15, 
#' #    err_eff = 0.1, 
#' #    err_tox = 0.1, 
#' #    err_all = 0.05, 
#' #    power_eff = 0.8, 
#' #    power_tox = 0.8, 
#' #    power_all = 0.8, 
#' #    pso_method = "default", 
#' #    nParallel = 3, 
#' #    seed = 5321, 
#' #    nSwarm = 64, 
#' #    maxIter = 100
#' #  )
#' # stop_cluster()  # Only if init_cluster() was used
#' #  
#' message("Run GBOP2_minSS_singleE() manually for real optimization.")
#' }
#' 
#' 
#' @details
#' Parallel computing is only used when the user explicitly sets nCore > 1. No more than 2 cores should be used
#' unless the user is aware and permits it. The function defaults to sequential execution. If multiple analyses
#' are planned, consider using `init_cluster(nCore)` and `stop_cluster()` manually to control the backend.
#'  
#' @export
#' @import globpso R6 RcppArmadillo 
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @importFrom foreach %dopar% foreach %do%
#' @importFrom doParallel registerDoParallel
#' @importFrom Rcpp sourceCpp cppFunction
#' @importFrom tidyr pivot_wider
#' @importFrom utils txtProgressBar setTxtProgressBar


GBOP2_minSS_TE <- function(design = "optimal", ## choose from "optimal", "minimax" and "unified"
                            unified.u = 1, ## specify when design = "unified", u in [0, 1]
                            nlooks = 1, ## number of interim looks (R)
                            skip_efficacy = NULL, ## A vector with a length equal to the number of stages. The default is NULL (no skip for efficacy). If skipping is enabled, set the corresponding stage to 1; otherwise, set it to 0. For example, c(1,0) indicates that futility monitoring is skipped in the first stage but applied in the second stage.
                            skip_toxicity = NULL, ## A vector with a length equal to the number of stages. The default is NULL (no skip for toxicity). If skipping is enabled, set the corresponding stage to 1; otherwise, set it to 0. For example, c(1,0) indicates that safety monitoring is skipped in the first stage but applied in the second stage.
                            maxPatients = 26, ## maximum number of patients
                            Nmin_cohort1 = 13, ## minimum number of first cohort
                            Nmin_increase = 13, ## minimum number of increase in each cohort
                            p01 = 0.3, ## efficacy under null
                            p02 = 0.4, ## toxicity under null
                            p03 = 0.2, ## efficacy and toxicity under null, quantifying the correlation between efficacy and toxicity
                            p11 = 0.6, ## efficacy under alternative
                            p12 = 0.2, ## toxicity under alternative
                            p13 = 0.15,## efficacy and toxicity under alternative, quantifying the correlation between efficacy and toxicity
                            err_eff = 0.1, ## Type I error for futile
                            err_tox = 0.1, ## Type I error for toxic
                            err_all = 0.05, ## Type I futile and toxic,
                            power_eff = 0.8, ## power for futile
                            power_tox = 0.8, ## power for toxic
                            power_all = 0.8, ## power for futile and toxic
                            pso_method = "all", ## choose from "all", ""default", "quantum" or "dexp"
                            nParallel = 3, ## number of PSO_ensemble
                            seed = 1324, ## seed for PSO
                            nSwarm = 32, ## nSwarm in PSO 
                            maxIter = 100 ## maxIter in PSO
                          ) {

  
  e1n <- p01  # H0 for Eff
  e2n <- p02  # H0 for Tox
  e3n <- p03 # H0 for Eff and Tox
  e1a <- p11  # Ha for Eff
  e2a <- p12  # Ha for Tox
  e3a <- p13 # Ha for Eff and Tox
  
  
##################################
  ## estimated total time
  message("\nGBOP2 process has started...\n")
  start_time <- Sys.time()  # Start timing
  
  one_task <- PSO_design_TE(design = design, unified.u = unified.u, method = "default", nlooks = nlooks, skip_efficacy = skip_efficacy,
                            skip_toxicity = skip_toxicity, maxPatients = maxPatients, Nmin_cohort1 = Nmin_cohort1,
                            Nmin_increase = Nmin_increase, e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,
                            e3a = e3a, err_eff = err_eff, err_tox = err_tox, err_all = err_all, power_eff = power_eff,
                            power_tox = power_tox, power_all = power_all, seed = seed,nSwarm = nSwarm, maxIter = 1)
  
  end_time <- Sys.time()  # End timing
  execution_time1T <- as.numeric(end_time - start_time)  # Convert to numeric (seconds)
  
  # Step 2: Estimate total execution time
  N_PSO <- nParallel * 3  # Total number of PSO_design calls
  nCore_used <- if (!is.null(get_cluster())) length(get_cluster()) else 1L
  total_time <- (N_PSO * execution_time1T * maxIter) / nCore_used
  
  
  # Step 3: Display estimated total execution time
  message("\nEstimated total execution time:", round(total_time, 2), "seconds\n")
  message("Or approximately:", round(total_time / 60, 2), "minutes\n")
  

  # Fake progress bar function to 99%
  fake_progress_bar <- function(total_time) {
    .GBOP2_env$pb <- txtProgressBar(min = 0, max = 101, style = 3)
    for (i in 1:99) {
      Sys.sleep(total_time / 100)
      setTxtProgressBar(.GBOP2_env$pb, i)
    }
  }
  fake_progress_bar(total_time + 30)
  
  
  #####################################################################
  
  # Set up parallel computing
  # Default to sequential unless cluster was manually started
  if (is.null(get_cluster())) {
    message("Running sequentially (single core). To use parallel computing, manually call init_cluster(nCore) before running this function.")
    foreach::registerDoSEQ()
  }
  

  ################################################
  # Define the seed list

  input <- list("seed" = seed)
  set.seed(input$seed)
  seeds_list <- round(runif(1000) * 1e4)
  
  
  `%operator%` <- if (!is.null(get_cluster())) {
    foreach::`%dopar%`
  } else {
    foreach::`%do%`
  }
  
  

  if (pso_method == "all") {
  # Perform parallel computation using foreach
  res <- foreach(i = 1:nParallel, .packages = c("dplyr", "globpso", "R6", "Rcpp", "RcppArmadillo"),
                 .combine = rbind) %operator% {

                   # source("BOP2_functions_EffTox.R")
                   # source("BOP2_TE_function.R")
                   # source("boundcode.R")
                   # Rcpp::sourceCpp(file = "Calculation2_original.cpp")
                   # source('PSO_design_TE.gbop2.R')
                   current_seed <- seeds_list[i]

                  
                     # Run PSO with three methods
                     r1 <- PSO_design_TE(design = design, unified.u = unified.u, method = "default", nlooks = nlooks, skip_efficacy = skip_efficacy,
                                        skip_toxicity = skip_toxicity, maxPatients = maxPatients, Nmin_cohort1 = Nmin_cohort1,
                                        Nmin_increase = Nmin_increase, e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,
                                        e3a = e3a, err_eff = err_eff, err_tox = err_tox, err_all = err_all, power_eff = power_eff,
                                        power_tox = power_tox, power_all = power_all, seed = current_seed, nSwarm = nSwarm, maxIter = maxIter)

                     r2 <- PSO_design_TE(design = design, unified.u = unified.u, method = "quantum", nlooks = nlooks, skip_efficacy = skip_efficacy,
                                        skip_toxicity = skip_toxicity, maxPatients = maxPatients, Nmin_cohort1 = Nmin_cohort1,
                                        Nmin_increase = Nmin_increase, e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,
                                        e3a = e3a, err_eff = err_eff, err_tox = err_tox, err_all = err_all, power_eff = power_eff,
                                        power_tox = power_tox, power_all = power_all, seed = current_seed, nSwarm = nSwarm, maxIter = maxIter)

                     r3 <- PSO_design_TE(design = design, unified.u = unified.u, method = "dexp", nlooks = nlooks, skip_efficacy = skip_efficacy,
                                        skip_toxicity = skip_toxicity, maxPatients = maxPatients, Nmin_cohort1 = Nmin_cohort1,
                                        Nmin_increase = Nmin_increase, e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,
                                        e3a = e3a, err_eff = err_eff, err_tox = err_tox, err_all = err_all, power_eff = power_eff,
                                        power_tox = power_tox, power_all = power_all, seed = current_seed, nSwarm = nSwarm, maxIter = maxIter)

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
                                   boudary_name, "expected_sample", "typeI_01", "typeI_10"       
                                   ,"typeI_00", "Power", "Utility" )
                     colnames(r1) <- listname
                     colnames(r2) <- listname
                     colnames(r3) <- listname
                     r_ensemble <- rbind(r1, r2,r3)
                     
                     r_ensemble <- r_ensemble |>
                       filter(Utility == min(Utility)) |>
                       filter(Power == max(Power))  
                     

                     
                     results <- r_ensemble
                   
                   return(results)
                 }
  
    res_final <- res |>
    distinct(Utility, .keep_all = TRUE) |>
    filter(Utility == min(Utility)) |>
      filter(Power == max(Power))
                   
  } else {
                     r <- PSO_design_TE(design = design, unified.u = unified.u, method = pso_method, nlooks = nlooks, skip_efficacy = skip_efficacy,
                                       skip_toxicity = skip_toxicity, maxPatients = maxPatients, Nmin_cohort1 = Nmin_cohort1,
                                       Nmin_increase = Nmin_increase, e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,
                                       e3a = e3a, err_eff = err_eff, err_tox = err_tox, err_all = err_all, power_eff = power_eff,
                                       power_tox = power_tox, power_all = power_all, nSwarm = nSwarm, maxIter = maxIter)

    
                     res_final <- r
                   }

                 
  # Update progress bar to 100% when computation finishes
  if (!is.null(.GBOP2_env$pb)) {
    setTxtProgressBar(.GBOP2_env$pb, 101)
    close(.GBOP2_env$pb)
  }
  
  
  
  if (pso_method == "all"){
    # Return the final result as a list
    res_final <- as.list(res_final)
    res_final[[1]] <- "GBOP2_minSS_TE"
  } else{
    res_final[[1]] <- "GBOP2_minSS_TE" 
  }
  
  class(res_final)<-"gbop2"
  
  on.exit(stop_cluster(), add = TRUE)
  
  return(res_final)
}












