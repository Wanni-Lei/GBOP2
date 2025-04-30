#' PSOGO: Optimal/Minimax design with dual boundaries
#'
#' This function implements PSOGO to find an optimal or minimax design with dual boundaries.
#'
#' @param design choose from "optimal", "minimax", or "unified"
#' @param unified.u specify when design = "unified", u in zero to one
#' @param nlooks number of interim looks
#' @param p0 Null hypothesis response rate
#' @param p1 Alternative hypothesis response rate
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
#'
#' @return A list on design parameters and operating characteristics
#' @examples
#' \donttest{
#' # init_cluster(2)
#' #  GBOP2_minSS_dualE(
#' #    design = "optimal", 
#' #    unified.u = unified.u, 
#' #    nlooks = 1, 
#' #    p0 = 0.2, 
#' #    p1 = 0.4, 
#' #    err1 = 0.05, 
#' #    minPower = 0.8, 
#' #    weight = 1, 
#' #    maxPatients = 25, 
#' #    Nmin_cohort1 = 10, 
#' #    Nmin_increase = 5, 
#' #    pso_method = "default", 
#' #    nParallel = 3, 
#' #    seed = 123, 
#' #    nSwarm = 64, 
#' #    maxIter = 200
#' #  )
#' # stop_cluster()  # Only if init_cluster() was used
#' #  
#' message("Run GBOP2_minSS_dualE() manually for real optimization.")
#' }
#'
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
#' @importFrom tidyr pivot_wider
#' @importFrom Rcpp sourceCpp cppFunction
#' @importFrom utils txtProgressBar setTxtProgressBar

GBOP2_minSS_dualE <- function(
    design = "optimal",
    unified.u = unified.u,
    weight = 1,
    nlooks = 1,
    p0 = 0.2,
    p1 = 0.4,
    err1 = 0.05,
    minPower = 0.8,
    maxPatients = 5,
    Nmin_cohort1 = 1,
    Nmin_increase = 1,
    pso_method = "default",
    nParallel = 3,
    seed = 123,
    nSwarm = 64,
    maxIter = 200){

  
  
  b1n <- p0
  b1a <- p1
  
  ##################################
  ## estimated total time
  message("\nGBOP2 process has started...\n")
  start_time <- Sys.time()  # Start timing
  
  one_task <- PSO_design_dual(
    design = design,
    unified.u = unified.u,
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
    seed = seed,
    nSwarm = nSwarm,
    maxIter = 1
  )
  
  end_time <- Sys.time()  # End timing
  execution_time1T <- as.numeric(end_time - start_time)  # Convert to numeric (seconds)
  
  # Step 2: Estimate total execution time
  N_PSO <- nParallel * 3  # Total number of PSO_design calls
  nCore_used <- if (!is.null(get_cluster())) length(get_cluster()) else 1L
  total_time <- (N_PSO * execution_time1T * maxIter) / nCore_used
  
  
  # Step 3: Display estimated total execution time
  message("\nEstimated total execution time:", round(total_time, 2), "seconds\n")
  message("Or approximately:", round(total_time / 60, 2), "minutes\n")
  
  
  #fake progress bar 
  fake_progress_bar <- function(total_time) {
    .GBOP2_env$pb <- txtProgressBar(min = 0, max = 101, style = 3)
    for (i in 1:99) {
      Sys.sleep(total_time / 100)
      setTxtProgressBar(.GBOP2_env$pb, i)
    }
  }
  fake_progress_bar(total_time + 30)
  
  
  #####################################################################
  
  # Default to sequential unless cluster was manually started
  if (is.null(get_cluster())) {
    message("Running sequentially (single core). To use parallel computing, manually call init_cluster(nCore) before running this function.")
    foreach::registerDoSEQ()
  }
  
  
  ################################################

  # Define the seed list
  #set.seed(seed)
  input <- list("seed" = seed)
  set.seed(input$seed)
  
  seeds_list <- round(runif(1000) * 1e4)

  
  `%operator%` <- if (!is.null(get_cluster())) {
    foreach::`%dopar%`
  } else {
    foreach::`%do%`
  }
  
  

  if (pso_method == "all") {
    # Perform parallel computation using foreach and %dopar%
    res <- foreach(i = 1:nParallel,
                   .packages = c("dplyr", "globpso", "R6", "Rcpp", "RcppArmadillo", "tidyr"),
                   .combine = rbind) %operator%   {
                     
                     # # Load necessary Rcpp and R scripts
                     # source("boundcode_equalrand_jsm.R")
                     # Rcpp::sourceCpp(file = "Calculation_twoboundaries_jsm.cpp", cacheDir = "cache")
                     # source('PSO_design_dual.gbop2.R')
                     
                     # Extract the seed for the current iteration
                     current_seed <- seeds_list[i]
                     
                     # Call PSO_design_dual with different methods
                     r1 <- PSO_design_dual(
                       design = design,
                       unified.u = unified.u,
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
                       unified.u = unified.u,
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
                       unified.u = unified.u,
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
                     
                     
                    
                     
                     r1$Corhort <- cohort_name
                     r2$Corhort <- cohort_name
                     r3$Corhort <- cohort_name
                     
                     
                     r1_wide <- r1 |>
                       pivot_wider(
                         names_from = Corhort,  # Use COHORT to create new column names
                         values_from = c(boundary.1, boundary.2),
                         names_glue = "{.value}{substr(Corhort, 7, 7)}"  # Extract last number for naming
                       )
                     
                     r2_wide <- r2 |>
                       pivot_wider(
                         names_from = Corhort,  # Use COHORT to create new column names
                         values_from = c(boundary.1, boundary.2),
                         names_glue = "{.value}{substr(Corhort, 7, 7)}"  # Extract last number for naming
                       )
                     
                     
                     r3_wide <- r3 |>
                       pivot_wider(
                         names_from = Corhort,  # Use COHORT to create new column names
                         values_from = c(boundary.1, boundary.2),
                         names_glue = "{.value}{substr(Corhort, 7, 7)}"  # Extract last number for naming
                       )
                     
                     
                     B1name <- c()
                     B2name <- c()
                     for(i in 1:(nlooks+1)){
                     
                     B1name[i] <- paste0("boundaryF", i)
                     B2name[i] <- paste0("boundaryE", i)
                     }
                     
                     
                     listname <- c("function", "design", "weight"
                                   ,"method", "parameters.lambda1", "parameters.lambda_grad1",
                                   "parameters.lambda_grad2", "parameters.Gamma_1", "parameters.Gamma_2", "parameters.Gamma_3", "parameters.delta0",  "parameters.delta1",
                                   cohort_name,  "Type.I.Error",
                                   "Power", "Expected.Sample.Size", "Utility",  B1name,  B2name  )

                     colnames(r1_wide) <- listname
                     colnames(r2_wide) <- listname
                     colnames(r3_wide) <- listname
                     r_ensemble <- rbind(r1_wide, r2_wide, r3_wide) 
                     # r_ensemble[r_ensemble == 999] <- NA 
                     
                     r_ensemble1 <- r_ensemble |>
                       filter(Utility == min(Utility)) |>
                       filter(Power == max(Power))
                     
                     # boundary1 <- t(as.vector(r_ensemble1$boundary.1))
                     # colnames(boundary1) <-c("cohort1bd1", "cohort2bd1")
                     # boundary2 <- t(as.vector(r_ensemble1$boundary.2))
                     # colnames(boundary2) <-c("cohort1bd2", "cohort2bd2")
                     # 
                     # r_ensemble2 <- r_ensemble1 |>
                     #   select(-c("boundary.1", "boundary.2")) |>
                     #   distinct()
                     # 
                     # r_ensemble1_final <- cbind(r_ensemble2, boundary1, boundary2)
                     # 
                     # r_ensemble1_final[[1]] <- "GBOP2_maxP_dual"
                     
                     results <- r_ensemble1
                     
                     return(results)
                   }
    
    res_final <- res |>
      distinct(Utility, .keep_all = TRUE) |>
      filter(Utility == min(Utility)) |>
      filter(Power == max(Power))
    
    
  } else { # Single method
                     r <- PSO_design_dual(   design = design,
                                             unified.u = unified.u,
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

                
  # Update progress bar to 100% when computation finishes
  
  if (!is.null(.GBOP2_env$pb)) {
    setTxtProgressBar(.GBOP2_env$pb, 101)
    close(.GBOP2_env$pb)
  }
  
  
  if (pso_method == "all"){
  # Return the final result as a list
  res_final <- as.list(res_final)
  res_final[[1]] <- "GBOP2_minSS_dualE"
  } else{
  res_final[[1]] <- "GBOP2_minSS_dualE" 
  }
  
  class(res_final)<-"gbop2"
  on.exit(stop_cluster(), add = TRUE)
  
  return(res_final)
}











   