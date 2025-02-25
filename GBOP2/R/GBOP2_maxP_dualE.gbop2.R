#' PSOGO: Power maximizing design with dual boundaries
#'
#' This function implements PSOGO to find a power maximizing design with dual boundaries.
#'
#' @name GBOP2_maxP_dualE
#' @param design fixed as "optimal", which can not be modified by user
#' @param nlooks number of interim looks
#' @param b1n Null hypothesis response rate
#' @param b1a Alternative hypothesis response rate
#' @param err1 Type I error rate
#' @param nParallel number of pso ensemble
#' @param minPower power
#' @param totalPatients total patients
#' @param Nmin_cohort1 minimum number of first cohort
#' @param Nmin_increase minimum number of increase in each cohort
#' @param pso_method "all" for using three distinct pso, otherwise indicate single pso method
#' @param seed seed for pso
#' @param nSwarm nSwarm for pso
#' @param maxIter maxIter for pso
#' @param nCore number of core
#'
#' @return A list on design parameters and operating characteristics
#' @examples 
#' GBOP2_maxP_dualE(
#'   nlooks = 1, 
#'   b1n = 0.2, 
#'   b1a = 0.4, 
#'   err1 = 0.05, 
#'   nParallel = 3, 
#'   minPower = 0.8, 
#'   totalPatients = 26, 
#'   Nmin_cohort1 = 10, 
#'   Nmin_increase = 5, 
#'   pso_method = "all", 
#'   seed = 1024, 
#'   nSwarm = 64, 
#'   maxIter = 200, 
#'   nCore = 8
#' )
#'
#' 
#' @export
#' @import globpso R6 Rcpp RcppArmadillo doParallel tidyr
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @importFrom foreach %dopar% foreach
#' @importFrom parallel makePSOCKcluster 
#' @importFrom doParallel registerDoParallel
#' @importFrom tidyr pivot_wider

utils::globalVariables(c("Power", "Corhort", "boundary.1", "boundary.2"))

GBOP2_maxP_dualE <- function(design = "optimal",
                           nlooks = 1,
                           b1n = 0.2,
                           b1a = 0.4,
                           err1 = 0.05,
                           nParallel = 3,
                           minPower = 0.8,
                           totalPatients = 50,
                           Nmin_cohort1 = 10,
                           Nmin_increase = 5,
                           pso_method = "all",
                           seed = 1024,
                           nSwarm = 64,
                           maxIter = 200,
                           nCore = 4) {




  ##########################################################################
  ## estimated total time
  cat("\nGBOP2 process has started...\n")
  start_time <- Sys.time()  # Start timing
  
  one_task <- PSO_power_dual(method = "default", totalPatients = totalPatients, nlooks = nlooks,
                             Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                             b1n = b1n, b1a = b1a, err1 = err1, minPower = minPower,
                             seed = seed, nSwarm = nSwarm, maxIter = 1)
  
  
  end_time <- Sys.time()  # End timing
  execution_time1T <- as.numeric(end_time - start_time)  # Convert to numeric (seconds)
  
  # Step 2: Estimate total execution time
  N_PSO <- nParallel * 3  # Total number of PSO_design calls
  total_time <- (N_PSO * execution_time1T * maxIter) / nCore  # Estimate based on parallel cores
  
  # Step 3: Display estimated total execution time
  cat("\nEstimated total execution time:", round(total_time, 2), "seconds\n")
  cat("Or approximately:", round(total_time / 60, 2), "minutes\n")
  
  #fake progress bar 
  # Fake progress bar function to 99%
  fake_progress_bar <- function(total_time) {
    pb <- txtProgressBar(min = 0, max = 101, style = 3)
    
    # Progress from 1% to 99%
    for (i in 1:99) {
      Sys.sleep(total_time / 100)
      setTxtProgressBar(pb, i)
    }
    
    # Do not wait; let the main function run
    #assign("progress_bar", pb, envir = .GlobalEnv)  # Store in global env to update later
    GBOP2_env <- new.env(parent = emptyenv())  # Create package-specific environment
    assign("progress_bar", pb, envir = GBOP2_env)
    
  }
  
  # Start the progress bar in the background
  total_time <- total_time +30  # Adjust this based on execution time
  fake_progress_bar(total_time)
  
  
  
  ##############################################################

  # Set up parallel computing
  cl <- makePSOCKcluster(nCore)
  registerDoParallel(cl)

  # Define the seed list
  #set.seed(123)
  input <- list("seed" = seed)
  set.seed(input$seed)
  seeds_list <- round(runif(1000) * 1e4)
 
 
if (pso_method == "all") {  # Perform parallel computation using foreach and %dopar%
  res <- foreach(i = 1:nParallel,
                 .packages = c("dplyr", "globpso", "R6", "Rcpp", "RcppArmadillo", "tidyr"),
                 .combine = rbind) %dopar%  {

                   # Load Rcpp source and custom functions for each worker
                   # source("boundcode_equalrand_jsm.R")
                   # Rcpp::sourceCpp(file = "Calculation_twoboundaries_jsm.cpp", cacheDir = "cache")
                   numOfSimForTiralSetting <- 10000   # Number of simulations
                   # source('PSO_power_dual.gbop2.R')

                   # Extract the seed for the current iteration
                   current_seed <- seeds_list[i]

                   
                     # Call PSO_power_dual with different methods
                     r1 <- PSO_power_dual(method = "default", totalPatients = totalPatients, nlooks = nlooks,
                                         Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                                         b1n = b1n, b1a = b1a, err1 = err1, minPower = minPower,
                                         seed = current_seed, nSwarm = nSwarm, maxIter = maxIter)

                     r2 <- PSO_power_dual(method = "quantum", totalPatients = totalPatients, nlooks = nlooks,
                                         Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                                         b1n = b1n, b1a = b1a, err1 = err1, minPower = minPower,
                                         seed = current_seed, nSwarm = nSwarm, maxIter = maxIter)

                     r3 <- PSO_power_dual(method = "dexp", totalPatients = totalPatients, nlooks = nlooks,
                                         Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                                         b1n = b1n, b1a = b1a, err1 = err1, minPower = minPower,
                                         seed = current_seed, nSwarm = nSwarm, maxIter = maxIter)

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
                     }


                     r1$Corhort <- cohort_name
                     r2$Corhort <- cohort_name
                     r3$Corhort <- cohort_name
                     
                     
                     r1_wide <- r1 %>%
                       pivot_wider(
                         names_from = Corhort,  # Use COHORT to create new column names
                         values_from = c(boundary.1, boundary.2),
                         names_glue = "{.value}{substr(Corhort, 7, 7)}"  # Extract last number for naming
                       )
                     
                     r2_wide <- r2 %>%
                       pivot_wider(
                         names_from = Corhort,  # Use COHORT to create new column names
                         values_from = c(boundary.1, boundary.2),
                         names_glue = "{.value}{substr(Corhort, 7, 7)}"  # Extract last number for naming
                       )
                     
                     
                     r3_wide <- r3 %>%
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
                     
                     
                     listname <- c("function", "design", "method", "parameters.lambda1",    
                                    "parameters.lambda_grad1", "parameters.lambda_grad2", "parameters.Gamma_1", "parameters.Gamma_2", "parameters.Gamma_3", "parameters.delta0",  "parameters.delta1",
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
                  
  } else{
    r <- PSO_power_dual(method = pso_method, totalPatients = totalPatients, nlooks = nlooks,
                        Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                        b1n = b1n, b1a = b1a, err1 = err1, minPower = minPower,
                        seed = seed, nSwarm = nSwarm, maxIter = maxIter)
    r <- unclass(r)
    res_final <- r
    
  }
  
  # Update progress bar to 100% when computation finishes
  if (exists("progress_bar", envir = .GlobalEnv)) {
    setTxtProgressBar(get("progress_bar", envir = .GlobalEnv), 101)
    close(get("progress_bar", envir = .GlobalEnv))
  }

  # Stop the cluster
  stopCluster(cl)
  registerDoSEQ()
  
  
  if (pso_method == "all"){
    # Return the final result as a list
    res_final <- as.list(res_final)
    res_final[[1]] <- "GBOP2_maxP_dualE"
  } else{
    res_final[[1]] <- "PSO_power_dual" 
  }
  
  class(res_final)<-"gbop2"
  return(res_final)
}




