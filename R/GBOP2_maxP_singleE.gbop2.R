#' PSOGO: Power maximizing design with single boundary
#'
#' This function implements PSOGO to find a power maximizing design with single boundary.
#'
#' @param nlooks number of interim looks
#' @param b1n Null hypothesis response rate
#' @param b1a Alternative hypothesis response rate
#' @param err1 Type I error rate
#' @param nParallel number of pso ensemble
#' @param minPower power
#' @param totalPatients total number of patients
#' @param Nmin_cohort1 minimum number of first cohort
#' @param Nmin_increase minimum number of increase in each cohort
#' @param pso_method "all" for using three distinct pso, otherwise indicate single pso method
#' @param nSwarm nSwarm for pso
#' @param maxIter maxIter for pso
#' @param seed  Random seed for reproducibility
#' @param nCore number of core
#'
#' @return A list on design parameters and operating characteristics
#' @examples 
#' GBOP2_maxP_singleE(
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
#' @import globpso R6 Rcpp RcppArmadillo doParallel magrittr
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @importFrom foreach %dopar% foreach registerDoSEQ
#' @importFrom parallel makePSOCKcluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom magrittr %>%

GBOP2_maxP_singleE <- function(
    nlooks = 1,
    b1n = 0.2,  # Null hypothesis response rate
    b1a = 0.4,  # Alternative hypothesis response rate
    err1 = 0.05,  # Type I error rate
    nParallel = 3,
    minPower = 0.8, ## power
    totalPatients = 50,  ## maximum number of patients
    Nmin_cohort1 = 10,
    Nmin_increase = 5,
    pso_method = "all", ## three different pso or three single pso
    seed = 1024,
    nSwarm = 64,
    maxIter = 200,
    nCore = 4){ ## how many cores to use
  ## option for which pso to use



  ##########################################################################
  ## estimated total time
  cat("\nGBOP2 process has started...\n")
  start_time <- Sys.time()  # Start timing
  
  one_task <- PSO_power(
    nlooks = nlooks,
    totalPatients = totalPatients,
    Nmin_cohort1 = Nmin_cohort1,
    Nmin_increase = Nmin_increase,
    method = "default",
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
  total_time <- (N_PSO * execution_time1T * maxIter) / nCore  # Estimate based on parallel cores
  
  # Step 3: Display estimated total execution time
  cat("\nEstimated total execution time:", round(total_time, 2), "seconds\n")
  cat("Or approximately:", round(total_time / 60, 2), "minutes\n")
  
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
  
  ####################################################################

  # Set up parallel computing
  cl <- makePSOCKcluster(nCore)  # Define cluster with 4 cores
  registerDoParallel(cl)

  # Define the seed list
  #set.seed(123)
  input <- list("seed" = seed)
  set.seed(input$seed)
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
    filter(Utility == min(Utility)) |>
      filter(Power == max(Power))
  
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
          filter(Utility == min(Utility)) |>
            filter(Power == max(Power))
      }

  
  # Update progress bar to 100% when computation finishes
  if (exists("progress_bar", envir = .GlobalEnv)) {
    setTxtProgressBar(get("progress_bar", envir = .GlobalEnv), 101)
    close(get("progress_bar", envir = .GlobalEnv))
  }
  
  # Stop the cluster
  stopCluster(cl)
  registerDoSEQ()

  # Return the final result as a list
  res_final <- as.list(res_final)
  res_final[[1]] <- "GBOP2_maxP_singleE"
  class(res_final)<-"gbop2"
  return(res_final)
}






