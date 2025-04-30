#' PSOGO: Power maximizing design with single boundary for futility
#'
#' This function implements PSOGO to find a power maximizing design with single boundary for futility.
#'
#' @param nlooks number of interim looks
#' @param p0 Null hypothesis response rate
#' @param p1 Alternative hypothesis response rate
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
#'
#' @return A list on design parameters and operating characteristics
#' @examples
#' \donttest{
#' # init_cluster(2)
#' #   GBOP2_maxP_singleE(
#' #   nlooks = 1, 
#' #   p0 = 0.2, 
#' #   p1 = 0.4, 
#' #   err1 = 0.05, 
#' #   minPower = 0.8, 
#' #   totalPatients = 26, 
#' #   Nmin_cohort1 = 10, 
#' #   Nmin_increase = 5, 
#' #   pso_method = "default", 
#' #   nParallel = 3, 
#' #   seed = 1024, 
#' #   nSwarm = 64, 
#' #   maxIter = 200
#' # )
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
#'     
#' @export
#' @import globpso R6 RcppArmadillo 
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @importFrom foreach %dopar% foreach registerDoSEQ %do%
#' @importFrom Rcpp sourceCpp cppFunction  
#' @importFrom tidyr pivot_wider
#' @importFrom utils txtProgressBar setTxtProgressBar

GBOP2_maxP_singleE <- function(
    nlooks = 1,
    p0 = 0.2,  # Null hypothesis response rate
    p1 = 0.4,  # Alternative hypothesis response rate
    err1 = 0.05,  # Type I error rate
    minPower = 0.8, ## power
    totalPatients = 5,  ## maximum number of patients
    Nmin_cohort1 = 1,
    Nmin_increase = 1,
    pso_method = "default", ## three different pso or three single pso
    nParallel = 3,
    seed = 1024,
    nSwarm = 64,
    maxIter = 200){ ## how many cores to use
  ## option for which pso to use
  

  b1n <- p0
  b1a <- p1

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
  nCore_used <- if (!is.null(get_cluster())) length(get_cluster()) else 1L
  total_time <- (N_PSO * execution_time1T * maxIter) / nCore_used
  
  
  # Step 3: Display estimated total execution time
  message("\nEstimated total execution time:", round(total_time, 2), "seconds\n")
  message("Or approximately:", round(total_time / 60, 2), "minutes\n")
  
  #################################################################
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
 

  
  #####################################################
  # Define the seed list
  #set.seed(123)
  input <- list("seed" = seed)
  set.seed(input$seed)
  seeds_list <- round(runif(1000) * 1e4)
  
  `%operator%` <- if (!is.null(get_cluster())) {
    foreach::`%dopar%`
  } else {
    foreach::`%do%`
  }
  
  
  
  # Perform parallel computation using foreach
  if (pso_method == "all") {
  res <- foreach(i = 1:nParallel,
                 .packages = c("dplyr", "globpso", "R6", "Rcpp", "RcppArmadillo"),
                 .combine = rbind) %operator%  {

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

                     r_ensemble1 <- r_ensemble |> distinct(Utility, .keep_all = TRUE)

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
  if (!is.null(.GBOP2_env$pb)) {
    setTxtProgressBar(.GBOP2_env$pb, 101)
    close(.GBOP2_env$pb)
  }
  
  # Return the final result as a list
  res_final <- as.list(res_final)
  res_final[[1]] <- "GBOP2_maxP_singleE"
  class(res_final)<-"gbop2"
  
  on.exit(stop_cluster(), add = TRUE)
  
  return(res_final)
}






