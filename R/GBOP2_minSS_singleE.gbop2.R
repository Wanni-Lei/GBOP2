#' PSOGO: Optimal/Minimax design with single boundary for futility
#'
#' This function implements PSOGO to find an optimal or minimax design with single boundary for futility.
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
#' @param nSwarm nSwarm for pso
#' @param maxIter maxIter for pso
#' @param seed  Random seed for reproducibility
#'
#' @return A list on design parameters and operating characteristics
#' @examples
#' \donttest{
#' # init_cluster(2)
#' #  GBOP2_minSS_singleE(
#' #   design = "optimal",
#' #    unified.u = 1,
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
#' #    seed = 1024,
#' #    nSwarm = 64,
#' #    maxIter = 200
#' #  )
#' # stop_cluster()  # Only if init_cluster() was used
#' #  
#' message("Run GBOP2_minSS_singleE() manually for real optimization.")
#' }
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


GBOP2_minSS_singleE <- function( design = "optimal", ## "optimal" or "minimax", "unified"
                                 unified.u = 1, ## specify when design = "unified", u in [0, 1]
                                 weight = 1, ## weight of sample size under null, w in [0, 1]
                                 nlooks = 2, ## number of interim looks. For 3-stage design, nlooks = 2.
                                 p0 = 0.2, ## response rate in null hypothesis 
                                 p1 = 0.4, ## response rate in alternative hypothesis
                                 err1 = 0.05, ## type I error
                                 minPower = 0.8, ## power
                                 maxPatients = 5, ## maximum number of patients
                                 Nmin_cohort1 = 1, ## minimum number of first cohort
                                 Nmin_increase = 1, ## minimum number of increase in each cohort
                                 pso_method = "default", ## choose from "all", "default", "quantum" or "dexp"
                                 nParallel = 3, ## number of PSO-ensemble, only effective when pso_method = "all"
                                 seed = 456,
                                 nSwarm = 64, ## nSwarm in PSO 
                                 maxIter = 200) {
  
  

  b1n <- p0
  b1a <- p1
 ##########################################################################
  ## estimated total time
  message("\nGBOP2 process has started...\n")
  start_time <- Sys.time()  # Start timing

  one_task <- PSO_design(
    design = design, unified.u = unified.u, nlooks = nlooks, b1n = b1n, b1a = b1a, err1 = err1,
    minPower = minPower, weight = weight, maxPatients = maxPatients,
    Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
    method = "default", seed = seed, nSwarm = nSwarm, maxIter = 1)

  end_time <- Sys.time()  # End timing
  execution_time1T <- as.numeric(end_time - start_time)  # Convert to numeric (seconds)

  # Step 2: Estimate total execution time
  N_PSO <- nParallel * 3  # Total number of PSO_design calls
  nCore_used <- if (!is.null(get_cluster())) length(get_cluster()) else 1L
  total_time <- (N_PSO * execution_time1T * maxIter) / nCore_used
  
  
  # Step 3: Display estimated total execution time
  message("\nEstimated total execution time:", round(total_time, 2), "seconds\n")
  message("Or approximately:", round(total_time / 60, 2), "minutes\n")

####################################################################
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
 
  
 
  #########################################################

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
  




  # Perform parallel computation using foreach and %dopar%
  if (pso_method == "all") {


    res <- foreach(i = 1:nParallel,
                   .packages = c("dplyr", "globpso", "R6", "Rcpp", "RcppArmadillo"),
                   .combine = rbind) %operator%  {


                     # # Load necessary Rcpp source and custom functions
                     # Rcpp::sourceCpp(file = "Calculation_minimizeN_twolambda_update.cpp", cacheDir = "cache")
                     # source('PSO_design.gbop2.R')

                     # Extract the seed for the current iteration
                     current_seed <- seeds_list[i]



                     # Call PSO_design with different methods
                     r1 <- PSO_design(
                       design = design, unified.u = unified.u, nlooks = nlooks, b1n = p0, b1a = p1, err1 = err1,
                       minPower = minPower, weight = weight, maxPatients = maxPatients,
                       Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                       method = "default", seed = current_seed, nSwarm = nSwarm, maxIter = maxIter
                     )

                     r2 <- PSO_design(
                       design = design, unified.u = unified.u, nlooks = nlooks, b1n = p0, b1a = p1, err1 = err1,
                       minPower = minPower, weight = weight, maxPatients = maxPatients,
                       Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                       method = "quantum", seed = current_seed, nSwarm = nSwarm, maxIter = maxIter
                     )

                     r3 <- PSO_design(
                       design = design, unified.u = unified.u, nlooks = nlooks, b1n = p0, b1a = p1, err1 = err1,
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

                     r_ensemble1 <- r_ensemble |>
                       distinct(Utility, .keep_all = TRUE)

                     index <- which(r_ensemble1$Utility == min(r_ensemble1$Utility))
                     results <- r_ensemble1[index, ]



                     return(results)
                   }
    res_final <- res |>
      distinct(Utility, .keep_all = TRUE) |>
      filter(Utility == min(Utility)) |>
      filter(Power == max(Power))

  } else{

    r <- PSO_design(
      design = design,
      unified.u = unified.u,
      nlooks = nlooks,
      b1n = p0,  # Null hypothesis response rate
      b1a = p1,  # Alternative hypothesis response rate
      err1 = err1,  # Type I error rate
      minPower = minPower,  # Power
      weight = weight,  # Weight of sample size under null
      maxPatients = maxPatients,  # Maximum number of patients
      Nmin_cohort1 = Nmin_cohort1,
      Nmin_increase = Nmin_increase,
      method = pso_method,  # PSO method
      seed = seed,  # Set seed to calculate OC
      nSwarm = nSwarm,
      maxIter = maxIter)



    r <- unclass(r)
    res_final <- as.data.frame(r) |>
      distinct(Utility, .keep_all = TRUE) |>
      filter(Utility == min(Utility) ) |>
      filter(Power == max(Power))
  } ## else


  # Update progress bar to 100% when computation finishes
  
  if (!is.null(.GBOP2_env$pb)) {
    setTxtProgressBar(.GBOP2_env$pb, 101)
    close(.GBOP2_env$pb)
  }
  
  
  res_final <- as.list( res_final)
  res_final[[1]] <- "GBOP2_minSS_singleE"
  class(res_final)<-"gbop2"
  
 
  on.exit(stop_cluster(), add = TRUE)
  return(res_final)
}














