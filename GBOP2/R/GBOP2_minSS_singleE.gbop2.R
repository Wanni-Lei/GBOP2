#' PSOGO: Optimal/Minimax design with single boundary
#'
#' This function implements PSOGO to find an optimal or minimax design with single boundary.
#'
#' @param design choose from "optimal", "minimax", or "unified"
#' @param unified.u specify when design = "unified", u in zero to one
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
#' @param seed  Random seed for reproducibility
#' @param nCore number of core
#'
#' @return A list on design parameters and operating characteristics
#' @examples 
#' GBOP2_minSS_singleE(
#'   design = "optimal", 
#'   unified.u = 1, 
#'   nlooks = 1, 
#'   b1n = 0.2, 
#'   b1a = 0.4, 
#'   err1 = 0.05, 
#'   nParallel = 3, 
#'   minPower = 0.8, 
#'   weight = 1, 
#'   maxPatients = 25, 
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
#' @import globpso R6 Rcpp RcppArmadillo doParallel doSNOW utils iterators snow
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @importFrom foreach %dopar% foreach registerDoSEQ
#' @importFrom parallel makePSOCKcluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW

GBOP2_minSS_singleE <- function(design = "optimal",
                                unified.u = 1,
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
                                seed = 1024,
                                nSwarm = 64,
                                maxIter = 200,
                                nCore = 4) {


 ##########################################################################
  ## estimated total time
  cat("\nGBOP2 process has started...\n")
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
  total_time <- (N_PSO * execution_time1T * maxIter) / nCore  # Estimate based on parallel cores

  # Step 3: Display estimated total execution time
  cat("\nEstimated total execution time:", round(total_time, 2), "seconds\n")
  cat("Or approximately:", round(total_time / 60, 2), "minutes\n")

####################################################################
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

  
  
#####################################################################
  # Set up parallel computing
  #cl <- parallel::makePSOCKcluster(nCore)  # Define cluster with specified number of cores
  cl <- makeCluster(nCore, type='SOCK')
  doSNOW::registerDoSNOW(cl)

  # Define the seed list
  #set.seed(123)
  input <- list("seed" = seed)
  set.seed(input$seed)

  seeds_list <- round(runif(1000) * 1e4)





  # Perform parallel computation using foreach and %dopar%
  if (pso_method == "all") {


    res <- foreach(i = 1:nParallel,
                   .packages = c("dplyr", "globpso", "R6", "Rcpp", "RcppArmadillo"),
                   ## progress bar
                   #.options.snow = opts,
                   .combine = rbind) %dopar%  {


                     # Notify only the first iteration
                     if (i == 1) {
                       message(paste0("Parallel computation is running: Task ", i, " is starting..."))
                     }




                     # # Load necessary Rcpp source and custom functions
                     # Rcpp::sourceCpp(file = "Calculation_minimizeN_twolambda_update.cpp", cacheDir = "cache")
                     # source('PSO_design.gbop2.R')

                     # Extract the seed for the current iteration
                     current_seed <- seeds_list[i]



                     # Call PSO_design with different methods
                     r1 <- PSO_design(
                       design = design, unified.u = unified.u, nlooks = nlooks, b1n = b1n, b1a = b1a, err1 = err1,
                       minPower = minPower, weight = weight, maxPatients = maxPatients,
                       Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                       method = "default", seed = current_seed, nSwarm = nSwarm, maxIter = maxIter
                     )

                     r2 <- PSO_design(
                       design = design, unified.u = unified.u, nlooks = nlooks, b1n = b1n, b1a = b1a, err1 = err1,
                       minPower = minPower, weight = weight, maxPatients = maxPatients,
                       Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                       method = "quantum", seed = current_seed, nSwarm = nSwarm, maxIter = maxIter
                     )

                     r3 <- PSO_design(
                       design = design, unified.u = unified.u, nlooks = nlooks, b1n = b1n, b1a = b1a, err1 = err1,
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
      b1n = b1n,  # Null hypothesis response rate
      b1a = b1a,  # Alternative hypothesis response rate
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
  if (exists("progress_bar", envir = .GlobalEnv)) {
    setTxtProgressBar(get("progress_bar", envir = .GlobalEnv), 101)
    close(get("progress_bar", envir = .GlobalEnv))
  }
  
  # Stop the cluster
  parallel::stopCluster(cl)
  foreach :: registerDoSEQ()

  res_final <- as.list( res_final)
  res_final[[1]] <- "GBOP2_minSS_singleE"
  class(res_final)<-"gbop2"
  
 
  
  return(res_final)
}














