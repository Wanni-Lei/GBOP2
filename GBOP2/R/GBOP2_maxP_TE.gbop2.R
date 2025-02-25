#' PSOGO: Power maximizing design with efficacy and toxicity boundaries
#'
#' This function implements PSOGO to find a power maximizing design with efficacy and toxicity boundaries.
#'
#' @name GBOP2_maxP_TE
#' @param design fixed as "optimal", cannot be modified by user
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
#' @param seed  Random seed for reproducibility
#' @param nCore number of core
#'
#' @return A list on design parameters and operating characteristics
#' @examples
#'   GBOP2_maxP_TE(
#'     design = "optimal",
#'     nlooks = 1,
#'     skip_efficacy = NULL,
#'     skip_toxicity = NULL,
#'     totalPatients = 20,
#'     Nmin_cohort1 = 10,
#'     Nmin_increase = 5,
#'     nParallel = 3,
#'     e1n = 0.15, e2n = 0.16, e3n = 0.024,
#'     e1a = 0.4, e2a = 0.08, e3a = 0.032,
#'     err_eff = 1, err_tox = 1, err_all = 0.1,
#'     power_eff = 0.8, power_tox = 0.8, power_all = 0.8,
#'     seed = 5321, pso_method = "all",
#'     nSwarm = 32, maxIter = 100, nCore = 8
#'   )
#'   
#' 
#' @export
#' @import globpso R6 Rcpp RcppArmadillo doParallel
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @importFrom foreach %dopar% foreach
#' @importFrom parallel makePSOCKcluster 
#' @importFrom doParallel registerDoParallel

utils::globalVariables(c("Power", "Corhort", "boundary.1", "boundary.2"))


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
                          seed = 5321,
                          pso_method = "all",
                          nSwarm = 32,
                          maxIter = 100,
                          nCore = 4) {




  ##########################################################################
  ## estimated total time
  cat("\nGBOP2 process has started...\n")
  start_time <- Sys.time()  # Start timing
  
  one_task <- PSO_power_TE(method = "default",nlooks = nlooks,
                           skip_efficacy = NULL, skip_toxicity = NULL ,totalPatients = totalPatients,
                           Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                           e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,  e3a = e3a, # Ha for Eff and Tox
                           err_eff = err_eff,  # Type I error rate: Efficacious but toxic
                           err_tox = err_tox,  err_all = err_all,  # Type I error rate: Futile and toxic
                           power_eff = power_eff,
                           power_tox = power_tox, power_all = power_all, seed = seed, nSwarm = nSwarm,maxIter = 1
  )
  
  end_time <- Sys.time()  # End timing
  execution_time1T <- as.numeric(end_time - start_time)  # Convert to numeric (seconds)
  
  # Step 2: Estimate total execution time
  N_PSO <- nParallel * 3  # Total number of PSO_design calls
  total_time <- (N_PSO * execution_time1T * maxIter) / nCore  # Estimate based on parallel cores
  
  # Step 3: Display estimated total execution time
  cat("\nEstimated total execution time:", round(total_time, 2), "seconds\n")
  cat("Or approximately:", round(total_time / 60, 2), "minutes\n")
  
  ####################################################################
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
  
  
  
  
  #####################################################################

  # Set up parallel computing
  cl <- makePSOCKcluster(nCore)
  registerDoParallel(cl)

  # Define the seed list
  #set.seed(123)
  input <- list("seed" = seed)
  set.seed(input$seed)
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
  r1 <- PSO_power_TE(method = "default",nlooks = nlooks,
                  skip_efficacy = NULL, skip_toxicity = NULL ,totalPatients = totalPatients,
                  Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                  e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,  e3a = e3a, # Ha for Eff and Tox
                  err_eff = err_eff,  # Type I error rate: Efficacious but toxic
                  err_tox = err_tox,  err_all = err_all,  # Type I error rate: Futile and toxic
                  power_eff = power_eff,
                  power_tox = power_tox, power_all = power_all, seed = seed, nSwarm = nSwarm,maxIter = maxIter
     )



                     r2 <-  PSO_power_TE(method = "quantum",nlooks = nlooks,
                                         skip_efficacy = NULL, skip_toxicity = NULL ,totalPatients = totalPatients,
                                         Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                                         e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,  e3a = e3a, # Ha for Eff and Tox
                                        err_eff = err_eff,  # Type I error rate: Efficacious but toxic
                                         err_tox = err_tox,  err_all = err_all,  # Type I error rate: Futile and toxic
                                         power_eff = power_eff,
                                         power_tox = power_tox, power_all = power_all, seed = current_seed, nSwarm = nSwarm,maxIter = maxIter
                     )

                     r3 <- PSO_power_TE(method = "dexp",nlooks = nlooks,
                                        skip_efficacy = NULL, skip_toxicity = skip_toxicity ,totalPatients = totalPatients,
                                        Nmin_cohort1 = Nmin_cohort1, Nmin_increase = Nmin_increase,
                                        e1n = e1n, e2n = e2n, e3n = e3n, e1a = e1a, e2a = e2a,  e3a = e3a, # Ha for Eff and Tox
                                        err_eff = err_eff,  # Type I error rate: Efficacious but toxic
                                        err_tox = err_tox,  err_all = err_all,  # Type I error rate: Futile and toxic
                                        power_eff = power_eff,
                                        power_tox = power_tox, power_all = power_all, seed = current_seed, nSwarm = nSwarm,maxIter = maxIter
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
                     
                     
                     listname <- c("function","method", "lambdae1",
                                   "lambdae2", "lambdat1", "lambdat2", "gamma" , cohort_name,      
                                   boudary_name,"typeI_01", "typeI_10"       
                                   ,"typeI_00", "Power", "Utility" )
                     colnames(r1) <- listname
                     colnames(r2) <- listname
                     colnames(r3) <- listname
                     r_ensemble <- rbind(r1, r2,r3)
                     
                     r_ensemble <- r_ensemble |>
                       filter(Utility == min(Utility)) |>
                       filter(Power == max(Power))
                     
                     r_ensemble[[1]] <- "GBOP2_minSS_TE"
                     results <- r_ensemble
                     return(results)
                 } 
  
  
  res_final <- res |>
    distinct(Utility, .keep_all = TRUE) |>
    filter(Utility == min(Utility)) |>
    filter(Power == max(Power))
  
} else {
                     r <- PSO_power_TE(
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
  res_final[[1]] <- "GBOP2_maxP_TE"
} else{
  res_final[[1]] <- "PSO_power_TE" 
}

class(res_final)<-"gbop2"
return(res_final)
}



