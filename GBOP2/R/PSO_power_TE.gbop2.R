#' PSO for maximizing power design with efficacy and toxicity boundaries
#'
#' @param design choose from "optimal" or "minimax"
#' @param method method for single PSO, choose from "default", "quantum" or "dexp"
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
#' @param seed seed for PSO
#' @param nSwarm nSwarm in PSO
#' @param maxIter maxIter in PSO
#'
#' @return A list on design parameters and operating characteristics
#' @export
#' @import globpso R6 Rcpp RcppArmadillo
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct



PSO_power_TE <- function(design = "optimal",
                         method = "default",
                         nlooks = 4,
                         skip_efficacy = NULL,
                         skip_toxicity = NULL,
                         totalPatients = 50,
                         Nmin_cohort1 = 10,
                         Nmin_increase = 5,
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
                         nSwarm = 32,
                         maxIter = 100
                         ){


  # library(globpso)
  # library(R6)
  # library(Rcpp)
  # library(RcppArmadillo)
  # source("BOP2_functions_EffTox.R")
  # source("BOP2_TE_function.R")
  # source("boundcode.R")
  # Rcpp::sourceCpp(file="Calculation2_original.cpp")
  # Rcpp::sourceCpp(file="Calculation_minimizeN.cpp",cacheDir="cache")
  numOfSimForTiralSetting = 10000

  ## Fixed parameters
   

  input <- list(
    e1n = 0.15,  # H0 for Eff
    e2n = 0.16,  # H0 for Tox
    e3n = 0.024, # H0 for Eff and Tox
    e1a = 0.4,  # Ha for Eff
    e2a = 0.08,  # Ha for Tox
    e3a = 0.032, # Ha for Eff and Tox
    ETprior1 = 0.45, # prior for Pr(Eff)
    ETprior2 = 0.3,  # prior for Pr(Tox)
    ETprior3 = 0.15, # prior for Pr(Eff&Eff)
    ETpriorSS = 1, # Prior effective sample size
    priorspec = TRUE,  # Indicates whether a prior specification is provided
    err_eff = 1,  # Type I error rate: Efficacious but toxic
    err_tox = 1,  # Type I error rate: Safe but futile
    err_all = 0.1,  # Type I error rate: Futile and toxic
    truet1e = TRUE, # strictly control type I error <= nominal value
    power_eff = 0.8,
    power_tox = 0.8,
    power_all = 0.8,
    numOfSimForTiralSetting = 10000,  # Number of simulations for trial settings
    seed = 1024,  
    # cohortInputs = "10 20 35 50",
    cohortInputs = "40 50",
    simSeed = 1024,  # Seed for simulations
    simunumber = 10000  # Number of simulations
  )

  # if(!is.null(skip_efficacy) && !is.null(skip_toxicity)){
  #   for(i in 1: nlooks){
  #     if (skip_efficacy[i] == 1 & skip_toxicity[i] ==1){
  #       stop("Error: Cannot skip both efficacy and toxicity at the same interim look.")
  #     }
  #   }
  # 
  # }


  miniPatients <- Nmin_cohort1 + nlooks*Nmin_increase
  
  if(totalPatients < miniPatients){
    stop(paste0("Error: Please increase totalPatients to more than ", miniPatients  ))
  }
  

  cohortSize = function(N, R, w, n_min_cohort1 = Nmin_cohort1, n_min_incre = Nmin_increase){
    # fix first stage to 10
    Nrest = N - R*n_min_incre - n_min_cohort1
    nobs = c(n_min_cohort1)
    extra = 0
    for ( i in 2:(R+1)){
      tmp = Nrest * w[i-1] + n_min_incre + nobs[i-1]
      extra = extra + round(Nrest * w[i-1])
      nobs = c(nobs, tmp)
    }
    # extra = extra - Nrest
    # w2 = w[-length(w)]
    # nobs[which.max(w2) + 1] = nobs[which.max(w2) + 1] - extra
    return(nobs)
  }
  



  r0 = input$e1n
  t0 = input$e2n
  r1 = input$e1a
  t1 = input$e2a
  
  scen = hypotheses_ind(r0, r1, t0, t1)
  input$scenario = scen
  
  inputlist <- input
  
  inputlist$R <- nlooks
  



  ## Build the utility function -----
  
  objf <- function(x, inputlist) {
    
    if (nlooks ==1){
      le = x[1]
      lt = x[2]
      g = x[3]
      # n = x[4]
      le2 = x[4]
      w1 = x[5]
      lt2 = x[6]
      w_list = c(w1, 1-w1)
      
    }else{
      le = x[1]
      lt = x[2]
      g = x[3]
      # n = x[4]
      le2 = x[4]
      
      lt2 = x[length(x)]
      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- x[5: (length(x)-1)] ## number of thetas
      w_list <- c()
      w_list[1] <- (cos(theta[1]))^2
      for( ii in 2: (n_cohort-1)){
        w_list[ii] <- (prod(sin(theta[1:(ii-1)]))*cos(theta[ii]))^2
        
      }
      w_list[n_cohort] <- (prod(sin(theta[1:(n_cohort-1)])))^2
    }
    
    
    if (round(sum(w_list)) != 1) {
      stop("Error: The sum of the elements in w_list must be equal to 1.")
    }
    
    nobs.seq <- cohortSize(N = totalPatients, R = inputlist$R, w = w_list, n_min_cohort1 = 10, n_min_incre = 8)
    nobs.seq[length(nobs.seq)] = totalPatients
    interm = interm.eff = interm.tox = nobs.seq
    
    boundary = get_boundarycpp_2lambda(interm.eff= interm.eff,interm.tox =interm.tox, 
                                       lambda_e=le, lambda_t=lt, 
                                       lambda_e2=le2, lambda_t2=lt2, gamma = g,
                                       prior=inputlist$scenario[1,], r0 = inputlist$e1n, t0 = input$e2n)
    
    interm = interm.eff = interm.tox = ceiling(cohortSize(N = totalPatients, R = inputlist$R, w = w_list))
    interm[length(interm)] = totalPatients
    nintm = length(interm)
    new_patient = c(0, interm)
    npt=rep(NA, nintm)
    for(j in seq(nintm)){
      npt[j] = new_patient[j+1] - new_patient[j]
    }
    
    bound_eff = boundary[[1]]
    bound_tox = boundary[[2]]
    r0 = inputlist$e1n
    t0 = inputlist$e2n
    r1 = inputlist$e1a
    t1 = inputlist$e2a
    temp_pe = exact_error_recursive2_Rcpp(interm.eff, bound_eff, r0, r1, 5)
    temp_pt = exact_error_recursive2_Rcpp(interm.tox, interm-bound_tox, 1-t0, 1-t1, 5)
    a00 = temp_pe$t1err * temp_pt$t1err
    a01 = temp_pe$t1err * temp_pt$power
    a10 = temp_pe$power * temp_pt$t1err
    a11 = temp_pe$power * temp_pt$power
    
    if (a00 > inputlist$err_all){
      result = 999
    } else {
      result = -a11
    }
    
    return(result)
    
  }
  
  
  # low_bound <- c(0.5, 0.5, 0, 0.5, 0.5, 0, 0, 0)
  # upp_bound <- c(0.99, 0.99, 1, 0.99, 0.99, pi/2, pi/2, pi/2)
  
  if(nlooks ==1){
    low_bound <- c(0.5, 0.5, 0, 0.5, 0, 0.5)
    upp_bound <- c(0.99, 0.99, 1,  0.99, 1, 0.99)
  } else{
    theta_L <- rep(0, nlooks) ## lower bound of theta
    theta_U <- rep(pi/2, nlooks) ## upper bound of theta
    low_bound <- c(0.5, 0.5, 0,  0.5, theta_L, 0.5)
    upp_bound <- c(0.99, 0.99, 1,  0.99, theta_U, 0.99)
  }
  
  
  
  n_sim = 1
  set.seed(123);
  seeds <- round(runif(10000)*10^8)
  
  ## PSO - comparison -----
  
  # alg_setting <- getPSOInfo(freeRun = 1, nSwarm = 32, maxIter = 100)
  if (method == "default"){
    ## default
    ## getPSOInfo:Create a list with PSO parameters for Minimization.
    alg_setting <- getPSOInfo(freeRun = 1, nSwarm = nSwarm, maxIter=maxIter) # default if "basic" Linearly Decreasing Weight PSO
  } else if (method == "quantum"){
    ## quantum:
    alg_setting <- getPSOInfo(psoType = "quantum", freeRun = 1, nSwarm = nSwarm, maxIter=maxIter)
  } else {
    alg_setting <- getPSOInfo(psoType = "dexp", freeRun = 1, nSwarm = nSwarm, maxIter = maxIter)
  }
  
  
  
  
  for ( i in 1){
    # print(paste("seeds:", seeds[i]))
    res <- globpso(objFunc = objf, lower = low_bound, upper = upp_bound,
                   fixed = NULL, PSO_INFO = alg_setting,
                   inputlist = inputlist, seed = seeds[i])
    # print(res$par)
    # print(res$val)
    pars = res$par
    
    if (nlooks ==1){
      le = pars[1]
      lt = pars[2]
      g = pars[3]
      # n = pars[4]
      le2 = pars[4]
      w1 = pars[5]
      lt2 = pars[6]
      w_list = c(w1, 1-w1)
      
    }else{
      le = pars[1]
      lt = pars[2]
      g = pars[3]
      # n = pars[4]
      le2 = pars[4]
      
      lt2 = pars[length(pars)]
      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- pars[5: (length(pars)-1)] ## number of thetas
      w_list <- c()
      w_list[1] <- (cos(theta[1]))^2
      for( ii in 2: (n_cohort-1)){
        w_list[ii] <- (prod(sin(theta[1:(ii-1)]))*cos(theta[ii]))^2
        
      }
      w_list[n_cohort] <- (prod(sin(theta[1:(n_cohort-1)])))^2
    }
    
    
    if (round(sum(w_list)) != 1) {
      stop("Error: The sum of the elements in w_list must be equal to 1.")
    }
    
    nobs.seq <- cohortSize(N = totalPatients, R = inputlist$R, w = w_list, n_min_cohort1 = 10, n_min_incre = 8)
    nobs.seq[length(nobs.seq)] = totalPatients
    interm = interm.eff = interm.tox = nobs.seq
    
    boundary = get_boundarycpp_2lambda(interm.eff= interm.eff,interm.tox =interm.tox, 
                                       lambda_e=le, lambda_t=lt, 
                                       lambda_e2=le2, lambda_t2=lt2, gamma = g,
                                       prior=inputlist$scenario[1,], r0 = inputlist$e1n, t0 = input$e2n)
    
    interm = interm.eff = interm.tox = ceiling(cohortSize(N = totalPatients, R = inputlist$R, w = w_list))
    interm[length(interm)] = totalPatients
    
    nintm = length(interm)
    new_patient = c(0, interm)
    npt=rep(NA, nintm)
    for(j in seq(nintm)){
      npt[j] = new_patient[j+1] - new_patient[j]
    }
    
    bound_eff = boundary[[1]]
    bound_tox = boundary[[2]]
    temp_pe = exact_error_recursive2_Rcpp(interm.eff, bound_eff, r0, r1, 5)
    temp_pt = exact_error_recursive2_Rcpp(interm.tox, interm-bound_tox, 1-t0, 1-t1, 5)
    a00 = temp_pe$t1err * temp_pt$t1err
    a01 = temp_pe$t1err * temp_pt$power
    a10 = temp_pe$power * temp_pt$t1err
    a11 = temp_pe$power * temp_pt$power
    
    Nt = temp_pt$pts
    Ne = temp_pe$pts
    Nte = temp_pe$pts+temp_pt$pts
    
    # default_tbl = rbind(default_tbl, c((res$cputime), res$par,
    #                                    interm, boundary[[1]],boundary[[2]], a00, a01, a10, a11,
    #                                    Nt,Ne, Nte, -res$val))
    
    
    
  }
  
  results_list <- list(
    "function" = "PSO_power_TE",
    "design" = design,"method" = method,
    "parameter" = list(
      "lambdae1" = le, "lambdae2" = le2,
      "lambdat1" = lt,
      "lambdat2" = lt2,
      "gamma" = g),
    "cohort" = as.list(interm),                  # Cohort sizes
    "boundary_effi" = as.list(boundary$boundary.eff),   # Boundary for efficacy
    "boundary_toxi" = as.list(boundary$boundary.tox),   # Boundary for toxicity
    # "expected_sample" =  expected_sample,
    "typeI_01" = a01,
    "typeI_10" = a10,
    "typeI_00" = a00,
    "power" = a11,
    "utility" = -res$val
  )
  
  
  
  results_list
  class(results_list)<-"gbop2"
  return(results_list)
}





