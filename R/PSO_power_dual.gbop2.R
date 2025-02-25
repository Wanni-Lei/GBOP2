#' PSO: Power maximizing design with dual boundaries
#'
#' This function implements Particle Swarm Optimization (PSO) to find an power maximizing design with dual boundaries.
#'
#' @param method choose from "optimal" or "minimax"
#' @param totalPatients number of total patients
#' @param nlooks number of interim looks
#' @param Nmin_cohort1 minimum number of first cohort
#' @param Nmin_increase minimum number of increase in each cohort
#' @param b1n Null hypothesis response rate
#' @param b1a Alternative hypothesis response rate
#' @param err1 Type I error rate
#' @param minPower power
#' @param seed seed for pso
#' @param nSwarm nSwarm for pso
#' @param maxIter maxIter for pso
#'
#' @return A list on design parameters and operating characteristics
#' @examples 
#' PSO_power_dual(
#'   method = "default", 
#'   totalPatients = 50, 
#'   nlooks = 1, 
#'   Nmin_cohort1 = 10, 
#'   Nmin_increase = 10, 
#'   b1n = 0.2, 
#'   b1a = 0.4, 
#'   err1 = 0.05, 
#'   minPower = 0.8, 
#'   seed = 1024, 
#'   nSwarm = 64, 
#'   maxIter = 200
#' )
#'
#' @export
#'
#' @import globpso R6 Rcpp RcppArmadillo doParallel
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @importFrom foreach %dopar% foreach
#'


PSO_power_dual <- function(
    method = "default", # "quantum", "dexp"
    totalPatients = 50,
    nlooks = 1,
    Nmin_cohort1 = 10,
    Nmin_increase = 10,
    b1n = 0.2 ,# Null hypothesis response rate
    b1a = 0.4 , # Alternative hypothesis response rate
    err1 = 0.05, # Type I error rate
    minPower = 0.8,
    seed = 1024,
    nSwarm = 64,
    maxIter = 200){

  # library(globpso)
  # library(R6)
  # library(Rcpp)
  # library(RcppArmadillo)
  # source("BOP2_functions_twoboundaries_v2.R")
  # source("boundcode_equalrand_jsm.R")
  # Rcpp::sourceCpp(file="Calculation_twoboundaries_jsm.cpp",cacheDir="cache")
  numOfSimForTiralSetting = 10000   # Number of simulations


  id <- 1
  ## Fixed parameters -----
  input <- list(
    "b1n" = b1n,  # Null hypothesis response rate
    "b1a" = b1a,  # Alternative hypothesis response rate
    "err1" = err1,  # Type I error rate
    "minPower" = minPower,
    "seed" = seed
  )

  miniPatients <- Nmin_cohort1 + nlooks*Nmin_increase

  if(totalPatients < miniPatients){
    stop(paste0("Error: Please increase maxPatients to more than ", miniPatients  ))
  }



  ## Set cohort size -----
  cohortSize = function(N, R , w, n_min_cohort1 = Nmin_cohort1, n_min_incre = Nmin_increase){
    Nrest = N - R*Nmin_increase - Nmin_cohort1
    nobs = c()
    extra = 0
    for ( i in 1:(R+1)){
      if (i == 1){
        tmp = Nrest * w[i] + Nmin_cohort1
      } else {
        tmp = Nrest * w[i] + Nmin_increase + nobs[i-1]
      }
      extra = extra + round(Nrest * w[i])
      nobs = c(nobs, tmp)
    }
    extra = extra - Nrest
    nobs[which.max(w)] = nobs[which.max(w)] - extra

    if (nobs[length(nobs)] > N){
      nobs[length(nobs)] = N
    }

    return(nobs)
  }

  ## Build the utility function -----
  n_sim = 1

  #set.seed(123);
  set.seed(input$seed)
  seeds <- round(runif(10000)*10^8)

  ## PSO - comparison -----
  ## three lambda, three gamma, two delta. 8 parameters
  objf <- function(x, inputlist, fcn) {

    if (nlooks ==1){ ## when nlooks = 1
      b = x[1] ## lambaf
      b_grad1 = x[2] ## lambdae 1
      b_grad2 = x[3] ## lambdae 2
      pow1 = x[4] ## gamma 1
      pow2 = x[5] ## gamma 2
      pow3 = x[6] ## gamma 3
      delta1 = x[7]
      delta0 = x[length(x)]
      #n = x[8]
      w1 = x[8]
      w_list = c(w1, 1-w1)

    } else{
      b = x[1] ## lambaf
      b_grad1 = x[2] ## lambdae 1
      b_grad2 = x[3] ## lambdae 2
      pow1 = x[4] ## gamma 1
      pow2 = x[5] ## gamma 2
      pow3 = x[6] ## gamma 3
      delta1 = x[7]
      delta0 = x[length(x)]
      # n = x[8]

      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- x[8: (length(x)-1)] ## number of thetas
      w_list <- c()
      w_list[1] <- (cos(theta[1]))^2
      for( ii in 2: (n_cohort-1)){
        w_list[ii] <- (prod(sin(theta[1:(ii-1)]))*cos(theta[ii]))^2

      }
      w_list[n_cohort] <- (prod(sin(theta[1:(n_cohort-1)])))^2
    }

    if (!all.equal(sum(w_list), 1, tolerance = 1e-6)) {
      stop("Error: The sum of the elements in w_list must be approximately equal to 1.")
    }


    nobs.seq <- cohortSize(N = totalPatients, R = nlooks, w = w_list)
    nobs.seq[length(nobs.seq)] = totalPatients

    temp=.Call("_GBOP2_GetocBiRcpp_dual",seed=inputlist$seed, nsim=numOfSimForTiralSetting,
                     contrast=inputlist$contrast, nobs=nobs.seq,
                     b = b, b_grad1 = b_grad1, b_grad2 = b_grad2, pow1 = pow1, pow2 = pow2, pow3=pow3,
                     dprior = inputlist$dprior, ptrue = inputlist$b1a, phi = inputlist$b1n,
                     delta0 = delta0,delta1 = delta1,fff=fcn);

    t1e = temp[[2]]; # t1e
    power = temp[[3]]; # power

    ## optimal power
    if (t1e>inputlist$err1){
      results = 999
    }else { ## optimal_power
      results = -1*temp[[3]]
    }

    return(results)

  }

  p.n = input$b1n
  p.a = input$b1a

  inputlist = input
  inputlist$contrast = as.matrix(1)
  inputlist$dprior = c(inputlist$b1n, 1-inputlist$b1n)



  #low_bound <- c(0.01, 0.01, 0.01, 0, 0, 0, 0, 25, 0, 0, 0, 0)
  #upp_bound <- c(0.99, 0.99, 0.99, 1, 1, 1, 1, 50, pi/2, pi/2, pi/2, 1) # delta range from 0 to 1


  if(nlooks ==1){
    low_bound <- c(0.01, 0.01, 0.01, 0, 0, 0, 0, 0, 0)
    upp_bound <- c(0.99, 0.99, 0.99, 1, 1, 1, 1, 1, 1)

  }else{
    theta_L <- rep(0, nlooks) ## lower bound of theta
    theta_U <- rep(pi/2, nlooks) ## upper bound of theta
    low_bound <- c(0.01, 0.01, 0.01, 0, 0, 0, 0, theta_L, 0)
    upp_bound <- c(0.99, 0.99, 0.99, 1, 1, 1, 1, theta_U, 1)
  }





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
    res <- globpso(objFunc = objf, lower = low_bound, upper = upp_bound,
                   fixed = NULL, PSO_INFO = alg_setting,
                   inputlist = inputlist, fcn = maxresp_dual, seed=seeds[i])
    # print(res$par)
    # print(res$val)

    pars = res$par

    if (nlooks ==1){ ## when nlooks = 1
      b = pars[1]
      b_grad1 = pars[2]
      b_grad2 = pars[3]
      pow1 = pars[4]
      pow2 = pars[5]
      pow3 = pars[6]
      delta1 = pars[7]
      # n = pars[8]
      w_list <- c(res$par[,8], 1-res$par[,8])
      delta0 = pars[length(pars)]
    } else{ ## when nlooks >=2
      b = pars[1]
      b_grad1 = pars[2]
      b_grad2 = pars[3]
      pow1 = pars[4]
      pow2 = pars[5]
      pow3 = pars[6]
      delta1 = pars[7]
      # n = pars[8]
      delta0 = pars[length(pars)]
      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- pars[8: (length(pars)-1)] ## number of thetas
      w_list <- c()
      w_list[1] <- (cos(theta[1]))^2
      for( ii in 2: n_cohort-1){
        w_list[ii] <- (prod(sin(theta[1:(ii-1)]))*cos(theta[ii]))^2

      }
      w_list[n_cohort] <- (prod(sin(theta[1:(n_cohort-1)])))^2
    }


    nobs2 <- cohortSize(N = totalPatients, R = nlooks, w = w_list)
    nobs2[length(nobs2)] = totalPatients

    bd = t(getboundary_dual(dprior=c(p.n, 1-p.n),contrast=as.matrix(1),
                       nobs=nobs2,b=b,b_grad1=b_grad1,b_grad2=b_grad2,
                       pow1=pow1, pow2 = pow2, pow3=pow3,
                       phi=input$b1n,delta0=delta0,delta1=delta1)); bd
    power_result = .Call("_GBOP2_GetocBiRcpp_dual",seed=input$seed, nsim=numOfSimForTiralSetting,
                               contrast=as.matrix(1), nobs=(nobs2),
                               b=b, b_grad1=b_grad1,b_grad2=b_grad2,
                               pow1=pow1, pow2 = pow2, pow3=pow3,
                               dprior= c(p.n,1-p.n), ptrue=p.a, phi=p.n,
                               delta0=delta0,delta1=delta1, fff=maxresp_dual)
    expected_sample = power_result[[4]]

  }

  design <- "optimal"
  results_list <- list("function" = "PSO_power_dual",
   "design"= design,"method"=method,
  "parameter" = list("lambda1" = b,"lambda_grad1" = b_grad1,"lambda_grad2" = b_grad2,
                    "gamma_1" = pow1, "gamma_2" =pow2, "gamma_3" =pow3,
                  "delta0" = delta0, "delta1" = delta1),
                  "cohort" = as.list(nobs2),"boundary"=list("1" = bd[, 2], "2" = bd[, 3]),
                  "Type I Error" = power_result[[2]],
                  "Power" = power_result[[3]],
                  "Expected Sample Size" = expected_sample,
                  "Utility" = res$val
  )


  class(results_list)<-"gbop2"

  return(results_list)


}




