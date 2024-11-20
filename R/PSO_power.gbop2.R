#' PSO maximizing power design with single boundary
#'
#' @param nlooks number of looks
#' @param totalPatients total number of patients
#' @param Nmin_cohort1 minimum number of first cohort
#' @param Nmin_increase minimum number of increase in each cohort
#' @param method method for single PSO, choose from "default", "quantum" or "dexp"
#' @param b1n Null hypothesis response rate
#' @param b1a Alternative hypothesis response rate
#' @param err1 Type I error rate
#' @param minPower power
#' @param seed seed for pso
#' @param nSwarm nSwarm for pso
#' @param maxIter maxIter for pso
#'
#' @return A list
#' @export
#' @import globpso R6 Rcpp RcppArmadillo
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @example 
#' PSO_power(
#'   nlooks = 1,
#'    totalPatients = 50, ## total patients in optimal_power
#'    Nmin_cohort1 = 10,
#'    Nmin_increase = 5,
#'    method = "default", # "quantum", "dexp"
#'    b1n = 0.2,  # Null hypothesis response rate
#'    b1a = 0.4,  # Alternative hypothesis response rate
#'    err1 = 0.05,  # Type I error rate
#'    minPower = 0.8, ## power
#'    seed = 1024,   ## set seed to calculate OC
#'    nSwarm = 64,
#'    maxIter = 200
#'    )
#'    b1a = 0.4,  # Alternative hypothesis response rate
#'    err1 = 0.05,  # Type I error rate
#'    minPower = 0.8, ## power
#'    seed = 1024,   ## set seed to calculate OC
#'    nSwarm = 64,
#'    maxIter = 200
#'    )

PSO_power <- function(
    nlooks = 1,
    totalPatients = 50, ## total patients in optimal_power
    Nmin_cohort1 = 10,
    Nmin_increase = 5,
    method = "default", # "quantum", "dexp"
    b1n = 0.2,  # Null hypothesis response rate
    b1a = 0.4,  # Alternative hypothesis response rate
    err1 = 0.05,  # Type I error rate
    minPower = 0.8, ## power
    seed = 1024,   ## set seed to calculate OC
    nSwarm = 64,
    maxIter = 200
){


  numOfSimForTiralSetting = 10000   # Number of simulations
  # library(globpso)
  # library(R6)
  # library(Rcpp)
  # library(RcppArmadillo)
  # source("boundcode_twolambda.R") # for two lambdas
  # Rcpp::sourceCpp(file="Calculation_minimizeN_twolambda_update.cpp",cacheDir="cache")


  input <- list(
    "b1n" = b1n,  # Null hypothesis response rate
    "b1a" = b1a,  # Alternative hypothesis response rate
    "err1" = err1,  # Type I error rate
    "minPower" = minPower, ## power
    "seed" = seed   ## set seed to calculate OC
  )


  miniPatients <- Nmin_cohort1 + nlooks*Nmin_increase

  if(totalPatients < miniPatients){
    stop(paste0("Error: Please increase maxPatients to more than ", miniPatients  ))
  }


  ## Set cohort size -----
  cohortSize = function(N, R, w, n_min_cohort1 = Nmin_cohort1, n_min_incre = Nmin_increase){
    Nrest = N - R*n_min_incre - n_min_cohort1 ## R interim looks
    nobs = c()
    extra = 0
    for ( i in 1:(R+1)){ ## R+1 cohorts
      if (i == 1){
        tmp = Nrest * w[i] + n_min_cohort1
      } else {
        tmp = Nrest * w[i] + n_min_incre + nobs[i-1]
      }
      extra = extra + round(Nrest * w[i])
      nobs = c(nobs, tmp)
    }
    extra = extra - Nrest
    nobs[which.max(w)] = nobs[which.max(w)] - extra
    return(nobs)
  }


  ## Build the utility function -----
  objf <- function(x, inputlist, fcn) { ## x is parameter set

    if (nlooks ==1){ ## when nlooks = 1
      b = x[1]
      pow = x[2]
      w1 = x[3] ## used to calculate cohortSize
      b2 = x[4] ## pars[1, 2,5] used to calculate boundary
      w_list = c(w1, 1-w1)

    } else{ ## when nlooks >=2
      b = x[1]
      b2 = x[length(x)]
      pow = x[2]
      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- x[3: (length(x)-1)] ## number of thetas
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


    nobs.seq <- cohortSize(N = totalPatients, R = nlooks, w = w_list)
    ## make sure the last cohort is the total sample size
    # nobs.seq[length(nobs.seq)] = totalPatients

    temp=.Call("_GBOP2_GetocBiRcpp", seed=inputlist$seed, nsim = n_sim,
                     contrast=inputlist$contrast, nobs=nobs.seq,
                     b = b, b2 = b2, pow2 = pow, dprior = inputlist$dprior, ptrue = inputlist$b1a, phi = inputlist$b1n, fff=fcn);

    t1e = temp[[2]]; # t1e
    power = temp[[3]]; # power


    ## optimal power
    if (t1e>inputlist$err1){
      results = 999
    }else { ## optimal_power
      results = -1*temp[[3]]
    }


    return(results)
  }  ## end of objt function

  ## lambda1, gamma, n, theta(# = nlooks), lambda2
  if(nlooks ==1){
    low_bound <- c(0.5, 0, 0, 0.5)
    upp_bound <- c(0.99, 1, 1, 0.99)

  }else{
    theta_L <- rep(0, nlooks) ## lower bound of theta
    theta_U <- rep(pi/2, nlooks) ## upper bound of theta
    low_bound <- c(0.5, 0, theta_L, 0.5)
    upp_bound <- c(0.99, 1, theta_U, 0.99)
  }



  p.n = input$b1n
  p.a = input$b1a


  inputlist = input
  inputlist$cutstart = 1
  inputlist$func = maxresp
  inputlist$contrast = as.matrix(1)
  inputlist$dprior = c(inputlist$b1n, 1-inputlist$b1n)


  n_sim = 1
  set.seed(123)
  seeds <- round(runif(10000)*10^8)

  ## PSO - comparison -----
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


  for ( i in 1:n_sim){
    res <- globpso(objFunc = objf, lower = low_bound, upper = upp_bound,
                   fixed = NULL, PSO_INFO = alg_setting,
                   inputlist = inputlist, fcn = maxresp, seed = seeds[i])

    pars = res$par

    if (nlooks ==1){ ## when nlooks = 1
      b = pars[1]
      pow = pars[2]
      w1 = pars[3] ## used to calculate cohortSize
      b2 = pars[4] ## pars[1, 2,5] used to calculate boundary
      w_list <- c(res$par[,3], 1-res$par[,3])
    } else{ ## when nlooks >=2
      b = pars[1]
      b2 = pars[length(pars)]
      pow = pars[2]
      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- pars[3: (length(pars)-1)] ## number of thetas
      w_list <- c()
      w_list[1] <- (cos(theta[1]))^2
      for( ii in 2: (n_cohort-1)){
        w_list[ii] <- (prod(sin(theta[1:(ii-1)]))*cos(theta[ii]))^2

      }
      w_list[n_cohort] <- (prod(sin(theta[1:(n_cohort-1)])))^2
    }
    ## cohort size
    nobs2 = cohortSize(N = totalPatients, R = nlooks, w = w_list)


    ## bd: boundary
    bd = t(getboundary(dprior=c(p.n, 1-p.n),contrast=as.matrix(1),
                       nobs=nobs2,b=res$par[,1],b2 = b2, pow=res$par[,2],phi=input$b1n))

    ##power_default[[2]] is type I error, [[3]] is power
    power_result = .Call("_GBOP2_GetocBiRcpp", seed=input$seed,  nsim = n_sim, contrast=as.matrix(1), nobs=(nobs2),
                               b=res$par[,1], b2= b2, pow2=res$par[,2],
                               dprior= c(p.n,1-p.n), ptrue=p.a, phi=p.n, maxresp)




    # power_result[[4]] is expected sample size under null
    design <- "Optimal"
    # results_tbl <- c(design, method,
    #                  res$cputime, round(b, 4), round(b2, 4), round(pow, 4),
    #                  nobs2, bd[,2],
    #                 round( power_result[[2]],4), round(power_result[[3]],4), round(power_result[[4]],4), round(res$val,4))

    result_list <- list("function" = "PSO_power", "design" = design, "method" = method,
                        "cputime" = res$cputime,
                        "lambda1" = round(b, 4), "lambda2" = round(b2, 4), "gamma" = round(pow, 4),
                        "cohort" = as.list(nobs2), "boundary" = as.list(bd[,2]),
                        "TypeI" = round( power_result[[2]],4), "Power" = round(power_result[[3]],4), "EN(P0)" = round(power_result[[4]],4), "Utility" = round(res$val,4))


  }



  class(result_list)<-"gbop2"
  return(result_list)



}




