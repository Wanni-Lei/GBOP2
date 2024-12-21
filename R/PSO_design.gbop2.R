#' PSO for optimal/minimax design with single boundary
#'
#' @param design choose from "optimal" or "minimax"
#' @param nlooks number of interim looks
#' @param b1n Null hypothesis response rate
#' @param b1a Alternative hypothesis response rate
#' @param err1 Type I error rate
#' @param minPower power
#' @param weight weight of sample size under null
#' @param maxPatients maximum number of patients
#' @param Nmin_cohort1 minimum number of first cohort
#' @param Nmin_increase minimum number of increase in each cohort
#' @param method method for single PSO, choose from "default", "quantum" or "dexp"
#' @param seed seed for pso
#' @param nSwarm nSwarm for pso
#' @param maxIter maxIter for pso
#'
#' @return A list
#' @import  globpso R6 Rcpp RcppArmadillo
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct
#' @useDynLib GBOP2, .registration = TRUE
#'
#' @export

 

PSO_design <- function(
    design = "optimal", #"minimax"
    nlooks = 1,
    b1n = 0.2,  # Null hypothesis response rate
    b1a = 0.4,  # Alternative hypothesis response rate
    err1 = 0.05,  # Type I error rate
    minPower = 0.8, ## power
    weight = 1, ## weight of sample size under null
    maxPatients = 50,  ## maximum number of patients
    Nmin_cohort1 = 10,
    Nmin_increase = 5,
    method = "default", # "quantum", "dexp"
    seed = 1024,   ## set seed to calculate OC
    nSwarm = 64,
    maxIter = 200){



  miniPatients <- Nmin_cohort1 + nlooks*Nmin_increase


  if(maxPatients < miniPatients){
    stop(paste0("Error: Please increase maxPatients to more than ", miniPatients  ))
  }

  numOfSimForTiralSetting = 10000   # Number of simulations
  # library(globpso)
  # library(R6)
  # library(Rcpp)
  # library(RcppArmadillo)
  # source("boundcode_twolambda.R") # for two lambdas
  # Rcpp::sourceCpp(file="Calculation_minimizeN_twolambda_update.cpp",cacheDir="cache")

  # Rcpp::sourceCpp(file="original_Calculation_minimizeN_twolambda.cpp",cacheDir="cache")
  #
  input <- list(
    "b1n" = b1n,  # Null hypothesis response rate
    "b1a" = b1a,  # Alternative hypothesis response rate
    "err1" = err1,  # Type I error rate
    "minPower" = minPower, ## power
    "seed" = seed   ## set seed to calculate OC
  )



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

    if (nobs[length(nobs)] > N){
      nobs[length(nobs)] = N
    }

    return(nobs)
  }


  ## Build the utility function -----
  objf <- function(x, inputlist, fcn) { ## x is parameter set

    if (nlooks ==1){ ## when nlooks = 1
      b = x[1]
      pow = x[2]
      n = x[3]
      w1 = x[4]
      b2 = x[5]
      w_list = c(w1, 1-w1)

    } else{ ## when nlooks >=2
      b = x[1]
      b2 = x[length(x)]
      pow = x[2]
      n = x[3]
      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- x[4: (length(x)-1)] ## number of thetas
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


    nobs.seq <- cohortSize(N = n, R = nlooks, w = w_list)
    ## should exact probability rather than simulation
    ## ge
    temp=GetocBiRcpp(seed=inputlist$seed, nsim = n_sim,
                     contrast=inputlist$contrast, nobs=nobs.seq,
                     b = b, b2 = b2, pow2 = pow, dprior = inputlist$dprior, ptrue = inputlist$b1a, phi = inputlist$b1n, fff=fcn); ## ptrue is under alternative, phi is under null

    t1e = temp[[2]]; # t1e
    power = temp[[3]]; # power


    if (t1e>inputlist$err1 | power<inputlist$minPower){
      results = 999
    } else{

      if(design =="optimal"){ ## optimal
        results = weight*temp[[4]] + (1-weight)*temp[[5]]

      }else{ ## minimax
        n_final = nobs.seq[length(nobs.seq)] ## n_final is total sample size
        results = n_final + (weight*temp[[4]] + (1-weight)*temp[[5]])/n_final
      }
    }


    # results = temp[[4]] + power_accept + t1e_accept
    return(results)
  }  ## end of objt function





  ## lambda1, gamma, n, theta(# = nlooks), lambda2
  if(nlooks ==1){
    low_bound <- c(0.5, 0,miniPatients, 0, 0.5)
    upp_bound <- c(0.99, 1, maxPatients, 1, 0.99)

  }else{
    theta_L <- rep(0, nlooks) ## lower bound of theta
    theta_U <- rep(pi/2, nlooks) ## upper bound of theta
    low_bound <- c(0.5, 0, miniPatients, theta_L, 0.5)
    upp_bound <- c(0.99, 1, maxPatients, theta_U, 0.99)
  }



  p.n = input$b1n
  p.a = input$b1a


  inputlist = input
  inputlist$cutstart = 1
  inputlist$func = maxresp
  inputlist$contrast = as.matrix(1)
  inputlist$dprior = c(inputlist$b1n, 1-inputlist$b1n)


  n_sim <- 1
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
                   inputlist = inputlist, fcn = maxresp, seed = seeds[i]
                   # ,init = c(0.75, 0.5, maxPatients, rep(0.5, nlooks), 0.75)
    )
    #upp_bound <- c(0.99, 1, maxPatients, theta_U, 0.99)
    ## init argument for sample size as maxPatients and other as null

    # print(res$par)
    # print(res$val)
    pars = res$par

    if (nlooks ==1){ ## when nlooks = 1
      b = pars[1]
      pow = pars[2]
      n = pars[3] ## total sample size
      w1 = pars[4] ## used to calculate cohortSize
      b2 = pars[5] ## pars[1, 2,5] used to calculate boundary
      w_list <- c(res$par[,4], 1-res$par[,4])
    } else{ ## when nlooks >=2
      b = pars[1]
      pow = pars[2]
      n = pars[3]
      b2 = pars[length(pars)]
      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- pars[4: (length(pars)-1)] ## number of thetas
      w_list <- c()
      w_list[1] <- (cos(theta[1]))^2
      for( ii in 2: n_cohort-1){
        w_list[ii] <- (prod(sin(theta[1:(ii-1)]))*cos(theta[ii]))^2

      }
      w_list[n_cohort] <- (prod(sin(theta[1:(n_cohort-1)])))^2
    }
    ## cohort size
    nobs2 = cohortSize(N = res$par[,3], R = nlooks, w = w_list)


    ## bd: boundary
    bd = t(getboundary(dprior=c(p.n, 1-p.n),contrast=as.matrix(1),
                       nobs=nobs2,b=res$par[,1],b2 = b2, pow=res$par[,2],phi=input$b1n))

    ##power_default[[2]] is type I error, [[3]] is power
    # power_result = GetocBiRcpp(seed=input$seed,  nsim = n_sim, contrast=as.matrix(1), nobs=(nobs2),
    #                            b=res$par[,1], b2= b2, pow2=res$par[,2],
    #                            dprior= c(p.n,1-p.n), ptrue=p.a, phi=p.n, maxresp)
    power_result = .Call("_GBOP2_GetocBiRcpp", seed=input$seed,  nsim = n_sim, contrast=as.matrix(1), nobs=(nobs2),
                         b=res$par[,1], b2= b2, pow2=res$par[,2],
                         dprior= c(p.n,1-p.n), ptrue=p.a, phi=p.n, maxresp)

    ## Calculate expected sample size
    expected_sample <- weight*power_result[[4]] + (1-weight)*power_result[[5]]



    results_list <- list("function" = "PSO_design",
                     "design" = design, "weight" = weight, "method" = method,
                     "cputime" = res$cputime, "parameter" = list("lambda1" = b, "lambda2" = b2, "gamma" = pow),
                     "cohort" = as.list(nobs2),"boundary"= as.list(bd[,2]),
                     "Type I Error" = power_result[[2]],
                     "Power" = power_result[[3]], "Expected Sample Size" = expected_sample,
                     "Utility" = res$val)

  }

  class(results_list)<-"gbop2"
  return(results_list)

}






