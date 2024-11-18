#' PSO for optimal/minimax design with efficacy and toxicity boundaries
#'
#' @param design choose from "optimal" or "minimax"
#' @param method method for single PSO, choose from "default", "quantum" or "dexp"
#' @param nlooks number of interim looks
#' @param skip_efficacy default is NULL, indicate skip efficacy as 1 and not skip as 0 in a vector
#' @param skip_toxicity default is NULL, indicate skip toxicity as 1 and not skip as 0 in a vector
#' @param maxPatients  maximum number of patients
#' @param Nmin_cohort1 minimum number of first cohort
#' @param Nmin_increase minimum number of increase in each cohort
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
#'
#' @return A list on design parameters and operating characteristics
#' @export
#' @import globpso R6 Rcpp RcppArmadillo
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct



PSODesign_TE <- function(design = "optimal",
                         method = "quantum",
                         nlooks = 2,
                         skip_efficacy = NULL, # NULL no skipping c(), 1 is skip
                         skip_toxicity = NULL , ## NULL c()
                         maxPatients = 50, ## maximum number of patients
                         Nmin_cohort1 = 10,
                         Nmin_increase = 5,
                         e1n = 0.3,  # H0 for Eff
                         e2n = 0.4,  # H0 for Tox
                         e3n = 0.2, # H0 for Eff and Tox
                         e1a = 0.6,  # Ha for Eff
                         e2a = 0.2,  # Ha for Tox
                         e3a = 0.15, # Ha for Eff and Tox
                         err_eff = 0.1,  # Type I error rate: Efficacious but toxic
                         err_tox = 0.1,  # Type I error rate: Safe but futile
                         err_all = 0.05,  # Type I error rate: Futile and toxic
                         power_eff = 0.8,
                         power_tox = 0.8,
                         power_all = 0.8,
                         nSwarm = 32,
                         maxIter = 100){


  if(!is.null(skip_efficacy) && !is.null(skip_toxicity)){
    for(i in 1: nlooks){
      if (skip_efficacy[i] == 1 & skip_toxicity[i] ==1){
        stop("Error: Cannot skip both efficacy and toxicity at the same interim look.")
      }
    }

  }

  # library(globpso)
  # library(R6)
  # library(Rcpp)
  # library(RcppArmadillo)
  # source("BOP2_functions_EffTox.R")
  # source("BOP2_TE_function.R")
  # source("boundcode.R")
  # Rcpp::sourceCpp(file="Calculation2_original.cpp")
  # Rcpp::sourceCpp(file="Calculation_minimizeN.cpp",cacheDir="cache")
  numOfSimForTiralSetting = 10000   # Number of simulations



  ## Fixed parameters
  input <- list(
    skip_efficacy = skip_efficacy, # if FALSE then skip tox
    skip_toxicity = skip_toxicity,
    e1n = 0.3,  # H0 for Eff
    e2n = 0.4,  # H0 for Tox
    e3n = 0.2, # H0 for Eff and Tox
    e1a = 0.6,  # Ha for Eff
    e2a = 0.2,  # Ha for Tox
    e3a = 0.15, # Ha for Eff and Tox
    err_eff = 0.1,  # Type I error rate: Efficacious but toxic
    err_tox = 0.1,  # Type I error rate: Safe but futile
    err_all = 0.05,  # Type I error rate: Futile and toxic
    power_eff = 0.8,
    power_tox = 0.8,
    power_all = 0.8
  )

  miniPatients <- Nmin_cohort1 + nlooks*Nmin_increase

  if(maxPatients < miniPatients){
    stop(paste0("Error: Please increase maxPatients to more than ", miniPatients  ))
  }

  cohortSize = function(N, R, w, n_min_cohort1 = Nmin_cohort1, n_min_incre = Nmin_increase){
    Nrest = N - R*n_min_incre - n_min_cohort1
    nobs = c()
    extra = 0
    for ( i in 1:(R+1)){
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

  r0 = input$e1n
  t0 = input$e2n
  r1 = input$e1a
  t1 = input$e2a
  t00 = input$e3n
  t11 = input$e3a

  # interm.eff = input$nobs.seq
  # interm.tox = input$nobs.seqTX


  ## toxicity and efficacy independent or correlated
  if(r0*t0 != t00 | r1*t1 != t11){ ## correlated
    scen = hypotheses_corr(r0,r1,t0,t1,PA_ET=input$e3a, PN_ET = input$e3n)

  } else{## independent
    scen = hypotheses_ind(r0, r1, t0, t1)
  }




  input$scenario = scen

  inputlist = input


  ## Build the utility function -----

  objf <- function(x, inputlist) {
    if (nlooks ==1){
      le = x[1]
      lt = x[2]
      g = x[3]
      n = x[4]
      le2 = x[5]
      w1 = x[6]
      lt2 = x[7]
      w_list = c(w1, 1-w1)

    }else{
      le = x[1]
      lt = x[2]
      g = x[3]
      n = x[4]
      le2 = x[5]

      lt2 = x[length(x)]
      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- x[6: (length(x)-1)] ## number of thetas
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


    ## interm is previous nobs.seq
    interm = interm.eff = interm.tox = (cohortSize(N = n, R = nlooks, w = w_list))

    boundary = get_boundarycpp_2lambda(interm.eff= interm.eff,interm.tox =interm.tox,
                                       lambda_e=le, lambda_t=lt,
                                       lambda_e2=le2, lambda_t2=lt2, gamma = g,
                                       prior=inputlist$scenario[1,], r0 = inputlist$e1n, t0 = input$e2n)

    interm = ceiling(cohortSize(N = n, R = nlooks, w = w_list))
    index_effi <-  which(skip_efficacy==1) ## skip which interim
    index_toxi <-  which(skip_toxicity==1) ## skip which interim
    if (is.null(skip_efficacy)){ ## skip efficacy
      boundary_eff <- boundary$boundary.eff
      boundary_eff[index_effi] <- -1
      # boundary_tox = boundary$boundary.tox
    } else if(is.null(skip_toxicity) ){ ## skip toxicity
      boundary_tox = boundary_tox
      boundary_tox[index_toxi] <- boundary_tox[index_toxi] + 1
      # boundary_eff = boundary$boundary.eff
    }


    ################################
    if(design =="minimax"){
      n_final = interm[length(interm)]
    }
    nintm = length(interm)
    new_patient = c(0, interm)
    npt=rep(NA, nintm)
    for(j in seq(nintm)){
      npt[j] = new_patient[j+1] - new_patient[j]
    }


    temp = .Call("_GBOP2_efftox_recursive_optim", interm=interm, npt=npt, p=inputlist$scenario[2,],
                                  bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
    N01 = temp$ptsa
    a01= temp$nonstop_prob

    if (a01 > inputlist$err_tox){
      # print("a01")
      result = 999

    } else {
      temp = .Call("_GBOP2_efftox_recursive_optim",interm=interm, npt=npt, p=inputlist$scenario[3,],
                                    bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
      N10 = temp$ptsa
      a10= temp$nonstop_prob

      if (a10 > inputlist$err_eff){
        # print("a10")
        result = 999
      } else {
        temp = .Call("_GBOP2_efftox_recursive_optim", interm=interm, npt=npt, p=inputlist$scenario[1,],
                                      bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
        N00 = temp$ptsa
        a00= temp$nonstop_prob

        if (a00 > inputlist$err_all){
          result = 999
        } else {
          temp = .Call("_GBOP2_efftox_recursive_optim", interm=interm, npt=npt, p=inputlist$scenario[4,],
                                        bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
          N11 = temp$ptsa
          a11= temp$nonstop_prob

          if (a11 < inputlist$power_all){
            result = 999
          } else {


            if(design =="optimal"){
              result = ((N00 + N01 + N10)/3 + N11)/2
            } else{
              result = n_final + (((N00 + N01 + N10)/3 + N11)/2)/n_final
            }
          }
        }
      }
    }


    return(result)

  }


  if(nlooks ==1){
    low_bound <- c(0.5, 0.5, 0, miniPatients,  0.5, 0, 0.5)
    upp_bound <- c(0.99, 0.99, 1, maxPatients, 0.99, 1, 0.99)
  } else{
    theta_L <- rep(0, nlooks) ## lower bound of theta
    theta_U <- rep(pi/2, nlooks) ## upper bound of theta
    low_bound <- c(0.5, 0.5, 0, miniPatients, 0.5, theta_L, 0.5)
    upp_bound <- c(0.99, 0.99, 1, maxPatients, 0.99, theta_U, 0.99)
  }

  n_sim = 1
  set.seed(123);
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


  for ( i in n_sim){

    res <- globpso(objFunc = objf, lower = low_bound, upper = upp_bound,
                   fixed = NULL, PSO_INFO = alg_setting,
                   inputlist = inputlist, seed = seeds[i])
    pars = res$par

    if (nlooks ==1){
      le = pars[1]
      lt = pars[2]
      g = pars[3]
      n = pars[4]
      le2 = pars[5]
      w1 = pars[6]
      lt2 = pars[7]
      w_list = c(w1, 1-w1)
    } else{ ## when nlooks >=2
      le = pars[1]
      lt = pars[2]
      g = pars[3]
      n = pars[4]
      le2 = pars[5]

      lt2 = pars[length(pars)]
      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- pars[6: (length(pars)-1)] ## number of thetas
      w_list <- c()
      w_list[1] <- (cos(theta[1]))^2
      for( ii in 2: n_cohort-1){
        w_list[ii] <- (prod(sin(theta[1:(ii-1)]))*cos(theta[ii]))^2

      }
      w_list[n_cohort] <- (prod(sin(theta[1:(n_cohort-1)])))^2
    }




    interm = interm.eff = interm.tox = (cohortSize(N = n, R = nlooks, w = w_list))

    boundary = get_boundarycpp_2lambda(interm.eff= interm.eff,interm.tox =interm.tox,
                                       lambda_e=le, lambda_t=lt,
                                       lambda_e2=le2, lambda_t2=lt2, gamma = g,
                                       prior=inputlist$scenario[1,], r0 = inputlist$e1n, t0 = input$e2n)

    interm = ceiling(cohortSize(N = n, R = nlooks, w = w_list))


    index_effi <-  which(skip_efficacy==1) ## skip which interim
    index_toxi <-  which(skip_toxicity==1) ## skip which interim
    if (is.null(skip_efficacy)){ ## skip efficacy
      boundary_eff <- boundary$boundary.eff
      boundary_eff[index_effi] <- -1
      # boundary_tox = boundary$boundary.tox
    } else if(is.null(skip_toxicity) ){ ## skip toxicity
      boundary_tox = boundary_tox
      boundary_tox[index_toxi] <- boundary_tox[index_toxi] + 1
      # boundary_eff = boundary$boundary.eff
    }


    if(design =="minimax"){
      n_final = interm[length(interm)]
    }

    ########################################################
    nintm = length(interm)
    new_patient = c(0, interm)
    npt=rep(NA, nintm)
    for(j in seq(nintm)){
      npt[j] = new_patient[j+1] - new_patient[j]
    }

    temp = .Call("_GBOP2_efftox_recursive_optim",interm=interm, npt=npt, p=inputlist$scenario[1,],
                                  bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
    N00 = temp$ptsa ## expected sample size under H00
    a00= temp$nonstop_prob

    temp = .Call("_GBOP2_efftox_recursive_optim",interm=interm, npt=npt, p=inputlist$scenario[4,],
                                  bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
    N11 = temp$ptsa
    a11= temp$nonstop_prob

    temp = .Call("_GBOP2_efftox_recursive_optim",interm=interm, npt=npt, p=inputlist$scenario[2,],
                                  bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
    N01 = temp$ptsa ## expected sample size under H01
    a01= temp$nonstop_prob ## probability of not stopping under H01

    temp = .Call("_GBOP2_efftox_recursive_optim", interm=interm, npt=npt, p=inputlist$scenario[3,],
                                  bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
    N10 = temp$ptsa ## expected sample size under H10
    a10= temp$nonstop_prob

    expected_sample <- ((N00 + N01 + N10)/3 + N11)/2
    if(design =="optimal"){
      utility = ((N00 + N01 + N10)/3 + N11)/2
    } else{
      utility = n_final + (((N00 + N01 + N10)/3 + N11)/2)/n_final
    }

  }


  # Concatenate values into a single vector, flattening sub-elements as needed
  results_list <- list(
    "function" = "PSODesign_TE",
    "design" = design,"method" = method,
    "parameter" = list(
      "lambdae1" = le, "lambdae2" = le2,
      "lambdat1" = lt,
      "lambdat2" = lt2,
      "gamma" = g),
    "cohort" = as.list(interm),                  # Cohort sizes
    "boundary_effi" = as.list(boundary$boundary.eff),   # Boundary for efficacy
    "boundary_toxi" = as.list(boundary$boundary.tox),   # Boundary for toxicity
    "expected_sample" =  expected_sample,
    "typeI_01" = a01,
    "typeI_10" = a10,
    "typeI_00" = a00,
    "power" = a11,
    "utility" = utility
  )

  # results_list
  class(results_list)<-"gbop2"
  return(results_list)
}




