#' PSO: Power maximization design with toxicity and efficacy boundaries
#'
#' This function implements Particle Swarm Optimization (PSO) to find an power maximizing design with toxicity and efficacy boundaries.
#'
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
#' @examples 
#' PSO_power_TE(
#'   method = "default", 
#'   nlooks = 1, 
#'   skip_efficacy = NULL, 
#'   skip_toxicity = NULL, 
#'   totalPatients = 26, 
#'   Nmin_cohort1 = 10, 
#'   Nmin_increase = 5, 
#'   e1n = 0.15, 
#'   e2n = 0.16, 
#'   e3n = 0.024, 
#'   e1a = 0.4, 
#'   e2a = 0.08, 
#'   e3a = 0.032, 
#'   err_eff = 1, 
#'   err_tox = 1, 
#'   err_all = 0.1, 
#'   power_eff = 0.8, 
#'   power_tox = 0.8, 
#'   power_all = 0.8, 
#'   seed = 1024, 
#'   nSwarm = 32, 
#'   maxIter = 100
#' )
#'
#'
#' @export
#' @import globpso R6 Rcpp RcppArmadillo
#' @importFrom stats dbinom na.omit pbeta pgamma rmultinom runif
#' @importFrom dplyr filter select distinct



PSO_power_TE <- function(method = "default",
                         nlooks = 3,
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

  if((!is.null(skip_efficacy) && (length(skip_efficacy) != (nlooks +1) ) )| (!is.null(skip_toxicity) && (length(skip_toxicity) != (nlooks +1)))){
    stop("skip_efficacy and skip_toxicity must be the length of (nlooks +1)")
  }
  
  
  
  if(!is.null(skip_efficacy) && !is.null(skip_toxicity)){
    for(i in 1: (nlooks+1)){
      if (skip_efficacy[i] == 1 && skip_toxicity[i] ==1){
        stop("Error: Cannot skip both efficacy and toxicity at the same stage")
      }
    }
    
  }
  
  boundary_tox <- rep(NA, nlooks + 1)  # Ensure it exists before modification
  boundary_eff <- rep(NA, nlooks + 1)
  
  # library(globpso)
  # library(R6)
  # library(Rcpp)
  # library(RcppArmadillo)
  # source("BOP2_functions_EffTox.R")
  # source("BOP2_TE_function.R")
  # source("boundcode.R")
  # Rcpp::sourceCpp(file="Calculation2_original.cpp")
  # Rcpp::sourceCpp(file="Calculation_minimizeN.cpp",cacheDir="cache")
  # numOfSimForTiralSetting = 10000

  ## Fixed parameters
   

  ## Fixed parameters
  input <- list(
    skip_efficacy = skip_efficacy, # if FALSE then skip tox
    skip_toxicity = skip_toxicity,
    e1n = e1n,  # H0 for Eff
    e2n = e2n,  # H0 for Tox
    e3n = e3n, # H0 for Eff and Tox
    e1a = e1a,  # Ha for Eff
    e2a = e2a,  # Ha for Tox
    e3a = e3a, # Ha for Eff and Tox
    err_eff = err_eff,  # Type I error rate: Efficacious but toxic
    err_tox = err_tox,  # Type I error rate: Safe but futile
    err_all = err_all,  # Type I error rate: Futile and toxic
    power_eff = power_eff,
    power_tox = power_tox,
    power_all = power_all,
    seed = seed
  )




  miniPatients <- Nmin_cohort1 + nlooks*Nmin_increase
  
  if(totalPatients < miniPatients){
    stop(paste0("Error: Please increase totalPatients to more than ", miniPatients  ))
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
    w2 = w[-length(w)]
    nobs[which.max(w)] = nobs[which.max(w)] - extra
    return(nobs)
  }
  

  
  r0 = input$e1n # H0 for Eff
  t0 = input$e2n # H0 for Tox
  t00 = input$e3n # H0 for Eff and Tox
  r1 = input$e1a
  t1 = input$e2a
  t11 = input$e3a
  
  ## toxicity and efficacy independent or correlated
  if(r0*t0 != t00 | r1*t1 != t11){ ## correlated
    scen = hypotheses_corr(r0,r1,t0,t1,PA_ET=input$e3a, PN_ET = input$e3n)
    
  } else{## independent
    scen = hypotheses_ind(r0, r1, t0, t1)
  }
  
  
  input$scenario = scen
  inputlist <- input
  inputlist$R <- nlooks
  input$R <- nlooks



  ## Build the utility function -----
  
  objf <- function(x, inputlist) {
    
    if (nlooks ==1){
      le = x[1]
      lt = x[2]
      g = x[3]
      #n = 26
      w1 = x[4]
      le2 = x[5]
      lt2 = x[6]

      w_list = c(w1, 1-w1)
    }else{
      le = x[1]
      lt = x[2]
      g = x[3]
      
      le2 = x[(length(x)-1)]
      lt2 = x[length(x)]
      
      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- x[4: (length(x)-2)] ## number of thetas
      w_list <- c()
      w_list[1] <- (cos(theta[1]))^2
      for( ii in 2: (n_cohort-1)){
        w_list[ii] <- (prod(sin(theta[1:(ii-1)]))*cos(theta[ii]))^2
        
      }
      w_list[n_cohort] <- (prod(sin(theta[1:(n_cohort-1)])))^2
    }
    
    #n = 26
    if (!all.equal(sum(w_list), 1, tolerance = 1e-6)) {
      stop("Error: The sum of the elements in w_list must be approximately equal to 1.")
    }
    
    nobs.seq <- cohortSize(N = totalPatients, R = nlooks, w = w_list, n_min_cohort1 = Nmin_cohort1, n_min_incre = Nmin_increase)
    nobs.seq[length(nobs.seq)] = totalPatients
    interm = interm.eff = interm.tox = nobs.seq
    
    boundary = get_boundarycpp_2lambda(interm.eff= interm.eff,interm.tox =interm.tox, 
                                       lambda_e=le, lambda_t=lt, 
                                       lambda_e2=le2, lambda_t2=lt2, gamma = g,
                                       prior=inputlist$scenario[1,], r0 = inputlist$e1n, t0 = input$e2n)
    
    interm = interm.eff = interm.tox = ceiling(cohortSize(N = totalPatients, R = nlooks, w = w_list, n_min_cohort1 = Nmin_cohort1, n_min_incre = Nmin_increase))
    interm[length(interm)] = totalPatients
    nintm = length(interm)
    new_patient = c(0, interm)
    npt=rep(NA, nintm)
    for(j in seq(nintm)){
      npt[j] = new_patient[j+1] - new_patient[j]
    }
    
    bound_eff = boundary[[1]]
    bound_tox = boundary[[2]]
    
    
    index_effi <-  which(skip_efficacy==1) ## skip which interim
    index_toxi <-  which(skip_toxicity==1) ## skip which interim
    if (!is.null(skip_efficacy)){ ## skip efficacy
      boundary_eff <- bound_eff
      boundary_eff[index_effi] <- -1 ## stop if <= cutoff
      # boundary_tox = boundary$boundary.tox
    } 
    
    if(!is.null(skip_toxicity) ){ ## skip toxicity
      boundary_tox <- bound_tox
      boundary_tox[index_toxi] <- 999 ## stop if >= cutoff
      # boundary_eff = boundary$boundary.eff
    }
    
    
    
    
    
    r0 = inputlist$e1n
    t0 = inputlist$e2n
    r1 = inputlist$e1a
    t1 = inputlist$e2a
    # temp_pe = exact_error_recursive2_Rcpp(interm.eff, bound_eff, r0, r1, 5)
    # temp_pt = exact_error_recursive2_Rcpp(interm.tox, interm-bound_tox, 1-t0, 1-t1, 5)
    
    
    if(is.null(skip_efficacy) && is.null(skip_toxicity)){
      # temp_pe = exact_error_recursive2_Rcpp(interm.eff, bound_eff, r0, r1, 5)
      # temp_pt = exact_error_recursive2_Rcpp(interm.tox, interm-bound_tox, 1-t0, 1-t1, 5)
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[2,], 
                                    bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
    }else if(!is.null(skip_efficacy) && is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[2,], 
                                    bound_eff=boundary_eff, bound_tox=boundary$boundary.tox)
    } else if(is.null(skip_efficacy) && !is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[2,], 
                                    bound_eff=boundary$boundary.eff, bound_tox=boundary_tox)
    } else{
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[2,], 
                                    bound_eff=boundary_eff, bound_tox=boundary_tox)
    }
    
    
    
    
    # a00 = temp_pe$t1err * temp_pt$t1err
    # a01 = temp_pe$t1err * temp_pt$power
    # a10 = temp_pe$power * temp_pt$t1err
    # a11 = temp_pe$power * temp_pt$power
    # if (a00 > inputlist$err_all){
    #   result = 999
    # } else {
    #   result = -a11
    # }
    
    
    N01 = temp$ptsa
    a01= temp$nonstop_prob
    
 
    
    if (a01 > inputlist$err_tox){
      # print("a01")
      result = 999
      
    } else {
     
      
      if(is.null(skip_efficacy) && is.null(skip_toxicity)){
  
        temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[3,], 
                                      bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
      }else if(!is.null(skip_efficacy) && is.null(skip_toxicity)){
        temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[3,], 
                                      bound_eff=boundary_eff, bound_tox=boundary$boundary.tox)
      } else if(is.null(skip_efficacy) && !is.null(skip_toxicity)){
        temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[3,], 
                                      bound_eff=boundary$boundary.eff, bound_tox=boundary_tox)
      } else{
        temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[3,], 
                                      bound_eff=boundary_eff, bound_tox=boundary_tox)
      }
      
      N10 = temp$ptsa
      a10= temp$nonstop_prob
      
      if (a10 > inputlist$err_eff){
        # print("a10")
        result = 999
      } else {
        
          if(is.null(skip_efficacy) && is.null(skip_toxicity)){
 
            temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[1,], 
                                          bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
          }else if(!is.null(skip_efficacy) && is.null(skip_toxicity)){
            temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[1,], 
                                          bound_eff=boundary_eff, bound_tox=boundary$boundary.tox)
          } else if(is.null(skip_efficacy) && !is.null(skip_toxicity)){
            temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[1,], 
                                          bound_eff=boundary$boundary.eff, bound_tox=boundary_tox)
          } else{
            temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[1,], 
                                          bound_eff=boundary_eff, bound_tox=boundary_tox)
          }
        
        N00 = temp$ptsa
        a00= temp$nonstop_prob
        
        if (a00 > inputlist$err_all){
          # print("a00")
          result = 999
        } else {
            if(is.null(skip_efficacy) && is.null(skip_toxicity)){
        
              temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[4,], 
                                            bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
            }else if(!is.null(skip_efficacy) && is.null(skip_toxicity)){
              temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[4,], 
                                            bound_eff=boundary_eff, bound_tox=boundary$boundary.tox)
            } else if(is.null(skip_efficacy) && !is.null(skip_toxicity)){
              temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[4,], 
                                            bound_eff=boundary$boundary.eff, bound_tox=boundary_tox)
            } else{
              temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[4,], 
                                            bound_eff=boundary_eff, bound_tox=boundary_tox)
            }
          N11 = temp$ptsa
          a11= temp$nonstop_prob
          
          result = -a11
          
        }
      } 
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
  set.seed(input$seed)
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
    
    pars = res$par
    
    if (nlooks ==1){
      le = pars[1]
      lt = pars[2]
      g = pars[3]
      #n = 26
      w1 = pars[4]
      le2 = pars[5]
      lt2 = pars[6]
      
      w_list = c(w1, 1-w1)
    }else{
      le = pars[1]
      lt = pars[2]
      g = pars[3]
      
      le2 = pars[(length(pars)-1)]
      lt2 = pars[length(pars)]
      
      n_cohort <- nlooks +1 ## n_cohort is number of cohort
      theta <- pars[4: (length(pars)-2)] ## number of thetas
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
    
    nobs.seq <- cohortSize(N = totalPatients, R = nlooks, w = w_list, n_min_cohort1 = Nmin_cohort1, n_min_incre = Nmin_increase)
    nobs.seq[length(nobs.seq)] = totalPatients
    interm = interm.eff = interm.tox = nobs.seq
    
    boundary = get_boundarycpp_2lambda(interm.eff= interm.eff,interm.tox =interm.tox, 
                                       lambda_e=le, lambda_t=lt, 
                                       lambda_e2=le2, lambda_t2=lt2, gamma = g,
                                       prior=inputlist$scenario[1,], r0 = inputlist$e1n, t0 = input$e2n)
    
    interm = interm.eff = interm.tox = ceiling(cohortSize(N = totalPatients, R = nlooks, w = w_list, n_min_cohort1 = Nmin_cohort1, n_min_incre = Nmin_increase))
    interm[length(interm)] = totalPatients
    
    nintm = length(interm)
    new_patient = c(0, interm)
    npt=rep(NA, nintm)
    for(j in seq(nintm)){
      npt[j] = new_patient[j+1] - new_patient[j]
    }
    
    bound_eff = boundary[[1]]
    bound_tox = boundary[[2]]
    index_effi <-  which(skip_efficacy==1) ## skip which interim
    index_toxi <-  which(skip_toxicity==1) ## skip which interim
    if (!is.null(skip_efficacy)){ ## skip efficacy
      boundary_eff <- bound_eff
      boundary_eff[index_effi] <- -1 ## stop if <= cutoff
   
    } 
    
    if(!is.null(skip_toxicity) ){ ## skip toxicity
      boundary_tox <- bound_tox
      boundary_tox[index_toxi] <- 999 ## stop if >= cutoff
      
    }
    
    # temp_pe = exact_error_recursive2_Rcpp(interm.eff, bound_eff, r0, r1, 5)
    # temp_pt = exact_error_recursive2_Rcpp(interm.tox, interm-bound_tox, 1-t0, 1-t1, 5)
    # boundary$boundary.eff <- c(4,10)
    # boundary$boundary.tox <- c(6,8)
   
    if(is.null(skip_efficacy) && is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[1,], 
                                    bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
     
    }else if(!is.null(skip_efficacy) && is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[1,], 
                                    bound_eff=boundary_eff, bound_tox=boundary$boundary.tox)
    } else if(is.null(skip_efficacy) && !is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[1,], 
                                    bound_eff=boundary$boundary.eff, bound_tox=boundary_tox)
    } else{
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[1,], 
                                    bound_eff=boundary_eff, bound_tox=boundary_tox)
    }
    
    N00 = temp$ptsa
    a00= temp$nonstop_prob
    
    if(is.null(skip_efficacy) && is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[4,], 
                                    bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
      
    }else if(!is.null(skip_efficacy) && is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[4,], 
                                    bound_eff=boundary_eff, bound_tox=boundary$boundary.tox)
    } else if(is.null(skip_efficacy) && !is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[4,], 
                                    bound_eff=boundary$boundary.eff, bound_tox=boundary_tox)
    } else{
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[4,], 
                                    bound_eff=boundary_eff, bound_tox=boundary_tox)
    }
    
    N11 = temp$ptsa
    a11= temp$nonstop_prob
    
    
    
    if(is.null(skip_efficacy) && is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[2,], 
                                    bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
      
    }else if(!is.null(skip_efficacy) && is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[2,], 
                                    bound_eff=boundary_eff, bound_tox=boundary$boundary.tox)
    } else if(is.null(skip_efficacy) && !is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[2,], 
                                    bound_eff=boundary$boundary.eff, bound_tox=boundary_tox)
    } else{
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[2,], 
                                    bound_eff=boundary_eff, bound_tox=boundary_tox)
    }
    
    N01 = temp$ptsa
    a01= temp$nonstop_prob
    
    
    
    if(is.null(skip_efficacy) && is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[3,], 
                                    bound_eff=boundary$boundary.eff, bound_tox=boundary$boundary.tox)
      
    }else if(!is.null(skip_efficacy) && is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[3,], 
                                    bound_eff=boundary_eff, bound_tox=boundary$boundary.tox)
    } else if(is.null(skip_efficacy) && !is.null(skip_toxicity)){
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[3,], 
                                    bound_eff=boundary$boundary.eff, bound_tox=boundary_tox)
    } else{
      temp = efftox_recursive_optim(interm=interm, npt=npt, p=inputlist$scenario[3,], 
                                    bound_eff=boundary_eff, bound_tox=boundary_tox)
    }
    
               
    N10 = temp$ptsa
    a10= temp$nonstop_prob
    
    # a00 = temp_pe$t1err * temp_pt$t1err
    # a01 = temp_pe$t1err * temp_pt$power
    # a10 = temp_pe$power * temp_pt$t1err
    # a11 = temp_pe$power * temp_pt$power
    # 
    # Nt = temp_pt$pts
    # Ne = temp_pe$pts
    # Nte = temp_pe$pts+temp_pt$pts
    
    # default_tbl <- c((res$cputime), res$par,
    #                                    interm, boundary[[1]],boundary[[2]], a00, a01, a10, a11,
    #                                    N00,N01,N10,N11, -res$val)
    
    
    
  }
  
  
  if(!is.null(skip_efficacy)){
    boundary$boundary.eff <- boundary_eff
  } 
  
  
  if(!is.null(skip_toxicity)){
    boundary$boundary.tox <- boundary_tox
  } 
  
  results_list <- list(
    "function" = "PSO_power_TE",
    "method" = method,
    "parameter" = list(
      "lambdae1" = le, "lambdae2" = le2,
      "lambdat1" = lt,
      "lambdat2" = lt2,
      "gamma" = g),
    "cohort" = as.list(interm),                  # Cohort sizes
    "boundary_effi" = as.list(boundary$boundary.eff),   # Boundary for efficacy
    "boundary_toxi" = as.list(boundary$boundary.tox),   # Boundary for toxicity
    # "expected_sample" =  expected_sample,
    "typeI_H01 (safe but futile)" = a01,
    "typeI_H10 (efficacious but toxic)" = a10,
    "typeI_H00 (futile and toxic)" = a00,
    "power" = a11,
    "utility" = -res$val
  )
  
  
  
  results_list
  class(results_list)<-"gbop2"
  return(results_list)
}





