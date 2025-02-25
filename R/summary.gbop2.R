#' Summary function
#' Summary function for gbop2 objects
#'@param object GBOP2_maxP_dualE GBOP2_maxP_singleE GBOP2_maxP_TE GBOP2_minSS_dualE GBOP2_minSS_singleE GBOP2_minSS_TE PSO_design PSO_design_dual PSO_design_TE PSO_power PSO_power_dual PSO_power_TE
#
#' @param ... ignored arguments
#'
#' @return A summary table
#' @examples
#' summary(GBOP2_maxP_singleE(
#'   nlooks = 1, b1n = 0.2, b1a = 0.4, err1 = 0.05, 
#'   nParallel = 3, minPower = 0.8, totalPatients = 50, 
#'   Nmin_cohort1 = 10, Nmin_increase = 5
#' ))
#'
#' @export

summary.gbop2 <- function(object, ...) {

  ########## PSO_design ########################################
  if (object$`function` == "PSO_design") {
    # Extract relevant fields from object
    design <- object$design
    weight <- object$weight
    method <- object$method
    cputime <- object$cputime
    parameter <- object$parameter
    cohort <- object$cohort
    boundary <- object$boundary
    typeI_error <- object$`Type I Error`
    power <- object$Power
    expected_sample <- object$`Expected Sample Size`
    utility <- object$Utility

    cat("------------------------------------------------\n")
    cat("PSO Optimal/Minimax design with single boundary\n")

    cat("design:           ", design, "\n")
    cat("weight:           ", weight, "\n")
    cat("method:           ", method, "\n")
    cat("CPU Time:         ", format(cputime, digits = 4), " seconds\n\n")

    cat("parameters:\n")
    cat("  lambda1:        ", format(parameter$lambda1, digits = 6), "\n")
    cat("  lambda2:        ", format(parameter$lambda2, digits = 6), "\n")
    cat("  gamma:          ", format(parameter$gamma, digits = 6), "\n\n")

    cat("cohort sizes:     ", paste(unlist(cohort), collapse = ", "), "\n")
    cat("boundary:  ", paste("<=", unlist(boundary), collapse = ", "), "\n\n")

    cat("type I error:     ", format(typeI_error, digits = 6), "\n")
    cat("power:            ", format(power, digits = 6), "\n")
    cat("expected sample size: ", format(expected_sample, digits = 6), "\n")
    cat("utility Value:    ", format(utility, digits = 6), "\n\n")

    cat("------------------------------------------------\n")
  }


  ########## PSO_power ########################################
  if (object$`function` == "PSO_power") {

    design <- object$design
    method <- object$method
    cputime <- object$cputime
    lambda1 <- object$lambda1
    lambda2 <- object$lambda2
    gamma <- object$gamma
    cohort <- object$cohort
    boundary <- object$boundary
    typeI_error <- object$TypeI
    power <- object$Power
    expected_sample <- object$`EN(P0)`
    utility <- object$Utility

    cat("------------------------------------------------\n")
    cat("PSO maximazing power with single boundary\n")

    cat("design:           ", design, "\n")
    cat("method:           ", method, "\n")
    cat("CPU Time:         ", format(cputime, digits = 4), " seconds\n\n")

    cat("parameters:\n")
    cat("  lambda1:        ", format(lambda1, digits = 6), "\n")
    cat("  lambda2:        ", format(lambda2, digits = 6), "\n")
    cat("  gamma:          ", format(gamma, digits = 6), "\n\n")

    cat("cohort sizes:     ", paste(unlist(cohort), collapse = ", "), "\n")
    cat("boundary:  ", paste("<=", unlist(boundary), collapse = ", "), "\n\n")
    
    cat("type I error:     ", format(typeI_error, digits = 6), "\n")
    cat("power:            ", format(power, digits = 6), "\n")
    cat("expected sample size (EN(P0)): ", format(expected_sample, digits = 6), "\n")
    cat("utility Value:    ", format(utility, digits = 6), "\n\n")

    cat("------------------------------------------------\n")
  }



  ########## GBOP2_minSS_single ########################################
  if (object$`function` == "GBOP2_minSS_singleE") {
    # Extract relevant fields from object
    design <- object$design
    weight <- object$weight
    method <- object$method
    cputime <- object$cputime
    lambda1 <- object$`parameter.lambda1`
    lambda2 <- object$`parameter.lambda2`
    gamma <- object$`parameter.gamma`
    typeI_error <- object$`Type.I.Error`
    power <- object$Power
    expected_sample <- object$`Expected.Sample.Size`
    utility <- object$Utility

    cohort <- grep("^cohort", names(object), value = TRUE)
    boundary <- grep("^boundary", names(object), value = TRUE)

    cat("------------------------------------------------\n")
    cat("PSOGO for optimal/minimax design with single boundary\n")

    cat("design:           ", design, "\n")
    cat("weight:           ", weight, "\n")
    cat("method:           ", method, "\n")
    cat("cpu time:         ", format(cputime, digits = 4), " seconds\n\n")

    cat("parameters:\n")
    cat("  lambda1:        ", format(lambda1, digits = 6), "\n")
    cat("  lambda2:        ", format(lambda2, digits = 6), "\n")
    cat("  gamma:          ", format(gamma, digits = 6), "\n\n")

    cat("cohort sizes:     ", paste(sapply(cohort, function(x) object[[x]]), collapse = ", "), "\n")
    cat("boundary values:  ", paste("stop if <=", sapply(boundary, function(x) object[[x]]), collapse = ", "), "\n\n")

    cat("type I error:     ", format(typeI_error, digits = 6), "\n")
    cat("power:            ", format(power, digits = 6), "\n")
    cat("expected sample size: ", format(expected_sample, digits = 6), "\n")
    cat("utility value:    ", format(utility, digits = 6), "\n\n")

    cat("------------------------------------------------\n")
  }


  ########## GBOP2_maxP_single ########################################
  if (object$`function` == "GBOP2_maxP_singleE") {
    # Extract relevant fields from object
    design <- object$design
    method <- object$method
    cputime <- object$cputime
    lambda1 <- object$lambda1
    lambda2 <- object$lambda2
    gamma <- object$gamma
    typeI_error <- object$TypeI
    power <- object$Power
    expected_sample <- object$EN.P0
    utility <- object$Utility

    # Extract cohort and boundary sizes
    cohort <- grep("^cohort", names(object), value = TRUE)
    boundary <- grep("^boundary", names(object), value = TRUE)

    cat("------------------------------------------------\n")
    cat("PSOGO for maximizing power with single boundary\n")
    cat("design:           ", design, "\n")
    cat("method:           ", method, "\n")
    cat("cpu time:         ", format(cputime, digits = 4), " seconds\n\n")

    cat("parameters:\n")
    cat("  lambda1:        ", format(lambda1, digits = 6), "\n")
    cat("  lambda2:        ", format(lambda2, digits = 6), "\n")
    cat("  gamma:          ", format(gamma, digits = 6), "\n\n")

    cat("cohort sizes:     ", paste(sapply(cohort, function(x) object[[x]]), collapse = ", "), "\n")
    # cat("boundary values:  ", paste(sapply(boundary, function(x) object[[x]]), collapse = ", "), "\n\n")
    cat("boundary values:  ", paste("stop if <=", sapply(boundary, function(x) object[[x]]), collapse = ", "), "\n\n")
    
    cat("type I error:     ", format(typeI_error, digits = 6), "\n")
    cat("power:            ", format(power, digits = 6), "\n")
    cat("expected sample size (en(p0)): ", format(expected_sample, digits = 6), "\n")
    cat("utility value:    ", format(utility, digits = 6), "\n\n")


    cat("------------------------------------------------\n")
  }


  ########## GBOP2_minSS_dual ########################################
  
  if (object$`function` == "GBOP2_minSS_dualE") {
    # Extract relevant fields from the object
    design <- object$design
    weight <- object$weight
    method <- object$method
    lambda1 <- object$parameters.lambda1
    lambda_grad1 <- object$parameters.lambda_grad1
    lambda_grad2 <- object$parameters.lambda_grad2
    gamma_1 <- object$parameters.Gamma_1
    gamma_2 <- object$parameters.Gamma_2
    gamma_3 <- object$parameters.Gamma_3
    delta0 <- object$parameters.delta0
    delta1 <- object$parameters.delta1
    
    # Extract cohort names and values
    cohort_name <- names(object)[grepl("cohort", names(object)) & !grepl("bd", names(object))]
    cohort_values <- unlist(object[cohort_name])  # Extract cohort values as a vector
    
    # Extract boundary values
    boundaryF <- unlist(object[grep("^boundaryF", names(object))])
    boundaryE <- unlist(object[grep("^boundaryE", names(object))])
    
    # Extract results
    typeI_error <- object$Type.I.Error
    power <- object$Power
    expected_sample <- object$Expected.Sample.Size
    utility <- object$Utility
    
    # Print summary
    cat("------------------------------------------------------\n")
    cat("PSOGO for minimizing sample size with dual-boundary\n")
    cat("Design:            ", design, "\n")
    cat("Weight:            ", weight, "\n")
    cat("Method:            ", method, "\n\n")
    
    cat("Cohort values:     ", paste(cohort_values, collapse = ", "), "\n")
    
    cat("Boundary values:\n")
    cat("  Boundary futile:      ", paste("stop if <=", boundaryF, collapse = ", "), "\n")
    cat("  Boundary effi:        ", paste("stop if >=", boundaryE, collapse = ", "), "\n\n")
    
    
    cat("Type I Error:      ", format(typeI_error, digits = 6), "\n")
    cat("Power:             ", format(power, digits = 6), "\n")
    cat("Expected Sample Size:", format(expected_sample, digits = 6), "\n")
    cat("Utility Value:     ", format(utility, digits = 6), "\n\n")
    
    cat("Parameters:\n")
    cat("  Lambda1:         ", format(lambda1, digits = 6), "\n")
    cat("  Lambda_grad1:    ", format(lambda_grad1, digits = 6), "\n")
    cat("  Lambda_grad2:    ", format(lambda_grad2, digits = 6), "\n")
    cat("  Gamma_1:         ", format(gamma_1, digits = 6), "\n")
    cat("  Gamma_2:         ", format(gamma_2, digits = 6), "\n")
    cat("  Gamma_3:         ", format(gamma_3, digits = 6), "\n")
    cat("  Delta0:          ", format(delta0, digits = 6), "\n")
    cat("  Delta1:          ", format(delta1, digits = 6), "\n\n")
    
    cat("------------------------------------------------\n")
  }


  ########## PSO_design_dual ########################################
  if (object$`function` == "PSO_design_dual") {
    # Extract relevant fields from the object
    design <- object$design
    weight <- object$weight
    method <- object$method
    parameters <- object$parameters
    lambda1 <- parameters$lambda1
    lambda_grad1 <-  parameters$lambda_grad1
    lambda_grad2 <-  parameters$lambda_grad2
    gamma_1 <-  parameters$Gamma_1
    gamma_2 <-  parameters$Gamma_2
    gamma_3 <-  parameters$Gamma_3
    delta0 <-  parameters$delta0
    delta1 <-  parameters$delta1
    
    # Extract cohort names and values
    cohort_name <- names(object)[grepl("cohort", names(object)) & !grepl("bd", names(object))]
    cohort_values <- unlist(object[cohort_name])  # Extract cohort values as a vector
    
    # Extract boundary values
    boundary1 <- object$boundary$`1`
    boundary2 <- object$boundary$`2`
    
    # Extract results
    typeI_error <- object$`Type I Error`
    power <- object$Power
    expected_sample <- object$`Expected Sample Size`
    utility <- object$Utility
    
    # Print summary
    cat("------------------------------------------------------\n")
    cat("PSO for minimizing sample size with dual-boundary\n")
    cat("Design:            ", design, "\n")
    cat("Weight:            ", weight, "\n")
    cat("Method:            ", method, "\n\n")
    
    cat("Cohort values:     ", paste(cohort_values, collapse = ", "), "\n")
    
    cat("Boundary values:\n")
    cat("  Boundary 1:      ", paste(boundary1, collapse = ", "), "\n")
    cat("  Boundary 2:      ", paste(boundary2, collapse = ", "), "\n\n")
    
    cat("Type I Error:      ", format(typeI_error, digits = 6), "\n")
    cat("Power:             ", format(power, digits = 6), "\n")
    cat("Expected Sample Size:", format(expected_sample, digits = 6), "\n")
    cat("Utility Value:     ", format(utility, digits = 6), "\n\n")
    
    cat("Parameters:\n")
    cat("  Lambda1:         ", format(lambda1, digits = 6), "\n")
    cat("  Lambda_grad1:    ", format(lambda_grad1, digits = 6), "\n")
    cat("  Lambda_grad2:    ", format(lambda_grad2, digits = 6), "\n")
    cat("  Gamma_1:         ", format(gamma_1, digits = 6), "\n")
    cat("  Gamma_2:         ", format(gamma_2, digits = 6), "\n")
    cat("  Gamma_3:         ", format(gamma_3, digits = 6), "\n")
    cat("  Delta0:          ", format(delta0, digits = 6), "\n")
    cat("  Delta1:          ", format(delta1, digits = 6), "\n\n")
    
    cat("------------------------------------------------\n")
  }
  
  
 ##################GBOP2_maxP_dual######################
  if (object$`function` == "GBOP2_maxP_dualE") {
    # Extract relevant fields from the object
    design <- object$design
    method <- object$method
    lambda1 <- object$parameters.lambda1
    lambda_grad1 <- object$parameters.lambda_grad1
    lambda_grad2 <- object$parameters.lambda_grad2
    gamma_1 <- object$parameters.Gamma_1
    gamma_2 <- object$parameters.Gamma_2
    gamma_3 <- object$parameters.Gamma_3
    delta0 <- object$parameters.delta0
    delta1 <- object$parameters.delta1
    
    # Extract cohort names and values
    cohort_name <- names(object)[grepl("cohort", names(object)) & !grepl("bd", names(object))]
    cohort_values <- unlist(object[cohort_name])  # Extract cohort values as a vector
    
    # Extract boundary values
    boundaryF <- unlist(object[grep("^boundaryF", names(object))])
    boundaryE <- unlist(object[grep("^boundaryE", names(object))])
    
    # Extract results
    typeI_error <- object$Type.I.Error
    power <- object$Power
    expected_sample <- object$Expected.Sample.Size
    utility <- object$Utility
    
    # Print summary
    cat("------------------------------------------------------\n")
    cat("PSOGO for power maximization with dual-boundary\n")
    cat("Design:            ", design, "\n")
    cat("Method:            ", method, "\n\n")
    
    cat("Cohort values:     ", paste(cohort_values, collapse = ", "), "\n")
    
    cat("Boundary values:\n")
    cat("  Boundary futile:    ", paste("stop if <=", boundaryF, collapse = ", "), "\n")
    cat("  Boundary effi:      ", paste("stop if >=", boundaryE, collapse = ", "), "\n\n")
    
    
    cat("Type I Error:      ", format(typeI_error, digits = 6), "\n")
    cat("Power:             ", format(power, digits = 6), "\n")
    cat("Expected Sample Size:", format(expected_sample, digits = 6), "\n")
    cat("Utility Value:     ", format(utility, digits = 6), "\n\n")
    
    cat("Parameters:\n")
    cat("  Lambda1:         ", format(lambda1, digits = 6), "\n")
    cat("  Lambda_grad1:    ", format(lambda_grad1, digits = 6), "\n")
    cat("  Lambda_grad2:    ", format(lambda_grad2, digits = 6), "\n")
    cat("  Gamma_1:         ", format(gamma_1, digits = 6), "\n")
    cat("  Gamma_2:         ", format(gamma_2, digits = 6), "\n")
    cat("  Gamma_3:         ", format(gamma_3, digits = 6), "\n")
    cat("  Delta0:          ", format(delta0, digits = 6), "\n")
    cat("  Delta1:          ", format(delta1, digits = 6), "\n\n")
    
    cat("------------------------------------------------\n")
  }
  
  
  
  ########## PSO_power_dual ########################################
  if (object$`function` == "PSO_power_dual") {
    # Extract relevant fields from the object
    design <- object$design
    weight <- object$weight
    method <- object$method
    parameters <- object$parameter
    lambda1 <- parameters$lambda1
    lambda_grad1 <-  parameters$lambda_grad1
    lambda_grad2 <-  parameters$lambda_grad2
    gamma_1 <-  parameters$gamma_1
    gamma_2 <-  parameters$gamma_2
    gamma_3 <-  parameters$gamma_3
    delta0 <-  parameters$delta0
    delta1 <-  parameters$delta1
    
    # Extract cohort names and values
    cohort_name <- names(object)[grepl("cohort", names(object)) & !grepl("bd", names(object))]
    cohort_values <- unlist(object[cohort_name])  # Extract cohort values as a vector
    
    # Extract boundary values
    boundary1 <- object$boundary$`1`
    boundary2 <- object$boundary$`2`
    
    # Extract results
    typeI_error <- object$`Type I Error`
    power <- object$Power
    expected_sample <- object$`Expected Sample Size`
    utility <- object$Utility
    
    # Print summary
    cat("------------------------------------------------------\n")
    cat("PSO for minimizing sample size with dual-boundary\n")
    cat("Design:            ", design, "\n")
    cat("Weight:            ", weight, "\n")
    cat("Method:            ", method, "\n\n")
    
    cat("Cohort values:     ", paste(cohort_values, collapse = ", "), "\n")
    
    cat("Boundary values:\n")
    cat("  Boundary 1:      ", paste(boundary1, collapse = ", "), "\n")
    cat("  Boundary 2:      ", paste(boundary2, collapse = ", "), "\n\n")
    
    cat("Type I Error:      ", format(typeI_error, digits = 6), "\n")
    cat("Power:             ", format(power, digits = 6), "\n")
    cat("Expected Sample Size:", format(expected_sample, digits = 6), "\n")
    cat("Utility Value:     ", format(utility, digits = 6), "\n\n")
    
    cat("Parameters:\n")
    cat("  Lambda1:         ", format(lambda1, digits = 6), "\n")
    cat("  Lambda_grad1:    ", format(lambda_grad1, digits = 6), "\n")
    cat("  Lambda_grad2:    ", format(lambda_grad2, digits = 6), "\n")
    cat("  Gamma_1:         ", format(gamma_1, digits = 6), "\n")
    cat("  Gamma_2:         ", format(gamma_2, digits = 6), "\n")
    cat("  Gamma_3:         ", format(gamma_3, digits = 6), "\n")
    cat("  Delta0:          ", format(delta0, digits = 6), "\n")
    cat("  Delta1:          ", format(delta1, digits = 6), "\n\n")
    
    cat("------------------------------------------------\n")
  }
  
#######################GBOP2_minSS_TE#################
  if (object$`function` == "GBOP2_minSS_TE") {
    # Extract general information
    design <- object$design
    method <- object$method
    lambdae1 <- object$lambdae1
    lambdae2 <- object$lambdae2
    lambdat1 <- object$lambdat1
    lambdat2 <- object$lambdat2
    gamma <- object$gamma
    expected_sample <- object$expected_sample
    typeI_01 <- object$typeI_01
    typeI_10 <- object$typeI_10
    typeI_00 <- object$typeI_00
    power <- object$Power
    utility <- object$Utility
    
    # Dynamically extract cohort sizes and stopping boundaries
    cohort <- unlist(object[grep("^cohort\\d+$", names(object))])  
    boundary_effi <- unlist(object[grep("^boundary_effi\\d+$", names(object))])  
    boundary_toxi <- unlist(object[grep("^boundary_toxi\\d+$", names(object))])  
    
    # Print the summary
    cat("------------------------------------------------\n")
    cat("PSOGO for minimizing sample size with TE boundary\n")
    cat("Design:            ", design, "\n")
    cat("Method:            ", method, "\n\n")
    
    cat("Cohort Sizes:      ", paste(cohort, collapse = ", "), "\n")
    cat("Boundary Values:\n")
    cat("  Efficacy:        ", paste("stop if <=", boundary_effi, collapse = ", "), "\n")
    cat("  Toxicity:        ", paste("stop if >=", boundary_toxi, collapse = ", "), "\n\n")
    
    cat("Type I Errors:\n")
    cat("  Type I (H01):    ", format(typeI_01, digits = 6), "\n")
    cat("  Type I (H10):    ", format(typeI_10, digits = 6), "\n")
    cat("  Type I (H00):    ", format(typeI_00, digits = 6), "\n\n")
    
    cat("Expected Sample Size:", format(expected_sample, digits = 6), "\n")
    cat("Power (H11):         ", format(power, digits = 6), "\n")
    cat("Utility Value:     ", format(utility, digits = 6), "\n\n")
    
    cat("Parameters:\n")
    cat("  Lambda_e1:       ", format(lambdae1, digits = 6), "\n")
    cat("  Lambda_e2:       ", format(lambdae2, digits = 6), "\n")
    cat("  Lambda_t1:       ", format(lambdat1, digits = 6), "\n")
    cat("  Lambda_t2:       ", format(lambdat2, digits = 6), "\n")
    cat("  Gamma:           ", format(gamma, digits = 6), "\n\n")
    
    cat("------------------------------------------------\n")
  }
  
  
  ########## PSO_design_TE ########################################
  if (object$`function` == "PSO_design_TE") {
    # Extract general information
    design <- object$design
    method <- object$method
    
    # Extract parameters from nested list
    lambdae1 <- object$parameter$lambdae1
    lambdae2 <- object$parameter$lambdae2
    lambdat1 <- object$parameter$lambdat1
    lambdat2 <- object$parameter$lambdat2
    gamma <- object$parameter$gamma
    
    # Extract cohorts and boundaries dynamically from nested lists
    cohort <- unlist(object$cohort)  
    boundary_effi <- unlist(object$boundary_effi)  
    boundary_toxi <- unlist(object$boundary_toxi)  
    
    # Extract expected sample, errors, power, and utility
    expected_sample <- object$expected_sample
    typeI_H01 <- object$`typeI_H01 (safe but futile)`
    typeI_H10 <- object$`typeI_H10 (efficacious but toxic)`
    typeI_H00 <- object$`typeI_H00 (futile and toxic)`
    power <- object$power
    utility <- object$utility
    
    # Print the summary
    cat("------------------------------------------------\n")
    cat("PSO optimal and minimax design with toxicity and efficacy boundaries\n")
    cat("------------------------------------------------\n")
    cat("Design:            ", design, "\n")
    cat("Method:            ", method, "\n\n")
    
    # Parameters
    cat("Parameters:\n")
    cat("  Lambda_e1:       ", format(lambdae1, digits = 6), "\n")
    cat("  Lambda_e2:       ", format(lambdae2, digits = 6), "\n")
    cat("  Lambda_t1:       ", format(lambdat1, digits = 6), "\n")
    cat("  Lambda_t2:       ", format(lambdat2, digits = 6), "\n")
    cat("  Gamma:           ", format(gamma, digits = 6), "\n\n")
    
    # Cohort Sizes
    cat("Cohort Sizes:      ", paste(cohort, collapse = ", "), "\n")
    
    # Boundary Values
    cat("Boundary Values:\n")
    cat("  Efficacy:        ", paste("stop if <=", boundary_effi, collapse = ", "), "\n")
    cat("  Toxicity:        ", paste("stop if >=", boundary_toxi, collapse = ", "), "\n\n")
    
    # Expected Sample, Errors, and Power
    cat("Expected Sample:   ", format(expected_sample, digits = 6), "\n")
    cat("Type I Errors:\n")
    cat("  Type I H01 (Safe but Futile): ", format(typeI_H01, digits = 6), "\n")
    cat("  Type I H10 (Efficacious but Toxic): ", format(typeI_H10, digits = 6), "\n")
    cat("  Type I H00 (Futile and Toxic): ", format(typeI_H00, digits = 6), "\n\n")
    
    cat("Power:             ", format(power, digits = 6), "\n")
    cat("Utility Value:     ", format(utility, digits = 6), "\n")
    cat("------------------------------------------------\n")
  }
  
  
  

  
  ########## PSO_power_TE ########################################
  if (object$`function` == "PSO_power_TE") {
    # Extract and format data
    design <- object$design
    method <- object$method
    parameters <- object$parameter
    cohort <- unlist(object$cohort)
    boundary_effi <- unlist(object$boundary_effi)
    boundary_toxi <- unlist(object$boundary_toxi)
    typeI_01 <- object$typeI_01
    typeI_10 <- object$typeI_10
    typeI_00 <- object$typeI_00
    power <- object$power
    utility <- object$utility
    
    # Display Summary
    cat("------------------------------------------------\n")
    cat("PSO for maximizing power design with toxicity and efficacy boundary\n")
    cat("Design:            ", design, "\n")
    cat("Method:            ", method, "\n\n")
    
    cat("Type I Errors:\n")
    cat("  H01 (Eff+Tox):   ", format(typeI_01, digits = 6), "\n")
    cat("  H10 (Eff only):  ", format(typeI_10, digits = 6), "\n")
    cat("  H00 (Tox only):  ", format(typeI_00, digits = 6), "\n")
    cat("Power:             ", format(power, digits = 6), "\n")
    cat("Utility Value:     ", format(utility, digits = 6), "\n")
    cat("------------------------------------------------\n")
    
    cat("Parameters:\n")
    cat("  Lambda_e1:       ", format(parameters$lambdae1, digits = 6), "\n")
    cat("  Lambda_e2:       ", format(parameters$lambdae2, digits = 6), "\n")
    cat("  Lambda_t1:       ", format(parameters$lambdat1, digits = 6), "\n")
    cat("  Lambda_t2:       ", format(parameters$lambdat2, digits = 6), "\n")
    cat("  Gamma:           ", format(parameters$gamma, digits = 6), "\n\n")
    
    cat("Cohort Sizes:      ", paste(cohort, collapse = ", "), "\n")
    cat("Boundary Values:\n")
    cat("  Efficacy:        ", paste("stop if <=", boundary_effi, collapse = ", "), "\n")
    cat("  Toxicity:        ", paste("stop if >=", boundary_toxi, collapse = ", "), "\n\n")
  }


  
  ########## GBOP2_maxP_TE  ########################################
  if (object$`function` == "GBOP2_maxP_TE") {
    # Extract and format data
    design <- object$design
    method <- object$method
    lambdae1 <- object$lambdae1
    lambdae2 <- object$lambdae2
    lambdat1 <- object$lambdat1
    lambdat2 <- object$lambdat2
    gamma <- object$gamma
    cohort <- c(object$cohort1, object$cohort2, object$cohort3, object$cohort4, object$cohort5)
    boundary_effi <- c(object$boundary_effi1, object$boundary_effi2, object$boundary_effi3, 
                       object$boundary_effi4, object$boundary_effi5)
    boundary_toxi <- c(object$boundary_toxi1, object$boundary_toxi2, object$boundary_toxi3, 
                       object$boundary_toxi4, object$boundary_toxi5)
    typeI_01 <- object$typeI_01
    typeI_10 <- object$typeI_10
    typeI_00 <- object$typeI_00
    power <- object$Power
    utility <- object$Utility
    
    # Print the summary
    cat("------------------------------------------------\n")
    cat("PSOGO for maximizing power with TE boundaries\n")
    cat("Design:             ", design, "\n")
    cat("Method:             ", method, "\n\n")
    
    cat("Cohort Sizes:       ", paste(cohort, collapse = ", "), "\n")
    cat("Boundary Values:\n")
    cat("  Efficacy:         ", paste("stop if <=", boundary_effi, collapse = ", "), "\n")
    cat("  Toxicity:         ", paste("stop if >=", boundary_toxi, collapse = ", "), "\n\n")
    
    cat("Type I Errors:\n")
    cat("  Type I (H01):     ", format(typeI_01, digits = 6), "\n")
    cat("  Type I (H10):     ", format(typeI_10, digits = 6), "\n")
    cat("  Type I (H00):     ", format(typeI_00, digits = 6), "\n\n")
    
    cat("Power (H11):        ", format(power, digits = 6), "\n")
    cat("Utility Value:      ", format(utility, digits = 6), "\n\n")
    
    cat("Parameters:\n")
    cat("  Lambda_e1:        ", format(lambdae1, digits = 6), "\n")
    cat("  Lambda_e2:        ", format(lambdae2, digits = 6), "\n")
    cat("  Lambda_t1:        ", format(lambdat1, digits = 6), "\n")
    cat("  Lambda_t2:        ", format(lambdat2, digits = 6), "\n")
    cat("  Gamma:            ", format(gamma, digits = 6), "\n")
    cat("------------------------------------------------\n")
  }
  
  
  
} ## end of function


