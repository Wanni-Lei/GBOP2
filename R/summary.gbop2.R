## summary function
#' summary function for PSOdesign, PSO_power,PSOGO_Design, PSOGO_power, PSODesign_dual,PSOPower_dual, PSOGDesign_Dual, PSOGPower_dual, PSODesign_TE, PSOGO_design_TE
#'
#' @param object PSOdesign, PSO_power,PSOGO_Design, PSOGO_power, PSODesign_dual,PSOPower_dual, PSOGDesign_Dual, PSOGPower_dual, PSODesign_TE, PSOGO_design_TE
#' @param ... ignored arguments
#'
#' @return A summary table
#' @export

summary.gbop2 <- function(object, ...) {

  ########## PSOdesign ########################################
  if (object$`function` == "PSOdesign") {
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
    cat("boundary:  ", paste(unlist(boundary), collapse = ", "), "\n\n")

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
    cat("Boundary:  ", paste(unlist(boundary), collapse = ", "), "\n\n")

    cat("type I error:     ", format(typeI_error, digits = 6), "\n")
    cat("power:            ", format(power, digits = 6), "\n")
    cat("expected sample size (EN(P0)): ", format(expected_sample, digits = 6), "\n")
    cat("utility Value:    ", format(utility, digits = 6), "\n\n")

    cat("------------------------------------------------\n")
  }



  ########## PSOGO_Design ########################################
  if (object$`function` == "PSOGO_Design") {
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
    cat("boundary values:  ", paste(sapply(boundary, function(x) object[[x]]), collapse = ", "), "\n\n")

    cat("type I error:     ", format(typeI_error, digits = 6), "\n")
    cat("power:            ", format(power, digits = 6), "\n")
    cat("expected sample size: ", format(expected_sample, digits = 6), "\n")
    cat("utility value:    ", format(utility, digits = 6), "\n\n")

    cat("------------------------------------------------\n")
  }


  ########## PSOGO_power ########################################
  if (object$`function` == "PSOGO_power") {
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
    cat("boundary values:  ", paste(sapply(boundary, function(x) object[[x]]), collapse = ", "), "\n\n")

    cat("type I error:     ", format(typeI_error, digits = 6), "\n")
    cat("power:            ", format(power, digits = 6), "\n")
    cat("expected sample size (en(p0)): ", format(expected_sample, digits = 6), "\n")
    cat("utility value:    ", format(utility, digits = 6), "\n\n")


    cat("------------------------------------------------\n")
  }


  ########## PSODesign_dual ########################################
  if (object$`function` == "PSODesign_dual") {
    # Extract relevant fields from object
    design <- object$design
    weight <- object$weight
    method <- object$method
    parameters <- object$parameters
    cohort <- object$cohort
    boundary <- object$boundary
    typeI_error <- object$`Type I Error`
    power <- object$Power
    expected_sample <- object$`Expected Sample Size`
    utility <- object$Utility

    cat("------------------------------------------------\n")
    cat("PSO for optimal/minimax design with dual-boundary \n")
    cat("design:           ", design, "\n")
    cat("weight:           ", weight, "\n")
    cat("method:           ", method, "\n\n")

    cat("parameters:\n")
    cat("  lambda1:        ", format(parameters$lambda1, digits = 6), "\n")
    cat("  lambda_grad1:   ", format(parameters$lambda_grad1, digits = 6), "\n")
    cat("  lambda_grad2:   ", format(parameters$lambda_grad2, digits = 6), "\n")
    cat("  gamma_1:        ", format(parameters$Gamma_1, digits = 6), "\n")
    cat("  gamma_2:        ", format(parameters$Gamma_2, digits = 6), "\n")
    cat("  gamma_3:        ", format(parameters$Gamma_3, digits = 6), "\n")
    cat("  delta0:         ", format(parameters$delta0, digits = 6), "\n")
    cat("  delta1:         ", format(parameters$delta1, digits = 6), "\n\n")

    cat("cohort sizes:     ", paste(unlist(cohort), collapse = ", "), "\n")
    cat("boundary values:\n")
    cat("  boundary 1:     ", paste(boundary[[1]], collapse = ", "), "\n")
    cat("  boundary 2:     ", paste(boundary[[2]], collapse = ", "), "\n\n")

    cat("type I error:     ", format(typeI_error, digits = 6), "\n")
    cat("power:            ", format(power, digits = 6), "\n")
    cat("expected sample size: ", format(expected_sample, digits = 6), "\n")
    cat("utility value:    ", format(utility, digits = 6), "\n\n")

    cat("------------------------------------------------\n")
  }



  ########## PSOPower_dual ########################################
  if (object$`function` == "PSOPower_dual") {
    # Extract relevant fields from object
    design <- object$design
    method <- object$method
    parameters <- object$parameter
    cohort <- object$cohort
    boundary <- object$boundary
    typeI_error <- object$`Type I Error`
    power <- object$Power
    expected_sample <- object$`Expected Sample Size`
    utility <- object$Utility

    cat("------------------------------------------------\n")
    cat("PSO for maximizing power with dual-boundary\n")
    cat("design:           ", design, "\n")
    cat("method:           ", method, "\n\n")

    cat("parameters:\n")
    cat("  lambda1:        ", format(parameters$lambda1, digits = 6), "\n")
    cat("  lambda_grad1:   ", format(parameters$lambda_grad1, digits = 6), "\n")
    cat("  lambda_grad2:   ", format(parameters$lambda_grad2, digits = 6), "\n")
    cat("  gamma_1:        ", format(parameters$gamma_1, digits = 6), "\n")
    cat("  gamma_2:        ", format(parameters$gamma_2, digits = 6), "\n")
    cat("  gamma_3:        ", format(parameters$gamma_3, digits = 6), "\n")
    cat("  delta0:         ", format(parameters$delta0, digits = 6), "\n")
    cat("  delta1:         ", format(parameters$delta1, digits = 6), "\n\n")

    cat("cohort sizes:     ", paste(unlist(cohort), collapse = ", "), "\n")
    cat("boundary values:\n")
    cat("  boundary 1:     ", paste(boundary[[1]], collapse = ", "), "\n")
    cat("  boundary 2:     ", paste(boundary[[2]], collapse = ", "), "\n\n")

    cat("type I error:     ", format(typeI_error, digits = 6), "\n")
    cat("power:            ", format(power, digits = 6), "\n")
    cat("expected sample size: ", format(expected_sample, digits = 6), "\n")
    cat("utility value:    ", format(utility, digits = 6), "\n\n")

    cat("------------------------------------------------\n")
  }


  ########## PSOGDesign_Dual ########################################
  if (object$`function` == "PSOGDesign_Dual") {
    # Extract relevant fields from object
    design <- object$design
    weight <- object$weight
    method <- object$method
    parameters <- object[grepl("parameter", names(object))]
    cohort <- c(object$cohort1, object$cohort2)
    boundary1 <- object[grepl("cohort.*bd1", names(object))]
    boundary2 <- object[grepl("cohort.*bd2", names(object))]
    typeI_error <- object$`Type.I.Error`
    power <- object$Power
    expected_sample <- object$`Expected.Sample.Size`
    utility <- object$Utility

    cat("------------------------------------------------\n")
    cat("PSOGO for optimal/minimax with dual-boundary\n")
    cat("design:           ", design, "\n")
    cat("weight:           ", weight, "\n")
    cat("method:           ", method, "\n\n")

    cat("parameters:\n")
    cat("  lambda1:        ", format(parameters$parameter.lambda1, digits = 6), "\n")
    cat("  lambda_grad1:   ", format(parameters$parameter.lambda_grad1, digits = 6), "\n")
    cat("  lambda_grad2:   ", format(parameters$parameter.lambda_grad2, digits = 6), "\n")
    cat("  gamma_1:        ", format(parameters$parameter.gamma_1, digits = 6), "\n")
    cat("  gamma_2:        ", format(parameters$parameter.gamma_2, digits = 6), "\n")
    cat("  gamma_3:        ", format(parameters$parameter.gamma_3, digits = 6), "\n")
    cat("  delta0:         ", format(parameters$parameter.delta0, digits = 6), "\n")
    cat("  delta1:         ", format(parameters$parameter.delta1, digits = 6), "\n\n")

    cat("cohort sizes:     ", paste(unlist(cohort), collapse = ", "), "\n")
    cat("boundary values:\n")
    cat("  boundary 1:     ", paste(unlist(boundary1), collapse = ", "), "\n")
    cat("  boundary 2:     ", paste(unlist(boundary2), collapse = ", "), "\n\n")

    cat("type I error:     ", format(typeI_error, digits = 6), "\n")
    cat("power:            ", format(power, digits = 6), "\n")
    cat("expected sample size: ", format(expected_sample, digits = 6), "\n")
    cat("utility value:    ", format(utility, digits = 6), "\n\n")


    cat("------------------------------------------------\n")
  }



  ########## PSOGPower_dual ########################################
  if (object$`function` == "PSOGPower_dual") {
    # Extract relevant fields from object
    design <- object$design
    method <- object$method
    parameters <- object[grepl("parameter", names(object))]
    cohort <- c(object$cohort1, object$cohort2)  # Extract cohort sizes
    boundary1 <- c(object$cohort1bd1, object$cohort2bd1)  # Extract boundary 1 values
    boundary2 <- c(object$cohort1bd2, object$cohort2bd2)  # Extract boundary 2 values
    typeI_error <- object$`Type.I.Error`
    power <- object$Power
    expected_sample <- object$`Expected.Sample.Size`
    utility <- object$Utility

    cat("------------------------------------------------\n")
    cat("PSOGO for maximizing power with dual-boundary\n")

    cat("design:           ", design, "\n")
    cat("method:           ", method, "\n\n")

    cat("parameters:\n")
    cat("  lambda1:        ", format(parameters$parameter.lambda1, digits = 6), "\n")
    cat("  lambda_grad1:   ", format(parameters$parameter.lambda_grad1, digits = 6), "\n")
    cat("  lambda_grad2:   ", format(parameters$parameter.lambda_grad2, digits = 6), "\n")
    cat("  gamma_1:        ", format(parameters$parameter.gamma_1, digits = 6), "\n")
    cat("  gamma_2:        ", format(parameters$parameter.gamma_2, digits = 6), "\n")
    cat("  gamma_3:        ", format(parameters$parameter.gamma_3, digits = 6), "\n")
    cat("  delta0:         ", format(parameters$parameter.delta0, digits = 6), "\n")
    cat("  delta1:         ", format(parameters$parameter.delta1, digits = 6), "\n\n")

    cat("cohort sizes:     ", paste(cohort, collapse = ", "), "\n")
    cat("boundary values:\n")
    cat("  boundary 1:     ", paste(boundary1, collapse = ", "), "\n")
    cat("  boundary 2:     ", paste(boundary2, collapse = ", "), "\n\n")

    cat("type I error:     ", format(typeI_error, digits = 6), "\n")
    cat("power:            ", format(power, digits = 6), "\n")
    cat("expected sample size: ", format(expected_sample, digits = 6), "\n")
    cat("utility:    ", format(utility, digits = 6), "\n\n")

    cat("------------------------------------------------\n")
  }


  ########## PSODesign_TE ########################################
  if (object$`function` == "PSODesign_TE") {
    # Extract relevant fields from object
    design <- object$design
    method <- object$method
    parameters <- object$parameter
    cohort <- unlist(object$cohort)  # Cohort sizes
    boundary_effi <- unlist(object$boundary_effi)  # Efficacy boundaries
    boundary_toxi <- unlist(object$boundary_toxi)  # Toxicity boundaries
    expected_sample <- object$expected_sample
    typeI_01 <- object$typeI_01
    typeI_10 <- object$typeI_10
    typeI_00 <- object$typeI_00
    power <- object$power
    utility <- object$utility

    cat("------------------------------------------------\n")
    cat("PSO for minimax/optimal design with toxicity and efficacy boundary\n")
    cat("design:           ", design, "\n")
    cat("method:           ", method, "\n\n")

    cat("parameters:\n")
    cat("  lambda_e1:      ", format(parameters$lambdae1, digits = 6), "\n")
    cat("  lambda_e2:      ", format(parameters$lambdae2, digits = 6), "\n")
    cat("  lambda_t1:      ", format(parameters$lambdat1, digits = 6), "\n")
    cat("  lambda_t2:      ", format(parameters$lambdat2, digits = 6), "\n")
    cat("  gamma:          ", format(parameters$gamma, digits = 6), "\n\n")

    cat("cohort sizes:     ", paste(cohort, collapse = ", "), "\n")
    cat("boundary values:\n")
    cat("  efficacy boundary: ", paste(boundary_effi, collapse = ", "), "\n")
    cat("  toxicity boundary: ", paste(boundary_toxi, collapse = ", "), "\n\n")

    cat("expected sample size: ", format(expected_sample, digits = 6), "\n")
    cat("type I error (H01):   ", format(typeI_01, digits = 6), "\n")
    cat("type I error (H10):   ", format(typeI_10, digits = 6), "\n")
    cat("type I error (H00):   ", format(typeI_00, digits = 6), "\n")
    cat("power:                ", format(power, digits = 6), "\n")
    cat("utility value:        ", format(utility, digits = 6), "\n\n")

    cat("------------------------------------------------\n")
  }


  ########## PSOGO_design_TE Summary ########################################
  if (object$`function` == "PSOGO_design_TE") {
    # Extract relevant fields from the object
    design <- object$design
    method <- object$method
    parameters <- object$parameter
    cohort <- unlist(object$cohort)  # Cohort sizes
    boundary_effi <- unlist(object$boundary_effi)  # Efficacy boundaries
    boundary_toxi <- unlist(object$boundary_toxi)  # Toxicity boundaries
    expected_sample <- object$expected_sample
    typeI_01 <- object$typeI_01
    typeI_10 <- object$typeI_10
    typeI_00 <- object$typeI_00
    power <- object$power
    utility <- object$utility

    cat("------------------------------------------------\n")
    cat("PSOGO for maximum/optimal design with toxicity and efficacy boundary\n")
    cat("design:           ", design, "\n")
    cat("method:           ", method, "\n\n")

    cat("parameters:\n")
    cat("  lambda_e1:      ", format(parameters$lambdae1, digits = 6), "\n")
    cat("  lambda_e2:      ", format(parameters$lambdae2, digits = 6), "\n")
    cat("  lambda_t1:      ", format(parameters$lambdat1, digits = 6), "\n")
    cat("  lambda_t2:      ", format(parameters$lambdat2, digits = 6), "\n")
    cat("  gamma:          ", format(parameters$gamma, digits = 6), "\n\n")

    cat("cohort sizes:     ", paste(cohort, collapse = ", "), "\n")
    cat("boundary values:\n")
    cat("  efficacy boundary: ", paste(boundary_effi, collapse = ", "), "\n")
    cat("  toxicity boundary: ", paste(boundary_toxi, collapse = ", "), "\n\n")

    cat("expected sample size: ", format(expected_sample, digits = 6), "\n")
    cat("type I error (H01):   ", format(typeI_01, digits = 6), "\n")
    cat("type I error (H10):   ", format(typeI_10, digits = 6), "\n")
    cat("type I error (H00):   ", format(typeI_00, digits = 6), "\n")
    cat("power:                ", format(power, digits = 6), "\n")
    cat("utility value:        ", format(utility, digits = 6), "\n\n")


    cat("------------------------------------------------\n")
  }












} ## end of summary function

## example
# summary.gbop2(design)








