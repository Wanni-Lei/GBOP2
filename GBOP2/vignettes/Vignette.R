## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
#install.packages("GBOP2")
library(GBOP2)

## ----echo=FALSE---------------------------------------------------------------
library(knitr)

# Create the data frame
data <- data.frame(
  Category = c("Sample Size Minimization", "Power Maximization", ""),
  `Efficacy - Single Boundary` = c("GBOP2_minSS_singleE", "GBOP2_maxP_singleE", ""),
  `Efficacy - Dual Boundary` = c("GBOP2_minSS_dualE", "GBOP2_maxP_dualE", ""),
  `Toxicity and Efficacy` = c("GBOP2_minSS_TE", "GBOP2_maxP_TE", "")
)

# Render the table
kable(data, align = "c", col.names = c("", "Efficacy - Single Boundary", "Efficacy - Dual Boundary", "Toxicity and Efficacy"))


## ----message=FALSE, warning=FALSE, eval = FALSE-------------------------------
#  ## PSO-GO
#  optimal_GO_single <- GBOP2_minSS_singleE(
#                    design = "optimal", ## "optimal" or "minimax", "unified"
#                    unified.u = 1, ## specify when design = "unified", u in [0, 1]
#                    nlooks = 2, ## number of interim looks
#                    b1n = 0.2, ## response rate in null hypothesis
#                    b1a = 0.4, ## response rate in alternative hypothesis
#                    err1 = 0.05, ## type I error
#                    nParallel = 3, ## number of PSO-ensemble
#                    minPower = 0.8, ## power
#                    weight = 1, ## weight of sample size under null
#                    maxPatients = 50, ## maximum number of patients
#                    Nmin_cohort1 = 5, ## minimum number of first cohort
#                    Nmin_increase = 5, ## minimum number of increase in each cohort
#                    pso_method = "all", ## choose from "all", ""default", "quantum" or "dexp"
#                    seed = 456,
#                    nSwarm = 64, ## nSwarm in PSO
#                    maxIter = 200, ## maxIter in PSO
#                    nCore = 8) ## number of core

## ----echo=FALSE, out.width="900px", warning=FALSE, message=FALSE--------------
link = system.file("Flowchart", "progressbar.png", package = "GBOP2")
knitr::include_graphics(link)


