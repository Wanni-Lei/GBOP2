% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PSO_power_dual.gbop2.R
\name{PSO_power_dual}
\alias{PSO_power_dual}
\title{PSO: Power maximizing design with dual boundaries}
\usage{
PSO_power_dual(
  method = "default",
  totalPatients = 50,
  nlooks = 1,
  Nmin_cohort1 = 10,
  Nmin_increase = 10,
  b1n = 0.2,
  b1a = 0.4,
  err1 = 0.05,
  minPower = 0.8,
  seed = 1024,
  nSwarm = 64,
  maxIter = 200
)
}
\arguments{
\item{method}{choose from "optimal" or "minimax"}

\item{totalPatients}{number of total patients}

\item{nlooks}{number of interim looks}

\item{Nmin_cohort1}{minimum number of first cohort}

\item{Nmin_increase}{minimum number of increase in each cohort}

\item{b1n}{Null hypothesis response rate}

\item{b1a}{Alternative hypothesis response rate}

\item{err1}{Type I error rate}

\item{minPower}{power}

\item{seed}{seed for pso}

\item{nSwarm}{nSwarm for pso}

\item{maxIter}{maxIter for pso}
}
\value{
A list on design parameters and operating characteristics
}
\description{
This function implements Particle Swarm Optimization (PSO) to find an power maximizing design with dual boundaries.
}
\examples{
PSO_power_dual(
  method = "default", 
  totalPatients = 50, 
  nlooks = 1, 
  Nmin_cohort1 = 10, 
  Nmin_increase = 10, 
  b1n = 0.2, 
  b1a = 0.4, 
  err1 = 0.05, 
  minPower = 0.8, 
  seed = 1024, 
  nSwarm = 64, 
  maxIter = 200
)

}
