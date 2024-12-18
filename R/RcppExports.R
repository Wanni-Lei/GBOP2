# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

create_grid <- function(n) {
    .Call(`_GBOP2_create_grid`, n)
}

filter_and_calculate <- function(grid, n) {
    .Call(`_GBOP2_filter_and_calculate`, grid, n)
}

calculate_and_add_probabilities <- function(nobs_grid, p) {
    .Call(`_GBOP2_calculate_and_add_probabilities`, nobs_grid, p)
}

den_cpp <- function(x, y, n, p) {
    .Call(`_GBOP2_den_cpp`, x, y, n, p)
}

create_tt <- function(s, npt, t_bound) {
    .Call(`_GBOP2_create_tt`, s, npt, t_bound)
}

create_rr <- function(s, npt, r_bound, r_interm) {
    .Call(`_GBOP2_create_rr`, s, npt, r_bound, r_interm)
}

efftox_recursive_optim <- function(interm, npt, p, bound_eff, bound_tox) {
    .Call(`_GBOP2_efftox_recursive_optim`, interm, npt, p, bound_eff, bound_tox)
}

set_seed <- function(seed) {
    invisible(.Call(`_GBOP2_set_seed`, seed))
}

Exacterror <- function(nobs, ncut, pnull, palter) {
    .Call(`_GBOP2_Exacterror`, nobs, ncut, pnull, palter)
}

exact_error_recursive2_Rcpp <- function(nobs, ncut, pnull, palter, ntotal) {
    .Call(`_GBOP2_exact_error_recursive2_Rcpp`, nobs, ncut, pnull, palter, ntotal)
}

my_dbinom <- function(x, size, prob) {
    .Call(`_GBOP2_my_dbinom`, x, size, prob)
}

my_pbinom <- function(q, size, prob) {
    .Call(`_GBOP2_my_pbinom`, q, size, prob)
}

dbinom_product <- function(vXs, nobs, prob) {
    .Call(`_GBOP2_dbinom_product`, vXs, nobs, prob)
}

GetocBiRcpp <- function(seed, nsim, contrast, nobs, b, b2, pow2, dprior, ptrue, phi, fff) {
    .Call(`_GBOP2_GetocBiRcpp`, seed, nsim, contrast, nobs, b, b2, pow2, dprior, ptrue, phi, fff)
}

GridSearchBiRcpp <- function(seed, contrast, nobs, dprior, b1, b2, pow, pn, pa, cutstart, nsim, err1, ffff) {
    .Call(`_GBOP2_GridSearchBiRcpp`, seed, contrast, nobs, dprior, b1, b2, pow, pn, pa, cutstart, nsim, err1, ffff)
}

set_seed_dual <- function(seed) {
    invisible(.Call(`_GBOP2_set_seed_dual`, seed))
}

exact_error_recursive_Rcpp <- function(nobs, ncut1, ncut2, pnull, palter, ntotal) {
    .Call(`_GBOP2_exact_error_recursive_Rcpp`, nobs, ncut1, ncut2, pnull, palter, ntotal)
}

GetocBiRcpp_dual <- function(seed, nsim, contrast, nobs, b, b_grad1, b_grad2, pow1, pow2, pow3, dprior, ptrue, phi, delta0, delta1, fff) {
    .Call(`_GBOP2_GetocBiRcpp_dual`, seed, nsim, contrast, nobs, b, b_grad1, b_grad2, pow1, pow2, pow3, dprior, ptrue, phi, delta0, delta1, fff)
}

