# Internal environment to store the cluster object
.GBOP2_env <- new.env(parent = emptyenv())

#' Initialize parallel cluster
#'
#' Creates and registers a parallel backend using the specified number of cores.
#' Falls back to sequential execution if nCore <= 1.
#'
#' @param nCore Number of cores to use (default is 2).
#' @export
init_cluster <- function(nCore = 2) {
  if (nCore > 2L) {
    warning("Using more than 2 cores may not be allowed on some systems. Proceed only if you are certain this is permitted.")
  }
  
  if (is.null(.GBOP2_env$cl) && nCore > 1L) {
    .GBOP2_env$cl <- parallel::makePSOCKcluster(nCore)
    doParallel::registerDoParallel(.GBOP2_env$cl)
    message("Cluster initialized with ", nCore, " cores.")
  } else {
    foreach::registerDoSEQ()
    message("Sequential mode (single core).")
  }
  
  invisible(.GBOP2_env$cl)
}

#' Get current cluster
#'
#' Returns the current parallel cluster object, if initialized.
#'
#' @return A cluster object or NULL.
#' @export
get_cluster <- function() {
  .GBOP2_env$cl
}

#' Stop and clean up the cluster
#'
#' Stops the currently running parallel cluster and reverts to sequential execution.
#'
#' @export
stop_cluster <- function() {
  if (!is.null(.GBOP2_env$cl)) {
    parallel::stopCluster(.GBOP2_env$cl)
    .GBOP2_env$cl <- NULL
    foreach::registerDoSEQ()
    message("Cluster stopped, returned to sequential backend.")
  }
  invisible(NULL)
}


