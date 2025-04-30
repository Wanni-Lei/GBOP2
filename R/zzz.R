# In a file like zzz.R

.onUnload <- function(libpath) {
  # Clean up any parallel backends if you used them
  stop_cluster()  # Only if you defined this helper
  library.dynam.unload("GBOP2", libpath)
}

#' @useDynLib GBOP2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
