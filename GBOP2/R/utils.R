#' @keywords internal
`%operator%` <- function(lhs, rhs) {
  if (!is.null(get_cluster())) {
    foreach::`%dopar%`(lhs, rhs)
  } else {
    foreach::`%do%`(lhs, rhs)
  }
}
