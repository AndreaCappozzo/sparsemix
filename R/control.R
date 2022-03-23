#' @export
control_EM <- function(itermax = 1000,
                       tol = 1e-08,
                       err = .Machine$double.xmax / 2,
                       tol_zero_var = sqrt(.Machine$double.eps)) {
  list(
    itermax = itermax,
    tol = tol,
    err = err,
    tol_zero_var = tol_zero_var
  )
}
