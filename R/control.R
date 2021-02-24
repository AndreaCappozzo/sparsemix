#' @export
control_EM <- function(itermax = 1000,
                       tol = 1e-08,
                       err = .Machine$double.xmax / 2) {
  list(itermax = itermax,
       tol = tol,
       err = err)
}
