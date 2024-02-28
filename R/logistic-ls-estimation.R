
#' @rdname fit_logis_ls
#' @export
get_X_logis <- function(n) {
  di <- digamma(1:n)

  X <- cbind(1,
             di - rev(di))

  return(X)
}

#' @rdname fit_logis_ls
#' @export
fit_logis_ls.quadprog <- function(y){
  n <- length(y)
  X <- get_X_logis(n)
  y <- sort(y)

  A <- rbind(0, 1)
  b <- 0

  out_solveQP <- quadprog::solve.QP(
    Dmat = crossprod(X, X),
    dvec = crossprod(X, y),
    Amat = A,
    bvec = b
  )

  return(out_solveQP$solution)
}

#' Least squares estimation for the logistic distribution
#'
#' @param y Data from which to estimate the parameters.
#' @param n Sample size.
#'
#' @return
#' `fit_logis_ls()` outputs a vector of length 2 containg the estimates of the parameters.
#'
#' `get_X_logis()` outputs a \eqn{n \times 2} matrix where the i-th row contains the
#' coefficients for expected value of the i-th order statistic of the logistic distribution.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' x <- rlogis(100, 3, 2)
#' fit_logis_ls(x)
#'
fit_logis_ls <- function(y){
  n <- length(y)
  X <- get_X_logis(n)
  y <- sort(y)

  theta <- qr.solve(X, y)

  if (theta[2] < 0) {
    theta <- c(qr.solve(X[,-2], y), 0)
  }

  return(theta)
}


#' @rdname fit_logis_ls
#' @export
fit_logis_ls_obj <- function(y){
  n <- length(y)
  X <- get_X_logis(n)
  y <- sort(y)

  qrX <- qr(X)
  theta <- qr.coef(qrX, y)

  tol_var <- 1e-2

  if (theta[2] <= tol_var) {
    qrX <- qr(X[,-2])
    theta <- c(qr.solve(qrX, y), tol_var)
  }

  res <- qr.resid(qrX, y)
  obj <- sum(res ^ 2)

  return(list("theta" = theta,
              "obj" = obj,
              "res" = res))
}
