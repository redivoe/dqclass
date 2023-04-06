
#' @rdname fit_fgldh_ls
#' @export
get_X_fgldh <- function(n){
  x3 <- exp(lgamma(n:1 - 0.5) + lgamma(n + 1) - lgamma(n:1) - lgamma(n - 0.5 + 1))
  x5 <- digamma(1:n) - digamma(n + 1)

  X <- cbind(1,
             (1:n) / (n + 1),
             x3,
             rev(-x3),
             x5,
             rev(-x5))

  return(X)
}


#' Least squares estimation for the \emph{fgldh}
#'
#' @description
#' `fit_fgldh_ls()` carries out the least squares estimation of the
#' \emph{fgld} (flattened generalized logistic distribution) with a linear
#' reparametrization.
#'
#' `get_X_fgldh()` returns the design matrix needed for the least squares estimation. The design matrix contains the coefficients of the
#' expected order statistics.
#'
#' @param y Data from which to estimate the parameters.
#' @param n Sample size.
#'
#' @return
#' `fit_fgldh_ls()` outputs a vector of length 6 containg the estimates of the parameters.
#'
#' `get_X_fgldh()` outputs a \eqn{n \times 6} matrix where the i-th row contains the
#' coefficients for expected value of the i-th order statistic of the fgld.
#'
#' @export
#'
#' @examples
#' x <- rnorm(10)
#' # get_X_fgldh(10)
#' fit_fgldh_ls(x)
#'
fit_fgldh_ls<-function(y){
  n <- length(y)
  X <- get_X_fgldh(n)
  y <- sort(y)

  A <- rbind(0, diag(1, 5))
  b <- rep(0, 5)

  out_solveQP <- quadprog::solve.QP(
    Dmat = crossprod(X, X),
    dvec = crossprod(X, y),
    Amat = A,
    bvec = b
  )

  return(out_solveQP$solution)
}



#' The flattened generalized logistic distribution with heavy tails (\emph{fgldh})
#'
#' Quantile and cumulative distribution functions for the \emph{fgldh}.
#'
#' @param u Vector of probabilities.
#' @param x,q Vector of quantiles.
#' @param theta Vector of parameters (of length 4) for the \emph{fgldh}.
#' @param log Logical indicating whether the density should be returned on the log scale.
#'
#' @return
#' `qfgld` gives the quantile function, `pfgld` gives the distribution functions.
#' The length of the output is the same as that of the input `x/u`.
#'
#'
#' @export
#'
#' @examples
#' theta <- c(1, 2, 0, 0, 1, 2)
#' u <- seq(0.1, 0.9, by = 0.1)
#' qfgldh(u, theta)
#'
qfgldh <- function(u, theta) {
  theta[1] +
    theta[2] * u +
    theta[3] / sqrt(1 - u) -
    theta[4] / sqrt(u) +
    theta[5] * log(u) -
    theta[6] * log(1 - u)
}


#' @rdname qfgldh
qfgldh_deriv <- function(u, theta) {
  ifelse(u >= 1 - 1e-10 | u <= 1e-10,
         1e10,
         theta[2] + theta[3] * 0.5 * (1 - u) ^ (-1.5) + theta[4] * 0.5 * u ^ (-1.5) +
           theta[5] / u + theta[6] / (1 - u)
  )
}

#' @rdname qfgldh
#' @export
pfgldh <- function(q, theta) {
  invQ(q, theta, Q = qfgldh)
}

#' @rdname qfgldh
#' @export
dfgldh <- function(x, theta, log = FALSE) {
  dQ(x, theta, log = log, Q = qfgldh, Qprime = qfgldh_deriv)
}
