
#' @rdname fit_fgld_ls
#' @export
get_X_fgld<-function(n){
  x2 <- digamma(1:n) - digamma(n + 1)

  X <- cbind(1,
             (1:n) / (n + 1),
             x2,
             rev(-x2))
  return(X)
}

#' @rdname fit_fgld_ls
#' @param n_ Sample size2.
#' @export
get_X_fgld_rob<-function(n, n_){
  x2 <- digamma(seq(1, n, len = n_)) - digamma(n + 1)

  X <- cbind(seq(1, n, len = n_) / (n + 1),
             x2,
             rev(-x2))
  return(X)
}

#' Least squares estimation for the \emph{fgld}
#'
#' @description
#' `fit_fgld_ls()` carries out the least squares estimation of the
#' \emph{fgld} (flattened generalized logistic distribution) with a linear
#' reparametrization.
#'
#' `get_X_fgld()` returns the design matrix needed for the least squares estimation. The design matrix contains the coefficients of the
#' expected order statistics.
#'
#' @param y Data from which to estimate the parameters.
#' @param n Sample size.
#'
#' @return
#' `fit_fgld_ls()` outputs a vector of length 4 containg the estimates of the parameters.
#'
#' `get_X_fgld()` outputs a \eqn{n \times 4} matrix where the i-th row contains the
#' coefficients for expected value of the i-th order statistic of the fgld.
#'
#' @export
#'
#' @examples
#' x <- rnorm(5)
#' get_X_fgld(5)
#' fit_fgld_ls(x)
#'
fit_fgld_ls<-function(y){
  n <- length(y)
  X <- get_X_fgld(n)
  y <- sort(y)
  theta <- qr.solve(X, y)

  if (!any(theta[-1] < 0)) {
    return(theta)
  } else{

    removed <- list(2, 3, 4, c(2, 3), c(2, 4), c(2, 3, 4))
    theta <- lapply(removed, function(i) {
      theta <- qr.solve(X[, -i], y)
      for (j in 1:length(i)) {
        theta <- append(theta, 0, i[j] - 1)
      }
      return(theta)
    })

    feasible <- sapply(theta, \(t) !any(t[-1] < 0))
    res <- sapply(theta, \(t) sum((y - X %*% t)^2))
    i <- which.min(res + ifelse(feasible, 0, Inf))

    return(theta[[i]])
  }
}

invQ <- function(x, theta, Q) {
  sapply(x, \(xi) tryCatch(stats::uniroot(function(u) Q(u, theta) - xi,
                                          interval = c(1 - 1e-10, 1e-10))$root,
                           error = function(cond)
                             ifelse(xi > theta[1], 1 - 1e-10, 1e-10)
                           )
         )
}

dQ <- function(x, theta, log = FALSE, Q, Qprime) {
  u <- invQ(x, theta, Q)
  if (isTRUE(log)) {
    return(- log(Qprime(u, theta)))
  } else{
    return(1 / Qprime(u, theta))
  }
}


#' The flattened generalized logistic distribution (\emph{fgld})
#'
#' Quantile and cumulative distribution functions for the \emph{fgld}.
#'
#' @param u Vector of probabilities.
#' @param x,q Vector of quantiles.
#' @param theta Vector of parameters (of length 4) for the \emph{fgld}.
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
#' theta <- c(1, 2, 0, 0)
#' u <- seq(0.1, 0.9, by = 0.1)
#' qfgld(u, theta)
#'
qfgld <- function(u, theta) {
  theta[1] + theta[2] * u + theta[3] * log(u) - theta[4] * log(1 - u)
}

#' @rdname qfgld
qfgld_deriv <- function(u, theta) {
  ifelse(u >= 1 - 1e-10 | u <= 1e-10,
         1e10,
         theta[2] + theta[3] / u + theta[4] / (1 - u))
}

#' @rdname qfgld
#' @export
pfgld <- function(q, theta) {
  invQ(q, theta, Q = qfgld)
}

#' @rdname qfgld
#' @export
dfgld <- function(x, theta, log = FALSE) {
  dQ(x, theta, log = log, Q = qfgld, Qprime = qfgld_deriv)
}


#' From original to linear parametrization of the \emph{fgld} and viceversa
#'
#' @param alpha Parameters in the original parametrization.
#' @param theta Parameters in the linear parametrization.
#'
#' @return
#' A vector containing parameter values for the other parametrization.
#'
#' @export
#'
#' @examples
#' alpha <- c(1, 2, 0.9, 2)
#' (theta <- origin2lin(alpha))
#' lin2origin(theta)
#'
origin2lin <- function(alpha) {
  theta <- numeric()
  theta[1] <- alpha[1]
  theta[2] <- alpha[2] * alpha[4]
  theta[3] <- alpha[2] * (1 - alpha[3])
  theta[4] <- alpha[2] * alpha[3]

  return(theta)
}

#' @rdname origin2lin
#' @export
lin2origin <- function(theta) {
  alpha <- numeric()
  alpha[1] <- theta[1]
  alpha[2] <- theta[3] + theta[4]
  alpha[3] <- theta[4] / alpha[2]
  alpha[4] <- theta[2] / alpha[2]

  return(alpha)
}


plot_qfgld <- function(theta, add = FALSE, col = 2, ...){
  lim <- 1e-3
  u <- seq(lim, 1-lim, len = 3e2)
  if(add){
    graphics::points(u, qfgld(u, theta), type = "l", col = col, lwd = 4, ylab = "Q(u)", ...)
  }else{
    plot(u, qfgld(u, theta), type = "l", col = col, lwd = 4, ylab = "Q(u)", ...)
  }
}

plot_dfgld <- function(theta, lims = c(0.05, 0.95), add = FALSE, col = 2, ...){
  xlims <- qfgld(lims, theta)
  x_seq <- seq(xlims[1], xlims[2], len = 3e2)
  if(add){
    graphics::points(x_seq, dfgld(x_seq, theta), type = "l", col = col, lwd = 4, ylab = "f(x)", ...)
  }else{
    plot(x_seq, dfgld(x_seq, theta), type = "l", col = col, lwd = 4, ylab = "f(x)", ...)
  }
}

