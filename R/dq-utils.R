##--------------
##  Directions
##--------------

#' Random directions
#'
#' @param B Number of random directions to generate.
#' @param dim Dimension for each direction.
#'
#' @return A matrix of dimension \eqn{dim \times B} where each column is a direction.
#' @export
#'
#' @examples
#' get_S(10, 3)
#'
get_S <- function(B, dim){
  S <- matrix(stats::rnorm(n = dim * B), nrow = B, ncol = dim)
  apply(X = S, 1, function(x) x / sqrt(sum(x ^ 2)), simplify = TRUE)
}

#' Equally spaced directions in two dimensions
#'
#' Function used for two dimensional samples that creates equally spaced directions by dividing the unit circle in the required number of angles.
#'
#' @param B Number of directions.
#'
#' @return A matrix of dimension \eqn{2 \times B} where each column is a direction in the 2D plane.
#' @export
#'
#' @examples
#' get_equi_S(5)
get_equi_S <- function(B){
  angle_seq <- seq(0, 2 * pi, length.out = B + 1)
  S <- rbind(cos(angle_seq), sin(angle_seq))
  return(S[, -1])
}


##-------------------
##  Depth functions
##-------------------

depth_from_u <- function(u) {
  1 - 2 * abs(u - 0.5)
}

get_depths <- function(X, theta_k, distr, n_dir){

  cdf <- switch(
    distr,
    fgld = pfgld,
    normal = function(x, theta) stats::pnorm(q = x, mean = theta[1], sd = theta[2]),
    # kde = function(x, theta) rowMeans(pnorm(q = outer(x, theta$x, "-"), mean = 0, sd = theta$bw))
    kde = function(x, theta) stats::approx(theta$x, theta$Fhat, xout = x, rule = 2, ties = "ordered")$y,
    ecdf = function(x, theta) theta(x)
  )

  depths <- sapply(1:n_dir, \(i) cdf(X[, i], theta_k[[i]])) |>
    depth_from_u()

  return(depths)
}

get_depth <- function(x, theta, distr){

  cdf <- switch(
    distr,
    fgld = pfgld,
    normal = function(x, theta) stats::pnorm(q = x, mean = theta[1], sd = theta[2]),
    # kde = function(x, theta) rowMeans(pnorm(q = outer(x, theta$x, "-"), mean = 0, sd = theta$bw))
    kde = function(x, theta) stats::approx(theta$x, theta$Fhat, xout = x, rule = 2, ties = "ordered")$y,
    ecdf = function(x, theta) theta(x)
  )

  depths <- cdf(x, theta) |>
    depth_from_u()

  return(depths)
}


##--------------------------
##  Sphering and centering
##--------------------------

#' Compute a sphering matrix
#'
#' Computes a sphering (or whitening) matrix for a data matrix \eqn{X}, that is a matrix \eqn{W}, such that \eqn{Z = X
#' W} has the identity matrix as sample covariance matrix. Uses the so-called ZCA (zero-phase component analysis) or
#' Mahalanobis whitening procedure. The whitening matrix is the inverse matrix
#' square root of the sample covariance matrix \eqn{\hat{\Sigma}}.
#' \eqn{W^{ZCA} = \hat{\Sigma}^{-\frac{1}{2}}}
#' The inverse square root is computed via the `svd` decomposition of the
#' centered data matrix.
#'
#' @param X A matrix of dimensions \eqn{n \times p}, where \eqn{n} is the number
#'   of points and \eqn{p} the number of dimensions.
#' @param centered A logical indicating whether the data matrix in argument `X`
#'   has already been centered.
#'
#' @return A sphering matrix of dimensions \eqn{p \times p} that premultiplied
#'   by data matrix `X` gives a sphered data matrix.
#' @export
#'
#' @examples
#' library(mvtnorm)
#' Sigma <- matrix(c(1, 0.9, 0.9, 2), 2, 2)
#' X <- mvtnorm::rmvnorm(n = 100, sigma = Sigma)
#' var(X)
#' W <- sphere_zca(X)
#' Z <- X %*% W
#' var(Z)
#'
sphere_zca <- function(X, centered = FALSE){
  if(centered){
    X_centered <- X
  }else{
    X_centered <- scale(X, TRUE, FALSE)
  }

  tol <- 1e-4
  n <- nrow(X_centered)
  out_svd <- svd(X_centered, nu = 0)
  rank <- sum(out_svd$d > tol * out_svd$d[1])
  return(out_svd$v[, 1:rank] %*% diag(1 / out_svd$d[1:rank]) %*% t(out_svd$v[, 1:rank]) * sqrt(n - 1))
}


scale_centered <- function(Xc, U, xbar){
  k <- ncol(U)
  Xcs <- Xc %*% U
  # xbars <- c(xbar %*% U)
  return(t(t(Xcs) + xbar[1:k]))
}
