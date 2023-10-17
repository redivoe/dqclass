get_quantiles <- function(u, theta, distr){

  qf <- switch(
    distr,
    fgld = qfgld,
    normal = function(u, theta) stats::qnorm(p = u, mean = theta[1], sd = theta[2]),
    kde = function(u, theta) stats::approx(theta$Fhat, theta$x, xout = u, rule = 2, ties = "ordered")$y
  )

  points <- sapply(theta, \(theta_b) qf(u, theta_b))

  return(points)
}


#' Integrated rank-weighted depth
#'
#' Computes the integrated rank-weighted depth based on a data matrix \eqn{X}.
#'
#' @param z Point(s) whose depth is to be computed. If omitted the depths are computed for the points in `X`.
#' @param X Data matrix with respect to which the depth is to be computed, with observations in rows and variables in columns.
#' @param S A matrix containing the directions in the columns.
#' @param n_dir Number of directions to consider. They are equally spaced in two dimensions and uniformly sampled for higher dimensions.
#' @param distr The distribution that is fit to each univariate directional projection.
#' @param weighted A logical indicating whether weights for the directions should be computed.
#' @param sphered A logical indicating whether the depth should be computed on the sphered data.
#'
#' @return A list containing:
#' * `theta`: a list of dimension `n_dir` containing the estimated parameters for each univariate directional projection.
#' * `depths`: a vector with the computed depths for the targe points in `z`.
#'
#' @export
#'
#' @examples
#' library(ggplot2)
#' library(mvtnorm)
#' set.seed(1)
#' X <- rmvnorm(n = 200, mean = c(10, 2), sigma = matrix(c(1, 0.9, 0.9, 2), 2, 2))
#' depths <- dqdepth(X = X, weighted = FALSE)
#'
#' data.frame("x1" = X[,1], "x2" = X[,2], "d" = depths) |>
#'   ggplot(aes(x1, x2, col = d))+
#'   geom_point()+
#'   scale_color_viridis_c()
#'
dqdepth <- function(z,
                    X,
                    S,
                    n_dir = 1e3,
                    distr = "normal",
                    sphered = FALSE,
                    weighted = FALSE) {

  out <- dqdepth_fit(X = X, S = S, n_dir = n_dir, distr = distr, sphered = sphered, weighted = weighted)

  if(missing(z)){
    z <- X
  }
  if(is.null(dim(z))){
    z <- matrix(z, ncol = length(z))
  }

  zc <- t(t(z) - out$mu)
  if(sphered){
    zc_proj <- zc %*% out$W %*% out$S
  }else{
    zc_proj <- zc %*% out$S
  }

  depths <- get_depths(zc_proj, out$theta, distr, n_dir)

  if(weighted){
    depths <- sweep(depths, MARGIN = 2, STATS = out$w, FUN = "*") |>
      rowSums()
    depths <- depths / sum(out$w)

  }else{
    if(is.null(dim(depths))){
      depths <- sum(depths)
    }else{
      depths <- rowMeans(depths)
    }
  }

  return(depths)
}

dqdepth_fit <- function(X, S, n_dir = 1e3, distr, sphered, weighted){

  fit_fun <- switch(
    distr,
    fgld = fit_fgld_ls,
    normal = function(x) c(mean(x), sd_floor(x)),
    kde = kdeF,
    ecdf = function(x) stats::ecdf(x = x)
  )

  if(missing(S)){
    if(ncol(X) == 2){
      S <- get_equi_S(B = n_dir)
    }else{
      S <- get_S(B = n_dir, dim = ncol(X))
    }
  }else{
    n_dir <- ncol(S)
  }

  mu <- colMeans(X)
  Xc <- t(t(X) - mu)

  if(sphered){
    W <- sphere_zca(Xc, centered = TRUE)
    Z <- Xc %*% W
    Z_proj <- Z %*% S
  }else{
    Z_proj <- Xc %*% S
  }

  theta <- apply(Z_proj, 2, fit_fun, simplify = FALSE)

  out <- list("theta" = theta,
              "S" = S,
              "mu" = mu,
              "distr" = distr)

  if(sphered){
    out$W <- W
  }

  if(weighted){
    depths <- get_depths(Z_proj, theta, distr, n_dir)
    delta <- colSums(depths)
    delta <- ifelse(delta < 0, 0, delta)
    w <- delta / sqrt(sum(delta ^ 2))
    out$w <- w
  }

  return(out)
}


#' Integrated rank-weighted depth contours
#'
#' Computes an approximation of the contour level, giving a set points with a predifined value of the depth.
#'
#' @inheritParams dqdepth
#' @param u Level of the contour that is required.
#'
#' @return
#' A matrix containing as rows the points that form the contour.
#'
#' @export
#'
#' @examples
#' library(mvtnorm)
#' set.seed(2)
#' X <- rmvnorm(n = 200, mean = c(10, 2), sigma = matrix(c(1, 0.9, 0.9, 2), 2, 2))
#' contour <- dqcontour(u = 0.9, X = X)
#'
#' plot(X, cex = 0.5)
#' points(rbind(contour, contour[1,]), col = 2, lwd = 3, type = "l")
#'
dqcontour <- function(u,
                      X,
                      S,
                      n_dir = 1e3,
                      distr = "normal",
                      sphered = FALSE,
                      weighted = FALSE){
  out <- dqdepth_fit(X = X, S = S, n_dir = n_dir, distr = distr, sphered = sphered, weighted = weighted)

  q <- get_quantiles(u = u,
                     theta = out$theta,
                     distr = out$distr)
  C <- sweep(out$S, 2, q, "*")
  if(sphered){
    C <- t(t(t(C) %*% solve(out$W)) + out$mu)
  }else{
    C <- t(C + out$mu)
  }

  return(C)
}
