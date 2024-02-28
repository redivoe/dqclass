# Function to check if an object is a vector
is_vector <- function(x) {
  is.numeric(x) && is.null(dim(x))
}

# Function to check if an object is a matrix
is_matrix <- function(x) {
  is.numeric(x) && !is.null(dim(x)) && length(dim(x)) == 2
}

#' Assignment problem with two vectors of different size
#'
#' @param x1 vector/matrix 1
#' @param x2 vector/matrix 2
#'
#' @return Indices of `x2` which minimize the sum of distances to `x1`.
#' @export
#'
#' @examples
#' x1 <- c(-10, 3, 4, 5)
#' x2 <- c(3, 4)
#' assign_diff_lengths(x1, x2)
#'
assign_diff_lengths <- function(x1, x2){

  if(is_vector(x1) & is_vector(x2)){
    n1 <- length(x1)
    n2 <- length(x2)

    distances <- array(dim = c(n1, n2))
    for(i in 1:n1) for(j in 1:n2) distances[i, j] <- (x1[i] - x2[j])^2

  }else if(is_matrix(x1) & is_matrix(x2)){
      if(ncol(x1) != ncol(x2)){
        stop("Points in x1 and x2 must have the same dimension")
      }
      n1 <- nrow(x1)
      n2 <- nrow(x2)

      distances <- array(dim = c(n1, n2))
      for(i in 1:n1) for(j in 1:n2) distances[i, j] <- sum((x1[i, ] - x2[j, ])^2)
  }

  lp_sol <- lpSolve::lp.transport(distances,
                                  row.signs = rep("=", n1),
                                  row.rhs = rep(1, n1),
                                  col.signs = rep("<=", n2),
                                  col.rhs = rep(1, n2))
  idx <- max.col(lp_sol$solution)

  return(idx)
}

#' Robust least squares estimation for the \emph{fgld}
#'
#' @param y Data vector.
#' @param alpha Proportion of contamination.
#'
#' @return List containing:
#'
#' * `theta`: estimated vector of parameters.
#' * `idx`: indices of the observations that have been identified as outliers.
#' * `obj`: vector of the objective evaluations at each iteration.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- c(rlogis(18), -10, 10)
#' fit_fgld_ls_robust(x, alpha = 0.1)
#'
fit_fgld_ls_robust <- function(y, alpha){
  n <- length(y)
  h <- floor(n * (1 - alpha))
  X <- get_X_fgld(h)

  tol <- 1e-4
  max_iter <- 10
  obj <- c(1e10, rep(NA, max_iter))

  iter <- 1

  # initial estimation
  idx <- sample(n, 5)

  repeat{
    ls_out <- fit_fgld_ls_obj(y[idx])
    theta <- ls_out$theta

    obj[iter + 1] <- ls_out$obj
    ratio <- abs(obj[iter + 1] - obj[iter]) / obj[iter]

    if(iter >= max_iter || ratio < tol) break

    yhat <- c(X %*% theta)
    idx <- assign_diff_lengths(yhat, y)
    iter <- iter + 1
  }

  out_idx <- setdiff(1:n, idx)

  return(list(
    "theta" = theta,
    "idx" = out_idx,
    "obj" = obj[3:(iter + 1)],
    "res" = ls_out$res
  ))
}
