
#' Directional quantile classifier
#'
#' @param X A data matrix with observations in rows and variables in columns.
#' @param y A vector of class labels for each sample of the training set.
#' @param distr The distribution that is fit to each univariate directional projection.
#' @param n_dir Number of directions to consider. They are equally spaced in two dimensions and uniformly sampled for higher dimensions.
#' @param S A matrix containing the directions in the columns.
#' @param weighted A logical indicating whether weights for the directions should be computed.
#' @param sphered A logical indicating whether groups should be sphered.
#' @param predicted A logical indicating whether the predicted values from the training sample should be returned.
#'
#' @return A list containing:
#' * `theta`: a list of dimension `n_dir` containing the estimated parameters for each univariate directional projection.
#' * `sphere`: a list of group specific sphering matrices.
#' @export
#'
#' @examples
#' data("iris")
#' out_train <- dqclass_train(iris[, -5], iris[, 5], predicted = TRUE)
#' table(iris[, 5], out_train$class)
#'
dqclass_train <- function(X, y,
                          distr = c("normal", "kde", "fgld"),
                          n_dir = 1e3,
                          S,
                          predicted = FALSE,
                          weighted = TRUE,
                          sphered = TRUE){

  distr <- match.arg(distr)

  fit_fun <- switch(
    distr,
    fgld = fit_fgld_ls,
    normal = function(x) c(mean(x), sd_floor(x)),
    kde = kdeF
  )

  if(missing(S)){
    if(ncol(X) == 2){
      S <- get_equi_S(B = n_dir)
    }else{
      S <- get_S(B = n_dir, dim = ncol(X))
    }
  }else{
    if(nrow(S) == ncol(X)){
      n_dir <- ncol(S)
    }else{
      stop("The matrix of directions S that has been provided does not have the
            appropriate dimensions. The number of rows must equal the dimension of the samples in X.")
    }
  }



  K <- length(unique(y))

  X_split <- lapply(by(X, y, identity), as.matrix)

  if(sphered){
    sphere <- lapply(X_split, sphere_zca)
    X_split <- mapply(\(X, W) X %*% W,
                      X_split, sphere,
                      SIMPLIFY = FALSE)
    X_split_proj <- lapply(X_split, \(X) X %*% S)
  }else{
    X_split_proj <- lapply(X_split, \(X) X %*% S)
  }

  theta <- lapply(X_split_proj, \(x) apply(x, 2, fit_fun, simplify = FALSE))

  cl <- match.call()
  out <- list(
    "theta" = theta,
    "S" = S,
    "n_dir" = n_dir,
    "distr" = distr,
    "weighted" = weighted,
    "sphered" = sphered,
    "call" = cl
  )

  if(weighted){
    dd_train <- lapply(X_split_proj,
                       \(X_k) lapply(theta,
                                     \(theta_k) get_depths(X_k, theta_k, distr, n_dir)))

    delta <- list()
    for(i in 1:K){
      delta[[i]] <- dd_train[[i]][[i]] - do.call(pmax, dd_train[[i]][-i])
    }
    delta <- do.call(rbind, delta) |>
      colSums()
    delta <- ifelse(delta < 0, 0, delta)

    out$w <- delta / sqrt(sum(delta ^ 2))
  }

  if(sphered){
    out$sphere <-  sphere
  }
  if(predicted){
    out <- c(out,
             predict_dqclass(X, out))
  }

  return(out)
}


predict_dqclass <- function(X, out_train){

  X <- as.matrix(X)

  if(out_train$sphered){
    Xs_split <- lapply(out_train$sphere, \(W_k) X %*% W_k)
    Xs_split_proj <- lapply(Xs_split, \(X_k) X_k %*% out_train$S)
    depths <- mapply(FUN = \(X, theta) get_depths(X, theta, out_train$distr, out_train$n_dir),
                     Xs_split_proj, out_train$theta,
                     SIMPLIFY = FALSE)
  }else{
    X_proj <- X %*% out_train$S
    depths <- lapply(out_train$theta,
                     \(theta_k) get_depths(X_proj, theta_k , out_train$distr, out_train$n_dir))
  }

  if(out_train$weighted){
    depths <- sapply(depths, \(x) sweep(x, 2, out_train$w, "*") |> rowSums())
    depths <- depths / sum(out_train$w)
  }else{

    depths <- sapply(depths, rowMeans)

  }

  class <- apply(depths, 1, which.max)

  return(list("class" = class,
              "depths" = depths))
}


dq_train_test <- function(train, test, n_dir = 500,
                          S,
                          distr = "fgld",
                          weighted = TRUE, sphered = TRUE){

  X_train <- train[, -ncol(train)]
  y_train <- train[, ncol(train)]
  X_test <- test[, -ncol(train)]
  y_test <- test[, ncol(train)]

  out_train <- dqclass_train(
    X = X_train,
    y = y_train,
    S = S,
    distr = distr,
    weighted = weighted,
    sphered = sphered
  )

  y_test_hat <- predict_dqclass(X = X_test, out_train)$class

  return(mean(y_test != y_test_hat))
}
