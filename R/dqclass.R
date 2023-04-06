
#' Directional quantile classifier
#'
#' @inheritParams dqdepth
#' @param X A data matrix with observations in rows and variables in columns.
#' @param y A vector of class labels for each sample of the training set.
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
                          weighted = FALSE,
                          sphered = TRUE){

  distr <- match.arg(distr)

  fit_fun <- switch(
    distr,
    fgld = fit_fgld_ls,
    normal = function(x) c(mean(x), sd_floor(x)),
    kde = kdeF
  )

  # univariate data
  if(is.vector(X) || ncol(X) == 1){
    X_split <- lapply(by(X, y, identity), as.matrix)
    K <- length(unique(y))
    theta <- lapply(X_split, fit_fun)

    cl <- match.call()
    out <- list(
      "theta" = theta,
      "distr" = distr,
      "univariate" = TRUE,
      "call" = cl
    )

    if(predicted){
      out <- c(out,
               predict_dqclass(out, X))
    }

    return(out)
  }

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
    "univariate" = FALSE,
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
             predict_dqclass(out, X))
  }

  return(out)
}


# predict.dq <- function(dq, ...){
#   predict()
# }

#' @rdname dqclass_train
#' @param out_train Ouput list from `dqclass_train`.
#' @export
predict_dqclass <- function(out_train, X){

  X <- as.matrix(X)

  if(out_train$univariate){

    depths <- sapply(out_train$theta, \(theta) get_depth(X, theta, out_train$distr))

  }else{

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
  }

  class <- apply(depths, 1, which.max)

  return(list("class" = class,
              "depths" = depths))
}

#' @rdname dqclass_train
#' @param train,test Matrices for training and testing data sets.
#'        The last column must contain the response class coded via integer numbers.
#' @export
dq_train_test <- function(train, test, n_dir = 500,
                          S,
                          distr = "fgld",
                          weighted = FALSE, sphered = TRUE){

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
