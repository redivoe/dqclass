#' @rdname fit_fgld_ls
#' @export
get_X_fgld_std <-function(n){
  x2 <- digamma(1:n) - digamma(n + 1)

  X <- cbind((1:n) / (n + 1),
             x2,
             rev(-x2))
  return(X)
}

#' @rdname fit_fgld_ls
#' @export
fit_fgld_std_ls<-function(y){
  n <- length(y)
  X <- get_X_fgld_std(n)
  y <- sort(y)

  A <- cbind(c(0, 1, 1), diag(rep(1, 3)))
  b <- c(1, rep(0, 3))

  out_solveQP <- quadprog::solve.QP(
    Dmat = crossprod(X, X),
    dvec = crossprod(X, y),
    Amat = A,
    bvec = b
  )

  return(list("theta" = out_solveQP$solution,
              "obj" = out_solveQP$value))
}
