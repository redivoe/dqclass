
get_X_fgld<-function(n){
  x2 <- digamma(1:n) - digamma(n + 1)

  X <- cbind(1,
             (1:n) / (n + 1),
             x2,
             rev(-x2))
  return(X)
}

fit_fgld_ls<-function(y){
  n <- length(y)

  X <- get_X_fgld(n)

  y <- sort(y)

  A <- sapply(seq(1e-10, 1 - 1e-10, len = 1e3),
              \(u) c(0, 1, 1 / u, 1 / (1 - u)))

  out_solveQP <- quadprog::solve.QP(
    Dmat = crossprod(X, X),
    dvec = crossprod(X, y),
    Amat = A,
    bvec = rep(1e-3, 1e3)
  )

  return(out_solveQP$solution)
}

invQ <- function(x, theta, Q) {
  sapply(
    x,
    ~ tryCatch(
      stats::uniroot(function(u) Q(u, theta) - .x,
                     interval = c(1 - 1e-10, 1e-10)
      )$root,
      error = function(cond) {
        ifelse(.x > theta[1], 1 - 1e-10, 1e-10)
      }
    )
  )
}

qfgld <- function(u, theta) {
  theta[1] + theta[2] * u + theta[3] * log(u) - theta[4] * log(1 - u)
}

qfgld_deriv <- function(u, theta) {
  ifelse(u >= 1 - 1e-10 | u <= 1e-10,
         1e10,
         theta[2] + theta[3] / u + theta[4] / (1 - u))
}

pfgld <- function(x, theta) {
  invQ(x, theta, Q = qfgld)
}
