
test_that("sphering works", {
  n <- 100
  p <- 20
  X <- rnorm(n * p) |> matrix(n, p)
  W <- sphere_zca(X)
  expect_equal(var(X %*% W),
               diag(1, p))
})
