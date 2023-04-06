# library(dqclass)
#
# data("iris")
# out_train <- dqclass_train(iris[, -5], iris[, 5], predicted = TRUE)
# table(out_train$class, iris[,5])
#
#
# library(ggplot2)
# library(mvtnorm)
# set.seed(1)
# X <- rmvnorm(n = 200, mean = c(10, 2), sigma = matrix(c(1, 0.9, 0.9, 2), 2, 2))
#
# depths <- dqdepth(X = X, n_dir = 1, weighted = FALSE)
#
# data.frame("x1" = X[,1], "x2" = X[,2], "d" = depths) |>
#   ggplot(aes(x1, x2, col = d))+
#   geom_point()+
#   scale_color_viridis_c()
#
#
# n <- c(100, 50)
# x <- c(rnorm(n[1], 10), rnorm(n[2], -10))
# y <- rep(c(1,2), times = n)
# plot(x, col = y)
# out_train <- dqclass_train(x, y, predicted = TRUE)
# table(out_train$class, y)
