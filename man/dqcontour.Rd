% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dqdepth.R
\name{dqcontour}
\alias{dqcontour}
\title{Integrated rank-weighted depth contours}
\usage{
dqcontour(
  u,
  X,
  S,
  n_dir = 1000,
  distr = "normal",
  sphered = FALSE,
  weighted = FALSE
)
}
\arguments{
\item{u}{Level of the contour that is required.}

\item{X}{Data matrix with respect to which the depth is to be computed, with observations in rows and variables in columns.}

\item{S}{A matrix containing the directions in the columns.}

\item{n_dir}{Number of directions to consider. They are equally spaced in two dimensions and uniformly sampled for higher dimensions.}

\item{distr}{The distribution that is fit to each univariate directional projection.}

\item{sphered}{A logical indicating whether the depth should be computed on the sphered data.}

\item{weighted}{A logical indicating whether weights for the directions should be computed.}
}
\value{
A matrix containing as rows the points that form the contour.
}
\description{
Computes an approximation of the contour level, giving a set points with a predifined value of the depth.
}
\examples{
library(mvtnorm)
set.seed(2)
X <- rmvnorm(n = 200, mean = c(10, 2), sigma = matrix(c(1, 0.9, 0.9, 2), 2, 2))
contour <- dqcontour(u = 0.9, X = X)

plot(X, cex = 0.5)
points(rbind(contour, contour[1,]), col = 2, lwd = 3, type = "l")

}
