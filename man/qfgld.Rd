% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fgld-ls-estimation.R
\name{qfgld}
\alias{qfgld}
\alias{qfgld_deriv}
\alias{pfgld}
\alias{dfgld}
\title{The flattened generalized logistic distribution (\emph{fgld})}
\usage{
qfgld(u, theta)

qfgld_deriv(u, theta)

pfgld(q, theta)

dfgld(x, theta, log = FALSE)
}
\arguments{
\item{u}{Vector of probabilities.}

\item{theta}{Vector of parameters (of length 4) for the \emph{fgld}.}

\item{x, q}{Vector of quantiles.}

\item{log}{Logical indicating whether the density should be returned on the log scale.}
}
\value{
\code{qfgld} gives the quantile function, \code{pfgld} gives the distribution functions.
The length of the output is the same as that of the input \code{x/u}.
}
\description{
Quantile and cumulative distribution functions for the \emph{fgld}.
}
\examples{
theta <- c(1, 2, 0, 0)
u <- seq(0.1, 0.9, by = 0.1)
qfgld(u, theta)

}
