% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fgld-ls-estimation.R
\name{origin2lin}
\alias{origin2lin}
\alias{lin2origin}
\title{From original to linear parametrization of the \emph{fgld} and viceversa}
\usage{
origin2lin(alpha)

lin2origin(theta)
}
\arguments{
\item{alpha}{Parameters in the original parametrization.}

\item{theta}{Parameters in the linear parametrization.}
}
\value{
A vector containing parameter values for the other parametrization.
}
\description{
From original to linear parametrization of the \emph{fgld} and viceversa
}
\examples{
alpha <- c(1, 2, 0.9, 2)
(theta <- origin2lin(alpha))
lin2origin(theta)

}
