% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fgld-ls-estimation.R,
%   R/fgld-std-ls-estimation.R
\name{get_X_fgld}
\alias{get_X_fgld}
\alias{fit_fgld_ls}
\alias{fit_fgld_ls_obj}
\alias{fit_fgld_ls.quadprog}
\alias{get_X_fgld_std}
\alias{fit_fgld_std_ls}
\alias{fit_fgldid_ls}
\title{Least squares estimation for the \emph{fgld}}
\usage{
get_X_fgld(n)

fit_fgld_ls(y)

fit_fgld_ls_obj(y)

fit_fgld_ls.quadprog(y)

get_X_fgld_std(n)

fit_fgld_std_ls(y)

fit_fgldid_ls(y)
}
\arguments{
\item{n}{Sample size.}

\item{y}{Data from which to estimate the parameters.}
}
\value{
\code{fit_fgld_ls()} outputs a vector of length 4 containg the estimates of the parameters.

\code{get_X_fgld()} outputs a \eqn{n \times 4} matrix where the i-th row contains the
coefficients for expected value of the i-th order statistic of the fgld.
}
\description{
\code{fit_fgld_ls()} carries out the least squares estimation of the
\emph{fgld} (flattened generalized logistic distribution) with a linear
reparametrization.

\code{get_X_fgld()} returns the design matrix needed for the least squares estimation. The design matrix contains the coefficients of the
expected order statistics.
}
\examples{
x <- rnorm(5)
get_X_fgld(5)
fit_fgld_ls(x)

}
