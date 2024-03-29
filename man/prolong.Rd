% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prolong.R
\name{prolong}
\alias{prolong}
\title{Fit prolong Model}
\usage{
prolong(
  x,
  y,
  lambda1 = NULL,
  lambda2 = NULL,
  lambdar = NULL,
  groups = TRUE,
  foldids = NULL,
  optimvals = c(1, 0.01)
)
}
\arguments{
\item{x}{Input covariate array, with n rows, p columns, and t slices}

\item{y}{Input response matrix, with n rows and t columns}

\item{lambda1}{Lasso/group lasso parameter, if left as \code{NULL} this parameter
will be chosen via a cross-validation that keeps each subject's observations
across time points together. It is recommended to save the lambda2 and
lambdar values from the first run of prolong for future runs, since the
optimization step takes the longest}

\item{lambda2}{Laplacian penalty parameter, if left as \code{NULL} will be chosen along
with lambdar via MLE, see supplementary material in
\insertCite{spatial-connectivity}{prolong}}

\item{lambdar}{Nuisance parameter, if left as \code{NULL} will be chosen along with
lambdar via MLE. Added to the the diagonal elements of the laplacian matrix
to get the invertibility required for the MLE}

\item{groups}{Optional pre-specified groups. If \code{NULL} or \code{FALSE}, lasso will be
used. If left as \code{TRUE}, then there will be p groups each containing the
observations across time points}

\item{foldids}{Optional pre-specified foldids for the cv. Should be of length
n. If left as \code{NULL}, subjects will be automatically split into 5 folds}

\item{optimvals}{Initial values of lambda2 and lambdar to be used by
\code{optim()} for the MLE of lambda2 and lambdar}
}
\value{
A list object with S3 class \code{prolong}
\tabular{ll}{
\code{beta} \tab \code{p*t*t-1/2 *length(lambda1)} matrix of coefficients. Currently, \code{lambda1} is chosen via cross-validation and is just a single value \cr
\code{selected} \tab the names of the variables with at least one non-zero coefficient \cr
\code{df} \tab number of non-zero coefficients (out of \code{p*t(t-1)/2}, not \code{p})\cr
\code{dim} \tab dimension of full coefficient matrix over \code{lambda1} values. Currently, \code{lambda1} is chosen via cross-validation and is just a single value \cr
\code{lambda1} \tab sequence of \code{lambda1} values used in final \code{gglasso}/\code{glmnet} call \cr
\code{lambda2} \tab \code{lambda2} value either passed by user or chosen via MLE, parameter of interest for the network penalty \cr
\code{lambdar} \tab \code{lambdar} value either passed by user or chosen via MLE, nuisance parameter needed to estimate \code{lambda2} via MLE \cr
\code{npasses} \tab total number of iterations summed over \code{lambda1} values for final \code{gglasso}/\code{glmnet} call \cr
\code{jerr} \tab error flag for final \code{gglasso}/\code{glmnet} call, 0 if no error \cr
\code{group} \tab vector of consecutive integers describing the grouping of coefficients \cr
\code{call} \tab the \code{gglasso}/\code{glmnet} call that produced this object \cr
}
}
\description{
\code{prolong()} is designed to take an \eqn{n \times t} outcome matrix and an
\eqn{n \times p \times t} array of covariates and automatically run the
entire processing and model fitting process, returning a list of coefficients
with corresponding variable names.
}
\details{
First, the data is reshaped into the first-differenced vectorized
\eqn{n(t-1)} length Y and first-differenced matricized \eqn{n(t-1) \times
pt(t-1)/2} X with a corresponding dependence matrix.
Next, hyperparameters for the graph laplacian-based network
penalty are found via MLE and the parameter for lasso/group lasso is found
via a careful implementation of cross-validation.
Lastly, a group lasso + laplacian or lasso + laplacian model is implemented,
and its bias-adjusted coefficients are returned.
}
\examples{
\dontrun{
promod <- prolong(Ymatrix, Xarray)
promod$beta
promod$selected

promod <- prolong(Ymatrix, Xarray, lambda2 = .001, lambdar = 10, groups = FALSE)
promod$beta
promod$selected
}

}
\references{
\insertRef{prolong}{prolong}
\insertAllCited{}
\insertRef{gglasso}{prolong}
\insertRef{group-lasso}{prolong}
\insertRef{glmnet}{prolong}
\insertRef{lasso}{prolong}
}
