% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotfuns.R
\name{delta_network}
\alias{delta_network}
\title{Experimental Network Plots for Delta-Scale Pairwise (Partial) Correlations}
\usage{
delta_network(
  x,
  timediff = "2-1",
  partial = TRUE,
  corr_thresh = 0.75,
  interactive = TRUE,
  method = c("pearson", "kendall", "spearman")
)
}
\arguments{
\item{x}{Input covariate array, with n rows, p columns, and t slices}

\item{timediff}{Pair of time points to use in heatmap. Should be in format \code{'t2-t1'}}

\item{partial}{If \code{TRUE}, partial correlations using the \code{ppcor} package will be used. If \code{FALSE}, correlations will be used}

\item{corr_thresh}{Correlation value threshold for edge inclusion in the graph}

\item{interactive}{If \code{TRUE}, an interactive \code{visNetwork} plot will open in the Viewer panel. If \code{FALSE}, a static \code{igraph} plot will be produced}

\item{method}{A character string indicating which correlation coefficient is to be computed, with 'pearson' as default}
}
\value{
Either an interactive \code{visNetwork} or static \code{igraph} plot
}
\description{
Experimental Network Plots for Delta-Scale Pairwise (Partial) Correlations
}
\examples{
\dontrun{
delta_network(Xarray)
delta_network(Xarray, timediff = '4-2', corr_thresh = .9, method = 'spearman')
}
}
