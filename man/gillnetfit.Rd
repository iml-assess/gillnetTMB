% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit.R
\name{gillnetfit}
\alias{gillnetfit}
\title{Fit gillnet selectivity model}
\usage{
gillnetfit(x, par)
}
\arguments{
\item{x}{input data (list). mesh needs to be in increasing order.}

\item{par}{start parameters (list)}
}
\description{
Fit gillnet selectivity model
}
\details{
Fits gillnet selecitivty model (obs = N * sel)
}
\examples{
fit <- gillnetfit(data, defpar(par))
}
