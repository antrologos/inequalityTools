% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantile.R
\name{make_quantileFunc}
\alias{make_quantileFunc}
\title{Makes a Lorenz Curve Function out of a income vector}
\usage{
make_quantileFunc(x, w = NULL)
}
\arguments{
\item{x}{A vector of incomes}

\item{w}{(optional) A vector of sample weights}
}
\value{
Returns a function which takes a vector of probabilities as inputs (p) and gives points at the Lorenz Curve as outputs
}
\description{
Makes a Lorenz Curve Function out of a income vector
}
