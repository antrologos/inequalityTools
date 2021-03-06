% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decomp_dattRavallion.R
\name{decomp_dattRavallion}
\alias{decomp_dattRavallion}
\title{Makes a Lorenz Curve Function out of a income vector}
\usage{
decomp_dattRavallion(z, mu0, mu1, lorenz0, lorenz1)
}
\arguments{
\item{z}{A numeric income value representing the poverty line}

\item{mu0}{A numeric value for the mean of the group/time 0}

\item{mu1}{A numeric value for the mean of the group/time 1}

\item{lorenz0}{A vector valued function which takes cumulative probabilities as input and returns points at the Lorenz Curve for group/time 0}

\item{lorenz1}{A vector valued function which takes cumulative probabilities as input and returns points at the Lorenz Curve for group/time 1}
}
\value{
Returns a 1x5 tibble with:
    - p0_t0: Poverty rate (P0) for group/time 0
    - p0_t1: Poverty rate (P0) for group/time 1
    - growth: The difference in poverty rates due to the Growth Effect
    - redistr: The difference in poverty rates due to the Redistribution Effect
    - povertyDiff: The total difference in poverty rates between groups/times 0 and 1 (Growth Effect + Redistribution Effect).
}
\description{
Makes a Lorenz Curve Function out of a income vector
}
