% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decomp_meanIncomeBySources.R
\name{dynamicDecomp_meanIncomeBySources}
\alias{dynamicDecomp_meanIncomeBySources}
\title{Dynamic decomposition of the mean income by sources}
\usage{
dynamicDecomp_meanIncomeBySources(
  decomp_meanIncomeBySources0,
  decomp_meanIncomeBySources1
)
}
\arguments{
\item{decomp_meanIncomeBySources0}{A data.frame with with results of a static decomposition the mean income (produced by decomp_meanIncomeBySources) for group/time 0}

\item{decomp_meanIncomeBySources1}{A data.frame with with results of a static decomposition the mean income (produced by decomp_meanIncomeBySources) for group/time 1}
}
\value{
Returns a tibble with the following coluns:
    - Income Source: name of the income source
    - mean_effect: change in mean total income due to the change in the mean of each income source
    - recipients_effect: change in mean total income due to the change in proportion of recipients of each income source
    - TOTAL: total change in mean total income due to each income source (mean_effect + recipients_effect)
}
\description{
Dynamic decomposition of the mean income by sources
}
