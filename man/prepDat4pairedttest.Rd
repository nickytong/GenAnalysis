% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TestAndCorr_p.R
\name{prepDat4pairedttest}
\alias{prepDat4pairedttest}
\title{given a vector of data (continuous) and class (binary) and sampleID, prepared matched data by ordering for paired t test}
\usage{
prepDat4pairedttest(data, class, sampleID = NULL, levels = NULL)
}
\arguments{
\item{data}{continuous response vector}

\item{class}{binary treatment vector}

\item{sampleID}{a vector of sample ID which can be used to automatically track the paired correspondance.}
}
\value{
a list
}
\description{
given a vector of data (continuous) and class (binary) and sampleID, prepared matched data by ordering for paired t test
}

