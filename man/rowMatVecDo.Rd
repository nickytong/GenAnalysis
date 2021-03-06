% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TestAndCorr_p.R
\name{rowMatVecDo}
\alias{rowMatVecDo}
\title{perform operations between rows of a mat and a vector}
\usage{
rowMatVecDo(mat, vec, fn, ...)
}
\arguments{
\item{mat}{expression matrix}

\item{vec}{a vector of response}

\item{fn}{a function applied to the data frame with columns gene, expr, vec}

\item{...}{additional parameters to fn}
}
\value{
a data frame
}
\description{
This should be a general approach for rowcorMatVec rowlmMatVec rowranktestMatVec rowttestMatVec
In principle, we want to apply arbitrary functions and return a data frame.
Currently this is implemented through ddply; later we might use dply to speed up computation
plyr can parallel; but dplyr cannot.
}

