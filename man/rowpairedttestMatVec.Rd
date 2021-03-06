% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TestAndCorr_p.R
\name{rowpairedttestMatVec}
\alias{rowpairedttestMatVec}
\title{Paired t test between rows of a matrix and a vector}
\usage{
rowpairedttestMatVec(mat, vec, sampleID = NULL, matIsX = TRUE,
  meta = list(tag = NULL, row = "expression", response = "response"),
  minN = 2, more = TRUE, levels = NULL, addData = TRUE, core = 1,
  base = 2, ...)
}
\arguments{
\item{mat}{a matrix, i.e. gene expression}

\item{vec}{a vector of phenotype, categorical variable}

\item{meta}{meta information, a list. This can specify what type of analysis (i.e. RPPA correlated with EMT score);
row, which represents what row variables are (rows of mat); response represents what does vec mean}
}
\value{
an object belonging to data frame and a class
}
\description{
Paired t test between rows of a matrix and a vector
}

