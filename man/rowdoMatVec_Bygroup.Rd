% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TestAndCorr_p.R
\name{rowdoMatVec_Bygroup}
\alias{rowdoMatVec_Bygroup}
\title{the driver of grouped version calculations, enpowering rowcorMatVec, rowttestMatVec}
\usage{
rowdoMatVec_Bygroup(mat, vec, group, fn, addAll = FALSE, addAllName = "all",
  parallel = FALSE, core = 1, returnArray = TRUE, returnEarly = FALSE,
  ...)
}
\arguments{
\item{mat}{a matrix, i.e. gene expression}

\item{vec}{a vector of phenotype, continuous variable}

\item{group}{a vector, categorical variable}

\item{fn}{function to apply, e.g. rowttestMatVec, rowcorMatVec}

\item{addAll}{whether to calculate correlation by using all data ignoring group; this is the overall correlation}

\item{parallel}{deprecated; parallel is automatically enabled when core>1}

\item{core}{number of cores to use; when core>1, parallel execution is enabled}

\item{returnEarly}{return early for the fitted list; useful for rowanovamatVec since array is not feasible due to different dimensions of fits}

\item{...}{additional parameters to fn}
}
\value{
a 3D array
}
\description{
This calculated result between rows of mat and vec but at each group at a time as specified by the group vector
It should be more flexible than rowcorMatMat which requires the same genomic mat; instead, here different group of analysis can
have different sample size in genomics
}

