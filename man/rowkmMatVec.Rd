% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TestAndCorr_p.R
\name{rowkmMatVec}
\alias{rowkmMatVec}
\title{Perform log rank test between rows of mat (categorical) and a vector of OS and OScensoring}
\usage{
rowkmMatVec(mat, OS, OScensoring, meta = list(tag = NULL, row = "group",
  response = "Overal survival"), more = TRUE, levels = NULL,
  parallel = FALSE, core = 1, addData = TRUE, ...)
}
\arguments{
\item{mat}{a matrix with rows as genes}

\item{OS}{a vector giving the survival time. This can be for example overall survival, PFS and others}

\item{OScensoring}{a vector giving the censor status. 1 is event and 0 is censored.}

\item{meta}{meta information, a list. This can specify what type of analysis (i.e. KRAS vs NF1 mutation);
row, which represents what row variables are (rows of mat); response represents what does OS mean}

\item{more}{not implemented}

\item{levels}{levels specified for each row of mat}

\item{core}{number of cores to use; this override parallel=TRUE if core=1 (no parallel)}

\item{addData}{logical indicating if data (mat) should be attached as an attribute (may be used for scatter plot)}

\item{...}{additional parameters to plotKM}
}
\value{
a rowkmMatVec class; meta information is stored as attributes
}
\description{
Notice that ggplot2 figures are attached; this may takes a long time to cache the objects in knitr
}

