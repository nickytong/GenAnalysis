% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aheatmat.R
\name{aheatmatcor}
\alias{aheatmatcor}
\title{plot correlation heatmap for rows from a matrix}
\usage{
aheatmatcor(mat, method = "pearson", CexRow = NULL, CexCol = NULL,
  oma = c(2, 1, 4, 3) + 0.1, dendrogram = "row", main = "", ykey = 1.05,
  xkey = 0.05, stretch = 0.08, colFn = bluered, ncol = 128,
  return = FALSE)
}
\arguments{
\item{mat}{a matrix where correlation will be computed from its rows}

\item{method}{correlation method}

\item{main}{main title}

\item{ykey}{y value for correlation key}
}
\description{
plot correlation heatmap for rows from a matrix
}

