% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatting_p.R
\name{augmentMatrixColumnByID}
\alias{augmentMatrixColumnByID}
\title{Augment a matrix by ID in the direction of columns}
\usage{
augmentMatrixColumnByID(mat, ID, computeNA = FALSE, pNAthr = 0.8,
  data.matrix = TRUE)
}
\arguments{
\item{mat}{a matrix to be augmented}

\item{ID}{IDs targeted for; should be simialr to colnames(mat)}

\item{computeNA}{whether to recompute columns with all NA}

\item{pNAthr}{percentage of NA to be called as NA columns}
}
\value{
an augMat object (also inherits matrix)
}
\description{
similar to fillVecByID, we might want to augment a matrix by column names specified by ID
}
\details{
It is possible to augment an already augmented matrix. In this case, we'd be carefull about filledColID and naColID
}

