% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatting_p.R
\name{swap}
\alias{swap}
\title{swap two variables in a given environment}
\usage{
swap(x, y, envir = parent.frame())
}
\arguments{
\item{x}{variable. in theory, any R object}

\item{y}{variable. in theory, any R object}
}
\description{
swap two variables in a given environment
}
\examples{
a <- 3
b <- 4
swap(a, b)
}

