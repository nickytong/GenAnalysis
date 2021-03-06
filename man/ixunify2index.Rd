% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatting_p.R
\name{ixunify2index}
\alias{ixunify2index}
\title{unify index/name indexing into global index given global name}
\usage{
ixunify2index(ix, name_global)
}
\arguments{
\item{ix}{either integer index or string vector}

\item{name_global}{global name index}
}
\value{
integer indexing (global reference)
}
\description{
usually people use integer index or name index. it is important to unify this.
this function maps either index or name into global integer index 
the algorithm is that when input is integer, it is not changed; when input is string, it is matched into global name
to identify the integer index
}

