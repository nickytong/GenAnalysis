% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qplot_survival_p.R
\name{plotKM}
\alias{plotKM}
\title{flexible KM curve wrapper}
\usage{
plotKM(x, legendTitle, tag, ...)
}
\arguments{
\item{x}{the stratification variable, either numeric or factor or character}

\item{...}{other parameters passed to either ggKMcat or ggKM}
}
\description{
plot KM curve for a cont/cat x.
}
\details{
the program detects x as categorical if it is a factor or it is a character variable.
}

