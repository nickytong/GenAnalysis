% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TestAndCorr_p.R
\name{lm_xy}
\alias{lm_xy}
\title{lm for cont x and y}
\usage{
lm_xy(x, y, collapse = TRUE)
}
\arguments{
\item{x}{continuous vector}

\item{y}{continuous vector}

\item{collapse}{when vec is non-numeric which will be treated as factor, there are multiple coefs for different levels. We thus need to
summarize them into a scaler (by mean). if collapse=FALSE, all coefs will be returned.}
}
\value{
a named vector
}
\description{
when x is categorical, lm p value is equivalent to ANOVA p value. Coef becomes the mean; Rsquare?
}

