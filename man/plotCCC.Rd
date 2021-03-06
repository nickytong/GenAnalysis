% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myplots_p.R
\name{plotCCC}
\alias{plotCCC}
\title{Calculates Concordance Correlation Coefficient (CCC) to access reproducibility}
\usage{
plotCCC(x, y, labels = NULL, cols = NULL, cexlab = 0.8,
  colab = "purple", poslab = 1, xlab = NA, ylab = NA, tag = "",
  col_legend = "red", pch = 20, main = NA, cex = 2, plot = TRUE,
  ylim = NA, xlim = NA, maskBeyondLim = FALSE, plotOutlier = FALSE,
  sampleName = NA, alpha_outlier = 0.01, show = c("identity", "h0", "v0"))
}
\arguments{
\item{x}{x vector}

\item{y}{y vector}

\item{cols}{color vector for the points}

\item{xlab}{xlab}

\item{ylab}{ylab}

\item{tag}{add a tag}

\item{pch}{pch of the points}

\item{main}{main title}

\item{cex}{for the dots}

\item{plot}{whether draw a figure}

\item{maskBeyondLim}{whether mask values beyond xlim (for x) or ylim (for y) when xlim and ylim is specified. Default is FALSE so that even xlim ylim specified, CCC will not be affected}

\item{plotOutlier}{logical whether add an outlier plot comparing Diff vs Mean plot and related statistics}

\item{sampleName}{symbols of the text to be shown in outlier plot}

\item{alpha_outlier}{alpha to call outliers based on the observed differences assuming from Normal(mean_diff, sd_diff)}

\item{cex.main}{cex for main title}

\item{cex.legend}{cex for legend text}
}
\value{
a vector of c('ccc', 's_shift', 'l_shift', 'ccc_lo', 'ccc_hi', 'Cb', 'corr'). s_shift is scale shift; l_shift is location shift;
		ccc_lo and ccc_hi represent 95% confidence interval for CCC. Cb for the bias correction term satisfying CCC=corr*Cb where corr is the
			Pearson correlation
}
\description{
This function calculates and plots concordance for paired data.
}
\details{
This function calls the epi.ccc() function from epiR package for the calculation.

infinite values in x and y are masked as NA so as to compute CCC and corr (otherwise, result becomes nan).
}

