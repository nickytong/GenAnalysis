% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TestAndCorr_p.R
\name{rowttestMatVec}
\alias{rowttestMatVec}
\title{Perform t test between rows of mat and a vector; Notice paired t test also works by setting paired=TRUE}
\usage{
rowttestMatVec(mat, vec, matIsX = TRUE, meta = list(tag = NULL, row =
  "expression", response = "response"), minN = 2, more = TRUE,
  var.equal = TRUE, levels = NULL, parallel = FALSE, base = 2,
  core = 1, addData = TRUE, symbol = NULL, sepSymProbe = ", ",
  paired = FALSE, FDR = NULL, pcutoff = 0.05, alphas = c(0.01, 0.05,
  0.1, 0.2, 0.3, 0.4, 0.5), ...)
}
\arguments{
\item{mat}{a matrix with rows as genes}

\item{vec}{a vector giving the group info. Internally, vec is converted to a vector by ttestP for computation}

\item{matIsX}{set matIsX=FALSE if each row of mat is binary (binary response vector)
  while set matIsX=TRUE if each row of mat is continuous. By default, mat rows are continuous (i.e. expr) and vec is binary; but for MUT data where mat is binary and the response might be continuous. In this case, we need to 
specify matIsX=FALSE which is passed to switchXY as TRUE in t test. By default, rows in mat is not response (matIsX=TRUE) and hence, switchXY is disabled.}

\item{meta}{meta information, a list. This can specify what type of analysis (i.e. KRAS vs NF1 mutation);
row, which represents what row variables are (rows of mat); response represents what does vec mean}

\item{base}{base parameter passed to ttestP}

\item{core}{number of cores to use; this override parallel=TRUE if core=1 (no parallel)}

\item{addData}{logical indicating if data (mat) should be attached as an attribute (may be used for scatter plot)}

\item{symbol}{in case this is probe level data, user can specify a vector of symbols matching}

\item{sepSymProbe}{separator for symbol-probe. That is, symbol+sepSymProbe+probe(rowname) is the new rowname}

\item{...}{additional parameters to ttestP}

\item{whether}{use paired t test (assumes data is special being paired)}
}
\value{
a rowttestMatVec class; meta information is stored as attributes
}
\description{
Perform t test between rows of mat and a vector; Notice paired t test also works by setting paired=TRUE
}

