% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_p.R
\name{plainRmd}
\alias{plainRmd}
\title{generate html with plain style}
\usage{
plainRmd(input, ...)
}
\arguments{
\item{input}{input Rmd file}

\item{...}{additional parameter to knit2html (actually markdownToHTML)}
}
\description{
Notice that envir is set to be globalenv() so that after this evaluation, all objects are still in the console.
Otherwise, some child script cannot run since it is not in the global.
}

