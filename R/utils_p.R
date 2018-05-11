#' Test if an object exist
#' @param object the object name, a variable
#' @return TRUE or FALSE
#' @export
testObject <- function(object)
{
   exists(as.character(substitute(object)))
}

###currenly I do not have the ability to write such a meta-function (not enough meta-programming skills). QUIT!----> solved; it is not that difficult
#' a wrapper for file overwrite
#' 
#' typical usage is to pass any function that has a file argument, i.e. save(), pdf() and so on
#'
#' this allows us to decide whether to overwrite existing file in a directory.
#' 
#' @param FUN a function to pass
#' @param file a file including path
#' @param overwrite TRUE or FALSE
#' @param ... additional parameters for FUN
#' @return Ord indicator if there is overwriting happen so that further code can be conditionally run, i.e. plot
#' @examples
#' pdfOR(file='test.pdf', overwrite=TRUE)
#' plot(rnorm(100))
#' dev.off()
#' overwritef(pdf, file='test.pdf', overwrite=TRUE)
#' plot(1:3)
#' dev.off()
#' @export
overwritef <- function(FUN, file, overwrite=TRUE, ...){
	#browser()
	if(overwrite){
		FUN(..., file=file)
		ORd <- TRUE
	} else {
		if(!file.exists(file)) {
		FUN(..., file=file)
		ORd <- TRUE
		} else {
			ORd <- FALSE
		}
	}
	return(ORd)
}
#' pdf() with overwrite option
#' @param return default to have return 
#' @export
pdfOR <- function(file, overwrite=TRUE, ..., return=TRUE){
	Ord <- overwritef(FUN=pdf, file=file, overwrite=overwrite, ...)
	if(return) return(Ord)
}
#' save() with overwrite option
#' @param return default to no return 
#' @export
#saveOR <- function(file, overwrite=TRUE, return=FALSE, ...){ # this is wrong: return will be the first argument, i.e. cor_RPPA_RPPA! make sure we get the right position for ...
saveOR <- function(file, overwrite=TRUE, ..., return=FALSE){
	Ord <- overwritef(FUN=save, file=file, overwrite=overwrite, ...)
	#browser()
	if(return==TRUE) return(Ord)
}
#saveOR(cor_RPPA_RPPA,cor_RPPA_MIR, cor_RPPA_MET, cor_RNA_MIR, cor_RNA_RNA, cor_RNA_MET, file=ff, overwrite=overwrite)

#' write csv with overwrite option
#' @param return default to no return 
#' @export
write.csvOR <- function(file, overwrite=TRUE, ..., return=FALSE){
	Ord <- overwritef(FUN=write.csv, file=file, overwrite=overwrite, ...)
	if(return) return(Ord)
}

# this enables interface with overwritef
WriteXLS_file <- function(x, file, ...) {
	WriteXLS(x=x, ExcelFileName=file, ...)
}
# this correct unaccepted chars in list name when writing a list as xls
# usually the list is from rowttestMatMat: drugs*statistics*genes; genes is the list dimension
# by default, genes are into different sheets and each sheet is drug*statistics 
# when there are too many genes (sheets), the output file is blank for unknown reason
# in this case, we might want to save different sheets as drugs 
#
WriteXLS_list <- function(x, file, ...){
	#browser()
	xList <- get(x, envir=.GlobalEnv)
	onames <- names(xList)
	# sheet names does not allow following chars
	# []:*?/\
	#browser()
	names(xList) <- str_replace_all(onames, '\\[|\\]|:|\\*|\\?|/|\\\\', '.')
	# use local environment to search the variable
	WriteXLS_file('xList', file=file, envir=environment(), ...)
}
WriteXLS_array <- function(x, file, margin=1, transpose=TRUE, ...){
	#browser()
	xVal <- get(x, envir=.GlobalEnv)
	tox <- function(x, margin){
		if(margin %in% c(1, 2)){
			if(transpose){
					res <- t(x) # transpose if in margin 1 or 2
			} else {
				res <- x
			}
		} else {
			res <- x
		}
		data.frame(res)
	}
	xList <- alply(xVal, margin, function(x) tox(x, margin))
	#browser()
	onames <- dimnames(xVal)[[margin]]
	# sheet names does not allow following chars
	# []:*?/\
	names(xList) <- str_replace_all(onames, '\\[|\\]|:|\\*|\\?|/|\\\\', '.')
	#attr(xList, 'split_type') <- NULL
	#attr(xList, 'split_labels') <- NULL
	# use local environment to search the variable
	WriteXLS_file(x='xList', file=file, envir=environment(), ...)
}
#WriteXLS_array(x=c("asso_RPPA_array"), file=file.path(DirOutput, "asso_RPPA.xls"), row.names=TRUE)

#
# this is not working: get(x) not found!
#
#' write xls (using WriteXLS) with overwrite option from an array
#'  default slice in margin 1, which means different rows are stored into different sheets
#' @export
WriteXLS_arrayOR  <- function(file, overwrite=TRUE, transpose=TRUE, ..., return=FALSE){
	#browser()
	Ord <- overwritef(FUN=WriteXLS_array, file=file, overwrite=overwrite, transpose=transpose, ...)
	if(return) return(Ord)
}
#WriteXLS_arrayOR(x=c("asso_RPPA_array"), file=file.path(DirOutput, "asso_RPPA.xls"), row.names=TRUE, overwrite=overwrite)


#' write xls (using WriteXLS) with overwrite option
#' @param return default to no return 
#' @export
WriteXLSOR <- function(file, overwrite=TRUE, ..., return=FALSE){
	#browser()
	Ord <- overwritef(FUN=WriteXLS_file, file=file, overwrite=overwrite, ...)
	if(return) return(Ord)
}
#' write xls (using WriteXLS) with overwrite option from a list
#' @param return default to no return 
#' @export
WriteXLS_listOR  <- function(file, overwrite=TRUE, ..., return=FALSE){
	#browser()
	Ord <- overwritef(FUN=WriteXLS_list, file=file, overwrite=overwrite, ...)
	if(return) return(Ord)
}
#WriteXLS_listOR(c("asso_RPPA"), file="asso_RPPA.xls", overwrite=TRUE, row.names=TRUE)
#WriteXLS_list(c("asso_RPPA"), file="asso_RPPA.xls", row.names=TRUE)

#' generate html from Rmd with docco-linear style as in \url{http://cran.r-project.org/web/packages/knitr/vignettes/docco-linear.html}
#'
#' This generates kind of good looking html
#' 
#' Notice that envir is set to be globalenv() so that after this evaluation, all objects are still in the console.
#' Otherwise, some child script cannot run since it is not in the global.
#' @param input input Rmd file
#' @param ... additional parameter to knit2html (actually markdownToHTML)
#' @export
doccoRmd <- function(input, ...){
	knit2html(input=input, ..., stylesheet = system.file("misc", 
        "docco-classic.css", package = "knitr"), template = system.file("misc", 
        "docco-template.html", package = "knitr"), options = c('toc', markdown::markdownHTMLOptions(TRUE)), envir = globalenv())
}
#' generate html with plain style
#'
#' 
#' Notice that envir is set to be globalenv() so that after this evaluation, all objects are still in the console.
#' Otherwise, some child script cannot run since it is not in the global.
#' @param input input Rmd file
#' @param ... additional parameter to knit2html (actually markdownToHTML)
#' @export
plainRmd <- function(input, ...){
	knit2html(input=input, ..., options = c('toc', markdown::markdownHTMLOptions(TRUE)), envir = globalenv())
}