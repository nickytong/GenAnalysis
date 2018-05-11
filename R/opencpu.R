# cd /home/ptong1/Backup/GitHub/
#' build an openCPU app with roxygenize+build+install+opencpu$browse
#'
#' this function wrappers roxygenize(), build(), install(), and browse()
#'
#' Importance: need to navigate to the mater folder first; cannot work on full path!
#'
#' @param package package name, a string
#' @param reopen whether to open the browser
#' @export
buildapp <- function(package, reopen=TRUE){
	library(roxygen2)	
	library(devtools)
	roxygenize(package)
	build(package)
	install(package)
	library(opencpu)
	if(reopen)
		opencpu$browse(sprintf("/library/%s/www", package))
}
# cannot specify a full path: needs to navigate to the directory to get the things right
#buildapp('appdemo')


#' build an openCPU app with build+install+opencpu$browse
#'
#' this function wrappers build(), install(), and browse()
#'
#' Importance: need to navigate to the mater folder first; cannot work on full path!
#'
#' @param package package name, a string
#' @param reopen whether to open the browser
#' @export
buildapp0 <- function(package, reopen=TRUE){
	#library(roxygen2)	
	library(devtools)
	#roxygenize(package)
	build(package)
	install(package)
	library(opencpu)
	if(reopen)
		opencpu$browse(sprintf("/library/%s/www", package))
}
# cannot specify a full path: needs to navigate to the directory to get the things right
#buildapp0('appdemo')
