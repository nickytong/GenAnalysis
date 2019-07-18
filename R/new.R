showDuplicateInfo <- function(vec){
	if(anyDuplicated(vec)==0){
		cat('No duplicates are observed\n')
		return(NULL)
	}
	#ind <- which(duplicated(vec)==TRUE)# there is a bug: when 3 times, it will give twice of the index!
	ind <- match(unique(vec[duplicated2(vec)]), vec) # this fix it!
	#browser()
	for(ii in ind){
		val <- vec[ii]
		nDup <- sum(vec==val)
		indDup <- which(vec==val)
		cat(sprintf('%s is duplicated in %d times at: %s\n', val, nDup, paste(indDup, collapse=', ')))
	}
}

#' return an indicator if any element is duplicated (first appearance is flagged duplicated)
#'
#' the duplicated function in R is odd when there are 3 or more replicated instances. This is due to the fromLast option.
#' now we make this right so that any duplicated element will be indicated as TRUE
#' @param x same as x in duplicated. a vector or a data frame or an array or NULL.
#' @return a logical vector
#' @export
#' @examples duplicated(c(1, 2, 2, 2, 3))
#' @examples duplicated2(c(1, 2, 2, 2, 3))
duplicated2 <- function(x) duplicated(x) | duplicated(x, fromLast=TRUE)

