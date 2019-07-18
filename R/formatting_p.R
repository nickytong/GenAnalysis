####
#### to build
####
if(FALSE){
# cd /home/ptong1/Backup/Package/; R0
## notice: @param @return @exprt all required for a successful export!
library(devtools)
library(roxygen2)
#library(roxygen) # not working
roxygenize("GenAnalysis")

library(devtools)
build('GenAnalysis')
install('GenAnalysis')

load_all() # load all R code
load_all('/home/ptong1/Backup/Package/GenAnalysis') # load all R code

# reload in a working R session
detach("package:GenAnalysis", unload=TRUE)
library(GenAnalysis)

}

#ixunify2index(ix=as.factor(imTargets), name_global=rownames(rna))
#' map local index to global index given reference
#'
#' @param local local index based on ref
#' @param ref ref index based on global index
#' @return index based on global as reference
#' @export
ixlocal2global <- function(local, ref){
	ref[local]
}
#' map global index to local index given reference
#'
#' @param global global index
#' @param ref ref index based on global index
#' @return index based on local reference
#' @export
ixglobal2local <- function(global, ref){
	match(global, ref)
}

#' unify index/name indexing into global index given global name
#'
#' usually people use integer index or name index. it is important to unify this.
#' this function maps either index or name into global integer index 
#' the algorithm is that when input is integer, it is not changed; when input is string, it is matched into global name
#' to identify the integer index
#' @param ix either integer index or string vector
#' @param name_global global name index
#' @return integer indexing (global reference)
#' @export
ixunify2index <- function(ix, name_global){
	#browser()
	if(inherits(ix, 'character') | inherits(ix, 'factor')){
		# input is character or factor (possibly converted from character)
		res <- match(as.character(ix), name_global)
	} else {
		res <- ix
	}
	res
}
#' general getter function for attributes in an object
#'
#' @param obj object with attributes
#' @param what specify what attribute to get 
#' @return attributes extracted 
#' @export
Getter <- function(obj, what='avail'){
	#browser()
	if(what %in% names(attributes(obj))){
		res <- attributes(obj)[[what]]
	} else if (what =='avail'){
		res <- names(attributes(obj))
	} else {
		res <- NULL
	}
	res
}

#' there may be mroe elegant ways to do this
emptydfFromDF <- function(df, nrow=NULL){
	#browser()
	if(is.null(nrow)) nrow <- nrow(df)
	res <- foreach(i=1:nrow, .combine='rbind') %do% {
		zz <- df[1, ]
		zz[zz!=3.1415926] <- NA
		zz
		#browser()
	}
	# modified on 06/22/2016: colnames(res) <- df
	colnames(res) <- colnames(df)
	res
}
#emptydfFromDF(clin_v2)
#' Augment a data frame by ID in the direction of rows
#'
#' similar to fillVecByID, we might want to augment a matrix by column names specified by ID
#' We cannot directly using augmentMatrixColumnByID with transposed input since data frame has side effects converting all numeric to string
#'
#' @param mat a matrix or data frame to be augmented
#' @param ID IDs targeted for; should be simialr to rownames(mat)
#' @param computeNA whether to recompute columns with all NA
#' @param pNAthr percentage of NA to be called as NA columns
#' @return an augMat object (also inherits matrix)
#' @export
augmentMatrixRowByID <- function(mat, ID, computeNA=FALSE, pNAthr=0.8){
	#res <- data.frame(matrix("", nrow=length(ID), ncol=ncol(mat)))
	#res <- emptydfFromDF(mat, nrow=length(ID)) # this is to inheric attribute of each column so as nto to convert to numeric (especially for factors!)
	res <- data.frame(matrix(NA, ncol=ncol(mat), nrow=length(ID)))
	rownames(res) <- ID
	colnames(res) <- colnames(mat)
	IDfill <- intersect(ID, rownames(mat))
	if(length(IDfill)<1) stop('ID has no shared elements with rownames(mat)!')
	# mat <- data.matrix(mat) # otherwise res is NULL
	res[IDfill, ] <- mat[IDfill, ]
	if(!computeNA) {
		naColID <- setdiff(ID, IDfill) # this might be inaccurate
	} else {
		pNA <- apply(res, 1, percentNA)
		naColID <- ID[pNA>pNAthr]
		IDfill <- setdiff(ID, naColID)
	}	
	#browser()
	filledColID <- IDfill
	naColInd <- match(naColID, ID)
	filledColInd <- match(IDfill, ID)
	attr(res, 'filledColID') <- filledColID
	attr(res, 'filledColInd') <- filledColInd
	attr(res, 'naColID') <- naColID
	attr(res, 'naColInd') <- naColInd
	attr(res, 'isNAcol') <- ID %in% naColID
	#browser()
	# bug: attr(res, 'effectiveMat') <- mat[, filledColInd]
	attr(res, 'N') <- length(filledColInd) # effective N
	if(length(filledColInd) < 1000)
		attr(res, 'effectiveMat') <- res[filledColInd, ]
	#attr(res, 'class') <- c('augMat', 'matrix')
	#browser()
	class(res) <- append(class(res),"augMat")
	# to check class: inherits(res, 'augMat')
	res
}
#clin_v2_aug <- augmentMatrixRowByID(clin_v2, ID=ID_aug)




#' Augment a matrix by ID in the direction of columns
#'
#' similar to fillVecByID, we might want to augment a matrix by column names specified by ID
#' 
#' It is possible to augment an already augmented matrix. In this case, we'd be carefull about filledColID and naColID
#'
#' @param mat a matrix to be augmented
#' @param ID IDs targeted for; should be simialr to colnames(mat)
#' @param computeNA whether to recompute columns with all NA
#' @param pNAthr percentage of NA to be called as NA columns
#' @return an augMat object (also inherits matrix)
#' @export
augmentMatrixColumnByID <- function(mat, ID, computeNA=FALSE, pNAthr=0.8, data.matrix=TRUE){
	res <- matrix(NA, nrow=nrow(mat), ncol=length(ID))
	rownames(res) <- rownames(mat)
	colnames(res) <- ID
	IDfill <- intersect(ID, colnames(mat)) # ID's that have value filled from mat, no NA
	if(length(IDfill)<1) stop('ID has no shared elements with colnames(mat)!')
	#browser()
	#if(data.matrix){
	#	mat <- data.matrix(mat) # otherwise res is NULL
	#} else {
	#	res <- data.frame(res) # cannot assign a data frame to matrix in: res[, IDfill] <- mat[, IDfill]
	#}	
	mat <- as.matrix(mat) # will not convert non-string to NA
	res[, IDfill] <- mat[, IDfill]
	if(!computeNA) {
		naColID <- setdiff(ID, IDfill) # this might be inaccurate
	} else {
		pNA <- apply(res, 2, percentNA)
		naColID <- ID[pNA>pNAthr]
		IDfill <- setdiff(ID, naColID)
	}	
	if(inherits(mat, 'augMat')){
		naColID0 <- Getter(mat, 'naColID')
		IDfill <- setdiff(IDfill, naColID0)
		naColID <- union(naColID, naColID0)
	}
	#browser()
	filledColID <- IDfill
	naColInd <- match(naColID, ID)
	filledColInd <- match(IDfill, ID)
	attr(res, 'filledColID') <- filledColID
	attr(res, 'filledColInd') <- filledColInd
	attr(res, 'naColID') <- naColID
	attr(res, 'naColInd') <- naColInd
	attr(res, 'isNAcol') <- ID %in% naColID
	#browser()
	# bug: attr(res, 'effectiveMat') <- mat[, filledColInd]
	attr(res, 'N') <- length(filledColInd) # effective N
	if(length(filledColInd) < 1000)
		attr(res, 'effectiveMat') <- res[, filledColInd]
	#attr(res, 'class') <- c('augMat', 'matrix')
	#browser()
	class(res) <- append(class(res),"augMat")
	# to check class: inherits(res, 'augMat')
	res
}
#rppaHNT_aug <- augmentMatrixColumnByID(rppaHNT, ID=colnames(rnaHNT))
#MIR_NT_aug <- augmentMatrixColumnByID(MIR_NT, ID=ID_all_NT)

#' update augMat after it is subsetted
#' 
#' @param ID IDs targeted for; should be simialr to colnames(mat)
#' @param computeNA whether to recompute columns with all NA
#' @return augMat
#' @export
augmentMatrixColumnByID_update <- function(mat, ID, pNAthr=0.8){
	augmentMatrixColumnByID(mat=mat, ID=ID, computeNA=TRUE, pNAthr=pNAthr)
}

#RPPA_aug <- augmentMatrixColumnByID(RPPA_ST, ID=colnames(RNA))

#' Generic get method
#'
#' @param x an object
#' @param ... other arguments
#' @export
getter <- function(x, ...)
{
    UseMethod("getter",x)
}

#' Generic testree method (for phm object produced by plotHeatmap function)
#'
#' @param x an object
#' @param ... other arguments
#' @export
testree <- function(x, ...)
{
    UseMethod("testree",x)
}

#' get method for class augMat
#' 
#' this function is named getter, not get to minimize conflict with the usual get function
#'
#' @param obj an augMat object
#' @param what specify what to get from augMat, i.e. effectiveMat, isNAcol
#' @method getter augMat
#' @return an object as requested
#' @export
getter.augMat <- function(obj, what='effectiveMat'){
	#browser()
	if(what %in% names(attributes(obj))){
		res <- attributes(obj)[[what]]
	} else if (what =='avail'){
		res <- names(attributes(obj))
	} else {
		res <- NULL
	}
	res
}
#zz=getter(rppa, 'naColInd')
#zz=getter(rppa, 'effectiveMat')


### logical indicator if the elements is top N
#' return a logical vector if elements in a vector is top N largest
#' @param vec a numeric vector
#' @param N top N 
#' @return logic vector
#' @export
isTopN <- function(vec, N=100){
	cutoff0 <-  sort(vec, decreasing=TRUE, na.last=TRUE)[N]
	if(!is.na(cutoff0))	cutoff <- cutoff0 # when too many NA, make sure the cutoff is not NA
	else cutoff <- min(vec, na.rm=T)
	vec >= cutoff	
	#browser()
}
#' get top N elements from a vector and return a vector of index that is also ordered based on values from largest to smallest 
#' @param vec a vector of values
#' @param N number of elemetns to select
#' @return indeces corresponding to top N values from largest to smallest
#' @export
getTopNindex <- function(vec, N=NULL){
	if(is.null(N)) N <- length(vec)
	if(N>length(vec)) N <- length(vec)
	order(vec, decreasing=T)[1:N]
}
#' get top N elements from a vector 
#' @param vec a vector of values
#' @param N number of elemetns to select
#' @return indeces corresponding to top N values from largest to smallest
#' @export
getTopNvalue <- function(vec, N=NULL){
	if(is.null(N)) N <- length(vec)
	sort(vec, decreasing=T)[1:N]
}


#' force to factor when non-numeric
#' 
#' for a given data frame, this function converts non-numeric columns into factor while the numeric columns are intact.
#' @param df it is better to be a data frame. if a vector or matrix, it will be internally converted to data.frame
#' @return a data frame
#' @export
nonNumeric2FactorByColumn <- function(df){
	df1 <- as.data.frame(df)
	isNumeric <- colwise(class)(df1)
	res <- df1
	for(i in seq_along(isNumeric)){
		if(isNumeric[i]!='numeric'){
			res[, i] <- as.factor(df1[, i])
		} else {
			res[, i] <- df1[, i]
		}
	}
	res
}

#' conver logic vector to sorted index vector
#'
#' given a logic vector and corresponding value for ordering, reuturn an index vector where logicVec is TRUE and sort the index by valVec
#'
#' @param logicVec logic vector
#' @param valVec value vector; should have 1-to-1 correspondence to logicVec
#' @param decreasing used in sorting
#' @return a vector of index
#' @export
logic2sortedIndex <- function(logicVec, valVec, decreasing=TRUE){
	#browser()
	ind <- which(logicVec)
	ind[order(valVec[ind], decreasing=decreasing)]
}

#' swap two variables in a given environment
#'
#' @param x variable. in theory, any R object
#' @param y variable. in theory, any R object
#' @examples
#' a <- 3
#' b <- 4
#' swap(a, b)
#' @export
swap <- function(x, y, envir=parent.frame()) {
	xvar <- substitute(x)
	yvar <- substitute(y)
	xold <- x; yold <- y
	#browser()
	assign(deparse(xvar), yold, envir=envir)
	assign(deparse(yvar), xold, envir=envir)
}

#' log transformation of matrix, i.e. count matrix
#' 
#' After specifying a minimum value, any value smaller than it would be filled with it and then applied with a log transform.
#'
#' The user can specify a logfun for the transformation. 
#' @param dat a vector or a matrix to be safe-log
#' @param min a value specifying values to be filled if <=0
#' @param logfun a log function to be applied, i.e. log, log10 and so on
#' @return the safe-logged result, a vector or a matrix depending on the dat input
#' @export
safelog <- function(dat, min=1e-5, logfun='log2') {
#browser()
	# works for both matrix and vector
	dat[dat<=min] <- min
	do.call(logfun, list(x=dat))
}

#' force selected columns of a data frame as factor
#'
#' this is just a general (coarse) approach; if need to specify levels, do it manually.
#'
#' @param dat a data frame
#' @param cols either an index vector or a string character for colnames to be applied with as.factor
#' @export
df_force_factor <- function(dat, cols){
	if(is.numeric(cols)){
		# specified column indeces
		columns <- cols
	} else {
		# specified vars as col names
		columns <- rmNAinVec(match(cols, colnames(dat)))
	}
	# remove invalid cols, i.e. too large
	clms <- columns[columns %in% 1:ncol(dat)]
	vars <- colnames(dat)[clms]
	# plyr:::mutate, plyr:::transform both does not work when given a variable name. quit!
	dat[, clms] <- lapply(dat[, clms], as.factor) # the trick of lapply is smart! 
	dat
}

#' remove NA cases for either a vector, matrix or data frame
#' 
#' obtain complete records from your data
#'
#' @param data a vector, data frame or matrix
#' @param returnIndex whether just to return index of rows
#' @return same format as dat
#' @export
noNA <- function(dat, returnIndex=FALSE){
	#browser()
	sel <- complete.cases(dat)
	if(returnIndex)
		return(sel)
	if(is.null(dim(dat))) { 
		res <- dat[sel]
	} else {
		res <- dat[sel, ]
	}
	res
}
#tt<- noNA(t(L$ICmat[c(i, j), ]))
#tpDF2 <- noNA(tpDF[, 3])
#tpDF2 <- noNA(tpDF)


#' return the number of unique values
#' @param x a vector
#' @param ignoreNA logical, whether to count NA as a separate unique value; default is count NA as one unique value
#' @return integer
#' @export
n_unique <- function(x, ignoreNA=FALSE){
	if(ignoreNA){
		# when NA is ignored, just extract complete values of vector x
		x <- noNA(x)
	}
	return(length(unique(x)))
}
#' unique values sorted
#' @export
s_unique <- function(x, decreasing=FALSE){
	sort(unique(x), decreasing=decreasing)
}

#' guess if a vector is continuous
#'
#' we treat character as categorial; double as continuous; for integer/numeric/double, we can specify a minimum categories to be continuous variable
#'
#' a vector of class character is for sure categorical and would be good to handled as a factor in R. The class numeric
#' can be either integer (is.integer) or double (is.double). If is.double, then we should treat it as continuous when it has enough categories (binary
#' can be stored as double mode!).
#' For is.integer, when there is less than 5 categories, i.e., it is better treated as categorical.
#'
#' @param vec a vector
#' @param minCat an integer specifying minimum categories required to be a continuous variable for integers
#' @return TRUE for continuous and FALSE for anything else
#' @export
is_continuous <- function(vec, minCat=10){
	#browser()
	if(is.character(vec) | is.logical(vec)) return(FALSE) # character always not continuous
	if(is.factor(vec)) return(FALSE) # factor always not continuous
	#if(is.double(vec)) return(TRUE) # double always continuous----> this is wrong: OScensoring can be double! it just specifies the storage mode!
	# for numeric
	if(is.numeric(vec)) {
		if(n_unique(vec)>=minCat) {
			return(TRUE) # integer & has many categories, treat it as continuous
		} else {
			return(FALSE)
		}		
	}
	#browser()
}	
#is_continuous(c(TRUE, TRUE, FALSE))
#colwise(is_continuous)(clin)
#is_continuous(clin$OScensoring)

#' format a data frame by imposing non-continuous variables as factor though guessing
#'
#' continuous variables are guessed by the is_continuous function
#'
#' End product: all character/logical variables are forced to be factor; some integers or numeric would be forced to factor when they do not
#' satisfy minimum categories
#'
#' @param dat a data frame
#' @param dat a data frame
#' @return a formated data frame
#' @export
df_format_df <- function(dat, minCat=10){
	cols <- which(!colwise(is_continuous, minCat=minCat)(dat))# cols not continuous are forced to factor
	df_force_factor(dat, cols=cols)
}
#df_format_df(clin)


#' dichotomize expression into Low/High group 
#'
#' It's not recommeded to do this. But biologists wants to do it...
#'
#' values exceeding Upper limit are called High; values smaller than Lower are called low. Values inbetween will be masked as NA.
#' Warning: When Lower and Upper are the same (when i.e. 40% not expressed, 1/3 and 2/3 quantiles are the same value!), only one group is made!
#'
#' @param vec expressionv ector
#' @param q1 quantile 1 (smaller)
#' @param q2 quantile 2 (Laerger) when this specified, it overrides Lower/Upper specification
#' @param labels what you want to call the two groups
#' @return a character vector with middle values as NA
#' @export
dichotomizeExpr_quantile <- function(vec, q1=NA, q2=NA, Lower=NULL, Upper=NULL, labels=c('Low', 'High')){
	if(!is.na(q1) & !is.na(q2)){
		# this overrides the Lower/Upper limit
		Lower <- quantile(vec, prob=q1, na.rm=TRUE)
		Upper <- quantile(vec, prob=q2, na.rm=TRUE)
	}
	#browser()
	res <- rep(NA, length(vec))
	res[vec<Lower] <- labels[1]
	res[vec>=Upper] <- labels[2]
	if(n_unique(res)==1) warning('Only one group is formed due to equal Upper and Lower cutoff!')
	res
}
#dichotomizeExpr(response, q1=1/3, q2=2/3)

#' dichotomize expression into binary based on BI or quantile
#'
#' This is a wrapper for dichotomizeExpr_quantile() and dichotomizeExpr_BI()
#'
#' @param vec expressionv ector
#' @param method either quantile or BI
#' @param q1 quantile 1 (smaller)
#' @param q2 quantile 2 (Laerger) when this specified, it overrides Lower/Upper specification
#' @param labels what you want to call the two groups
#' @return a character vector 
#' @export
dichotomizeExpr <- function(vec, method='quantile', q1=NA, q2=NA, Lower=NULL, Upper=NULL, labels=c('Low', 'High'), plot=TRUE, forceFactor=FALSE){
	if(method=='quantile'){
		res <- dichotomizeExpr_quantile(vec, q1=q1, q2=q2, Lower=Lower, Upper=Upper, labels=labels)
	} else {
		res <- dichotomizeExpr_BI(vec, labels=labels, plot=plot)
	}
	if(forceFactor){
		res <- orderFactor(res, new.order=labels)
	}
	res
}

#' dichotomize expression into binary based on BI
#'
#' This is usually needed for KM curve or building surragote resonse variable from continuous variable
#'
#' @param vec expressionv ector
#' @param labels what you want to call the two groups
#' @return a character vector 
#' @export
dichotomizeExpr_BI <- function(vec, labels=c('Low', 'High'), plot=FALSE){
	#browser()
	fit <- plotNormalMix(vec, plot=plot)
	thr <- fit['thr']
	as.character(cut(vec, breaks=c(-Inf, thr, Inf), labels=labels))
}
#dichotomizeExpr_BI(rna['RBM10', ])


#' convert R variable names (i.e. from read table with default make.names) to natural variable names
#'
#' this replaces . with ' '; removes leading X and ending .; 
#' @param vars string vector for variable names
#' @return names vector of formated names
#' @export
toNaturalVar <- function(vars){
	#browser()
	vars1 <- str_replace(vars, '^X|\\.$', '') # remove leading X and ending .
	vars1 <- str_replace_all(vars1, '\\.+', ' ') # replace . or .. with ' '
	vars1 <- str_replace_all(vars1, '^ +| $', '') # remove leading white space or ending white space
	names(vars1) <- vars
	vars1
}


#' extract string given a pattern using back reference implemented by sub
#'
#' this function is an improved version of str_extract where we want to match the desired substring with a pattern using back reference. 
#' Notice that str_extract can describe what we want, which may be difficult. Now str_extract2 describes what we don't want and what is left is the desired result.
#'
#' @param string a string vector
#' @param pattern a pattern that to be passed to sub
#' @return returned by sub
#' @export
str_extract2 <- function(string, pattern){
	#browser()
	sub(pattern, '\\1', string)
}
#str_extractPat(rg[, 1], '.*Sample_(.*)/tophat_out.*')


showWithRestrict <- function(vec, restrict, Short=TRUE){
	if(length(vec)>restrict) { 
		ind <- 1:restrict
		cat(sprintf('***A subset of elements is shown below (the first %d elements):***\n', restrict))
	} else {
		ind <- 1:length(vec)
	}
	### whether to break lines so that the whole vector can be shown in limited width 
	if(Short) {
	#browser()
	res <- longVec2short(paste(vec[ind], collapse=', '), print=FALSE, maxLen=70, sep='\n')
	} else {
	res <- vec[ind]
	}
	res
}
longVec2short <- function(vec, maxLen=80, sep='\n\t', print=TRUE){
	Len <- nchar(vec)
	nSeg <- Len %/% maxLen
	if(nSeg<1)
		res <- vec
	else {
		st <- (0:nSeg)*maxLen+1
		pieces <- sapply(st, function(x) {
			ed <- x+maxLen-1
			if(ed>Len) 
				ed <- Len
			substr(vec, x, ed)	
		})
		res <- paste(pieces, collapse=sep)
	}
	#browser()
	if(print)
		cat(res, '\n')
	else {
		res
	}
}
#' Compare two vectors
#' @param name1 vector 1 of strings
#' @param name2 vector 2 of strings
#' @param more whether to show detailed info (e.g. examples)
#' @param restrict number of examples to show
#' @export
compareNames <- function(name1, name2, more=TRUE, restrict=50){
	vector1str <- deparse(substitute(name1))
	vector2str <- deparse(substitute(name2))
	ind1 <- name1 %in% name2
	ind2 <- name2 %in% name1
	#browser()
	if(all(ind1)) { ## name1 is a subset of name2
		if(all(ind2)) { # name1=name2
		cat("-->The two vectors are identical sets\n")
		if(identical(as.character(name1), as.character(name2))){
			cat('-->The two vectors are identical when converted to string\n')
		}
		if(identical(name1, name2)){
			cat('-->Further, the two vectors are identical\n')
			} else {
			cat('-->However, the two vectors are not identical, due to different ordering or duplicates or different storage mode\n')
			}
		} else { # name1 is a subset of name2
		cat(sprintf("-->Elements of %s are all in %s\n", vector1str, vector2str))
		cat(sprintf("-->These shared elements are: \n%s\n", showWithRestrict(name1, restrict)))
		#browser()
		cat(sprintf('-->%s has %d additional elements\n', vector2str, length(setdiff(name2, name1))))
		if(more) {
			#browser()
			cat(sprintf("-->Elements only in %s: \n%s\n", vector2str, showWithRestrict(sort(setdiff(name2, name1)), restrict)))
			}
		}
	} else {
		if(all(ind1==FALSE)) cat("-->None of the elements in vector 1 is in vector 2\n") ## not any of it in!
#browser()
		cat(sprintf('-->Elements only in %s: (N=%d)\n', vector1str, length(setdiff(name1, name2))))
		cat('------------------------------------------------\n')
		if(more)
			cat(showWithRestrict(sort(setdiff(name1, name2)), restrict), '\n')
		
		cat('------------------------------------------------\n')
		cat(sprintf('-->Shared elements: (N=%d)\n', length(sort(intersect(name1, name2)))))		
		if(more){
			if(length(intersect(name1, name2))>0)
				cat(showWithRestrict(intersect(name1, name2), restrict), '\n')
			else
				cat('-->No shared elements is observed!', '\n')
		}
		cat('------------------------------------------------\n')
		cat(sprintf('-->Elements only in %s: (N=%d)\n', vector2str, length(setdiff(name2, name1))))
		if(more){
			if(length(setdiff(name2, name1))>0)
				#browser()
				cat(showWithRestrict(sort(setdiff(name2, name1)), restrict), '\n')
			else
				{
					#browser()
					cat(sprintf('-->%s is a subset of %s. No additional elements!\n', vector2str, vector1str))
				}	
				
		}
		cat('------------------------------------------------\n')
	}
}

# sometimes the colbar has outlier that distorts the data. To save color, we truncate it.
## this also applies for a heatmap matrix where outlier data might dominate
### q1: lower quantile to truncate
### q1: upper quantile to truncate


#' truncate a matrix by quantile
#'
#' sometimes we explicitely truncate with upper or lower quantile
#'
#' @param q1 lower quantile. 
#' @param q2 upper quantile.  
#' @return formated matrix
#' @export
truncByQuantile <- function(x, q1, q2) {
	if(missing(q1)) q1 <- 0.01
	if(missing(q2)) q2 <- 0.99
	v <- quantile(x, c(q1, q2), na.rm=T)
	x[x<v[1]] <- v[1]
	x[x>v[2]] <- v[2]
	x
}
#' truncate a matrix by limit
#'
#' sometimes we explicitely truncate with upper or lower limit
#'
#' @param Lower lower cutoff value
#' @param Upper upper cutoff value
#' @return formated matrix
#' @export
truncByLimit <- function(x, Lower, Upper) {
	x[x<Lower] <- Lower
	x[x>Upper] <- Upper
	x
}
### we need to remove samples with IC50=8 in the icMat. This means we can just mask unwanted ICs with NA
#' mask values as NA in a matrix
#' 
#' @param x a matrix
#' @param Lower lower cutoff value
#' @param Upper upper cutoff value
#' @return formated matrix
#' @export
maskByLimit <- function(x, Lower, Upper) {
	#browser()
	if(any(x<Lower, na.rm=TRUE)){
	x[x<Lower] <- NA
	}
	if(any(x>Upper, na.rm=TRUE)){
	x[x>Upper] <- NA
	}
	x
}

## as.numeric will give warning for NA
	#t2FU <-  maskByLimit(as.numeric(td1[, 'last_contact_days_to']), Lower=0, Upper=Inf)/30 # to month
	#t2death <- maskByLimit(as.numeric(td1[, 'death_days_to']), Lower=0, Upper=Inf)/30 # to month

 
#' truncate matrix used to generate heatmap; add option to scale
#'
#' this is especially useful for heatmap
#'
#' @param mat the matrix to be formated
#' @param scale whether to scale the data before formatting
#' @param Lower lower cutoff value
#' @param Upper upper cutoff value
#' @param q1 lower quantile. It is effective only if Lower is not specified
#' @param q2 upper quantile. It is effective only if Lower is not specified
#' @return formated matrix
#' @export
truncMat <- function(mat, scale=TRUE, Lower, Upper, q1=0.01, q2=0.99){
	if(scale){
		mat <- t(scale(t(mat))) # scale within row
	}
	# if missing Lower, use quantile truncation
	if(missing(Lower)) {
		v <- quantile(mat, c(q1, q2), na.rm=T)
		Lower <- v[1]
		Upper <- v[2]
	}
	res <- mat
	res[mat > Upper] <- Upper
	res[mat < Lower] <- Lower
	#browser()
	res
}

#' guess alpha
#' @param G number of genes
#' @return a scalar between 0 and 1
#' @export
guessAlpha <- function(G) {
	if(G<300) {
		return(1)
	} else if(G<1000) {	
		return(0.8)
	} else if(G<2000) {	
		return(0.7)
	} else if (G<5000) {	
		return(0.5)
	} else if (G<10000) {
		return(0.4)
	} else if (G<30000) {
		return(0.3)
	} else {
		return(0.2)
	}
}

#' guess alpha in a line plot
#' @param G number of genes or lines
#' @return a scalar between 0 and 1
#' @export
guessAlpha_line <- function(G) {
	if(G<20) {
		return(0.9)
	} else if(G<100){
		return(0.5)
	} else if(G<200) {	
		return(0.3)
	} else if(G<500) {	
		return(0.2)
	} else if (G<1000) {	
		return(0.1)
	} else if (G<10000) {
		return(0.05)
	} else {
		return(0.01)
	}
}

#' given number of genes G in a heatmap, guess rowCEX to be passed to cexRow
#' @param G number of genes
#' @return a scalar between 0 and 1
#' @export
guessCEX <- function(G) {
	if(G<5) {
		rowCEX <- 2
	} else if(G<10) {	
		rowCEX <- 1.5
	} else if(G<20) {	
		rowCEX <- 1.3
	} else if(G<30) {	
		rowCEX <- 1.2
	} else if(G<50) {	
		rowCEX <- 1.0
	} else if(G<80) {	
		rowCEX <- 0.75
	} else if (G<100) {	
		rowCEX <- 0.5
	} else if (G<200) {
		rowCEX <- 0.25
	} else if (G<300) {
		rowCEX <- 0.17
	} else if (G<400) {
		rowCEX <- 0.12
	} else if (G<500) {
		rowCEX <- 0.10
	} else if (G<600) {
		rowCEX <- 0.08
	} else if (G<700) {
		rowCEX <- 0.05
	} else {
		rowCEX <- 0.02
	}
	rowCEX
}

