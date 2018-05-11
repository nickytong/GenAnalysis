#' given a p value vector and a vector of cutoffs, find number of genes as well as indeces for genes passing these cutoffs
#'
#' @param P a vector of P values
#' @param thr a vector of thersholds
#' @return a list with elements smryTab, selList, cutoff
#' @export
getPvecSmry <- function(P, thr=c(1e-3, 1e-2, 0.05)){
	names(thr) <- thr
	nSel <- sapply(thr, function(x) sum(P<x, na.rm=TRUE))
	lSel <- lapply(thr, function(x) {
		ind <- which(P<x)
		rr <- ind[order(P[ind], decreasing=FALSE)]# now order the index by P from smallest to largest
		rr
		})
	smryTab <- data.frame(cutoff=thr, nSel=nSel)
	#browser()
	res <- list(smryTab=smryTab, selList=lSel, cutoff=unname(thr))
	res
}
#printPvecSmry(ttest_MUT[, 'pval'])

#getPvecSmry(ttest_MUT_RNA[, 'pval'], thr=c(1e-3))
#getPvecSmry(ttest_MUT_RNA[, 'pval'], thr=c(1e-3, 1e-2, 0.05))


#' given a p value vector and a vector of cutoffs, find number of genes as well as indeces for genes passing these cutoffs
#'
#' this function just print the summary, a companion of getPvecSmry
#' 
#' @param P a vector of P values
#' @param thr a vector of thersholds
#' @return a list with elements smryTab, selList, cutoff
#' @export
printPvecSmry <- function(P, thr=c(1e-3, 1e-2, 0.05)){
	getPvecSmry(P=P, thr=thr)$smryTab
}




################# WARNING ###################
## for small sample size, cor.test is more accurate than rcorr when use spearman: the fortran library in rcorr is not good!!!
## i.e. x <- [1] -5.76420  1.59860 -3.86690 -3.67550  0.83585
##      y <- [1] -0.48691601  4.00699993 -0.05349526  0.46392579  1.43104670
### correlation between two vectors
################# WARNING ###################

## computes correlation between two vectors: more accurate
## (1) cor.test: getAnywhere('cor.test.default') shows this is written in R. More accurate but slower.
#' correlation betweent two vector
#' 
#' @param x x
#' @param y y
#' @param method method use
#' @return a vector
#' @export
corTest <- function(x, y, method='spearman', minN=5, ...) {
	# check if minimum sample size is satisified
	df <- data.frame(x, y)
	df1 <- df[complete.cases(df), ]
	#browser()
	if(!missing(minN) & !is.na(minN)) { #  !is.na(minN) added for backward compatibility. note analysis1.Rnw might not work now
		#browser()
		if(nrow(df1)<minN) {
			res <- c(NA, NA)
			return(res)
		}	
	}
	### this continues if sample size > minN
	tm <- cor.test(x, y, method=method, ...)
	res <- c(tm$estimate, tm$p.value)
	names(res) <- c('cor', 'pval')
	res
}

## computes correlation between two vectors
## (2) rcorr for speed
#### added on 2013/09/05: use rcorr for speed
#### correlation between two vectors
## modified on 2013/12/10 to make minN strictly applied (previous version fails when all IC50 are NA as in db2.Rnw)
#' computes correlation between two vectors
#' @export
rcorrTest <- function(x, y, method='spearman', minN=3) {
	# check if minimum sample size is satisified
	df <- data.frame(x, y)
	df1 <- df[complete.cases(df), ]
	#browser()
	if(!missing(minN) & !is.na(minN)) { #  !is.na(minN) added for backward compatibility. note analysis1.Rnw might not work now
		#browser()
		if(nrow(df1)<minN) {
			res <- c(NA, NA)
			return(res)
		}	
	}
	### this continues if sample size > minN
	tp <- rcorr(x=df1$x, y=df1$y, type=method)
	res <- c(tp$r[1, 2], tp$P[1, 2])
	#browser()
	res
}
#rcorrTest(ic, score, method='spearman', minN=5)


## computes pair-wise correlation among columns use rcorr
## typical application might be searching for co-mutation, mutual exclusivity
rcorrMat <- function(mat, method='pearson'){
	#browser()
	tp <- rcorr(mat, type=method)
	res <- list(corMat=tp$r, pvalMat=tp$P)
	res
}
#res <- rcorrMat(t(mutCCLE_glioma))

#rcorrTest(c(rnorm(50), NA), c(rnorm(50), NA), minN=100)
### correlation between each row in a matrix and a vector
### typical application is to compute for each gene, cor between IC50 and mRNA for a drug
### use rcorr() from Hmisc package: (1) NAs deleted in pairs rather than in all rows (2) fast implementation in Fortran (3) provides p values.
### usercorr: added on 2013/09/18: whether to use rcorrTest (fast) or 

updateRowname <- function(mat, symbol=NULL, sepSymProbe=', '){
	res <- mat
	## when symbol is not specified, it is assumed to be gene level data; nothing needs to be done
	if(!is.null(symbol)){
		if(length(symbol)!=nrow(mat)){
			stop('symbol is not NULL thus assuming probe level data; it has unequal length as nrow(mat)!')
		} else {
			#indNA <- which(is.na(symbol))
			#browser()
			symbol[is.na(symbol)] <- 'N.A.'
			rownames(res) <- str_c(symbol, rownames(mat), sep=sepSymProbe)
		}
	}
	res
}
#' Calculate correlation between rows of a matrix and a vector
#' 
#' @param mat a matrix, i.e. gene expression
#' @param vec a vector of phenotype, continuous variable
#' @param method method passed to correlation calculation. It can be spearman or pearson
#' @param minN minimum sample size for correlation calculation; If not achieved, correlation is returned as NA
#' @param usercorr whether to use rcorr package for speeded calculation (using Fortran)
#' @param meta meta information, a list. This can specify what type of analysis (i.e. RPPA correlated with EMT score);
#'  row, which represents what row variables are (rows of mat); response represents what does vec mean
#' @param symbol in case this is probe level data, user can specify a vector of symbols matching 
#' @param sepSymProbe separator for symbol-probe. That is, symbol+sepSymProbe+probe(rowname) is the new rowname
#' @param FDR FDR to create bum for the p values, passed to create_bum()
#' @param pcutoff pcutoff passed to create_bum()
#' @param alphas alphas passed to create_bum()
#' @return a matrix with two columns: cor and pval
#' @export
rowcorMatVec <- function(mat, vec, meta=list(tag='correlation', row='expression', response='response'), 
	method='spearman', minN=5, usercorr=TRUE, addData=TRUE, symbol=NULL, sepSymProbe=', ',
	FDR=NULL, pcutoff = 0.05, alphas = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
	addBum=NULL) {
	#browser()
	## update rownames to deal with probe level data
	mat <- updateRowname(mat=mat, symbol=symbol, sepSymProbe=sepSymProbe)
	if(ncol(mat)!=length(vec)) stop(sprintf("mat ncol: %d\nvec length: %d\n", nrow(mat), length(vec)))
	if(usercorr) {
		# for speed
		#browser()
		tmp <- t(mapply(rcorrTest, data.frame(t(mat)), data.frame(vec), MoreArgs=list(method=method, minN=minN)))
	} else {	
		# for accuracy
		tmp <- t(mapply(corTest, data.frame(t(mat)), data.frame(vec), MoreArgs=list(method=method, minN=minN)))
	}
	#browser()
	if(is.null(addBum)){
		if(nrow(mat)>100){
			addBum <- TRUE
		} else {
			addBum <- FALSE
		}
	}
	res <- data.frame(mat[, 1:2]) # hope to preserve the row names
	res[, 1:2] <- tmp
	rownames(res) <- rownames(mat)
	colnames(res) <- c('cor', 'pval')
	if(addBum)
		obj_bum <- try(create_bum(res$pval, FDR=FDR, pcutoff=pcutoff, alphas=alphas), silent=TRUE)
	attr(res, 'meta') <- meta
	attr(res, 'method') <- capitalize(method)
	attr(res, 'minN') <- minN
	attr(res, 'addData') <- addData
	if(addData){
		attr(res, 'mat') <- mat
		attr(res, 'vec') <- vec
	}
	attr(res, 'hybrid') <- FALSE
	if(addBum)
		attr(res, 'obj_bum') <- obj_bum
	attr(res, 'class') <- c('rowcorMatVec', 'data.frame')
	#browser()
	res
}

#' Heatmap plot for rowcorMatVec class
#'
#' @method heatm rowcorMatVec
#' @param x a rowcorMatVec class (a data frame)
#' @param gSel string specifying what genes are selected in the elements of 'obj_bum'
#' @param revlog10 sometimes vec is log10 IC50; it is better to show IC50 in original scale in the heatmap; this will be possible by specifying revlog10=TRUE
#' @export
heatm.rowcorMatVec <- function(x, plot=TRUE, return=TRUE, gSelForce=NULL, gSel='gSelectInd_P05', colbarname='vec', 
	revlog10=FALSE, plotLegend='TRUE', ...){
	#browser()
	obj_bum <- Getter(x, 'obj_bum')
	vec <- Getter(x, 'vec')
	if(revlog10) {
		vec <- 10^vec
		code_vec <- sprintf("vec <- 10^Getter(x, 'vec')")
	} else {
		code_vec <- sprintf("vec <- Getter(x, 'vec')")
	}
	mat <- Getter(x, 'mat')
	code0 <- sprintf("\tx <- %s\n\tgSel <- '%s'", as.character(substitute(x)), as.character(substitute(gSel)))
	code_gSel1 <- "if(!is.null(gSelForce)){
		gSel1 <- gSelForce
	} else {
		obj_bum <- Getter(x, 'obj_bum')
		gSel1 <- obj_bum[[gSel]]
	}
	"
	#browser()
	code1=sprintf("
	obj_bum <- Getter(x, 'obj_bum')
	%s
	mat <- Getter(x, 'mat')
	%s
	pheno <- data.frame(vec=vec)
	colnames(pheno) <- '%s'
	rownames(pheno) <- colnames(mat)
	colpal <- list(%s=colorpalette2colvec('greenred'))
	phm <- aheatmat(mat=mat, sSel1=which(!is.na(pheno$%s)), gSel=gSel1, 
    	pheno=pheno, q1=0.05, q2=0.95, clustering_distance_rows = 'pearson', 
		cluster_cols=F, colOrderIndex=order(pheno$%s, decreasing=F), 
		cluster_rows=T, 
   		clustering_method = 'ward.D', yd=0.2, colpal=colpal, plotLegend=%s, ...)
	\n", code_vec, code_gSel1, colbarname, colbarname, colbarname, colbarname, plotLegend)
	code <- sprintf('%s%s', code0, code1)
	res <- list(code=code)
	if(plot){
		#pander::evals(code)
		#browser()
		phm <- eval(parse(text=code))
		#browser()
		res$phm <- phm
		mapL <- phm
		stretchL <- llply(mapL, get_stretch)
		res$stretchL <- stretchL
		# obj_bum <- Getter(x, 'obj_bum')
		# vec <- Getter(x, 'vec')
		# mat <- Getter(x, 'mat')
		# gSel1 <- obj_bum[[gSel]]
		# pheno <- data.frame(vec=vec)
		# rownames(pheno) <- colnames(mat)
		# colpal <- list(vec=colorpalette2colvec('greenred'))
		# phm <- aheatmat(mat=mat, sSel1=which(!is.na(pheno$vec)), gSel=gSel1, 
  # 		  	pheno=pheno, q1=0.05, q2=0.95, clustering_distance_rows = 'pearson', 
		# 	cluster_cols=F, colOrderIndex=order(pheno$vec, decreasing=F), 
		# 	cluster_rows=T, 
  # 	 		clustering_method = 'ward.D', yd=0.2, colpal=colpal, plotLegend=T)	
	}
	if(return){
		cat(code)
		return(res)
	}
}


#cor_EMTscore <- rowcorMatVec(mat=ic50_aug, vec=pData$EMTscore, method='spearman', meta=list(tag='correlation', row='IC50', response='EMT score'))

#' when pearson and spearman values both are desired, we can fit them separately and then combine
#'
#' @param fit1 first fit from rowcorMatVec
#' @param fit2 second fit from rowcorMatVec
#' @return rowcorMatVec object, slightly modified
#' @export
combine_cormat <- function(fit1, fit2, tag=c('spearman', 'pearson')){
	fit <- moveColumnToRowName(merge(fit1, fit2, by.x=0, by.y=0, suffixes=c('_spearman', '_pearson')), column=1)
	fit <- fit[rownames(fit1), ]
	# first tag as primary
	#fit$cor <- fit1$cor
	#fit$pval <- fit1$pval
	#browser()
	addData <- Getter(fit1, 'addData')
	attr(fit, 'class') <- c('rowcorMatVec', 'data.frame')
	# following attributes are wrong; but just to get started
	meta <- Getter(fit1, 'meta')
	meta$tag <- sprintf('Correlation between %s and %s', meta$row, meta$response)
	attr(fit, 'meta') <- meta
	attr(fit, 'addData') <- addData
	if(addData){
		attr(fit, 'mat') <- Getter(fit1, 'mat')
		attr(fit, 'vec') <- Getter(fit1, 'vec')
	}
	attr(fit, 'hybrid') <- TRUE
	attr(fit, 'method') <- tag
	attr(fit, 'minN') <- Getter(fit1, 'minN')
	fit
}
#fit <- combine_cormat(fit1=cor_Cisplatin, fit2=cor_Cisplatin_pearson, tag=c('spearman', 'pearson'))
#' scatter plot for rowcorMatVec
#'
#' this function needs to be polished and put into TestAndCorr_p.R
#'
#' @param mat expression matrix used in rowcorMatVec 
#' @param vec vec as used in rowcorMatVec 
#' @param FUN plotting function
#' @param geneSel index for plotting; default is plotting all sorted by P
#' @param rowcorMatVec rowcorMatVec object
#' @param coord_flip whether to flip x-y coordinate
#' @param plot.margin when coord_flip is true, margin may be too small; here we can adjust it
#' @param ... additonal parameters to plotScatter
#' @return multiple ggplot figures printed
#' @export
scatter_rowcorMatVec <- function(mat, vec, FUN, geneSel, rowcorMatVec, smooth=FALSE, lm=TRUE, labels=NULL, coord_flip=FALSE, plot.margin = NULL, ...){
	if(missing(geneSel)){
		#geneSel <- order(rowcorMatVec$pval, decreasing=FALSE)
		geneSel <- setdiff(order(rowcorMatVec$pval, decreasing=FALSE), which(is.na(rowcorMatVec$pval))) # remove NA records where no data points
	}
	method <- getter(rowcorMatVec, 'method')
	if(is.null(smooth)) smooth <- ifelse(method=='pearson', FALSE, TRUE)
	if(is.null(lm)) lm <- ifelse(method=='pearson', TRUE, FALSE)
	if(missing(FUN)){
		# notice that xlab and ylab are from meta info
		# meta$tag is used as the overall header, e.g. BATTLE 2, KRAS mutants
		# xlab is gene
		FUN <- function(expr, vec, gene, meta, ...) {
			#browser()
			if(length(method)==1){
				tag <- ifelse(is.null(meta$tag), gene, meta$tag)
				main <- sprintf('%s\n%s corr=%.3f\np value=%.3e', tag, method, rowcorMatVec[gene, 'cor'], rowcorMatVec[gene, 'pval'])
				p <- plotScatter(expr, vec, xlab=gene, ylab=meta$response, main=main, 
					plot=F, smooth=smooth, lm=lm, labels=labels, ...)$p
			
			} else {
				#browser()
				tag <- ifelse(is.null(meta$tag), gene, meta$tag)
				main <- sprintf('%s\n%s corr=%.3f, p value=%.3g\n%s corr=%.3f, p value=%.3g', tag, 'Spearman', rowcorMatVec[gene, 'cor_spearman'], rowcorMatVec[gene, 'pval_spearman'],
																 'Pearson', rowcorMatVec[gene, 'cor_pearson'], rowcorMatVec[gene, 'pval_pearson'])
				p <- plotScatter(expr, vec, xlab=gene, ylab=meta$response, main=main, 
					plot=F, smooth=smooth, lm=lm, labels=labels, ...)$p
			}
			if(coord_flip){
				p <- p+ coord_flip()
			}
			if(!is.null(plot.margin)){
				p <- p+theme(plot.margin = plot.margin)
			}
			p
		}
	}
	meta <- getter(rowcorMatVec, 'meta') # meta info
	for(g in geneSel){
		#browser()
		gene <- rownames(mat)[g]
		#browser()
		tt=try(FUN(expr=mat[g, ], vec=vec, gene=gene, meta, ...), silent=T)
		if(!inherits(tt, 'try-error')) {
			print(tt)
		}	
	}
}
#scatter_rowcorMatVec(mat=ic50_aug, vec=pData$EMTscore, rowcorMatVec=cor_EMTscore)

#' Scatter plot for rowcorMatVec class: each row produces a printted figure
#'
#' @method scatter rowcorMatVec
#' @param x a rowcorMatVec class (a data frame)
#' @param ... additional parameters to scatter_rowcorMatVec
#' @param addLabel whether to add sample names as text in the scatter plot
#' @export
scatter.rowcorMatVec <- function(x, maxRow=1000, main=NULL, addLabel=FALSE, ...){
	hybrid <- Getter(x, 'hybrid')
	if(!hybrid){
		pvec <- x[, 'pval']
	} else {
		tag <- Getter(x, 'method')
		# choose primary method, the first one
		pvec <- x[, str_detect(colnames(x), tag[1])&str_detect(colnames(x), 'pval')]
	}
	# when too many rows, select top rows based on P
	#browser()
	if(nrow(x)>maxRow){
		ind <- which(isTopN(-pvec, N=maxRow))[1:maxRow] # [1:N] in case too many 0 p values. 
		ind <- ind[order(pvec[ind], decreasing=FALSE)] # ordering
		main <- sprintf('%s\nTop %d genes selected based on p value', main, maxRow)
	} else {
		ind <- order(pvec, decreasing=FALSE)
	}
	# removes NA records
	ind <- setdiff(ind, which(is.na(pvec)))
	# requires addData=TRUE
	if(is.null(getter(x, 'mat'))) 
		stop('No mat (e.g. expression matrix) found! Please specify addData=TRUE\n')
	#browser()
	mat <- getter(x, 'mat')
	vec <- getter(x, 'vec')
	#
	# fix bug: need to also subset mat !
	#
	if(!identical(rownames(x), rownames(mat))){
		mat <- mat[rownames(x), ] # subsetting
	}
	if(addLabel) {
		labels <- colnames(mat)
	} else {
		labels <- NULL
	}
	scatter_rowcorMatVec(mat=mat, vec=vec, geneSel=ind, rowcorMatVec=x, labels=labels, ...)
} 
#scatter.rowcorMatVec(cor_EMTscore, lm=T)

#scatter_rowcorMatVec(cor_EMTscore)
#' Bar plot for rowcorMatVec class
#'
#' @method plot rowcorMatVec
#' @param x a rowcorMatVec class (a data frame)
#' @param ... additional parameters to plot4rowttestMatVec
#' @export
plot.rowcorMatVec <- function(x, maxRow=400, main=NULL, tag='', ...){
	# when too many rows, select top rows based on P
	#browser()
	if(is.null(main)) main <- sprintf('%s%s', getAnalysisName(x), tag)
	hybrid <- Getter(x, 'hybrid')
	if(!hybrid){
		pvec <- x[, 'pval']
	} else {
		tag <- Getter(x, 'method')
		# choose primary method, the first one
		pvec <- x[, str_detect(colnames(x), tag[1])&str_detect(colnames(x), 'pval')]
	}
	if(nrow(x)>maxRow){
		ind <- which(isTopN(-pvec, N=maxRow))
		main <- sprintf('%s\nTop %d hits selected based on p value', main, maxRow)
	} else {
		ind <- 1:nrow(x)
	}
	plot4rowcorMatVec(x[ind, ], main=main, ...)
} 

#' Bum plot for rowcorMatVec class
#'
#' @method hist rowcorMatVec
#' @param x a rowcorMatVec class (a data frame)
#' @param ... additional parameters to plotBumFDR
#' @export
hist.rowcorMatVec <- function(x, main=NULL, tag='', plot=TRUE, return=FALSE, alphas = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
    FDR=NULL, maskPthr = 1, pcutoff = 0.05, ...){
	obj_bum <- create_bum(pval=x[, 'pval'], alphas=alphas, FDR=FDR, maskPthr=maskPthr, pcutoff=pcutoff)
	#lev <- attributes(x)$levels
	if(is.null(main)){
		main <- sprintf('%s%s', getAnalysisName(x), tag)
	}
	#browser()
	if(plot){
		plotBumFDR(obj_bum, main=main, ...)
	}
	if(return){
		return(obj_bum)
	}
}

#' get method for class rowcorMatVec
#' 
#' this function is named getter, not get to minimize conflict with the usual get function
#'
#' @param obj an rowcorMatVec object
#' @param what specify what to get from augMat, i.e. base, var.equal
#' @method getter rowcorMatVec
#' @return an object as requested
#' @export
getter.rowcorMatVec <- function(obj, what='meta'){
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

#' ordered bar plot for rowcorMatVec class
#'
#'
#' @param rowInfo alternative rownames to show on the figure
#' @param decreasing logic, decreasing order of bars or not.
#' @param dataMode a string used to label the x-axis. Default is 'IC50' for IC50 t-test; Please use meaning string relevant to the data. 
#' @param showP do not show P value on graph (as part of gene name) ('None'); show actual P ('Value'); show P category ('Category')
#' @param pthr filtering by p value. Rows with larger P value will be deleted. Default is to use all rows (pthr=1)
#' @param legend.position position for the legend if P value is show as category
#' @param sizex x axis label size
#' @param sizey y axis label size
#' @param legend.title title of the legend
#' @param truncCorlow cutoff to truncate correlation lower than this value. Defalt is -1, which means no truncation on small corr (negative) value
#' @param truncCorhigh cutoff to truncate correlation higher than this value. Defalt is 1, which means no truncation on large corr (usually positive) value
#' @param barWidth relative barwidth (0~1); if this is too small, some bars might be displayed as missing; so set it larger
#' @param corlim correlation limit; default is not specified which will use observed correlation limit; can specify c(-1, 1) to see this range.
#' @param coord_flip whether to flip coordinate; default sets y as the correlation; when coord_flip=TRUE, x is plotted as correlation
#' @return a ggplot2 object
#' @export
plot4rowcorMatVec <- function(res, rowInfo=NULL, decreasing=TRUE, dataMode='IC50', ylab=NULL, xlab=NULL, sizex, sizey=12, main=NULL, tag='', 
	showP='None', pthr=1, colorSignif=TRUE, legend.position=c(.85,0.85), legend.title='Significance', showLegend=TRUE, 
	truncCorlow=-1, truncCorhigh=1, barWidth=NULL, corlim, coord_flip=FALSE, plot=TRUE){
	if(is.null(xlab)) xlab <- sprintf('%s correlation', getter(res, 'method'))
	if(is.null(ylab)) ylab <- ''
	if(is.null(main)) main <- sprintf('%s', tag)
	res <- data.frame(res)
	### sometimes NA for pval: this leads to error
	res <- subset(res, !is.na(pval))
	#browser()
	# truncate FC 
	res <- mutate(res, cor=truncByLimit(cor, Lower=truncCorlow, Upper=truncCorhigh))
	# add P category
	#res$pcat <- categorizePvec(res[, 'pval'])
	res$pcat <- droplevels(categorizePvec(res[, 'pval']))
	#
	if(is.null(rowInfo)) rowInfo <- rownames(res)
	# showP: add P into drug info
	if(showP=='Value'){
		rowInfo <- mapply(function(x, y) sprintf('%s, P=%.3f', x, y), rowInfo, res$pval)
	} else if(showP=='Category'){
		rowInfo <- mapply(function(x, y) sprintf('%s, %s', x, y), rowInfo, res$pcat)
	}
	#browser()
	# format res
	res <- data.frame(res, rowInfo=rowInfo) # add rowInfo
	res <- subset(res, !is.na(cor) & pval<pthr) # remove rows with FC NA (not enough sample size)
	res$rowOrdered <- orderVecByVal(vec=res$rowInfo, val=res$cor, decreasing=decreasing) # levels added
	if(missing(corlim)) corlim <- range(res$cor)
	# when all cor are positive, specify corlim[1]>0 will makes the line disappear
	if(corlim[1]>0) corlim[1] <- 0
	if(corlim[2]<0) corlim[2] <- 0
	if(!coord_flip) {
		# y is cor: need to rotate gene name
		text_angle <- 90
	} else {
		text_angle <- 0
	}
	if(missing(sizex)) sizex <- guessSize(nrow(res))
	scale_axis_size <- theme(axis.text.x=element_text(size=sizex, colour = "black", face ='bold', angle=text_angle, hjust=0.5, vjust=0.5), 
		axis.text.y=element_text(colour = "black", face ='bold', size=sizey))
	#browser()
	#fillCol <- c('red', 'magenta', 'purple', 'blue')
	fillCol <- c('#ED281E', '#2834A4', '#EDB51E', 'gray')
	names(fillCol) <- c('P<0.001', 'P<0.01', 'P<0.05', 'P>0.05')
	fillColVec <- fillCol[sort(unique(as.character(res$pcat)), decreasing=FALSE)] # sorted order: first smaller p value
	if(is.null(barWidth)){
		barWidth <- ifelse(nrow(res)<50, 0.2, 0.5)
	}
	#browser()
	if(!colorSignif) {
	p <- ggplot(res, aes(y=cor, x=rowOrdered))+geom_bar(stat = "identity", position="dodge", width=barWidth)+ 
		ylab(xlab)+xlab(ylab)+ggtitle(main)+scale_axis_size+ylim(corlim)+
		theme(legend.position=legend.position)
	} else {
	p <- ggplot(res, aes(y=cor, x=rowOrdered, fill=pcat))+geom_bar(stat = "identity", position="dodge", width=barWidth)+ 
		ylab(xlab)+xlab(ylab)+ggtitle(main)+ylim(corlim)+
		scale_fill_manual(name=legend.title, breaks=names(fillColVec), values=fillColVec)+
		scale_axis_size+theme(legend.position=legend.position)		
	}	
	if(coord_flip) p <- p+coord_flip()
	if(showLegend==FALSE) {
		#browser()
		p <- p + theme(legend.position = "none")
	}
	#browser()
	if(plot==FALSE){
		return(list(p=p, dat=res))
	}
	p
	
}

#plot4rowcorMatVec(cor_EMTscore)

#lmMirWithREST_GBM <- rowcorMatVec(mat=mat_mir2[, pID_used_MIR], vec=restScore_GBMarray[pID_used_MIR])
#mat <- matrix(rnorm(200), ncol=50)
#vec <- rnorm(50)
#rowcorMatVec(mat, vec)

#' extension of rowcorMatVec
#'
#' @param mat a matrix, i.e. gene expression
#' @param matR a matrix for response; each row is a response; ncol should match with ncol(mat)
#' @param method method passed to correlation calculation. It can be spearman or pearson
#' @param minN minimum sample size for correlation calculation; If not achieved, correlation is returned as NA
#' @export
rowcorMatMat_L <- function(mat, matR, method='spearman', minN=5, usercorr=TRUE, addData=TRUE, symbol=NULL, sepSymProbe=', ', meta=list(responseSuffix='', row='', tag='correlation')){
	#browser()
	if(ncol(mat)!=ncol(matR)) stop(sprintf('ncol(mat)!=ncol(matR)!'))
	fitL <- foreach(i=1:nrow(matR))	%do% {
		rnR <- rownames(matR)[i]
		#browser()
		rowcorMatVec(mat=mat, vec=matR[i, ], method=method, meta=list(response=sprintf('%s%s', rnR, meta$responseSuffix), row=meta$row, tag=meta$tag), 
			minN=minN, usercorr=usercorr, addData=addData, symbol=symbol, sepSymProbe=sepSymProbe)
	}
	names(fitL) <- rownames(matR)
	suffix=sapply(fitL, function(x) Getter(x, 'meta')$response)
	#browser()
	table <- cbind2dfwithSuffix_fromList(fitL, suffix=str_c('_', suffix))
	array <- foreach(i=1:length(fitL), .combine='combIntoArray') %do% {
			fitL[[i]]
	}
	dimnames(array)[[3]] <- rownames(matR)
	list(fitL=fitL, table=table, array=array, meta=meta)
}
#corL_ets_rna_UTSW_NSCLC <- rowcorMatMat(mat=expr_NSCLC[geneAllison_ETS, ], matR=ic50utsw_aug_NSCLC[drugCisplatin, ], method='spearman', meta=list(responseSuffix=' IC50', row=' expression', tag='correlation'))

#' extension of rowttestMatVec
#'
#' @param mat a matrix, i.e. gene expression
#' @param matR a matrix for response; each row is a response; ncol should match with ncol(mat)
#' @param method method passed to correlation calculation. It can be spearman or pearson
#' @param minN minimum sample size for correlation calculation; If not achieved, correlation is returned as NA
#' @export
rowttestMatMat_L <- function(mat, matR, method='spearman', minN=5, usercorr=TRUE, addData=TRUE, symbol=NULL, sepSymProbe=', ', meta=list(responseSuffix='', row='', tag='correlation')){
	#browser()
	if(ncol(mat)!=ncol(matR)) stop(sprintf('ncol(mat)!=ncol(matR)!'))
	fitL <- foreach(i=1:nrow(matR))	%do% {
		rnR <- rownames(matR)[i]
		#browser()
		rowttestMatVec(mat=mat, vec=matR[i, ], method=method, meta=list(response=sprintf('%s%s', rnR, meta$responseSuffix), row=meta$row, tag=meta$tag), 
			minN=minN, usercorr=usercorr, addData=addData, symbol=symbol, sepSymProbe=sepSymProbe)
	}
	suffix=sapply(fitL, function(x) Getter(x, 'meta')$response)
	#browser()
	table <- cbind2dfwithSuffix_fromList(fitL, suffix=str_c('_', suffix))
	list(fitL=fitL, table=table)
}

#' the driver of grouped version calculations, enpowering rowcorMatVec, rowttestMatVec
#'
#' This calculated result between rows of mat and vec but at each group at a time as specified by the group vector
#' It should be more flexible than rowcorMatMat which requires the same genomic mat; instead, here different group of analysis can
#' have different sample size in genomics
#' @param mat a matrix, i.e. gene expression
#' @param vec a vector of phenotype, continuous variable
#' @param group a vector, categorical variable
#' @param addAll whether to calculate correlation by using all data ignoring group; this is the overall correlation
#' @param fn function to apply, e.g. rowttestMatVec, rowcorMatVec
#' @param parallel deprecated; parallel is automatically enabled when core>1
#' @param core number of cores to use; when core>1, parallel execution is enabled
#' @param returnEarly return early for the fitted list; useful for rowanovamatVec since array is not feasible due to different dimensions of fits
#' @param ... additional parameters to fn
#' @return a 3D array
#' @export
rowdoMatVec_Bygroup <- function(mat, vec, group, fn, addAll=FALSE, addAllName='all', parallel=FALSE, core=1, returnArray=TRUE, returnEarly=FALSE, ...){
	#browser()
	df <- data.frame(Group=group, 
		Response=vec, t(mat), check.names=FALSE)
	# addAll
	if(addAll){
		df <- dfaugWithAll(df, groupVar='Group', addAllName=addAllName)
	}
	if(missing(core)){
		core <- 1
		#core <- length(unique(df$Group))
	}
	#browser()
	if(core>1) parallel <- TRUE
	if(parallel){
		# required variables (including functions) and packages are exported here
		cl <- createCluster(core=core, logfile = "/dev/null", export = c('rowcorMatVec', 'rcorrTest'), lib = c('plyr', 'Hmisc'))
		on.exit(stopCluster(cl))
	}
	#browser()
	## the working horse: function applied to each piece
	rowdoMatVec_piece <- function(x){
		ttmat <- t(x[, -(1:2)])
		vec <- x$Response
		#browser()
		eval(fn(ttmat, vec, ...), envir = parent.frame())
	}
	resL <- dlply(df, .(Group), function(x) rowdoMatVec_piece(x), .parallel=parallel)
	if(returnEarly){
		res <- list(resL=resL, mat=mat, vec=vec, group=group)
		return(res)
	}
	resArr <- laply(resL, function(x) data.matrix(x))
	dimnames(resArr)[[2]] <- rownames(mat)
	resArr <- swithArrayDim12(resArr) # so that genes * group * [cor, pval]
	#browser()
	if(returnArray) {
		res <- resArr
	} else {
		# return list
		res <- list(resL=resL, resArr=resArr, mat=mat, vec=vec, group=group)
	}	
	res
}
#tp <- rowdoMatVec_Bygroup(mat=t(icMatCCLE), vec=AXL_mRNA_CCLE, fn=rowcorMatVec, symbol=targetMap_CCLE_csvformat$targetOriginal, sepSymProbe=' | ', group=primarySite_CCLEGroupedBy20, meta=list(response='AXL mRNA', row='IC50', tag='Spearman correlation between IC50 and AXL mRNA'), addAll=TRUE, returnArray=FALSE, method='spearman')


#' the grouped version for rowcorMatVec
#'
#' This calculated the correlation between rows of mat and vec but at each group at a time as specified by the group vector
#' It should be more flexible than rowcorMatMat which requires the same genomic mat; instead, here different group of analysis can
#' have different sample size in genomics
#' @param mat a matrix, i.e. gene expression
#' @param vec a vector of phenotype, continuous variable
#' @param group a vector, categorical variable
#' @param addAll whether to calculate correlation by using all data ignoring group; this is the overall correlation
#' @param parallel deprecated; parallel is automatically enabled when core>1
#' @param addData logical indicating if data (mat) should be attached as an attribute (may be used for scatter plot)
#' @param core number of cores to use; when core>1, parallel execution is enabled
#' @param ... additional parameters to rowcorMatVec
#' @return a 3D array
#' @export
rowcorMatVec_Bygroup <- function(mat, vec, group, addAll=FALSE, addAllName='all', parallel=FALSE, core=1, 
	method='spearman', minN=5, usercorr=TRUE, returnArray=TRUE, symbol=NULL, sepSymProbe=', ', ...){
	#browser()
	df <- data.frame(Group=group, 
		Response=vec, t(mat), check.names=FALSE)
	# addAll
	if(addAll){
		df <- dfaugWithAll(df, groupVar='Group', addAllName='all')
	}
	if(missing(core)){
		core <- 1
		#core <- length(unique(df$Group))
	}
	#browser()
	if(core>1) parallel <- TRUE
	if(parallel){
		# required variables (including functions) and packages are exported here
		cl <- createCluster(core=core, logfile = "/dev/null", export = c('rowcorMatVec', 'rcorrTest'), lib = c('plyr', 'Hmisc'))
		on.exit(stopCluster(cl))
	}
	#browser()
	## the working horse: function applied to each piece
	rowcorMatVec_piece <- function(x){
		ttmat <- t(x[, -(1:2)])
		vec <- x$Response
		#browser()
		eval(rowcorMatVec(ttmat, vec, method=method, minN=minN, usercorr=usercorr, symbol=symbol, sepSymProbe=sepSymProbe, ...), envir = parent.frame())
	}
	resL <- dlply(df, .(Group), function(x) rowcorMatVec_piece(x), .parallel=parallel)
	resArr <- laply(resL, function(x) data.matrix(x))
	dimnames(resArr)[[2]] <- rownames(mat)
	resArr <- swithArrayDim12(resArr) # so that genes * group * [cor, pval]
	#browser()
	if(returnArray) {
		res <- resArr
	} else {
		# return list
		res <- list(resL=resL, resArr=resArr, mat=mat, vec=vec, group=group)
	}	
	res
}
#corWithREST_GarnettGroupBy20 <- rowcorMatVec_Bygroup(mat=enExprGarnett, vec=enExprGarnett['REST', ], group=primarySite_garnettGroupedBy20, addAll=TRUE, method='pearson')
if(FALSE){
# for a test
mat <- matrix(rnorm(48), nrow=3)
group <- rep(letters[1:3], each=16)
r <- rnorm(48)
tt <- rowcorMatVec_Bygroup(mat=mat, vec=r, group=group, addAll=TRUE, method='pearson', minN=10)
}


# SIBER has a bug: when a vector all equal, i.e. all are 0, it takes infinite time!

#' Compute BI for a gene expression matrix
#'
#' @param y a vector
#' @param ... additional parameters passed to SIBER
#' @return a vector
#' @export
SIBER2 <- function(y, ...){
	#if(all(y==y[1])){ # y[1] might be NA
	if(n_unique(y)<4){
		res <- rep(NA, 7)
		names(res) <- c("mu1", "mu2", "sigma1", "sigma2", "pi1", "delta", "BI")
	} else {
		tt <- system.time(res <- SIBER(y=y, ...))
        #if(tt[1]>300) stop('Here takes too long!\n') # check what happens for slow samples
        #browser()
	}
	res
}




#' fit SIBER on a matrix, one row at a time
#'
#' this function enables parallel through plyr package
#' 
#' @param mat matrix of expression, methylation and so on
#' @param model model as in SIBER
#' @param prune logical, whether to prune sigma1, sigma2 for with a specified percentile; this is to obtain stable BI estimate due to extremely small
#' sigma estimate, especially in the V model
#' @param parallel logical
#' @param core number of cores to register for parallel computing
#' @return a matrix
#' @export
fitSIBERonMat <- function(mat, model='NL', prune=TRUE, q=0.03, parallel=TRUE, core=10){
    require(SIBER)
    if(parallel){
		# required variables (including functions) and packages are exported here
		cl <- createCluster(core=core, logfile = "/dev/null", export='SIBER2', lib = c('SIBER'))
		on.exit(stopCluster(cl))
	}
	tmpRes <- adply(mat, 1, SIBER2, model=model, .parallel=parallel)
    res <- moveColumnToRowName(tmpRes)
	#browser()
	if(prune){
		# truncate sigma
		# res[, 'sigma1'] <- truncByQuantile(res[, 'sigma1'], q1=q, q2=1)
		# res[, 'sigma2'] <- truncByQuantile(res[, 'sigma2'], q1=q, q2=1)
		tpmin <- quantile(c(res[, 'sigma1'], res[, 'sigma2']), prob=q, na.rm=T)
		res[, 'sigma1'] <- truncByLimit(res[, 'sigma1'], Lower=tpmin, Upper=Inf)
		res[, 'sigma2'] <- truncByLimit(res[, 'sigma2'], Lower=tpmin, Upper=Inf)
		#browser()
		# update delta and BI
		res <- updateBImat(res)
	}
    res
}

# tt <- fitSIBERonMat(RNA2[1:200, ])

#' fit SIBER on matrix but split by group
#' 
#' this is written to facilitate pancancer analysis where BI is computed across tumors, thus, a group variable
#' 
#' @param mat matrix of expression, methylation and so on
#' @param group a vector indicating the grouping membership
#' @param addAllName the name for the all group
#' @param ... additional parameters to fitSIBERonMat such as model='NL', prune=TRUE, q=0.01, parallel=TRUE, core=10
#' @return a data frame of BI
#' @export
fitSIBERonMat_Bygroup <- function(mat, group, addAll=TRUE, addAllName='all', ...){
	df <- data.frame(Group=group, t(mat))
	# addAll
	if(addAll){
		df <- dfaugWithAll(df, groupVar='Group', addAllName='all')
	}
	## the working horse: function applied to each piece
	# paralell is within fitSIBERonMat, so sequencially for the groups but parallel on each task which is more efficient usage of the cpus
	f_piece <- function(x){
		#browser()
		data.matrix(moveColumnToRowName(fitSIBERonMat(t(x[, -1]), ...), column=1)) # make sure this is a matrix
	}
	#browser()
	res <- daply(df, .(Group), f_piece, .parallel=FALSE)
	### I disable the reshape since the BI info columns becomes alphabetical! Ugly!
	# reshape so that gene ~ BI info ~ group
	#tt <- melt(res)
	#res <- acast(tt, Var.2~Var.3~Group)
	#browser()
	### therefore, res is group ~ gene ~ BI info
	res <- unname_dimnames(res)
	res
}
#tt <- fitSIBERonMat_Bygroup(mat=RNA2[1:100, 1:200], group=ANNO$Tumor[1:200], addAll=TRUE)




#' lm for cont x and y
#'
#' when x is categorical, lm p value is equivalent to ANOVA p value. Coef becomes the mean; Rsquare?
#' @param x continuous vector
#' @param y continuous vector
#' @param collapse when vec is non-numeric which will be treated as factor, there are multiple coefs for different levels. We thus need to
#'		  summarize them into a scaler (by mean). if collapse=FALSE, all coefs will be returned. 
#' @return a named vector
lm_xy <- function(x, y, collapse=TRUE){
	temp <- data.frame(x=x, y=y)
	tpDF <- nonNumeric2FactorByColumn(temp) # so that categorical variable also works
	tpDF <- tpDF[complete.cases(tpDF),] # na.rm very important; otherwise, cannot use annova as sample size might differ
	fit_0 <- try(lm(y~1, data=tpDF), silent=TRUE) ## this may also fail, i.e. data is all NA, tpDF is 0*0
	fit_x <- try(lm(y~x, data=tpDF), silent=TRUE)
	if(class(fit_x)!='try-error' & class(fit_0)!='try-error' & length(unique(tpDF$x))>1){
		coef2 <- coef(summary(fit_x))
		#####
		##### WARNING: here is possible places for incompatibility! #####
		#####
		# only works for continuous xx: pval_raw <- coef2['xx', 'Pr(>|t|)'] 
		# for non-numeric (factor) x, records the mean of coefs.
		pval_raw <- anova(fit_0,fit_x)[2, 'Pr(>F)']
		adjRsq <- summary(fit_x)$adj.r.squared
		if(collapse){
			coef_raw <- mean(coef2[-1, 'Estimate'])
		} else {
			coef_raw <- coef2[-1, 'Estimate']
		}
	} else {
		pval_raw <- NA
		adjRsq <- NA
		if(class(tpDF$x)!='numeric' & class(fit_x)!='try-error'){
			coef_raw <- rep(NA, length(levels(tpDF$x))-1)
		} else {
			coef_raw <- NA # numeric x
		}
	}	
	#browser()
	res <- c(coef_raw, pval_raw, adjRsq)
	names(res) <- c('coef', 'pval', 'adjRsquare')
	res 
	#browser()
}
#with(df_info, lm_xy(x=packYear, y=NmutationInFM))

### rowlmMatVec may not work if NA is included and categorical
### no adjustment needed: bug solved; no needed
#quick_rowlmMatVec <- function(mat, vec, parallel=FALSE, core=12){
#
#}

#### we frequently need to fit lm for each row (continuous or categorical) of a mat with a fixed outcome (continuous or categorical), i.e. EMT score. 
#### We may also need to adjust for covariates (continuous or categorical).
#### I hate to write a loop each time I need to run this and examine the lm output. Thus I write this wrapper.
####
#### currently parallel is not implemented since I got error for variables not present------------> technical difficulty solved. now implemented on 11/10/2013
####
#### to deal with continuous, categorical (int) and character values uniformly, we use class(x) to detect if the value is numeric.
#### a) for mat, use class(mat[1,]) to decide if it is continuous.
#### b) for vec, use class(vec)
#### c) for adj, use colwise(class)(adj)
#### what p val?
#### y~1 ----
####        |
####     ----- pval_raw, coef_raw for x
####        |
#### y~x ----
#### ---------------------for mat=mutation, no adjutment for covariates, the result is equivalent to t test (ttestP same result) 
#### y~adj    ----
####             |
####             ----- pval_adj, coef_adj for x
####             |
#### y~adj+x  ----
#### if class is not numeric, the program will push the value into a factor. This is done by nonNumeric2FactorByColumn()
#### WARNING: adj in theory can be a multi-column data frame. However, it is not straightforward to decide what result we want to keep.
####		  thus, currently the result is prepared for adjustment with 0 or 1 variable. There is no guarantee if multiple variables will work.
#' Fit a linear model with rows of mat and vec, possibly adding an adjustment variable (adj)
#' 
#' @param mat a matrix
#' @param vec a vector, length of ncol(mat)
#' @param adj a vector for adjustment variable
#' @param matIsX by default, mat rows are continuous (i.e. expr) and vec is categorical response; but for CN data where mat is categorical and the response might be continuous or not.
#' 	in this case, we want lm(y~x) treats rows of mat as y and vec as x. This can be achieved by specifying matIsX=FALSE that makes rows of mat as the response variable
#' @param adj this can be a vector or a data frame that contains multiple variables to be adjusted. It needs to have nrow=length(vec).----> now demands a 1 column for easy handling
#' @param collapse when vec is non-numeric which will be treated as factor, there are multiple coefs for different levels. We thus need to
#'		  summarize them into a scaler (by mean). if collapse=FALSE, all coefs will be returned. 
#' @param core number of cores to use for parallel
#' @return a data frame
#' @export
rowlmMatVec <- function(mat, vec, adj, matIsX=TRUE, collapse=TRUE, parallel=FALSE, core=12) {
	if(parallel){
		cl <- makeCluster(core)
		registerDoParallel(cl)
		on.exit(stopCluster(cl))
	}
	if(ncol(mat)!=length(vec)) stop(sprintf("mat ncol: %d\nvec length: %d\n", nrow(mat), length(vec)))
	needAdj <- !missing(adj) # TRUE if adjustment is specified
	if(needAdj){
		adj <- data.frame(adj) # pushed into data frame
		if(nrow(adj)!=length(vec)) stop(sprintf("adj nrow: %d\nvec length: %d\n", nrow(adj), length(vec)))
	} else {
		#Error in get(s, env, inherits = FALSE) : 
		#	argument "adj" is missing, with no default
		# add this so that when paralell, this error will not show.
		# this might be a bug in foreach that tries to import every variable in the frame
		adj <- NULL # 
	}
	#browser()
	res <- foreach(r=1:nrow(mat), .combine='rbind', .packages=c('plyr'), .export=c('nonNumeric2FactorByColumn')) %dopar% { # nrow(mat)
		if(needAdj){
			#temp <- data.frame(y=vec, x=mat[r, ], adj=adj[,1]) # removed on 2014/02/13: we do expr~factor(group), not the other way!
			temp <- data.frame(x=vec, y=mat[r, ], adj=adj[,1])
		} else {
			# no adjusted variable
			#temp <- data.frame(y=vec, x=mat[r, ])
			temp <- data.frame(x=vec, y=mat[r, ])
		}
		#browser()
		tpDF <- nonNumeric2FactorByColumn(temp)
		if(!matIsX){
			# when mat is y, rename x as y and y as x
			colnames(tpDF) <- mapNames(colnames(tpDF), rbind(c('x', 'y'), c('y', 'x')))
		}
		#browser()
		tpDF <- tpDF[complete.cases(tpDF),] # na.rm very important; otherwise, cannot use annova as sample size might differ
		fit_0 <- try(lm(y~1, data=tpDF), silent=TRUE) ## this may also fail, i.e. data is all NA, tpDF is 0*0
		fit_x <- try(lm(y~x, data=tpDF), silent=TRUE)
		###
		if(class(fit_x)!='try-error' & class(fit_0)!='try-error' & length(unique(tpDF$x))>1){
		coef2 <- coef(summary(fit_x))
		#####
		##### WARNING: here is possible places for incompatibility! #####
		#####
		# only works for continuous xx: pval_raw <- coef2['xx', 'Pr(>|t|)'] 
		# for non-numeric (factor) x, records the mean of coefs.
		pval_raw <- anova(fit_0,fit_x)[2, 'Pr(>F)']
		if(collapse){
			coef_raw <- mean(coef2[-1, 'Estimate'])
		} else {
			coef_raw <- coef2[-1, 'Estimate']
		}
		} else {
		pval_raw <- NA
		if(class(tpDF$x)!='numeric' & class(fit_x)!='try-error'){
			coef_raw <- rep(NA, length(levels(tpDF$x))-1)
		} else {
			coef_raw <- NA # numeric x
		}
		}
		#browser()
		if(needAdj){
			fit_adj <- try(lm(y~adj, data=tpDF), silent=TRUE)
			fit_full <- try(lm(y~., data=tpDF), silent=TRUE)
			if(class(fit_adj)!='try-error' & class(fit_full)!='try-error' & length(unique(tpDF$x))>1){
			coef1 <- coef(summary(fit_full))
			#####
			##### WARNING: here is possible places for incompatibility! #####
			#####
			#---> only works for continuous: pval_adj <- coef1['xx', 'Pr(>|t|)'] 
			pval_adj <- anova(fit_adj,fit_full)[2, 'Pr(>F)'] ## H0: beta_CN1=beta_CN2=...=0
			# coef_adj becomes complicated to record, and with limited benefit. Thus remove it.
			#---> only works for continuous: for ordinal, multiple xx! coef_adj <- coef1['xx', 'Estimate'] 
			if(collapse){
			# x related coef is found by grep
			coef_adj <- mean(coef1[grep('^x',rownames(coef1)), 'Estimate'])
			} else {
			coef_adj <- coef1[coef1[grep('^x',rownames(coef1)), 'Estimate'], 'Estimate']
			}
			} else {
			pval_adj <- NA
			if(class(tpDF$adj)!='numeric'){
				coef_adj <- rep(NA, length(levels(tpDF$adj))-1)
				} else {
				coef_adj <- NA # numeric x
			}
			}
		}
		if(needAdj){
		rr <- c(pval_raw, pval_adj, coef_raw, coef_adj)
		names(rr) <- c('pval_raw', 'pval_adj', paste('coef_raw', 1:length(coef_raw), sep=''), paste('coef_adj', 1:length(coef_adj), sep=''))
		} else {
		rr <- c(pval_raw, coef_raw)
		names(rr) <- c('pval_raw', paste('coef_raw', 1:length(coef_raw), sep=''))
		}
		rr
		#browser()
	}
	#browser()
	rownames(res) <- rownames(mat)
	res
}
#lmCNWithREST_GBM <- rowlmMatVec(mat=mat_CN[, pID_used_CN], vec=restScore_GBMarray[pID_used_CN], parallel=TRUE, core=12) # ANOVA P
#lm_dr <- rowlmMatVec(mat=DataFile_JH[, indSel], vec=ClinicalFile_JH[indSel, 'X.treatment.'], parallel=FALSE)
# ttestMutWithREST_GBM <- rowlmMatVec(mat=mutGBM2[, pID_GBM_mutWithRNA2], vec=restScore_GBMarray[pID_GBM_mutWithRNA2])
# tt <- rowlmMatVec(mat=MIR, vec=EMTscore, adj=data.frame(ANNO[, 'histo_groupBy10']))
#ttestMutWithREST_GBM <- rowlmMatVec(mat=pancan12mutBinary[, pID_GBM_mutWithRNA], vec=restScore_GBMarray[pID_GBM_mutWithRNA])




##### Wilcoxon rank sum test 
##### different from t test, there is no ordering / sign here. 
#' wilcox P value for a vector of response and a vector of group
#' 
#' If both x and y are given and paired is FALSE, a Wilcoxon rank sum test (equivalent to the Mann-Whitney test: see the Note) is carried out.
#'
#' @param data a vector of y
#' @param class a vector of group (binary categorical)
#' @param alternative types of test
#' @param switchXY sometimes we have data as binary and class as continuous (for some reason, we do not want to change this in calling ttestP) and we need to swtich
#' @param ... additional parameters passed to wilcox.test
#' @return the P value
#' @export
wilcoxP <- function(data, class, switchXY=FALSE, alternative='two.sided', ...)
{
	if(switchXY){ # in case x and y needs to be switched
		tmp <- data
		data <- class
		class <- tmp
	}
	C <- unique(class[!is.na(class)])
	x <- data[which(class==C[1])]
	y <- data[which(class==C[2])]
	if(all(is.na(x), T) | all(is.na(y), T)) 
	{
		p <- NA
	} else {
		fit <- try(wilcox.test(x, y, alternative=alternative, ...), silent=TRUE)
		if(class(fit)!='try-error') {
		p <- fit $p.value
		} else {
			p <- NA
		}
	}
	#browser()
	p
}

### 09/09/2013: added na.action=na.omit to deal with missing
# deal with NA
# var.equal=TRUE: this is handy: when one group N=1, would still be able to compute P value!
# minN: minimum sample size in each group. Default is 2 (when not specified)
### 09/19/2013: add option more: when it is set to true, will return c(pval=pval, tstat=unname(tt$statistic), mDiff=unname(mDiff))

#' given a vector of data (continuous) and class (binary) and sampleID, prepared matched data by ordering for paired t test
#' @param data continuous response vector
#' @param class binary treatment vector
#' @param sampleID a vector of sample ID which can be used to automatically track the paired correspondance.
#' @return a list
prepDat4pairedttest <- function(data, class, sampleID=NULL, levels=NULL) {
	data <- as.numeric(data) # in case it is a data frame
	if(is.null(sampleID)) {
		# assume the data is properly ordered if sampleID is not specified
		res <- list(data=data, class=class, sampleID=sampleID, index=1:length(class))
	} else {
		if(length(sampleID)!=length(class))
			stop(sprintf('%d elemetns in class != %d elemetns in sampleID!', length(class), length(sampleID)))
		if(n_unique(sampleID)!=length(sampleID)/2) 
			stop(sprintf('sampleID is not perfect duplication: %d elements != 2* %d unique elements!\n', length(class), n_unique(sampleID)))
		if(n_unique(class)!=2) 
			stop(sprintf('%d groups observed in class: only 2 allowed!\n', n_unique(class)))
		tp <- data.frame(data=data, class=class, sampleID=sampleID, index=1:length(sampleID))		
		tp <- tp[order(tp$sampleID), ] # order by sampleID
		# the clean version: updating
		data <- tp$data
		class <- tp$class	
		sampleID <- tp$sampleID	
		index <- tp$index	
		res <- list(data=data, class=class, sampleID=sampleID, index=index)
	}
	levels <- with(res, 
		if(is.null(levels)){
		if(is.factor(class)){
			levels(class)
		} else {
			makestring(sort(unique(class), decreasing=FALSE))
		}
	})
	C <- levels
	#browser()
	tp <- with(res, data.frame(x=data[class==C[1]], y=data[class==C[2]]))
	colnames(tp) <- C
	res$df4track <- tp
	res
}
#datL <- prepDat4pairedttest(data=dat[1, indSel_GvD], class=pData$condition[indSel_GvD], sampleID=pData$CL[indSel_GvD])

# convert data frame used in trackPlot to the format as returned by prepDat4pairedttest which is easy for paired t test
df2trackL <- function(dat, getP=TRUE, levels=NULL, base=2){
	#browser()
	data <- as.numeric(c(dat[, 1], dat[, 2]))
	class <- rep(colnames(dat), each=nrow(dat))
	sampleID <- rep(rownames(dat), 2)
	index <- 1:length(class)
	res <- list(data=data, class=class, sampleID=sampleID, index=index)
	#browser()
	if(getP)
		res$res_pairedttest <- with(res, pairedttestP(data, class, sampleID, more=TRUE, base=base, levels=levels))
	res
}
#df2trackL(datL$df4track)

#' paired t-test
#'
#' this calculates paired t test p and other statistics
#'
#' We need to ensure samples are properly matched for paired data.
#' This bottles down to the binary treatment variable (class) and sample id to track the samples.
#' By default, if sampleID is not specified, it is assumed the samples are correctly matched (the clean case). This
#' actually means class starts with C[1] and then C[2] and no mixture allowed and the sample ID are matched.
#' A more flexible way is to specify class and sampleID simultaneously (in which case
#' class can be mixed) and the code will match class by sampleID through
#' order function. Of course, class and sampleID should have 1-to-1 correspondance. 
#'
#' var.equal has no effect on paired t test;
#' @param data continuous response vector
#' @param class binary treatment vector
#' @param sampleID a vector of sample ID which can be used to automatically track the paired correspondance.
#' @param levels the user can specify the required ordering c(C1, C2): mean diff is C2-C1. 
#' @return a vector
#' @export
pairedttestP <- function(data, class, sampleID=NULL, switchXY=FALSE, minN, more=TRUE, levels=NULL, base=2){
	if(switchXY){ # in case x and y needs to be switched
		tmp <- data
		data <- class
		class <- tmp
	}
	data <- as.numeric(data)
	# properly align data and class
	datL <- prepDat4pairedttest(data, class, sampleID)
	data <- datL$data
	class <- datL$class
	# unifying levels of the factor
	if(is.null(levels)){
		if(is.factor(class)){
			levels <- levels(class)
		} else {
			levels <- makestring(sort(unique(class), decreasing=FALSE))
		}
	}
	class <- gdata::reorder.factor(factor(class), new.order=levels)
	#C <- rev(levels(class)) 
	C <- levels(class) # so that C2-C1 is actually level2-level1
	if(length(C)!=2) 
		stop(sprintf('%d groups found; only exactly 2 groups allowed!\n', length(C)))
	x <- as.numeric(data[class==C[1]]) # group 1 values
	y <- as.numeric(data[class==C[2]]) # group 2 values
	#browser()
	if(missing(minN)) minN <- 2
	### when pval, tstat, mDiff are requested
	# if equal variance (var.equal=TRUE), N=1 in one group is also OK!
	#browser()
	tt <- try(t.test(x, y, na.action=na.omit, var.equal=var.equal, paired=TRUE), silent=TRUE)
	N1 <- sum(!is.na(x))
	N2 <- sum(!is.na(y))
	if(more) {
		## return NA if sample size is not as required or fit failure (i.e. exact same data in a group, sd=0)
		if(N1<minN | N2<minN | class(tt)=='try-error') {
			# tt can be try-error simply because group1 has same IC50 (sensitive) and group2 has same IC50 (resistant): this is not meaningful!
			return(c(pval=NA, tstat=NA, mDiff=NA, FC=NA, N1=NA, N2=NA))
		} else {
		# mean_C2-mean_C1
		# for mutation 0, 1 it is mean_mutated-mean_WT---> sometimes it C=1, 0
		# the reason is that unique(x) extracts unique elements in the order of their emergence: unique(c('d', 'a', 'b'))='d', 'a', 'b'!!!
		#mDiff <- diff(tt$estimate)
		mDiff <- -tt$estimate # default is C[1]-C[2]; we need to switch it
		pval <- tt$p.value
		#browser() # why mDiff has a wrong sign? it is because unique() records elements in the order of their emergence!!!
		res <- c(pval=pval, tstat=-unname(tt$statistic), mDiff=unname(mDiff), FC=toFC(unname(mDiff), base=base), N1=as.integer(N1), N2=as.integer(N2))
		## add info as what diff is computed
		#attr(res, 'info') <- paste('diff by class label: ', C[1], '-', C[2], sep=''): this is wrong! found on 2014/05/28
		#browser()
		attr(res, 'info') <- sprintf('%s\nBase used in calculating fold change is %.2f\n', paste('diff by class label: ', C[2], '-', C[1], sep=''), base)
		#browser()
		#return(res)
		}
	} else {
		if(N1<minN | N2<minN | class(tt)=='try-error') 
		{
		#return(NA)
		res <- NA
		} else {
		#return(tt$p.value)
		res <- tt$p.value
		}
	}
	attr(res, 'levels') <- C
	attr(res, 'base') <- base
	attr(res, 'minN') <- minN
	attr(res, 'switchXY') <- switchXY
	attr(res, 'test') <- 'paired t-test'
	res	
}
#pairedttestP(data=dat[1, indSel], class=pData$condition[indSel], sampleID=pData$CL[indSel])
#pairedttestP(data=dat[1, c(sid_paired_DMSO, sid_paired_G24)], class=rep(c('DMSO', 'G24'), each=length(sid_paired_DMSO)))
#ttestP(data=dat[1, c(sid_paired_DMSO, sid_paired_G24)], class=rep(c('DMSO', 'G24'), each=length(sid_paired_DMSO)), levels=c('DMSO', 'G24'), more=T)

#' Paired t test between rows of a matrix and a vector
#' 
#' @param mat a matrix, i.e. gene expression
#' @param vec a vector of phenotype, categorical variable
#' @param meta meta information, a list. This can specify what type of analysis (i.e. RPPA correlated with EMT score);
#'  row, which represents what row variables are (rows of mat); response represents what does vec mean
#' @return an object belonging to data frame and a class
#' @export
# interface to pairedttestP and rowpairedttestMatVec
#pairedttestP1 <- function(x) pairedttestP(data=x$expr, class=x$vec, ...)
rowpairedttestMatVec <- function(mat, vec, sampleID=NULL, matIsX=TRUE, 
	meta=list(tag=NULL, row='expression', response='response'), 
	minN=2, more=TRUE, levels=NULL, addData=TRUE, core=1, base=2, ...){
	#browser()
	if(core>1){
		# FORK: zombies created!
		cl <- makeCluster(core)
		registerDoParallel(cl)
		on.exit(stopCluster(cl))
	}
	if(ncol(mat)!=length(vec)) stop(sprintf('%d columns in mat, not equal to length of vec %d', ncol(mat), length(vec)))
	`%doforeach%` <- if(core>1) `%dopar%` else `%do%` # a concensus operation in case cannot find connection
	res <- foreach(i=1:nrow(mat), .combine='rbind', .export=c('pairedttestP', 'toFC'), .packages=c('gdata')) %doforeach% {
		res_rank <- wilcoxP(data=mat[i, ], class=vec, switchXY=!matIsX, paired=TRUE)
		rr <- pairedttestP(data=mat[i, ], class=vec, sampleID=sampleID, switchXY=!matIsX, minN=minN, more=more, levels=levels, base=base, ...)
		rr[7] <- res_rank
		names(rr)[7] <- 'pval_rankTest'
		rr
	}
	res <- data.frame(res)
	if(!is.null(rownames(mat)))
		rownames(res) <- rownames(mat)
	attr(res, 'meta') <- meta
	attr(res, 'levels') <- levels
	attr(res, 'base') <- base
	attr(res, 'minN') <- minN
	attr(res, 'matIsX') <- matIsX
	if(addData){
		attr(res, 'mat') <- mat
		attr(res, 'vec') <- vec
		attr(res, 'sampleID') <- sampleID
	}
	attr(res, 'class') <- c('rowpairedttestMatVec', 'data.frame')
	# genes * (group, FC, pval)
	res
}

#' Scatter plot for rowpairedttestMatVec class: each row produces a printted figure
#'
#' @method scatter rowpairedttestMatVec
#' @param x a rowpairedttestMatVec class (a data frame)
#' @param ... additional parameters to scatter_rowpairedttestMatVec
#' @export
scatter.rowpairedttestMatVec <- function(x, maxRow=Inf, main=NULL, ...){
	# when too many rows, select top rows based on P
	if(nrow(x)>maxRow){
		ind <- which(isTopN(-x[, 'pval'], N=maxRow))
		#ind <- orderVecByVal(vec=1:nrow(x), val=-x[, 'pval'])
		main <- sprintf('%s\nTop %d genes selected based on p value', main, maxRow)
	} else {
		ind <- order(x$pval, decreasing=FALSE)
	}
	# browser()
	# requires addData=TRUE
	if(is.null(Getter(x, 'mat'))) 
		stop('No mat (e.g. expression matrix) found! Please specify addData=TRUE\n')
	#browser()
	# this is dangerous: it is getting the global var; subsetting on ttest would not change mat
	mat <- getter(x, 'mat') 
	vec <- getter(x, 'vec')
	#
	# fix bug: need to also subset mat !
	#
	#browser()
	if(!identical(rownames(x), rownames(mat))){
		# only subsetting with shared genes!
		if(!all(rownames(x) %in% rownames(mat))){
			warning(sprintf('%d genes not in mat but in rowttestMatVec object!\n', length(setdiff(rownames(x), rownames(mat)))))
		}
		gg <- intersect(rownames(x), rownames(mat))
		mat <- mat[gg, ] # subsetting
	}
	scatter_rowttestMatVec(mat=mat, vec=vec, geneSel=ind, rowttestMatVec=x, ...)
} 

# scatter plot for paired t test
# sampleID and others are extracted
scatter_rowpairedttestMatVec <- function(mat, vec, FUN=NULL, geneSel,  
		rowpairedttestMatVec, anglex=NULL) {
	#browser()
	if(missing(geneSel)){
		geneSel <- order(rowpairedttestMatVec$pval, decreasing=FALSE)
	}
	if(is.null(anglex)) anglex <- guess_angle(catVecU=Getter(rowpairedttestMatVec, 'levels'))
	if(is.null(FUN)){
		# notice that xlab and ylab are from meta info
		FUN <- function(datL) {
			#browser()
			trackPlot(datL$df4track, ylab=str_c(datL$gene, datL$meta$row, sep=' '), xlab=datL$meta$response, 
				main=sprintf('%s\nPaired t test p: %.4e\nFC=%.3f', datL$gene, datL$testRes['pval'], datL$testRes['FC']), anglex=anglex)
		}
	}
	#browser()
	# vec is assumed to be a factor (rowttestMatVec); thus forced it here
	vec <- factor(vec, levels=Getter(rowpairedttestMatVec, 'levels'))
	meta <- Getter(rowpairedttestMatVec, 'meta') # meta info
	sampleID <- Getter(rowpairedttestMatVec, 'sampleID') 
	#browser()
	for(g in geneSel){
		#browser()
		gene <- rownames(mat)[g]
		testRes <- rowpairedttestMatVec[gene, ]
		datL <- prepDat4pairedttest(data=mat[g, ], class=vec, sampleID=sampleID)
		datL$meta <- meta
		datL$gene <- gene
		datL$testRes <- testRes
		tt=try(FUN(datL), silent=T)
		if(!inherits(tt, 'try-error')) {
			print(tt)
		}	
	}
}
#scatter_rowpairedttestMatVec(mat=dat[, indSel_GvD], vec=pData$condition[indSel_GvD], rowpairedttestMatVec=pairedttest_GvD)

#' Scatter plot for rowpairedttestMatVec class: each row produces a printted figure
#'
#' @method scatter rowpairedttestMatVec
#' @param x a rowpairedttestMatVec class (a data frame)
#' @param ... additional parameters to scatter_rowttestMatVec
#' @export
scatter.rowpairedttestMatVec <- function(x, ind=NULL, maxRow=Inf, ...){
	# order by p value
	if(is.null(ind)) ind <- order(x$pval, decreasing=FALSE)
	# browser()
	# requires addData=TRUE
	if(is.null(Getter(x, 'mat'))) 
		stop('No mat (e.g. expression matrix) found! Please specify addData=TRUE\n')
	#browser()
	# this is dangerous: it is getting the global var; subsetting on ttest would not change mat
	mat <- Getter(x, 'mat') 
	vec <- Getter(x, 'vec')
	#
	# fix bug: need to also subset mat !
	#
	#browser()
	if(!identical(rownames(x), rownames(mat))){
		# only subsetting with shared genes!
		if(!all(rownames(x) %in% rownames(mat))){
			warning(sprintf('%d genes not in mat but in rowttestMatVec object!\n', length(setdiff(rownames(x), rownames(mat)))))
		}
		gg <- intersect(rownames(x), rownames(mat))
		mat <- mat[gg, ] # subsetting mat
		ind <- ind[rownames(x)[ind] %in% gg] # subsetting ind
	}
	scatter_rowpairedttestMatVec(mat=mat, vec=vec, geneSel=ind, rowpairedttestMatVec=x, ...)
} 
#scatter.rowpairedttestMatVec(x=pairedttest_GvD, ind=1:3)

#' Bum plot for rowpairedttestMatVec class
#'
#' @method hist rowpairedttestMatVec
#' @param x a rowpairedttestMatVec class (a data frame)
#' @param ... additional parameters to plotBumFDR
#' @export
hist.rowpairedttestMatVec <- function(x, main=NULL, plot=TRUE, return=FALSE, alphas = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
    FDR=NULL, maskPthr = 1, pcutoff = 0.05, ...){
	obj_bum <- create_bum(pval=x[, 'pval'], alphas=alphas, FDR=FDR, maskPthr=maskPthr, pcutoff=pcutoff)
	if(is.null(main)){
		main <- sprintf('%s\n', getAnalysisName_ttest(x))
	}
	if(plot){
		plotBumFDR(obj_bum, main=main, ...)
	}
	if(return){
		return(obj_bum)
	}
}
#hist.rowpairedttestMatVec(pairedttest_GvD)

#' Compute t test P value
#'
#' @param data a vector of response
#' @param class a vector of membership, binary variable. Internally, class is converted to a vector during computation
#' @param switchXY sometimes we have data as binary and class as continuous (for some reason, we do not want to change this in calling ttestP) and we need to swtich
#' the role of x and y; switchXY when specified as TRUE, it will treat data as class and class as data. 
#' @param levels the user can specify the required ordering c(C1, C2): mean diff is C2-C1. 
#'  If not specified, (1) the class is a factor, levels would be extracted; (2) otherwise, use alphabetical order from class.
#' Notice that levels should be a vector of character strings. This means if
#'  true level is 0 and 1, levels should be specified as levels=c('0', '1') rather than levels=c(0, 1). This is due to the bug
#'  in gdata::reorder.factor that does not recognize new.order=c(0, 1). 
#' @param  more whether to include additional info i.e. FC, mean difference
#' @param minN minimum sample size (in each group) to carry out calculation
#' @param var.equal same as in t.test
#' @param base base used to calculate fold change using FC(); default is 2, which is typical for microarray and RNAseq data. 
#' @param whether use paired t test (assumes data is special being paired)
#' @param FCbyRatio sometimes the values are in original scale all positive, it is possible to compute the ratio and corresponding FC. 
#' 		  Notice, this will override the FC calculation already computed through specified base.
#' @return a vector
#' @export
ttestP <- function(data, class, switchXY=FALSE, minN, more=FALSE, var.equal=TRUE, levels=NULL, base=2, paired=FALSE, FCbyRatio=FALSE){
	### modified on 10/01/2013 to fix mutation data that has mDiff wrong sign when mutation data has 1 emerged first.
	###### this is not efficient!
	#C <- unique(class[!is.na(class)])
	#if(is.na(levels)){ # the user can specify the required ordering: C2-C1
	#	if(is.factor(class)){
	#		C <- reverse(levels(class)) # so that C1-C2 is actually level2-level1
	#	} else {
	#		C <- sort(unique(class[!is.na(class)]), decreasing=FALSE) # alphabetical
	#	}
	#}	
	## the following first examine if user specified levels is available; if yes, force the levels; otherwise, use alphabetical
	#browser()
	if(switchXY){ # in case x and y needs to be switched
		tmp <- data
		data <- class
		class <- tmp
	}
	# fix bug when levels is not specified
	# change levels of the factor
	#if(is.na(levels[1])){
	#	class <- factor(class)
	#} else {
	#	class <- gdata::reorder.factor(factor(class), new.order=levels)
	#}
	# unifying levels of the factor
	#browser()
	if(is.null(levels)){
		#class <- factor(class)
		if(is.factor(class)){
			levels <- levels(class)
		} else {
			levels <- makestring(sort(unique(class), decreasing=FALSE))
		}		
	}
	data <- as.numeric(data) # force to vector in case data is one row data frame!
	if(is.numeric(levels)) levels <- as.character(levels) # when 0, 1, reorder.factor will be problematic; so makes it string
	class <- gdata::reorder.factor(factor(class), new.order=levels)
	#C <- rev(levels(class)) 
	C <- levels(class) # so that C1-C2 is actually level2-level1
	if(length(C)!=2) {
		warning('Only one unique value observed in class! No computation will be performed!')
		return(NA)
	}
	x <- data[class==C[1]] # group 1 values
	y <- data[class==C[2]] # group 2 values
	if(paired & sum(class==C[1]) != sum(class==C[2])){
		stop('Specified paired t test but number of samples in each class is unequal!')
	}
	#browser()
	if(missing(minN)) minN <- 2
	### when pval, tstat, mDiff are requested
	# if equal variance (var.equal=TRUE), N=1 in one group is also OK!
	tt <- try(t.test(x, y, na.action=na.omit, var.equal=var.equal, paired=paired), silent=TRUE)
	# pooled sd for effect size (mDiff/sp)
	x_ <- noNA(x); y_ <- noNA(y)
	sp <- sqrt((var(x_)*(length(x_)-1)+var(y_)*(length(y_)-1))/(length(x_)+length(y_)-2)) # pooled standard deviation
	#
	N1 <- sum(!is.na(x))
	N2 <- sum(!is.na(y))
	#browser()
	if(more) {
		## return NA if sample size is not as required or fit failure (i.e. exact same data in a group, sd=0)
		if(N1<minN | N2<minN | class(tt)=='try-error') {
			# tt can be try-error simply because group1 has same IC50 (sensitive) and group2 has same IC50 (resistant): this is not meaningful!
			#res <- c(pval=pval, tstat=-unname(tt$statistic), means, mDiff=unname(mDiff), FC=toFC(unname(mDiff), base=base), N1=as.integer(N1), N2=as.integer(N2))
			return(c(pval=NA, tstat=NA, mean1=NA, mean2=NA, mDiff=NA, FC=NA, EffSize=NA, N1=as.integer(N1), N2=as.integer(N2)))
		} else {
		# mean_C2-mean_C1
		# for mutation 0, 1 it is mean_mutated-mean_WT---> sometimes it C=1, 0
		# the reason is that unique(x) extracts unique elements in the order of their emergence: unique(c('d', 'a', 'b'))='d', 'a', 'b'!!!
		#browser()
		if(paired){
			mDiff <- -tt$estimate # only has the difference now for paired data; this is due to inconsistency of t test function
			means <- c(mean(x, na.rm=TRUE), mean(y, na.rm=TRUE))
			names(means) <- c('mean1', 'mean2')
		} else {
			mDiff <- diff(tt$estimate) # this is C[2] - C[1]
			means <- tt$estimate
			names(means) <- c('mean1', 'mean2')
		}
		pval <- tt$p.value
		EffSize <- unname(mDiff/sp)
		#browser() # why mDiff has a wrong sign? it is because unique() records elements in the order of their emergence!!!
		#res <- c(pval=pval, tstat=unname(tt$statistic), mDiff=unname(mDiff), FC=toFC(unname(mDiff), base=base), N1=as.integer(N1), N2=as.integer(N2))
		# bug fix: since we switch mDiff as C[2]-C[1], we need to swith the sign of test statistic
		res <- c(pval=pval, tstat=-unname(tt$statistic), means, mDiff=unname(mDiff), FC=toFC(unname(mDiff), base=base), EffSize=EffSize, N1=as.integer(N1), N2=as.integer(N2))
		# force FC calculation by ratio
		if(FCbyRatio){
			res['FC'] <- toFC(log2(res['mean2']/res['mean1']), base=2)
		}
		## add info as what diff is computed
		#attr(res, 'info') <- paste('diff by class label: ', C[1], '-', C[2], sep=''): this is wrong! found on 2014/05/28
		attr(res, 'info') <- sprintf('%s\nBase used in calculating fold change is %.2f\n', paste('diff by class label: ', C[2], '-', C[1], sep=''), base)
		#browser()
		#return(res)
		}
	} else {
		if(N1<minN | N2<minN | class(tt)=='try-error') 
		{
		#return(NA)
		res <- NA
		} else {
		#return(tt$p.value)
		res <- tt$p.value
		}
	}
	#browser()
	attr(res, 'levels') <- C
	attr(res, 'base') <- base
	attr(res, 'var.equal') <- var.equal
	attr(res, 'minN') <- minN
	attr(res, 'switchXY') <- switchXY
	attr(res, 'FCbyRatio') <- FCbyRatio
	attr(res, 'test') <- 't-test'
	res	
}
#ttestP(data=EMTscore_pancancore, class=mutpancanCore[5, ], minN=2, more=TRUE)
#with(df, ttestP(RNA, KRAS, more=T))
#ttestP(RPPA[1, ind2RPPA], pheno_0$R_MUT[ind2core], levels=c('0', '1'), more=TRUE)

#ttestP(data=ic50_comb[56, ], class=BRAFstat_factor, more=TRUE)
#ttestP(data=EMTscore_pancancore, class=mutpancanCore[gSel2percent[1], ], more=TRUE)


#' compute rank test P as well as a t test result
#'
#' this uses the wilcoxP() and ttestP(); the interface is similar to rowranktestMatVec.
#' @param mat a matrix with rows as genes
#' @param vec a vector giving the group info
#' @param core number of cores to use
#' @param ... additional pars for rowranktestMatVec
#' @return a data frame
#' @export
rowranktestMatVec <- function(mat, vec, core=12, ...){
	parallel <- ifelse(core>1, TRUE, FALSE) # this overwrides the parallel option
	if(parallel){
		# FORK: zombies created!
		cl <- makeCluster(core) #ttestP not found: PSOCK, MPI
		registerDoParallel(cl)
		on.exit(stopCluster(cl))
	}
	#browser()
	`%doforeach%` <- if(core>1) `%dopar%` else `%do%` # a concensus operation in case cannot find connection
	res <- foreach(i=1:nrow(mat), .combine='rbind', .export=c('ttestP', 'toFC', 'wilcoxP'), .packages=c('gdata')) %doforeach% {
		res_t <- ttestP(data=mat[i, ], class=vec, more=TRUE, ...)
		res_rank <- wilcoxP(data=mat[i, ], class=vec)
		rr <- c(res_rank, res_t[c(1, 3)])
		names(rr) <- c('pval', 'pval_ttest', 'mDiff')
		rr
	}
	if(!is.null(rownames(mat)))
		rownames(res) <- rownames(mat)
	res
}
#ranktest_MUT_MET <- rowranktestMatVec(mat=MET[1:10, ], vec=pheno_0$R_MUT, levels=c('0', '1'), core=1) 
#

#' given expression mean difference, transforms back to fold change: negative means down-regulate
#'
#' @param vec a vector of difference
#' @param base base for transformation
#' @return a vector
#' @export
toFC <- function(vec, base=2) {
	if(base!=1) {
		val <- base^abs(vec)
		res <- val*sign(vec)
	} else {
		# if base==1, this means mDiff is good; there is no ratio since values can be negative; not meaningful to comptute FC
		res <- rep(NA, length(vec))
	}	
	res
}

#' formatting a given vector according to given levels
#' @param vec a vector
#' @param levels a vector of levels specified for vec
#' @return a factor
#' @export
formatLevel <- function(vec, levels=NULL){
	# compute levels
	if(is.null(levels)){
		if(is.factor(vec)){
			levels <- levels(vec)
		} else {
			levels <- makestring(sort(unique(vec), decreasing=FALSE)) # alphabetical
		}		
	} else {
		# in case some level is not observed in vec
		levels <- levels[levels %in% makestring(unique(vec))]
	}
	#browser()
	res <- gdata::reorder.factor(factor(vec), new.order=as.character(levels)) # reorder.factor requries character input
	res
}
#' Perform log rank test between rows of mat (categorical) and a vector of OS and OScensoring
#' 
#' Notice that ggplot2 figures are attached; this may takes a long time to cache the objects in knitr
#'
#' @param mat a matrix with rows as genes
#' @param OS a vector giving the survival time. This can be for example overall survival, PFS and others
#' @param OScensoring a vector giving the censor status. 1 is event and 0 is censored.
#' @param core number of cores to use; this override parallel=TRUE if core=1 (no parallel)
#' @param addData logical indicating if data (mat) should be attached as an attribute (may be used for scatter plot)
#' @param ... additional parameters to plotKM
#' @param meta meta information, a list. This can specify what type of analysis (i.e. KRAS vs NF1 mutation);
#'  row, which represents what row variables are (rows of mat); response represents what does OS mean
#' @param more not implemented
#' @param levels levels specified for each row of mat
#' @return a rowkmMatVec class; meta information is stored as attributes
#' @export
rowkmMatVec <- function(mat, OS, OScensoring,  
	meta=list(tag=NULL, row='group', response='Overal survival'), 
	more=TRUE, levels=NULL, parallel=FALSE, core=1, addData=TRUE, ...){
	parallel <- ifelse(core>1, TRUE, FALSE) # this overwrides the parallel option
	#
	# parallel rowttestMatVec seems not helpful at all!
	#
	if(parallel){
		# FORK: zombies created!
		cl <- makeCluster(core)
		registerDoParallel(cl)
		on.exit(stopCluster(cl))
	}
	if(ncol(mat)!=length(OS)) stop(sprintf('%d columns in mat, not equal to length of OS %d', ncol(mat), length(OS)))
	# deprecated: ttestP can handle all these
	# change levels of the factor
	#if(is.na(levels[1])){
	#	vec <- factor(vec)
	#} else {
	#	vec <- gdata::reorder.factor(factor(vec), new.order=levels)
	#}
	#C <- rev(levels(vec)) 
	#levels <- levels(vec) # levels[2]-levels[1] for mean diff
	#browser()
	`%doforeach%` <- if(core>1) `%dopar%` else `%do%` # a concensus operation in case cannot find connection
	resL <- foreach(i=1:nrow(mat), .export=c('ttestP', 'toFC'), .packages=c('gdata')) %doforeach% {
		#browser()
		# as.numeric is very important!
		#browser()
		gene <- rownames(mat)[i]
		ptag <- ifelse(is.null(meta$tag), gene, sprintf('%s: %s', gene, meta$tag))
		xx <- formatLevel(vec=mat[i, ], levels=levels)
		pRes <- try(plotKM(x=xx, OS=OS, OScensoring=OScensoring, xlab=meta$response, tag=ptag, legendTitle=gene, plot=FALSE), silent=TRUE)
		if(!inherits(pRes, 'try-error')){
			coxfitGrouped <- pRes$coxfitGrouped	
			zzL <- list(Pval=coxfitGrouped$P, OR=coxfitGrouped$coef_exp, CI_L=coxfitGrouped$ciL, CI_U=coxfitGrouped$ciU, N1=pRes$N1, N2=pRes$N2, p=pRes$p)
		} else {
			xxtab <- table(xx)
			zzL <- list(Pval=NA, OR=NA, CI_L=NA, CI_U=NA, N1=xxtab[1], N2=xxtab[2], p=NULL)
		}
		# this is a list
		zzL
	}
	names(resL) <- rownames(mat)
	res <- ldply(resL, .fun=function(x) c(pval=x$Pval, HR=x$OR, CI_L=x$CI_L, CI_U=x$CI_U, N1=x$N1, N2=x$N2))[, -1]
	pList <- lapply(resL, function(x) x$p)
	if(!is.null(levels)){
		Nnames <- str_c('N_', levels, sep='')
	} else {
		Nnames <- c('N1', 'N2')
	}
	colnames(res) <- c('pval', 'HR', 'CI_L', 'CI_U', Nnames)
	rownames(res) <- rownames(mat)
	#browser()
	# add meta info
	attr(res, 'meta') <- meta
	attr(res, 'levels') <- levels
	attr(res, 'pList') <- pList
	if(addData){
		attr(res, 'mat') <- mat
		attr(res, 'OS') <- OS
		attr(res, 'OScensoring') <- OScensoring
	}
	attr(res, 'class') <- c('rowkmMatVec', 'data.frame')
	res	
}
#res_km <- rowkmMatVec(mat=mutmat_binary, OS=clin$OS, OScensoring=clin$OScensoring, meta=list(tag=NULL, row='group', response='Overal survival'), more=TRUE, levels=c('0', '1'), parallel=FALSE, core=1, addData=TRUE)


#' Scatter plot for rowkmMatVec class: each row produces a printted figure
#'
#' @method scatter rowkmMatVec
#' @param x a rowkmMatVec class (a data frame)
#' @export
scatter.rowkmMatVec <- function(x, maxRow=500, main=NULL, ...){
	# when too many rows, select top rows based on P
	if(nrow(x)>maxRow){
		ind <- which(isTopN(-x[, 'pval'], N=maxRow))
		main <- sprintf('%s\nTop %d genes selected based on p value', main, maxRow)
	} else {
		ind <- order(x$pval, decreasing=FALSE)
	}
	#browser()
	pList <- Getter(x, 'pList')
	for(i in ind){
		print(pList[[i]])
	}	
} 

# based on plotHeatmapTable to extract fisher or chisq p value

#' compute fisher exact test (chisq test when failure) p value for two categorical variable
#' @param x categorical var 1
#' @param y categorical var 2
#' @param isBinary if both are binary variable 0, 1, we can compute OR to determine direction of association
#' @export
getfisherP <- function(x, y, isBinary=FALSE){
	tab <- table(y, x)
	mtab <- melt(tab)
	colnames(mtab) <- c("x", "y", "count") # for melt data, x is x from table(x, y)
	mtab <- mutate(mtab, x=as.factor(x), y=as.factor(y)) # in case 0/1 will be interpretated as continous, giving 0.5 in the axis
	fitP <- try(fisher.test(tab)$p.value, silent=TRUE)
	#browser()
	if(class(fitP)=='try-error'){ 
		# add: min(dim(tab)); otherwise, 2-by-1 table would be given a value very significant!
		if(min(dim(tab))>1){
			fitP <- try(chisq.test(tab)$p.value, silent=TRUE)
			test <- 'Chi-squre'
			if(class(fitP)=='try-error'){
				fitP <- NA
			}	
		} else { # not enough data to perform the calculation
			fitP <- NA
			test <- 'Fisher/chisq exact' # tried both
		}
		
	} else {
		test <- 'Fisher exact'
	}
	res <- list(pval=fitP, test=test)
	if(isBinary) {
		#browser()
		if(tab['1', '0']*tab['0', '1']!=0){
			OR <- tab['1', '1']*tab['0', '0']/(tab['1', '0']*tab['0', '1'])
		} else {
			OR <- NA
		}
		res$OR <- OR
	}
	return(res)
}

#' Perform Fisher Exact test between rows of mat (categorical) and a factor
#' 
#' This will perform fisher exact test; if this fails (e.g. due to large dimention or counts), chisq test will be used to get the result.
#' Originally, ggplot2 figures are attached; this may takes a long time to cache the objects in knitr as well as saving it. Now ggplot2 objects are not
#' computed and thus not returned to improve efficiency
#'
#' @param mat a matrix with rows as genes
#' @param vec a factor giving the response
#' @param core number of cores to use; this override parallel=TRUE if core=1 (no parallel)
#' @param addData logical indicating if data (mat) should be attached as an attribute (may be used for scatter plot)
#' @param ... additional parameters to plotKM
#' @param meta meta information, a list. This can specify what type of analysis (i.e. KRAS vs NF1 mutation);
#'  row, which represents what row variables are (rows of mat); response represents what does OS mean
#' @param more not implemented
#' @param levels levels specified for the response vec; rows of mat are formatted with as.factor currently; (this is not tested if it controls levels for plotting although it does not affect p value calculation)
#' @return a rowfisherMatVec class; meta information is stored as attributes
#' @export
rowfisherMatVec <- function(mat, vec,  
	meta=list(tag=NULL, row='group', response='response'), 
	more=TRUE, levels=NULL, parallel=FALSE, core=1, addData=TRUE, ...){
	parallel <- ifelse(core>1, TRUE, FALSE) # this overwrides the parallel option
	#
	# parallel rowttestMatVec seems not helpful at all!
	#
	if(parallel){
		# FORK: zombies created!
		cl <- makeCluster(core)
		registerDoParallel(cl)
		on.exit(stopCluster(cl))
	}
	if(ncol(mat)!=length(vec)) stop(sprintf('%d columns in mat, not equal to length of vec %d', ncol(mat), length(vec)))
	`%doforeach%` <- if(core>1) `%dopar%` else `%do%` # a concensus operation in case cannot find connection
	resL <- foreach(i=1:nrow(mat), .export=c('ttestP', 'toFC'), .packages=c('gdata', 'GenAnalysis')) %doforeach% { # 1:nrow(mat)
		#browser()
		# as.numeric is very important!
		#browser()
		gene <- rownames(mat)[i]
		ptag <- ifelse(is.null(meta$tag), gene, sprintf('%s\n%s', gene, meta$tag))
		yy <- formatLevel(vec=vec, levels=levels)
		resfit <- getfisherP(x=as.factor(mat[i, ]), y=yy)
		#browser()
		# the resulting ggplot2 objects are too expensive to save
		#pRes <- try(plotScatter(x=as.factor(mat[i, ]), y=yy, xlab=meta$row, ylab=meta$response, tag=ptag, plot=FALSE), silent=TRUE)
		#if(!inherits(pRes, 'try-error')){
		#	zzL <- list(pval=pRes$Pval, test=pRes$test, p=pRes$p)
		#} else {
		#	zzL <- list(pval=NA, test=NA, p=NULL)
		#}
		zzL <- list(pval=resfit$pval, test=resfit$test)
		# this is a list
		zzL
	}
	names(resL) <- rownames(mat)
	res <- ldply(resL, .fun=function(x) data.frame(pval=x$pval, testName=x$test))[, -1]
	#pList <- lapply(resL, function(x) x$p)
	colnames(res) <- c('pval', 'testName')
	rownames(res) <- rownames(mat)
	#browser()
	# add meta info
	attr(res, 'meta') <- meta
	attr(res, 'levels') <- levels
	#attr(res, 'pList') <- pList
	#attr(res, 'class') <- 'rowttestMatVec'
	if(addData){
		attr(res, 'mat') <- mat
		attr(res, 'vec') <- vec
	}
	attr(res, 'class') <- c('rowfisherMatVec', 'data.frame')
	res	
}
#anovaRes_cn_KRASmut <- rowfisherMatVec(mat=cn[sample(1:10000, 100), ], vec=mut['KRAS', ], meta=list(tag='association between CN and KRAS mutation', row='Copy number', response='KRAS mutation'), more=TRUE, levels=NULL, parallel=FALSE, core=1, addData=TRUE)

#' Scatter plot for rowfisherMatVec class: each row produces a printted figure
#'
#' @method scatter rowfisherMatVec
#' @param x a rowfisherMatVec class (a data frame)
#' @export
scatter.rowfisherMatVec <- function(x, maxRow=500, main=NULL, ...){
	# when too many rows, select top rows based on P
	if(nrow(x)>maxRow){
		ind0 <- which(isTopN(-x[, 'pval'], N=maxRow))
		ind <- ind0[order(x[ind0, 'pval'], decreasing=FALSE)]
		main <- sprintf('%s\nTop %d genes selected based on p value', main, maxRow)
	} else {
		ind <- order(x$pval, decreasing=FALSE)
	}
	#browser()
	vec <- Getter(x, 'vec')
	levels <- Getter(x, 'levels')
	meta <- Getter(x, 'meta')
	mat0 <- Getter(x, 'mat')
	mat <- mat0[rownames(x), ]
	for(i in ind){
		#browser()
		yy <- formatLevel(vec=vec, levels=levels)
		gene <- rownames(mat)[i]
		ptag <- ifelse(is.null(meta$tag), gene, sprintf('%s\n%s', gene, meta$tag))
		pRes <- try(plotScatter(x=as.factor(mat[i, ]), y=yy, xlab=meta$row, ylab=meta$response, tag=ptag, plot=FALSE), silent=TRUE)
		#browser()
		if(!inherits(pRes, 'try-error')) print(pRes$p)
	}	
} 

#' Bum plot for rowfisherMatVec class
#'
#' @method hist rowfisherMatVec
#' @param x a rowfisherMatVec class (a data frame)
#' @param ... additional parameters to plotBumFDR
#' @export
hist.rowfisherMatVec <- function(x, main=NULL, plot=TRUE, return=FALSE, alphas = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
    FDR=NULL, maskPthr = 1, pcutoff = 0.05, ...){
	obj_bum <- create_bum(pval=x[, 'pval'], alphas=alphas, FDR=FDR, maskPthr=maskPthr, pcutoff=pcutoff)
	if(is.null(main)){
		main <- sprintf('%s\n', getAnalysisName(x))
	}
	if(plot){
		plotBumFDR(obj_bum, main=main, ...)
	}
	if(return){
		return(obj_bum)
	}
}


#' Perform t test between rows of mat and a vector; Notice paired t test also works by setting paired=TRUE
#' 
#' @param mat a matrix with rows as genes
#' @param vec a vector giving the group info. Internally, vec is converted to a vector by ttestP for computation
#' @param matIsX set matIsX=FALSE if each row of mat is binary (binary response vector)
#'   while set matIsX=TRUE if each row of mat is continuous. By default, mat rows are continuous (i.e. expr) and vec is binary; but for MUT data where mat is binary and the response might be continuous. In this case, we need to 
#' specify matIsX=FALSE which is passed to switchXY as TRUE in t test. By default, rows in mat is not response (matIsX=TRUE) and hence, switchXY is disabled. 
#' @param core number of cores to use; this override parallel=TRUE if core=1 (no parallel)
#' @param base base parameter passed to ttestP
#' @param addData logical indicating if data (mat) should be attached as an attribute (may be used for scatter plot)
#' @param ... additional parameters to ttestP
#' @param meta meta information, a list. This can specify what type of analysis (i.e. KRAS vs NF1 mutation);
#'  row, which represents what row variables are (rows of mat); response represents what does vec mean
#' @param symbol in case this is probe level data, user can specify a vector of symbols matching 
#' @param sepSymProbe separator for symbol-probe. That is, symbol+sepSymProbe+probe(rowname) is the new rowname
#' @param whether use paired t test (assumes data is special being paired)
#' @return a rowttestMatVec class; meta information is stored as attributes
#' @export
rowttestMatVec <- function(mat, vec, matIsX=TRUE, 
	meta=list(tag=NULL, row='expression', response='response'), 
	minN=2, more=TRUE, var.equal=TRUE, levels=NULL, parallel=FALSE, base=2, core=1, addData=TRUE, 
	symbol=NULL, sepSymProbe=', ', paired=FALSE, 
	FDR=NULL, pcutoff = 0.05, alphas = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), ...){
	parallel <- ifelse(core>1, TRUE, FALSE) # this overwrides the parallel option
	#
	# parallel rowttestMatVec seems not helpful at all!
	#
	## update rownames to deal with probe level data
	mat <- updateRowname(mat=mat, symbol=symbol, sepSymProbe=sepSymProbe)
	if(parallel){
		# FORK: zombies created!
		cl <- makeCluster(core)
		registerDoParallel(cl)
		on.exit(stopCluster(cl))
	}
	if(ncol(mat)!=length(vec)) stop(sprintf('%d columns in mat, not equal to length of vec %d', ncol(mat), length(vec)))
	#if(is.null(levels)) levels <- sort(unique(vec), decreasing=FALSE)
	if(is.null(levels)){
		if(is.factor(vec)){
			levels <- levels(vec)
		} else {
			levels <- makestring(sort(unique(vec), decreasing=FALSE))
		}		
	}
	# force levels in case vec has degenerated levels (in its level but not present)
	# this cannot run if matIsX=FALSE
	if(matIsX) vec <- reorder(factor(vec), new.order=levels)
	# deprecated: ttestP can handle all these
	# change levels of the factor
	#if(is.na(levels[1])){
	#	vec <- factor(vec)
	#} else {
	#	vec <- gdata::reorder.factor(factor(vec), new.order=levels)
	#}
	#C <- rev(levels(vec)) 
	#levels <- levels(vec) # levels[2]-levels[1] for mean diff
	`%doforeach%` <- if(core>1) `%dopar%` else `%do%` # a concensus operation in case cannot find connection
	#browser()
	res <- foreach(i=1:nrow(mat), .combine='rbind', .export=c('ttestP', 'toFC'), .packages=c('gdata')) %doforeach% {
		#browser()
		# as.numeric is very important!
		res_rank <- wilcoxP(data=mat[i, ], class=vec, switchXY=!matIsX, paired=paired)
		rr <- ttestP(data=mat[i, ], class=vec, switchXY=!matIsX, minN=minN, more=more, var.equal=var.equal, levels=levels, base=base, paired=paired, ...)
		rr['pval_rankTest'] <- res_rank
		#browser()
		rr
	}
	#browser()
	res <- data.frame(res)
	if(!is.null(rownames(mat)))
		rownames(res) <- rownames(mat)
	# add levels in case it is NULL
	if(is.null(levels)){
		levels <- Getter(ttestP(data=as.numeric(mat[i, ]), class=vec, switchXY=!matIsX, minN=minN, more=more, var.equal=var.equal, levels=levels, base=base, ...), 'levels')
	}
	# add meta info
	test <- ifelse(paired, 'paired t test', 't test')
	obj_bum <- try(create_bum(res$pval, FDR=FDR, pcutoff=pcutoff, alphas=alphas), silent=TRUE)
	attr(res, 'meta') <- meta
	attr(res, 'levels') <- levels
	attr(res, 'paired') <- paired
	attr(res, 'test') <- test
	attr(res, 'base') <- base
	attr(res, 'var.equal') <- var.equal
	attr(res, 'minN') <- minN
	attr(res, 'matIsX') <- matIsX
	attr(res, 'obj_bum') <- obj_bum
	#attr(res, 'class') <- 'rowttestMatVec'
	if(addData){
		attr(res, 'mat') <- mat
		attr(res, 'vec') <- vec
	}
	attr(res, 'class') <- c('rowttestMatVec', 'data.frame')
	#class(res) <- c('rowttestMatVec', 'data.frame')
	#browser()
	res	
}
#' Heatmap plot for rowttestMatVec class
#'
#' @method heatm rowttestMatVec
#' @param x a rowttestMatVec class (a data frame)
#' @param gSel string specifying what genes are selected in the elements of 'obj_bum'
#' @param revlog10 sometimes vec is log10 IC50; it is better to show IC50 in original scale in the heatmap; this will be possible by specifying revlog10=TRUE
#' @export
heatm.rowttestMatVec <- function(x, colpalvec=NULL, plot=TRUE, return=TRUE, gSel='gSelectInd_P05', colbarname='vec', 
	revlog10=FALSE, plotLegend='TRUE', ...){
	#browser()
	obj_bum <- Getter(x, 'obj_bum')
	vec <- Getter(x, 'vec')
	if(revlog10) {
		vec <- 10^vec
		code_vec <- sprintf("vec <- 10^Getter(x, 'vec')")
	} else {
		code_vec <- sprintf("vec <- Getter(x, 'vec')")
	}
	mat <- Getter(x, 'mat')
	code0 <- sprintf("\tx <- %s\n\tgSel <- '%s'", as.character(substitute(x)), as.character(substitute(gSel)))
	#browser()
	code1=sprintf("
	obj_bum <- Getter(x, 'obj_bum')
	%s
	mat <- Getter(x, 'mat')
	gSel1 <- obj_bum[[gSel]]
	pheno <- data.frame(vec=vec)
	colnames(pheno) <- '%s'
	rownames(pheno) <- colnames(mat)
	colpal <- list(%s=%s)
	phm <- aheatmat(mat=mat, sSel1=which(!is.na(pheno$%s)), gSel=gSel1, 
    	pheno=pheno, q1=0.05, q2=0.95, clustering_distance_rows = 'pearson', 
		cluster_cols=F, colOrderIndex=order(pheno$%s, decreasing=F), 
		cluster_rows=T, 
   		clustering_method = 'ward.D', yd=0.2, colpal=colpal, plotLegend=%s, ...)
	\n", code_vec, colbarname, colbarname, colpalvec, colbarname, colbarname, plotLegend)
	code <- sprintf('%s%s', code0, code1)
	res <- list(code=code)
	#browser()
	if(plot){
		#pander::evals(code)
		#browser()
		phm <- eval(parse(text=code))
		res$phm <- phm
		mapL <- phm
		stretchL <- llply(mapL, get_stretch)
		res$stretchL <- stretchL
		# obj_bum <- Getter(x, 'obj_bum')
		# vec <- Getter(x, 'vec')
		# mat <- Getter(x, 'mat')
		# gSel1 <- obj_bum[[gSel]]
		# pheno <- data.frame(vec=vec)
		# rownames(pheno) <- colnames(mat)
		# colpal <- list(vec=colorpalette2colvec('greenred'))
		# phm <- aheatmat(mat=mat, sSel1=which(!is.na(pheno$vec)), gSel=gSel1, 
  # 		  	pheno=pheno, q1=0.05, q2=0.95, clustering_distance_rows = 'pearson', 
		# 	cluster_cols=F, colOrderIndex=order(pheno$vec, decreasing=F), 
		# 	cluster_rows=T, 
  # 	 		clustering_method = 'ward.D', yd=0.2, colpal=colpal, plotLegend=T)	
	}
	if(return){
		cat(code)
		return(res)
	}
}

#' get method for class rowttestMatVec
#' 
#' this function is named getter, not get to minimize conflict with the usual get function
#'
#' @param obj an rowttestMatVec object
#' @param what specify what to get from augMat, i.e. base, var.equal
#' @method getter rowttestMatVec
#' @return an object as requested
#' @export
getter.rowttestMatVec <- function(obj, what='meta'){
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

#' ordered horizontal bar plot for rowttestMatVec
#'
#'
#' @param rowInfo alternative rownames to show on the figure
#' @param decreasing logic, decreasing order of bars or not.
#' @param dataMode a string used to label the x-axis. Default is 'IC50' for IC50 t-test; Please use meaning string relevant to the data. 
#' @param showP do not show P value on graph ('None'); show actual P ('Value'); show P category ('Category')
#' @param pthr filtering by p value. Rows with larger P value will be deleted. Default is to use all rows (pthr=1)
#' @param legend.position position for the legend if P value is show as category
#' @param sizex x axis tick text size
#' @param sizey y axis tick text size
#' @param size.title.x x label size
#' @param size.title.y y label size
#' @param legend.title title of the legend
#' @param truncFClow cutoff to truncate FC lower than this value. Defalt is -Inf, which means no truncation on small FC (negative) value
#' @param truncFChigh cutoff to truncate FC higher than this value. Defalt is Inf, which means no truncation on large FC (usually positive) value
#' @param barWidth relative barwidth (0~1)
#' @param useFC whether to use FC or mDiff
#' @return a ggplot2 object
#' @export
plot4rowttestMatVec <- function(res, useFC=FALSE, rowInfo=NULL, decreasing=TRUE, dataMode='IC50', ylab=NULL, xlab=NULL, sizex=12, sizey, sizeyadj=1, size.title.x=14, size.title.y=14, main=NULL, 
	showP='None', pthr=1, colorSignif=TRUE, legend.position=c(.85,0.85), legend.title='Significance',
	truncFClow=-4, truncFChigh=4, barWidth=0.8, verbose=FALSE, fillCol=c('#ED281E', '#2834A4', '#EDB51E', 'gray'),
	coord_flip=TRUE, anglex=NULL, angley=NULL, plot=TRUE){
	levels <- attributes(res)$levels
	if(is.null(ylab)) ylab <- ''
	if(is.null(main)) main <- ''
	res <- data.frame(res)
	#browser()
	### sometimes NA for pval: this leads to error
	res <- subset(res, !is.na(pval))
	if(useFC){
		res$plotValue <- res$FC
		if(is.null(xlab)) xlab <- sprintf('Fold change of mean %s (%s/%s)', dataMode, levels[2], levels[1])
	} else {
		res$plotValue <- res$mDiff
		if(is.null(xlab)) xlab <- sprintf('Difference of mean %s (%s - %s)', dataMode, levels[2], levels[1])
	}
	toplot <- ifelse(nrow(res)>2, TRUE, FALSE) # may be all NA's and thus no figure
	if(toplot){
	#browser()
	# truncate FC 
	res <- mutate(res, plotValue=truncByLimit(plotValue, Lower=truncFClow, Upper=truncFChigh))
	# add P category
	res$pcat <- categorizePvec(res[, 'pval'])
	#
	if(is.null(rowInfo)) rowInfo <- rownames(res)
	# showP: add P into drug info
	if(showP=='Value'){
		rowInfo <- mapply(function(x, y) sprintf('%s, P=%.3g', x, y), rowInfo, res$pval)
	} else if(showP=='Category'){
		rowInfo <- mapply(function(x, y) sprintf('%s, %s', x, y), rowInfo, res$pcat)
	}
	#browser()
	# format res
	res <- data.frame(res, rowInfo=rowInfo) # add rowInfo
	res <- subset(res, !is.na(plotValue) & pval<pthr) # remove rows with FC NA (not enough sample size)
	res$rowOrdered <- orderVecByVal(vec=res$rowInfo, val=res$plotValue, decreasing=decreasing) # levels added
	if(missing(sizey)) sizey <- guessSize(nrow(res))
	sizey <- sizey*sizeyadj
	scale_axis_size <- theme(axis.text.x=element_text(size=sizex, colour = "black", face ='bold', angle=anglex), 
		axis.text.y=element_text(colour = "black", face ='bold', angle=angley, size=sizey),
		axis.title.x=element_text(size=size.title.x, colour = "black"), axis.title.y=element_text(size=size.title.y, colour = "black"))
	#browser()
	if(verbose) cat(sprintf('sizex=%.2f; sizey=%.2f\n', sizex, sizey))
	#fillCol <- c('red', 'magenta', 'purple', 'blue')
	names(fillCol) <- c('P<0.001', 'P<0.01', 'P<0.05', 'P>0.05')
	fillColVec <- fillCol[sort(unique(as.character(res$pcat)), decreasing=FALSE)]
	#browser()
	if(!colorSignif) {
	p <- ggplot(res, aes(y=plotValue, x=rowOrdered))+geom_bar(stat = "identity", position="dodge", width=barWidth)+ 
		ylab(xlab)+xlab(ylab)+ggtitle(main)+scale_axis_size+
		theme(legend.position=legend.position)
	} else {
	p <- ggplot(res, aes(y=plotValue, x=rowOrdered, fill=pcat))+geom_bar(stat = "identity", position="dodge", width=barWidth)+ 
		scale_fill_manual(name=legend.title, breaks=names(fillColVec), values=fillColVec)+
		ylab(xlab)+xlab(ylab)+ggtitle(main)+
		scale_axis_size+theme(legend.position=legend.position)
	}	
	if(coord_flip)
		p <- p+coord_flip()
	#browser()
	#return(p)
	} else {
		plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', main='No figure')
	}
	if(!plot){
		return(list(p=p, data=res[order(res$rowOrdered, decreasing=TRUE), ]))
	} else {
		return(p)
	}
	
}
#plot4rowttestMatVec(ttestRes_NF1mut_B2, dataMode='IC50', main='NF1 MT vs WT', truncFClow=-4, truncFChigh=4)

#plot4rowttestMatVec(ttestRes_NF1, rowInfo=drugNameShown, showP='None', main='NF1 mutation', legend.position=c(.2,0.15))

#' Generic scatter method
#'
#' @param x an object
#' @param ... other arguments
#' @export
scatter <- function(x, ...)
{
    UseMethod("scatter",x)
}

#' Generic heatm method
#'
#' @param x an object
#' @param ... other arguments
#' @export
heatm <- function(x, ...)
{
    UseMethod("heatm",x)
}

#' Generic venndiagram method
#'
#' @param x an object
#' @param ... other arguments
#' @export
venndiagram <- function(x, ...)
{
    UseMethod("venndiagram",x)
}


#' Scatter plot for rowttestMatVec
#' 
#' 
#' @param mat expression matrix used in rowttestMatVec 
#' @param vec vec as used in rowttestMatVec 
#' @param FUN plotting function
#' @param geneSel index for plotting; default is plotting all sorted by P
#' @param rowttestMatVec rowttestMatVec object
#' @param type type to be passed to plotContCat
#' @param maintagfunc a function to compute additional information for main 
#' @export
scatter_rowttestMatVec <- function(mat, vec, FUN, geneSel, rowttestMatVec, sizeForce=NULL, anglex=NULL, ylim=NULL, 
	matIsX=TRUE, paired=FALSE, addRankP=FALSE, type='jitter', col=NULL, shape=18, maintagfunc=NULL, errorDotPlot=FALSE){
	if(missing(geneSel)){
		geneSel <- order(rowttestMatVec$pval, decreasing=FALSE)
	}
	#browser()
	if(is.null(anglex)) anglex <- guess_angle(catVecU=Getter(rowttestMatVec, 'levels'))
	#browser()
	if(missing(FUN)){
		# notice that xlab and ylab are from meta info
		FUN <- function(expr, vec, gene, meta, main, matIsX) {
			#browser()
			ylab <- ifelse(matIsX, str_c(gene, meta$row, sep=', '), meta$response)
			xlab <- ifelse(matIsX, meta$response, str_c(gene, meta$row, sep=', '))
			plotScatter(vec, as.numeric(expr), type=type, col=col, shape=shape, labels=colnames(mat), ylab=ylab, xlab=xlab, main=main, ylim=ylim, anglex=anglex, sizeForce=sizeForce)
		}
		FUNpaired <- function(dat, vec, gene, meta, main, matIsX) {
			# when paired, expr and vec are aligned by samples (e.g. cell lines) 1-to-1
			#browser()
			ylab <- ifelse(matIsX, str_c(gene, meta$row, sep=', '), meta$response)
			xlab <- ifelse(matIsX, meta$response, str_c(gene, meta$row, sep=', '))
			trackPlot(dat, xlab=xlab, ylab=ylab, main=main)
		}
		# plot SEM and mean + jitter as KLKP paper from Figure 7D
		FUN_errorDotPlot <- function(expr, vec, gene, meta, main, matIsX){
			ylab <- ifelse(matIsX, str_c(gene, meta$row, sep=', '), meta$response)
			xlab <- ifelse(matIsX, meta$response, str_c(gene, meta$row, sep=', '))
			errorDotPlot(x=vec, y=expr, ylab=ylab, xlab=xlab, main=main)
		}
		if(errorDotPlot){
			FUN <- FUN_errorDotPlot # override boxplot with errorDotPlot
		}
	}
	#browser()
	meta <- getter(rowttestMatVec, 'meta') # meta info
	test <- getter(rowttestMatVec, 'test') # t test or paired t test
	levels <- getter(rowttestMatVec, 'levels') 
	#paired <- getter(rowttestMatVec, 'paired') # meta info
	#browser()
	for(g in geneSel){
		#browser()
		gene <- rownames(mat)[g]
		tp <- rowttestMatVec[gene, ]
		if(!is.na(tp['FC'])){
			main0 <- sprintf('%s\n%s p=%.3g; fold change=%.3g\n', gene, test, tp['pval'], tp['FC'])
		} else {
			main0 <- sprintf('%s\n%s p=%.3g; mean difference=%.3g\n', gene, test, tp['pval'], tp['mDiff'])
		}
		if(addRankP){
			main0 <- sprintf('%sMann-Whitney rank p=%.3g', main0, tp['pval_rankTest'])
		}
		#browser()
		if(!is.null(maintagfunc)) {
			tag0 <- maintagfunc(gene)
			main0 <- sprintf('%s%s', main0, tag0)
		}
		if(matIsX){
			# vec is group; mat row is continuous
			# vec is assumed to be a factor (rowttestMatVec); thus forced it here
			gg <- factor(vec, levels=getter(rowttestMatVec, 'levels'))	
			expr <- as.numeric(mat[g, ])
		} else {
			# vec is coninuous; mat row is group
			expr <- as.numeric(vec)
			gg <- factor(mat[g, ], levels=getter(rowttestMatVec, 'levels'))	
		}
		#browser()
		if(paired==FALSE){
			tt=try(FUN(expr=expr, vec=gg, gene=gene, meta=meta, main=main0, matIsX=matIsX), silent=T)
		} else {
			dat <- data.frame(cbind(expr[vec %in% levels[1]], expr[vec %in% levels[2]]), check.names=F)
			colnames(dat) <- levels
			tt=try(FUNpaired(dat, vec=gg, gene=gene, meta=meta, main=main0, matIsX=matIsX), silent=T)
		}
		if(!inherits(tt, 'try-error')) {
			print(tt)
		}	
	}
}
#
# scatter_rowttestMatVec(mat=ic50_aug, vec=pData[, 'TTF1_halfcut'], rowttestMatVec=ttest_TTF1)


#' Scatter plot for rowttestMatVec class: each row produces a printted figure
#'
#' @method scatter rowttestMatVec
#' @param x a rowttestMatVec class (a data frame)
#' @param addRankP whether to add rank test p value in the scatter annotation
#' @param maintagfunc a function to compute additional information for main
#' @param ... additional parameters to scatter_rowttestMatVec
#' @export
scatter.rowttestMatVec <- function(x, maxRow=Inf, main=NULL, maintagfunc=NULL, addRankP=FALSE, type='jitter', col=NULL, shape=18, ...){
	# when too many rows, select top rows based on P
	if(nrow(x)>maxRow){
		ind <- which(isTopN(-x[, 'pval'], N=maxRow))
		#ind <- orderVecByVal(vec=1:nrow(x), val=-x[, 'pval'])
		ind <- ind[order(x$pval[ind], decreasing=FALSE)]
		main <- sprintf('%s\nTop %d genes selected based on p value', main, maxRow)
	} else {
		ind <- order(x$pval, decreasing=FALSE)
	}
	#browser()
	# requires addData=TRUE
	if(is.null(getter(x, 'mat'))) 
		stop('No mat (e.g. expression matrix) found! Please specify addData=TRUE\n')
	#browser()
	# this is dangerous: it is getting the global var; subsetting on ttest would not change mat
	mat <- getter(x, 'mat') 
	vec <- getter(x, 'vec')
	paired <- getter(x, 'paired')
	matIsX <- getter(x, 'matIsX')
	#
	# fix bug: need to also subset mat !
	#
	#browser()
	if(!identical(rownames(x), rownames(mat))){
		# only subsetting with shared genes!
		if(!all(rownames(x) %in% rownames(mat))){
			warning(sprintf('%d genes not in mat but in rowttestMatVec object!\n', length(setdiff(rownames(x), rownames(mat)))))
		}
		gg <- intersect(rownames(x), rownames(mat))
		mat <- mat[gg, ] # subsetting
	}
	#browser()
	scatter_rowttestMatVec(mat=mat, vec=vec, geneSel=ind, rowttestMatVec=x, matIsX=matIsX, paired=paired, addRankP=addRankP, type=type, col=col, shape=shape, maintagfunc=maintagfunc, ...)
} 
#scatter.rowttestMatVec(obj, maxRow=Inf)
#scatter.rowttestMatVec(ttest_TTF1)

# get analysis name for rowttest and rowcor
getAnalysisName <- function(obj){
	meta <- Getter(obj, 'meta')
	matIsX <- Getter(obj, 'matIsX')
	if(is.null(matIsX)) matIsX <- TRUE # no matIsX
	#browser()
	response <- ifelse(matIsX, meta$response, meta$row) # mat row might be the binary response
	row <- ifelse(matIsX, meta$row, meta$response) # mat row might be the binary response
	if(!is.null(meta$tag)){
		res <- capitalize(meta$tag)
	} else {
		#res <- sprintf('Association between %s and %s', meta$response, meta$row)
		res <- capitalize(sprintf('%s-%s', response, row))
		#browser()
	}
	res
}
# get levels info
getAnalysisName_ttest <- function(obj){
	meta <- Getter(obj, 'meta')
	levels <- Getter(obj, 'levels')
	matIsX <- Getter(obj, 'matIsX')
	response <- ifelse(matIsX, meta$response, meta$row) # mat row might be the binary response
	#browser()
	if(!is.null(meta$tag)){
		#res <- capitalize(sprintf('%s for %s: %s vs %s', meta$tag, meta$row, levels[2], levels[1]))
		res <- capitalize(sprintf('%s\n %s: %s vs %s', meta$tag, response, levels[2], levels[1]))
		#res <- capitalize(sprintf('%s', meta$tag))
	} else {
		res <- capitalize(sprintf('%s', meta$tag))
	}
	res
}

#' Bar plot for rowttestMatVec class
#' @method plot rowttestMatVec
#' @param x a rowttestMatVec class (a data frame)
#' @param fillCol color for bar filling in the order of c('P<0.001', 'P<0.01', 'P<0.05', 'P>0.05') 
#' @param ... additional parameters to plot4rowttestMatVec
#' @export
plot.rowttestMatVec <- function(x, maxRow=400, main=NULL, tag='', dataMode=NULL, fillCol=c('#ED281E', '#2834A4', '#EDB51E', 'gray'), ...){
	# when too many rows, select top rows based on P
	#browser()
	mainSpecified <- !is.null(main)
	if(!mainSpecified) main <- sprintf('%s%s', getAnalysisName(x), tag)
	if(nrow(x)>maxRow){
		ind <- which(isTopN(-x[, 'pval'], N=maxRow))
		if(!mainSpecified)
			main <- sprintf('%s%s\nTop %d genes selected based on p value', main, tag, maxRow)
	} else {
		ind <- 1:nrow(x)
	}
	#browser()
	base <- Getter(x, 'base')
	matIsX <- Getter(x, 'matIsX')
	useFC <- ifelse(base==1, FALSE, TRUE) # base==1, use mDiff
	if(is.null(dataMode)) dataMode <- ifelse(matIsX, Getter(x, 'meta')$row, Getter(x, 'meta')$response)
	#browser()
	plot4rowttestMatVec(x[ind, ], useFC=useFC, main=main, dataMode=dataMode, fillCol=fillCol, ...)
} 

#' Bum plot for rowttestMatVec class
#'
#' @method hist rowttestMatVec
#' @param x a rowttestMatVec class (a data frame)
#' @param ... additional parameters to plotBumFDR
#' @export
hist.rowttestMatVec <- function(x, main=NULL, tag='', plot=TRUE, return=FALSE, alphas = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
    FDR=NULL, maskPthr = 1, pcutoff = 0.05, ...){
	obj_bum <- create_bum(pval=x[, 'pval'], alphas=alphas, FDR=FDR, maskPthr=maskPthr, pcutoff=pcutoff)
	#browser()
	if(is.null(main)){
		main <- sprintf('%s%s', getAnalysisName_ttest(x), tag)
	}
	if(plot){
		plotBumFDR(obj_bum, main=main, ...)
	}
	if(return){
		return(obj_bum)
	}
}
#hist.rowttestMatVec(ttestlist_mutRNA[[1]])
#ttest_MUT_RPPA <- rowttestMatVec(mat=RPPA[, ind2RPPA], vec=pheno_0$R_MUT[ind2core], levels=c('0', '1')) 

#tt <- rowttestMatVec(mat=ic50_moonshot, vec=BRAFstat_factor, minN=2, more=TRUE, var.equal=TRUE, levels=NA, parallel=TRUE, core=12)


####################################################################################################################
############## rowMatVecDo ANOVA
####################################################################################################################

#' perform operations between rows of a mat and a vector
#'
#' This should be a general approach for rowcorMatVec rowlmMatVec rowranktestMatVec rowttestMatVec
#' In principle, we want to apply arbitrary functions and return a data frame.
#' Currently this is implemented through ddply; later we might use dply to speed up computation
#' plyr can parallel; but dplyr cannot.
#' @param mat expression matrix
#' @param vec a vector of response
#' @param fn a function applied to the data frame with columns gene, expr, vec
#' @param ... additional parameters to fn
#' @return a data frame
#' @export
rowMatVecDo <- function(mat, vec, fn, ...) {
	dat <- rowMatVec2df(mat=mat, vec=vec)#[1:33*2, ]
	# anovaP may return different length due to level degeneration!
	#ll <- dlply(dat, .(gene), function(x)x)
	#z <- foreach(x= ll) %do% {
	#	 with(x, anovaP(expr=expr, label=vec, more=TRUE))
	#}
	#browser()
	zz <- dlply(dat, .(gene), fn, ...)
	# in case group degeneration, augment all groups
	IDaug0 <- noNA(unionID(llply(zz, function(x) names(x))))
	# this may include padj from 2 level t test res; thus remove this to avoid error; also remove FC and NA due to base=1
	IDaug <- IDaug0[str_detect(IDaug0, 'anovaP|_')]
	if(length(IDaug)==0){
		# fn_fisherP_comutation only has Pval and OR
		IDaug <- IDaug0	
	} 
	#browser()
	res <- foreach(z=zz, .combine='rbind') %do% {
		#browser()
		if(is.null(names(z))){
			# when z is all NA, no names!
			ID2 <- IDaug[1:length(z)]
		} else {
			ID2 <- names(z)
		}
		indwname <- !is.na(names(z)) # make sure name is not NA
		rr <- fillVecByID(ID1=IDaug, ID2=ID2[indwname], value=z[indwname])
		rr
	}
	res <- data.frame(res)
	rownames(res) <- names(zz)
	colnames(res) <- IDaug # make sure names are like:  FC_Squamous-NSCLC, not  FC_Squamous.NSCLC since this - means minus
	#browser()
	res
	# dplyr not easy to use
	#res <- dat %>%
	#	group_by(gene) %>% 
	#		do(polr=polrP(expr=.$expr, label=.$vec))
	#do(function(x) polrP(expr=x$expr, label=x$vec))
}
#anovaRes_pathology <- rowanovaMatVec(mat=ic50_aug, vec=pData$pathology, meta=list(tag='ANOVA', row='IC50', response='pathology'))

#PolrP
fn_polrP <- function(x) polrP(expr=x$expr, label=x$vec)
fn_anovaP <- function(x, matIsX=TRUE, ...) {
	if(matIsX){
		res <- anovaP(expr=x$expr, label=x$vec, ...)
	} else {
		# expr is group ; ve is response
		res <- anovaP(expr=x$vec, label=as.factor(x$expr), ...)
	}	
	#browser()
	res
}


#' two way anova test for a given response and two vectors
#'
#' @param varLabel the label for the 3 effects: main effect of primary variable, main effect of variable not of interest (confounding)
#' 	and interaction effect 
#' @param interaction whether to add interaction term
#' @return a fn_TWanova object with plot method
#' @export
fn_TWanova <- function(y, mainVec, confounderVec, interaction=TRUE, varLabel=c('Treatment', 'Cell line', 'Interaction')){
	Data <- data.frame(y=y, mainVec=mainVec, confounderVec=confounderVec)
	lData <- dlply(Data, .(confounderVec), function(x) x) # the cell line are ordered alphabatically; 
	# full model result containing: main effect, modifier effect, interaction
	res_full <- c('pval_main'=NA, 'pval_confounder'=NA, 'pval_interaction'=NA)
	# for each confounderGroup, anovaP (or t test) value for main effect: p value already adjusted
	#anovaP_stratified <- ldply(lData, function(x) anovaP(x$y, x$mainVec)['anovaP'])
	#res_tratified <- anovaP_stratified[, 2]
	#names(res_tratified) <- str_c('anovaP|', anovaP_stratified[, 1], sep='')
	# 20170418: modified to include not just anovaP, but all other FC, mdiff values
	anova_stratified <- ldply(lData, function(x) anovaP(x$y, x$mainVec))
	mm <- melt(anova_stratified)
	res_tratified <- mm$value
	names(res_tratified) <- with(mm, str_c(variable, confounderVec, sep='|'))
	#browser()
	# may add Hotelling p value as combining effects across cell lines
	if(interaction){
		model <- try(lm(y ~ mainVec + confounderVec + mainVec:confounderVec, data=Data), silent=TRUE)
		#model2 <- try(lm(y ~ mainVec * confounderVec, data=Data), silent=TRUE) # identical
		if(!inherits(model, 'try-error')){
			aov <- anova(model)
			res_full['pval_main'] <- aov['mainVec', 'Pr(>F)']
			res_full['pval_confounder'] <- aov['confounderVec', 'Pr(>F)']
			res_full['pval_interaction'] <- aov['mainVec:confounderVec', 'Pr(>F)']
		}
		
	} else {
		model <- try(lm(y ~ mainVec + confounderVec, data=Data), silent=TRUE)
		if(!inherits(model, 'try-error')){
			aov <- anova(model)
			res_full['pval_main'] <- aov['mainVec', 'Pr(>F)']
			res_full['pval_confounder'] <- aov['confounderVec', 'Pr(>F)']
		}
	}
	res <- c(res_full, res_tratified)
	#browser()
	attr(res, 'pval_full') <- res_full
	attr(res, 'res_tratified') <- res_tratified
	pval_stratified <- res_tratified[str_detect(names(res_tratified), 'anovaP')]
	#browser()
	attr(res, 'pval_stratified') <- pval_stratified
	#attr(res, 'groupsConditionedOn') <- anovaP_stratified[, 1]
	attr(res, 'groupsConditionedOn') <- anova_stratified[, 1]
	attr(res, 'Data') <- Data
	attr(res, 'varLabel') <- varLabel
	attr(res, 'class') <- 'fn_TWanova'
	res
}
#zz <- fn_TWanova(y=mymat1[i, ], mainVec=myanno1$Treatment2, confounderVec=myanno1$CellLine, interaction=F)

#' Plot method for fits from two-way data by fn_TWanova
#'
#' This in the title shows the p values for main effects (two variables namely treatment and confounder) and interaction
#' The figure is boxplot colored by confounder; at the bottome, stratified anova p value is showing (for each color)
#'
#' @param legendLabelSize legend label size
#' @param cols a named color vector where names are the unique group categories. Default color is specified with sorted group category.
#' @export
plot.fn_TWanova <- function(x, pshape=17, psize=3.5, tag='', cols=NULL, main=NULL, xlab=NULL, ylab='Expression', 
	legendLabelSize=11, anglex=NULL, sizex=12, sizey=12, angley=0, plot=T){
	Data <- Getter(x, 'Data')
	varLabel <- Getter(x, 'varLabel')
	pval_stratified <- Getter(x, 'pval_stratified')
	groupsConditionedOn <- Getter(x, 'groupsConditionedOn')
	groupsLabel <- foreach(i=1:length(pval_stratified), .combine='c') %do% {
		sprintf('%s, p=%.3g', groupsConditionedOn[i], pval_stratified[i])
	}
	groupsLabelOrder <- order(pval_stratified, decreasing=FALSE)
	ugs <- sort(unique(groupsConditionedOn), decreasing=F) # unique groups sorted
	if(is.null(cols)){
		cols <- colorpalette2colvec('set1')[1:n_unique(ugs)]
		names(cols) <- ugs
	}
	groupsShown <- groupsConditionedOn[groupsLabelOrder] # the order of legend labels sorted by p
	cols_shown <- cols[groupsShown] # corresponding color so that the same group would have the same color across different figures (usually in rowTWanovaMatVec )
	#ggplot(Data, aes(x=mainVec, y=y, col=confounderVec))+geom_point(size=psize, shape=pshape)
	Data$confounderVec <- factor(Data$confounderVec, levels=groupsConditionedOn[groupsLabelOrder]) # reorder level so legend is also ordered
	if(is.null(main)) {
		gtitle <- sprintf('%s%s: p=%.3g\n%s: p=%.3g, %s: p=%.3g\n', tag, varLabel[1], x['pval_main'], varLabel[2], x['pval_confounder'], varLabel[3], x['pval_interaction'])
	} else {
		gtitle <- main
	}
	if(is.null(xlab)) xlab <- varLabel[1]
	if(is.null(anglex)) {
		#browser()
		anglex <- ifelse(sum(nchar(as.character(unique(Data$mainVec))))>16, 45, 0)
	}
	#browser()
	scale_axis_size <- theme(axis.text.x=element_text(size=sizex, colour = "black", face ='bold', angle=anglex, hjust=0.5, vjust=0.5), 
		axis.text.y=element_text(colour = "black", face ='bold', size=sizey, angle=angley))
	set.seed(1000) # due to jitter, needs to be consistent
	p <- ggplot(Data, aes(x=mainVec, y=y, col=confounderVec))+
		geom_jitter(position=position_jitter(w=0.2, h=0), size=psize, shape=pshape)+
		scale_colour_manual(name="",values=cols_shown, labels=groupsLabel[groupsLabelOrder])+
		guides(color = guide_legend(nrow = ceiling(length(groupsLabel) / 2), byrow=TRUE))+
		theme(legend.position="bottom", legend.text = element_text(size = legendLabelSize))+ggtitle(gtitle)+
		xlab(xlab)+ylab(ylab)+scale_axis_size
	res <- list(p=p, fit=x)
	if(plot){
		return(p)
	} else {
		return (res)
	}

}
#scatter(res_twanova)
#plot(zz)

#' Two-way ANOVA between rows of a matrix and two vectors: two main effects and then interaction effect
#'  
#' Each row is to be fitted by fn_TWanova()
#'
#' @param mat a matrix, i.e. gene expression
#' @param vec a vector of phenotype, categorical variable
#' @param matIsX not implemented at this moment
#' @param meta meta information, a list. This can specify what type of analysis (i.e. RPPA correlated with EMT score);
#'  row, which represents what row variables are (rows of mat); response represents what does vec mean
#' @return an object belonging to data frame and a class
#' @export
rowTWanovaMatVec <- function(mat, vec1, vec2, matIsX=TRUE,  interaction=TRUE,
	meta=list(tag='ANOVA', row='expression', response='response'), 
	base=2, more=TRUE, addData=TRUE, ...){
	#browser()
	#res <- rowMatVecDo(mat=mat, vec=vec, fn=fn_anovaP, base=base, more=more)
	# no rowMatVecDo as vec is not a parameter
	resL <- foreach(i=1:nrow(mat)) %do% {
		fn_TWanova(y=mat[i, ], mainVec=vec1, confounderVec=vec2, interaction=interaction)
	}
	# caution: only works if no missing data so that no collapse of groups 
	if(n_unique(sapply(resL, length))>1) 
		stop(sprintf('fn_TWanova() returns different number of elements for various rows!'))
	res <- ldply(resL, function(x) x) # are the stratified anovaP consistent, e.g. cell liens matched?; would be a problem if missing data presents. 
	#browser()
	rownames(res) <- names(resL) <- rownames(mat)
	#res <- rowMatVecDo(mat=mat, vec=vec, fn=fn_TWanova, mainVec=vec1, confounderVec=vec2, interaction=interaction, base=base, more=more)
	attr(res, 'meta') <- meta
	attr(res, 'base') <- base
	if(addData){
		attr(res, 'mat') <- mat
		attr(res, 'vec1') <- vec1
		attr(res, 'vec2') <- vec2
	}
	attr(res, 'class') <- c('rowTWanovaMatVec', 'data.frame')
	attr(res, 'matIsX') <- matIsX 
	attr(res, 'resL') <- resL 
	res
}
#zz = rowTWanovaMatVec(mat=mymat1, vec1=myanno1$Treatment2, vec2=myanno1$CellLine, interaction=TRUE)

hist_rowTWanovaMatVec0 <- function(x, var='pval_main', varLabel='Treatment effect p value', main=NULL, plot=TRUE, return=FALSE, alphas = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
    FDR=NULL, maskPthr = 1, pcutoff = 0.05, ...){
	obj_bum <- create_bum(pval=x[, var], alphas=alphas, FDR=FDR, maskPthr=maskPthr, pcutoff=pcutoff)
	#browser()
	if(is.null(main)){
		main <- sprintf('%s\n', varLabel)
	}
	if(plot){
		plotBumFDR(obj_bum, main=main, ...)
	}
	if(return){
		return(obj_bum)
	}
}

#' Bum plot for rowTWanovaMatVec class
#'
#' @method hist rowTWanovaMatVec
#' @param x a rowTWanovaMatVec class (a data frame)
#' @param var what variable (pval_main, pval_confounder, pval_interaction, anovaP|HCC1806, ...) for bum plot; if NULL, all are plotted
#' @param varLabel label to be used for annotation purpose describing var
#' @param ... additional parameters to plotBumFDR
#' @export
hist.rowTWanovaMatVec <- function(x, var='pval_main', varLabel='Treatment effect p value', main=NULL, plot=TRUE, return=FALSE, alphas = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
    FDR=NULL, maskPthr = 1, pcutoff = 0.05){
	if(!is.null(var)){
		res <- hist_rowTWanovaMatVec0(x=x, var=var, varLabel=varLabel, main=main, plot=plot, return=return, alphas=alphas, FDR=FDR, maskPthr = maskPthr, pcutoff = pcutoff)
	} else {
		# plot all
		#browser()
		pvalName_all <- colnames(x)[str_detect(colnames(x), 'pval|anovaP')]
		res <- foreach(var=pvalName_all) %do% {
			rr <- hist_rowTWanovaMatVec0(x=x, var=var, varLabel=var, main=main, plot=plot, return=return, alphas=alphas, FDR=FDR, maskPthr = maskPthr, pcutoff = pcutoff)
			rr	
		}
	}
	res
}
#zz = hist(res_twanova, var=NULL)

#' venndiagram method for rowTWanovaMatVec class
#'
#' When number of cell lines is smaller than 5, venn2 would be used; otherwise, imageVenn would be used to
#' display genes selected at a given method (as selected using create_bum)
#' @method venndiagram rowTWanovaMatVec
#' @param x a rowTWanovaMatVec class (a data frame)
#' @param cat.pos Vector giving the position (in degrees) of each category name along the circle, with 0 at 12 o'clock, being on top
#' @param cat.dist Vector giving the distance (in npc units) of each category name from the edge of the circle (can be negative)
#' @param use either 'venn2' or 'imageVenn'; default is NULL, in which case 'venn2' would be selected if less or equal to 5 cell lines
#' @param ... additional parameters to plot function
#' @export
venndiagram.rowTWanovaMatVec <- function(x, FDR=NULL, selVar='gSelectInd', pcutoff=0.05, oma=c(6,2,1,2), use=NULL, cat.pos = 0, cat.dist = 0.03, main='', ...){
	pvalName_byCL <- colnames(x)[str_detect(colnames(x), 'anovaP')]
	CL <- str_replace(pvalName_byCL, 'anovaP\\|', '')
	genes <- rownames(x)
	selL <- foreach(pn=pvalName_byCL) %do% {
		#browser()
		indSel <- create_bum(x[, pn], FDR=FDR, pcutoff=pcutoff)[[selVar]]
		genes[indSel]
	}
	names(selL) <- CL
	if(is.null(use)) use <- ifelse(length(pvalName_byCL)<=5, 'venn2', 'imageVenn')
	measureMat <- imageVenn(selL, oma=oma, plot=F)$measureMat	
	if(use=='venn2'){
		if(length(selL)>5) stop(sprintf('%d lists are present; only able to do venn diagram in at most 5 elements', length(selL)))
		#browser()
		venn2(selL, cat.pos = cat.pos, cat.dist = cat.dist, ...)
	} else {
		imageVenn(selL, oma=oma, main=main, colFn=redscale, cexColAdj=1.1)	
	}
	# add sorting and meta info:
	#browser()
	myd <- data.frame(measureMat, check.names=F)
	myd$nShare <- apply(myd, 1, sum, na.rm=T)
	myd <- myd[order(myd$nShare, decreasing=TRUE), ]
	list(measureMat0=measureMat, measureMat=myd)
	#measureMat
}
#measureMat <- venndiagram(res_twanova, FDR=0.1, use='venn2')
#measureMat <- venndiagram(res_twanova, FDR=0.1, use='imageVenn')

#' Scatter plot for rowTWanovaMatVec class: each row produces a printted figure
#'
#' @method scatter rowTWanovaMatVec
#' @param x a rowTWanovaMatVec class (a data frame)
#' @param var what variable to be used to order the genes (pval_main, pval_confounder, pval_interaction, anovaP|HCC1806, ...)
#' @param ylabsuffix ylab suffix to be added after gene name
#' @param ... additional parameters to plot.fn_TWanova
#' @export
scatter.rowTWanovaMatVec <- function(x, var='pval_main', ylabsuffix=', expression', maxRow=300, main=NULL, ...){
	# when too many rows, select top rows based on P
	if(nrow(x)>maxRow){
		ind <- which(isTopN(-x[, var], N=maxRow))
		main <- sprintf('%s\nTop %d genes selected based on p value', main, maxRow)
	} else {
		ind <- order(x[, var], decreasing=FALSE)
	}
	#browser()
	resL <- Getter(x, 'resL')
	resL_trim <- resL[rownames(x)] # this is how the index is relative to
	for(i in ind){
		ylab <- sprintf('%s%s', rownames(x)[i], ylabsuffix)
		print(plot(resL_trim[[i]], ylab=ylab, ...))
	}
} 
#scatter.rowTWanovaMatVec(res_twanova[1:3, ])

#' scatter plot for rowanovaMatVec
#'
#' this function needs to be polished and put into TestAndCorr_p.R
#'
#' @param mat expression matrix used in rowanovaMatVec 
#' @param vec vec as used in rowanovaMatVec 
#' @param FUN plotting function
#' @param geneSel index for plotting; default is plotting all sorted by P
#' @param rowanovaMatVec rowanovaMatVec object
#' @param ... additional parameters to plotScatter
#' @return multiple ggplot figures printed
#' @export
scatter_rowTWanovaMatVec <- function(mat, vec, FUN, geneSel, rowanovaMatVec, matIsX=TRUE, ...){
	if(missing(geneSel)){
		anP <- getter(rowanovaMatVec, 'anovaP')
		geneSel <- order(anP, decreasing=FALSE)
	}
	if(missing(FUN)){
		# notice that xlab and ylab are from meta info
		FUN <- function(expr, vec, gene, meta, matIsX, ...) {
			ylab <- ifelse(matIsX, meta$row, meta$response)
			xlab <- ifelse(matIsX, meta$response, meta$row)
			#browser()
			plotContCat(vec, expr, ylab=ylab, xlab=xlab, tag=gene, plot=F, ...)$p
		}
	}
	# vec is assumed to be a factor (rowanovaMatVec); thus forced it here
	vec <- getter(rowanovaMatVec, 'vec')
	meta <- getter(rowanovaMatVec, 'meta') # meta info
	for(g in geneSel){
		gene <- rownames(mat)[g]
		if(matIsX){
			# vec is group; mat row is continuous
			# vec is assumed to be a factor (rowttestMatVec); thus forced it here
			gg <- factor(vec)	
			expr <- as.numeric(mat[g, ])
		} else {
			# vec is coninuous; mat row is group
			expr <- as.numeric(vec)
			gg <- factor(mat[g, ])	
		}
		#browser()
		tt=try(FUN(expr=expr, vec=gg, gene=gene, meta, matIsX=matIsX, ...), silent=T)
		if(!inherits(tt, 'try-error')) {
			print(tt)
		}	
	}
}


#' ANOVA between rows of a matrix and a vector
#' 
#' @param mat a matrix, i.e. gene expression
#' @param vec a vector of phenotype, categorical variable
#' @param levels levels of the vec; when vec is not a factor, levels would be alphabetical order which will determine the scatter plot. 
#' @param meta meta information, a list. This can specify what type of analysis (i.e. RPPA correlated with EMT score);
#'  row, which represents what row variables are (rows of mat); response represents what does vec mean
#' @return an object belonging to data frame and a class
#' @export
rowanovaMatVec <- function(mat, vec, matIsX=TRUE, levels=NULL, 
	meta=list(tag='ANOVA', row='expression', response='response'), 
	base=2, more=TRUE, addData=TRUE, ...){
	#browser()
	if(is.null(levels)) {
		if(is.factor(vec))
			levels <- levels(vec)
		else 
			levels <- sort(unique(vec))
	}
	if(length(levels)!=n_unique(noNA(vec)) | !all(levels %in% unique(vec))) {
		#browser()
		stop('levels and unique values in vec in consistent!')
	}
	vec <- reorder(factor(vec), new.order=levels)
	#res <- rowMatVecDo(mat=mat, vec=vec, fn=fn_anovaP, base=base, more=more)
	res <- rowMatVecDo(mat=mat, vec=vec, fn=fn_anovaP, matIsX=matIsX, base=base, more=more)
	attr(res, 'meta') <- meta
	attr(res, 'base') <- base
	if(addData){
		attr(res, 'mat') <- mat
		attr(res, 'vec') <- vec
	}
	attr(res, 'class') <- c('rowanovaMatVec', 'data.frame')
	# genes * (group, FC, pval)
	formatted <- format_rowanovaMatVec(res) 
	attr(res, 'formatted') <- formatted 
	attr(res, 'matIsX') <- matIsX 
	attr(res, 'levels') <- levels 
	attr(res, 'anovaP') <- attributes(formatted)$anovaP
	res
}
#tp <- rowanovaMatVec(mat=ic50_aug, vec=pData$pathology, fn=fn_anovaP)

#' get method for class rowanovaMatVec
#' 
#' this function is named getter, not get to minimize conflict with the usual get function
#'
#' @param obj an rowanovaMatVec object
#' @param what specify what to get from augMat, i.e. base, var.equal
#' @method getter rowanovaMatVec
#' @return an object as requested
#' @export
getter.rowanovaMatVec <- function(obj, what='avail'){
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

# this seems not necesary: only need to rename anovaP as pval
format_rowanovaMatVec <- function(x){
	#browser()
	base <- Getter(x, 'base')
	coln <- colnames(x)
	if(base!=1){
		fcn <- coln[str_detect(coln, '^FC_')]
		fcn2 <- str_replace(fcn, '^FC_', '')
	} else {
		fcn <- coln[str_detect(coln, '^diff_')]
		fcn2 <- str_replace(fcn, '^diff_', '')
	}
	names(fcn2) <- fcn
	pn <- coln[str_detect(coln, '^padj_')]
	pn2 <- str_replace(pn, '^padj_', '')
	names(pn2) <- pn		
	if(length(fcn)>0){
		FC <- melt(data.matrix(rename(x[, fcn], fcn2))) # matrix keeps both row and col names
		pval <- melt(data.matrix(rename(x[, pn], pn2)))
		#browser()
		res <- merge(FC, pval, by.x=1:2, by.y=1:2)
		colnames(res) <- c('row', 'group', 'FC', 'pval') # pval is adjusted p; FC is FC if base!=1; FC=mDiff if base==1
	} else {
		# in case FC is not suitable
		res <- x
	}
	anP <- x[, 'anovaP']
	names(anP) <- rownames(x)
	attr(res, 'anovaP') <- anP
	#browser()
	res
}
#format_rowanovaMatVec(anovaRes_pathology)

#' Bum plot for rowanovaMatVec class
#'
#' @method hist rowanovaMatVec
#' @param x a rowanovaMatVec class (a data frame)
#' @param ... additional parameters to plotBumFDR
#' @export
hist.rowanovaMatVec <- function(x, main=NULL, plot=TRUE, return=FALSE, alphas = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), 
    FDR=NULL, maskPthr = 1, pcutoff = 0.05, ...){
	obj_bum <- create_bum(pval=x[, 'anovaP'], alphas=alphas, FDR=FDR, maskPthr=maskPthr, pcutoff=pcutoff)
	if(is.null(main)){
		main <- sprintf('%s\n', getAnalysisName(x))
	}
	if(plot){
		plotBumFDR(obj_bum, main=main, ...)
	}
	if(return){
		return(obj_bum)
	}
}

#' ordered horizontal bar plot for rowanovaMatVec
#'
#'
#' @param rowInfo alternative rownames to show on the figure
#' @param decreasing logic, decreasing order of p value used to sort the rows (genes) based on anovaP.
#' @param dataMode a string used to label the x-axis. Default is 'IC50' for IC50 t-test; Please use meaning string relevant to the data. 
#' @param showP do not show P value on graph ('None'); show actual P ('Value'); show P category ('Category')
#' @param pthr filtering by p value. Rows with larger P value will be deleted. Default is to use all rows (pthr=1)
#' @param legend.position position for the legend if P value is show as category
#' @param sizex x axis label size
#' @param sizey y axis label size
#' @param size_group label size for facets
#' @param legend.title title of the legend
#' @param truncFClow cutoff to truncate FC lower than this value. Defalt is -Inf, which means no truncation on small FC (negative) value
#' @param truncFChigh cutoff to truncate FC higher than this value. Defalt is Inf, which means no truncation on large FC (usually positive) value
#' @param barWidth relative barwidth (0~1)
#' @param groupShown string for group to shown, leading a t test type (FC vs adj p value). e.g. 'HCC4006_ER1-HCC4006'
#' @return a ggplot2 object
#' @export
plot4rowanovaMatVec <- function(x, rowInfo=NULL, decreasing=FALSE, dataMode='IC50', groupShown=NULL, ylab=NULL, xlab=NULL, sizex=12, sizey, size_group=10, main=NULL, 
	showP='None', pthr=1, colorSignif=TRUE, legend.position=c(.85,0.2), legend.title='Significance',
	truncFClow=-Inf, truncFChigh=Inf, barWidth=0.8, verbose=FALSE){
	#dat <- format_rowanovaMatVec(anovaRes_pathology)
	#browser()
	if(is.null(ylab)) ylab <- ''
	if(is.null(main)) main <- ''
	#brower()
	#dat0 <- getter(x, 'formatted') # this does not deal with subset selection
	dat0 <- format_rowanovaMatVec(x)  # when a subset of rows selected, this will still work;
	anoP <- getter(x, 'anovaP')
	base <- getter(x, 'base')
	#rows_orderedByP <- names(anoP)[order(anoP, decreasing=decreasing)] # sort genes by anovaP
	### sometimes NA for pval: this leads to error
	dat <- subset(dat0, !is.na(pval))
	#browser()
	# truncate FC 
	if(base!=1){
		# with FC
		if(is.null(xlab)) xlab <- ifelse(is.null(groupShown), sprintf('Fold change of mean %s', dataMode), sprintf('Fold change of mean %s\n(%s)', dataMode, str_replace(groupShown, '-', '/')))
	} else {
		# no FC, use mDiff, FC is mDiff due to format_rowanovaMatVec; thus only need to change xlab
		if(is.null(xlab)) xlab <- sprintf('Difference of mean %s', dataMode)
	}
	
	# add P category
	dat$pcat <- categorizePvec(dat[, 'pval'])
	#browser()
	#
	if(is.null(rowInfo)) rowInfo <- as.character(dat$row)
	# showP: add P into drug info
	if(showP=='Value'){
		rowInfo <- mapply(function(x, y) sprintf('%s, P=%.3f', x, y), rowInfo, dat$pval)
	} else if(showP=='Category'){
		rowInfo <- mapply(function(x, y) sprintf('%s, %s', x, y), rowInfo, dat$pcat)
	}
	#browser()
	# format dat
	dat <- data.frame(dat, rowInfo=rowInfo) # add rowInfo
	dat <- subset(dat, !is.na(FC) & pval<pthr) # remove rows with FC NA (not enough sample size)
	dat <- mutate(dat, FC=truncByLimit(FC, Lower=truncFClow, Upper=truncFChigh))
	#browser()
	dat$anovaP <- anoP[match(dat$row, names(anoP))]
	dat$rowOrdered <- orderVecByVal(vec=dat$rowInfo, val=-dat$anovaP, decreasing=decreasing) # levels added
	dat <- mutate(dat, group2=str_replace(group, '-', ' VS\n'))
	if(missing(sizey)) sizey <- guessSize(n_unique(dat$row))/1.8
	scale_axis_size <- theme(axis.text.x=element_text(size=sizex, colour = "black", face ='bold'), 
		axis.text.y=element_text(colour = "black", face ='bold', size=sizey))
	#browser()
	fillCol <- c('#ED281E', '#2834A4', '#EDB51E', 'gray')
	names(fillCol) <- c('P<0.001', 'P<0.01', 'P<0.05', 'P>0.05')
	fillColVec <- fillCol[sort(unique(as.character(dat$pcat)), decreasing=FALSE)]
	#browser()
	theme_striptext <- theme(strip.text.x = element_text(size = size_group, colour = "black"))
	if(verbose) cat(sprintf('sizex=%.2f; sizey=%.2f\n', sizex, sizey))
	if(!colorSignif) {
	p <- ggplot(dat, aes(y=FC, x=rowOrdered))+geom_bar(stat = "identity", position="dodge", width=barWidth)+ 
		coord_flip()+
		ylab(xlab)+xlab(ylab)+ggtitle(main)+scale_axis_size+
		theme(legend.position=legend.position)
	} else {
	#browser()
		if(is.null(groupShown)){
			# show anova type: all groups
			p <- ggplot(dat, aes(y=FC, x=rowOrdered, fill=pcat))+geom_bar(stat = "identity", position="dodge", width=barWidth)+ 
			coord_flip()+scale_fill_manual(name=legend.title, breaks=names(fillColVec), values=fillColVec)+
			ylab(xlab)+xlab(ylab)+ggtitle(main)+
			scale_axis_size+theme(legend.position=legend.position)+facet_wrap(~ group2)+theme_striptext
		} else {
			# show t test type: one group
			# need to reorder by FC
			#browser()
			#if(is.null(xlab)) xlab <- sprintf('Fold change of mean %s\n(%s)', dataMode, str_replace(groupShown, '-', '/'))
			datShown <- subset(dat, group==groupShown)	
			datShown$rowOrdered <- orderVecByVal(vec=datShown$rowInfo, val=datShown$FC, decreasing=decreasing)
			p <- ggplot(datShown, aes(y=FC, x=rowOrdered, fill=pcat))+geom_bar(stat = "identity", position="dodge", width=barWidth)+ 
			coord_flip()+scale_fill_manual(name=legend.title, breaks=names(fillColVec), values=fillColVec)+
			ylab(xlab)+xlab(ylab)+ggtitle(main)+
			scale_axis_size+theme(legend.position=legend.position)+theme_striptext
		}

	}	
	#browser()
	p
}
#plot4rowanovaMatVec(anovaRes_pathology)


#' Bar plot for rowanovaMatVec class
#'
#' @method plot rowanovaMatVec
#' @param x a rowanovaMatVec class (a data frame)
#' @param ... additional parameters to plot4rowanovaMatVec
#' @export
plot.rowanovaMatVec <- function(x, maxRow=400, main=NULL, dataMode=NULL, groupShown=NULL, ...){
	# when too many rows, select top rows based on P
	if(is.null(main)) {
		if(is.null(groupShown)){
			main <- getAnalysisName(x)
		} else {
			main <- sprintf('%s', str_replace(groupShown,'-', ' vs '))
		}	
	}
	if(is.null(dataMode)) dataMode <- getter(x, 'meta')$row
	if(nrow(x)>maxRow){
		ind <- which(isTopN(-x[, 'anovaP'], N=maxRow))
		main <- sprintf('%s\nTop %d genes selected based on p value', main, maxRow)
	} else {
		ind <- 1:nrow(x)
	}
	#browser()
	plot4rowanovaMatVec(x[ind, ], main=main, dataMode=dataMode, groupShown=groupShown, ...)
} 


#' scatter plot for rowanovaMatVec
#'
#' this function needs to be polished and put into TestAndCorr_p.R
#'
#' @param mat expression matrix used in rowanovaMatVec 
#' @param vec vec as used in rowanovaMatVec 
#' @param FUN plotting function
#' @param geneSel index for plotting; default is plotting all sorted by P
#' @param rowanovaMatVec rowanovaMatVec object
#' @param maintagfunc a function to compute additional information for main 
#' @param ... additional parameters to plotScatter
#' @return multiple ggplot figures printed
#' @export
scatter_rowanovaMatVec <- function(mat, vec, FUN, geneSel, rowanovaMatVec, matIsX=TRUE, maintagfunc=NULL, groupShown=NULL, ...){
	if(missing(geneSel)){
		anP <- getter(rowanovaMatVec, 'anovaP')
		geneSel <- order(anP, decreasing=FALSE)
	}
	if(missing(FUN)){
		# notice that xlab and ylab are from meta info
		FUN <- function(expr, vec, gene, meta, matIsX, main, ...) {
			ylab <- ifelse(matIsX, meta$row, meta$response)
			xlab <- ifelse(matIsX, meta$response, meta$row)
			# if(!is.null(maintagfunc)) {
			# 	tag <- maintagfunc(gene)
			# } else {
			# 	tag <- gene
			# }
			# #browser()
			plotContCat(vec, expr, ylab=ylab, xlab=xlab, main=main, plot=F, ...)$p
		}
	}
	# vec is assumed to be a factor (rowanovaMatVec); thus forced it here
	vec <- getter(rowanovaMatVec, 'vec')
	meta <- getter(rowanovaMatVec, 'meta') # meta info
	for(g in geneSel){
		gene <- rownames(mat)[g]
		#browser()
		if(matIsX){
			# vec is group; mat row is continuous
			# vec is assumed to be a factor (rowttestMatVec); thus forced it here
			gg <- factor(vec)	
			expr <- as.numeric(mat[g, ])
		} else {
			# vec is coninuous; mat row is group
			expr <- as.numeric(vec)
			gg <- factor(mat[g, ])	
		}
		myd <- data.frame(expr=expr, gg=gg)
		if(!is.null(maintagfunc)) {
			tag <- maintagfunc(gene)
		} else {
			tag <- gene
		}
		# subsetting
		if(!is.null(groupShown)){
			# t test type
			groups <- rev(str_split(groupShown, '-')[[1]]) # just two groups; reversed as level
			myd <- subset(myd, gg %in% groups)
			# ordering
			myd <- mutate(myd, gg=orderFactor(gg, new.order=groups))
			if(!is.na(rowanovaMatVec[g, str_c('FC_', groupShown, sep='')])){
				text_test <- sprintf('p=%.3g; fold change=%.3g\n', rowanovaMatVec[g, str_c('padj_', groupShown, sep='')], rowanovaMatVec[g, str_c('FC_', groupShown, sep='')])
			} else {
				text_test <- sprintf('p=%.3g; mean difference=%.3g\n', rowanovaMatVec[g, str_c('padj_', groupShown, sep='')], rowanovaMatVec[g, str_c('diff_', groupShown, sep='')])
			}
			main <- sprintf('%s\n%s', tag, text_test)
			tt=try(FUN(expr=myd$expr, vec=myd$gg, gene=gene, meta, main=main, matIsX=matIsX, ...), silent=T)
		} else {
			# ANOVA type
			#browser()
			text_test <- sprintf('ANOVA test P value: %.3g', rowanovaMatVec[g, 'anovaP'])
			main <- sprintf('%s\n%s', tag, text_test)
			tt=try(FUN(expr=myd$expr, vec=myd$gg, gene=gene, meta, matIsX=matIsX, main=main, ...), silent=T)
		}
		#browser()
		if(!inherits(tt, 'try-error')) {
			print(tt)
		}	
	}
}
# scatter_rowanovaMatVec(mat=ic50_aug, vec=pData$pathology, rowanovaMatVec=anovaRes_pathology, anglex=45)
#plotContCat(x=factor(pData$pathology), y=ic50_aug[1, ], ignoreNA=TRUE)
#anovaRes_pathology <- rowanovaMatVec(mat=ic50_aug, vec=pData$pathology, meta=list(tag='ANOVA', row='IC50', response='pathology'))


#' Scatter plot for rowanovaMatVec class: each row produces a printted figure
#'
#' @method scatter rowanovaMatVec
#' @param x a rowanovaMatVec class (a data frame)
#' @param ... additional parameters to scatter_rowttestMatVec
#' @export
scatter.rowanovaMatVec <- function(x, maxRow=300, main=NULL, maintagfunc=NULL, groupShown=NULL, ...){
	# when too many rows, select top rows based on P
	if(nrow(x)>maxRow){
		ind <- which(isTopN(-x[, 'anovaP'], N=maxRow))
		main <- sprintf('%s\nTop %d genes selected based on p value', main, maxRow)
	} else {
		ind <- order(x$anovaP, decreasing=FALSE)
	}
	# requires addData=TRUE
	if(is.null(getter(x, 'mat'))) 
		stop('No mat (e.g. expression matrix) found! Please specify addData=TRUE\n')
	#browser()
	mat <- getter(x, 'mat')
	vec <- getter(x, 'vec')
	matIsX <- getter(x, 'matIsX')
	#
	# fix bug: need to also subset mat !
	#
	if(!identical(rownames(x), rownames(mat))){
		mat <- mat[rownames(x), ] # subsetting
	}
	#browser()
	# scatter_rowanovaMatVec(mat=ic50_aug, vec=pData$pathology, rowanovaMatVec=anovaRes_pathology, anglex=45)
	scatter_rowanovaMatVec(mat=mat, vec=vec, groupShown=groupShown, geneSel=ind, rowanovaMatVec=x, matIsX=matIsX, maintagfunc=maintagfunc, ...)
} 
#scatter.rowanovaMatVec(anovaRes_pathology)


#' convert P value into categories
#'
#' @param Pvec a vector of P value
#' @return result by cut function
#' @export
categorizePvec <- function(Pvec){
	cut(Pvec, breaks=c(0-0.01, 0.001, 0.01, 0.05, 1+0.01), labels=c('P<0.001', 'P<0.01', 'P<0.05', 'P>0.05'))
}


#### sometimes we have deliberate treatment and control assignment. Therefore, the difference is treatment-control. We do not need a label now as it is in the trt and ctrl data.
## mDiff: it is treatment mean - control mean
## tstat: it is adjusted such that it is treatment-control, the same sign as mDiff
## add paired t test option; this is included in the t.test()
ttestTrtVsCtr <- function(ctr, trt, minN=2, more=FALSE, paired=FALSE, var.equal=TRUE){
	if(paired){
		if(length(ctr)!=length(trt)) stop('ctr and trt should have equal length for paired t-test!\n')
	}
	if(more) {
		if(sum(!is.na(ctr))<minN | sum(!is.na(trt))<minN) {
			res <- c(pval=NA, tstat=NA, mDiff=NA)
		} else {
		# mean_C2-mean_C1
		# for mutation 0, 1 it is mean_mutated-mean_WT
		tt <- t.test(ctr, trt, na.action=na.omit, var.equal=var.equal, paired=paired)
		mDiff <- diff(tt$estimate) # of course, it is trt-ctr
		pval <- tt$p.value
		res <- c(pval=pval, tstat=-unname(tt$statistic), mDiff=unname(mDiff))
		#browser()
		## add info as what diff is computed
		#attr(res, 'info') <- paste('diff by class label: ', C[1], '-', C[2], sep='')
		}
		#browser()
		return(res)
	} else {
		if(sum(!is.na(ctr))<=minN | sum(!is.na(trt))<=minN) {
		res <- NA
		} else {
		res <- t.test(ctr, trt, na.action=na.omit, var.equal=var.equal, paired=paired)$p.value
		}
		return(res)
	}
}
#ttestTrtVsCtr(ctr=DataFile[1, indCtr], trt=DataFile[1, indTrt], minN=2, more=TRUE)

#' ordinal regression for association between an ordered categorical variable with a conitnuous variable
#' 
#' @param expr a vector of continuous values
#' @param label a categorical vector or factor
#' @param addANOVA whether to perform an ANOVA analysis
#' @param base additional parameters passed to anovaP; only effective if addANOVA is TRUE.
#' @param more additional parameters passed to anovaP; only effective if addANOVA is TRUE.
#' @return a vector
#' @export
polrP <- function(expr, label, addANOVA=TRUE, base=2, more=TRUE){
	#browser()
	require(MASS)
	m <- try(polr(label~expr, na.action=na.omit, Hess=TRUE), silent=TRUE)
	if(class(m)[1]!='try-error'){
		ctable <- coef(summary(m))
		p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
		ctable <- cbind(ctable, `p value` = p)
		res <- ctable[1, c('Value', 'p value')]
	} else {
		res <- rep(NA, 2)
	}
	names(res) <- c('OrdinalRegCoef', 'OrdinalRegP')
	res <- rev(res)
	if(addANOVA) {
		res <- c(res, anovaP(expr=expr, label=label, base=base, more=more))
	}
	res
}
#polrP(expr=ll_RPPA[[1]][1, ], label=ll_ICvecBin[[1]])





# anova P
## note that lm(expr~label) rather than lm(label~expr): bug fixed on 2013/09/09

#' ANOVA test P value
#'
#' This function calculates F test P value and Tukey HSD test. The adjusted P values, difference of means as well as FC are also calculated
#' @param expr a vector of continuous values
#' @param label a categorical vector or factor
#' @param base base for calculating Fold Change
#' @param more whether to supply just P value or add Tukey HSD info
#' @return a vector 
#' @export
anovaP <- function(expr, label, base=NULL, more=TRUE)
{
	if(is.null(base)) base <- 2
	# cann tolerate NA in y or x
	#fit <- try(lm(expr~label, na.action=na.omit), silent=TRUE)
	aovFit <- try(aov(expr~label, na.action=na.omit), silent=TRUE)
	#browser()
	# aov has a bug: when expr has identical value (e.g. are all 0.9828544), P value might be significant!
	# this is likely due to numerical precision issue
	#if(class(aovFit)[1]!='try-error') {
	#browser()
	if(class(aovFit)[1]!='try-error' & n_unique(expr, ignoreNA=TRUE)>2) {
		#browser()
		aovP <- summary(aovFit)[[1]][1, 5] # ANOVA P
		names(aovP) <- 'anovaP'
		tukey <- TukeyHSD(aovFit)[[1]]
		dif0 <- tukey[, 'diff']
		dif <- addNamePrefix(dif0, 'diff')
		FC <- addNamePrefix(toFC(dif0, base=base), 'FC')
		padj <- addNamePrefix(tukey[, 'p adj'], 'padj')
		res <- c(aovP, dif, padj, FC)
	} else {
		res <- rep(NA, length=1+3*n_unique(label, ignoreNA=TRUE))
		## problem: impossible to add names to res, leading a problem if we fillVecByID.
		#names(res)<- c('anovaP', 'diff_Mid-Low', 'diff_High-Low', 'diff_High-Mid', 'padj_Mid-Low', 'padj_High-Low', 'padj_High-Mid', 'FC_Mid-Low', 'FC_High-Low', 'FC_High-Mid')
		res
	}
	if(more){
		return(res)
	} else {
		return(res[1])
	}	
}
#anovaP(expr=ic50_aug[1, ], label=pData$pathology)
#anovaRes_pathology <- rowanovaMatVec(mat=ic50_aug, vec=pData$pathology, meta=list(tag='ANOVA', row='IC50', response='pathology'))


##### cox PH regression
coxphP <- function(expr, OS, OScensoring)
{
	fit <- try(coxph(Surv(OS, OScensoring)~expr, na.action=na.omit), silent=TRUE)
	# in case model fails, i.e. no one dies, returns NA
	if(class(fit)=='try-error')
		return (NA)
	res <- summary(fit)$logtest[3]
	res
}

## add OR (categorical) or HR (continuous)
# param name a name for eacg of the categories (continuous, just 1) in expr
# param useLRT whether to use LRT p value or wald test P; default is LRT since it is more powerful/suitable for small sample size although asymtotically they are the same
#' cox PH P value
#' @export
coxphP2 <- function(expr, OS, OScensoring, name, useLRT=TRUE){
	fit <- try(coxph(Surv(OS, OScensoring)~expr, na.action=na.omit), silent=TRUE)
	# in case model fails, i.e. no one dies, returns NA
	if(class(fit)=='try-error') {
		res <- list(P=NA, coef_exp=NA, coef_exp_P=NA, ciL=NA, ciU=NA, fit=fit)
		return (res)
	}
	#pval <- summary(fit)$logtest[3] # this is LRT version
	smry <- summary(fit)
	#browser()
	if(useLRT){
		pval <- smry$logtest[3]  # likelihood ratio test: more powerful at small sample size
	} else {
		pval <- smry$sctest[3]  # to make it consistent with coef_exp_P; also this is the logrank test version: wald test, performance not good; asymtotically, logrank test is the LRT for cox PH comparing two groups
	}
	cc <- smry$coefficients
	ci <- smry$conf.int
	ciL <- ci[, 'lower .95']
	ciU <- ci[, 'upper .95']
	coef_exp <- cc[, 2]
	coef_exp_P <- cc[, 5]
	newName <- str_replace(rownames(smry$coef), 'expr', '')
	# it turns out smry$coef is always a matrix; when a cont var is specified,  rownames(smry$coef)='expr', which we cannot make it '' (uninformative)
	if(newName[1]==''){
		newName <- deparse(substitute(expr))[1] # here coef_exp has no name, thus extract from input
		#newName <- str_replace(rownames(smry$coef), 'expr', '')
	}
	# override the names when a user has specified a proper name for the OR/CI
	if(!missing(name)) {
		if(length(name)==length(newName))
			newName <- name
	}
	#browser()
	names(coef_exp) <- names(coef_exp_P) <- names(ciL) <- names(ciU) <- newName # use  more informative names
	res <- list(P=pval, coef_exp=coef_exp, coef_exp_P=coef_exp_P, ciL=ciL, ciU=ciU, fit=fit)
	res
}
#tt <- coxphP2(expr=myy, OS=tdf$OS, OScensoring=tdf$OScensoring, name='AXL')

#tt <- coxphP2(expr=cut(myy, breaks=c(-Inf, quantile(myy, prob=c(0.5), na.rm=TRUE), Inf)), OS=tdf$OS, OScensoring=tdf$OScensoring)
#tt <- coxphP2(expr=cut(myy, breaks=c(-Inf, quantile(myy, prob=c(1/3, 2/3), na.rm=TRUE), Inf)), OS=tdf$OS, OScensoring=tdf$OScensoring)

#tt <- coxphP2(expr=myy, OS=tdf$OS, OScensoring=tdf$OScensoring)

#' parse the coxphP2() returned value to text so as to show OR, P and so on
#'
#' @param x results returned by coxphP2()
#' @param type either OR from categorical variable or HR from continuous variable
#' @param digits digits to keep in OR and CI
#' @return a list of text string. i.e.
#'
#' $coefExp_text
#'
#' [1] "OR: 1.41"
#'
#' $CI_text
#'
#' [1] "CI: [0.88, 2.27]"
#'
#' parse coxphP2 result
#' @export
parse_coxphP2 <- function(x, type='HR', digits=2){
	coefExp <- x$coef_exp
	fit <- x$fit
	ciL <- x$ciL
	ciU <- x$ciU
	#browser()
	coefExp_text <- sprintf('%s: %s', type, str_c(names(coefExp), round(coefExp, digits), sep=': ', collapse=' '))
	CI_text <- sprintf('CI: %s', str_c(names(coefExp), ': [',  round(ciL, digits), ', ', round(ciU, digits), ']', 
				sep='', collapse=' '))
	list(coefExp_text=coefExp_text, CI_text=CI_text)
}
#parse_coxphP2(tt)