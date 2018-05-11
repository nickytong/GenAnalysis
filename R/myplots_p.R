####
#### to build
####
if(FALSE){
# cd /home/ptong1/Backup/Package/; R0

library(roxygen2)
#library(roxygen) # not working
roxygenize("GenAnalysis")

library(devtools)
build('GenAnalysis', binary=TRUE)
install('GenAnalysis')

check('GenAnalysis')
check_cran('GenAnalysis')

}
#' Calculates Concordance Correlation Coefficient (CCC) to access reproducibility
#' 
#' This function calculates and plots concordance for paired data.
#'
#' This function calls the epi.ccc() function from epiR package for the calculation. 
#'
#' infinite values in x and y are masked as NA so as to compute CCC and corr (otherwise, result becomes nan).
#' @param x x vector
#' @param y y vector 
#' @param xlab xlab
#' @param ylab ylab
#' @param maskBeyondLim whether mask values beyond xlim (for x) or ylim (for y) when xlim and ylim is specified. Default is FALSE so that even xlim ylim specified, CCC will not be affected
#' @param tag add a tag
#' @param cols color vector for the points
#' @param pch pch of the points
#' @param main main title 
#' @param cex for the dots 
#' @param cex.main cex for main title
#' @param cex.legend cex for legend text
#' @param plot whether draw a figure
#' @param plotOutlier logical whether add an outlier plot comparing Diff vs Mean plot and related statistics
#' @param sampleName symbols of the text to be shown in outlier plot
#' @param alpha_outlier alpha to call outliers based on the observed differences assuming from Normal(mean_diff, sd_diff)
#' @return a vector of c('ccc', 's_shift', 'l_shift', 'ccc_lo', 'ccc_hi', 'Cb', 'corr'). s_shift is scale shift; l_shift is location shift;
#'			ccc_lo and ccc_hi represent 95% confidence interval for CCC. Cb for the bias correction term satisfying CCC=corr*Cb where corr is the
#' 			Pearson correlation
#' @export
plotCCC <- function(x, y, labels=NULL, cols=NULL, cexlab=0.8, colab='purple', poslab=1, xlab=NA, ylab=NA, tag='', col_legend='red', 
	pch=20, main=NA, cex=2, plot=TRUE, ylim=NA, xlim=NA, maskBeyondLim=FALSE, 
	plotOutlier=FALSE, sampleName=NA, alpha_outlier=0.01,
	show=c('identity', 'h0', 'v0')){
	# show=c('lm', 'identity', 'anno')
	require(epiR)
	require(calibrate) # text around points in base graphics
	if(is.na(xlab)) xlab <- deparse(substitute(x)) # should happen at very beginning
	if(is.na(ylab)) ylab <- deparse(substitute(y))
	if(!is.null(labels) & length(labels)!=length(x)) stop('labels should have same length as data points!')
	#browser()
	if(!is.vector(x)) x <- as.vector(x)
	if(!is.vector(y)) y <- as.vector(y)
	if(length(x)!=length(y)) stop('Length of x and y is unequal!')
	if(!is.na(ylim[1]) & maskBeyondLim) y <- maskByLimit(y, Lower=ylim[1], Upper=ylim[2])
	if(!is.na(xlim[1]) & maskBeyondLim) x <- maskByLimit(x, Lower=xlim[1], Upper=xlim[2])
	## inf value will lead to NaN in CCC and cor: for simplicity, mask as NA
	isInf_x <- is.infinite(x)
	isInf_y <- is.infinite(y)
	if(sum(isInf_x)>0) warning(sprintf('%d infinite values observed in x, masked as NA in order to compute CCC and Cor!', sum(isInf_x)))
	if(sum(isInf_y)>0) warning(sprintf('%d infinite values observed in y, masked as NA in order to compute CCC and Cor!', sum(isInf_y)))
	x[isInf_x] <- NA
	y[isInf_y] <- NA
	## remove NA: so that cor and CCC will not be NA
	if(is.null(cols)) cols <- rep('black', length(x)) # point colors
	tdf <- noNA(data.frame(x=x, y=y, cols=cols))
	sel <- noNA(data.frame(x=x, y=y), returnIndex=TRUE)
	x <- tdf$x
	y <- tdf$y
	cols <- tdf$cols
	## the default epi.ccc use sd(y)/sd(x) and y-x as the scale and location shift!
	ccc_fit <- epi.ccc(x, y, ci = "z-transform", conf.level = 0.95) 
	############### outlier calculation
	dat_outlier <- data.frame(avg=(x+y)*0.5, dif=y-x)
	### so inconsistent!: this is x-y!; use my own calc!
	#diff_mean <- mean(ccc_fit$blalt$delta)
	#diff_sd <- sqrt(var(ccc_fit$blalt$delta))
	diff_mean <- mean(dat_outlier$dif, na.rm=TRUE)
	diff_sd <- sd(dat_outlier$dif, na.rm=TRUE)
	### basic statistics that may be used for filtering
	mean_x <- mean(x, na.rm=TRUE)
	mean_y <- mean(y, na.rm=TRUE)
	sd_x <- sd(x, na.rm=TRUE)
	sd_y <- sd(y, na.rm=TRUE)
	#browser()
	extremes <- qnorm(c(alpha_outlier/2, 1-alpha_outlier/2), diff_mean, diff_sd) # lower and upper bound for outliers
	dat_outlier$isOutlier <- dat_outlier[, 'dif'] < extremes[1] | dat_outlier[, 'dif'] > extremes[2] 
	if(is.na(sampleName[1])) sampleName <- 1:nrow(dat_outlier)
	dat_outlier$sampleName <- sampleName
	#diff_mean-2*diff_sd
	#diff_mean+2*diff_sd
	############### CCC calculation
	tt <- unlist(ccc_fit$rho.c)
	s_shift <- ccc_fit$s.shift # scale shift , that is sigma1/sigma2
	l_shift <- ccc_fit$l.shift # location shift, that is (mu1-mu2)/sqrt(sigma1*sigma2)
	ccc <- tt['est'] # CCC estimate
	ccc_lo <- tt['lower'] # 95% CI lower
	ccc_hi <- tt['upper'] # 95% CI upper
	corr <- cor(x, y, use='complete.obs')
	Cb <- ccc_fit$C.b # bias correction factor; cor*Cb=CCC
	CCCinfo0 <- c(ccc, s_shift, l_shift, ccc_lo, ccc_hi, Cb, corr)
	names(CCCinfo0) <- c('ccc', 's_shift', 'l_shift', 'ccc_lo', 'ccc_hi', 'Cb', 'corr') 
	#res <- list(CCCinfo=res, mean_x=mean_x, mean_y=mean_y, sd_x=sd_x, sd_y=sd_y)
	aux <- c(mean_x=mean_x, mean_y=mean_y, sd_x=sd_x, sd_y=sd_y)
	CCCinfo <- c(CCCinfo0, aux)
	res <- CCCinfo
	#browser()
	if(plot){
		## when x is a mat, this xlab is the element of the mat and really long!
		if(length(xlab)>2) xlab <- 'x'
		if(length(ylab)>2) ylab <- 'y'
		#browser()
		if(is.na(main)) main <- sprintf('Concordance plot: CCC=%.3f, Corr=%.3f', ccc, corr)
		rr <- range(c(x, y), na.rm=TRUE) # using largest range possible so that we have a symmetric view
		if(!is.na(ylim[1])){
			ylim <- ylim
		} else {
			ylim <- rr
		}
		if(!is.na(xlim[1])){
			xlim <- xlim
		} else {
			xlim <- rr
		}
		plot(x, y, xlab=xlab, ylab=ylab, main=str_c(tag, '\n', main, collapse=''), col=cols, pch=pch, cex=cex, xlim=xlim, ylim=ylim)
		if(!is.null(labels)){
			#browser()
			text(x, y, labels=labels[sel], cex=cexlab, col=colab, pos=poslab)
		}
		#abline(0, 1, col=80, lty=2, lwd=2)
		if('identity' %in% show) abline(a = 0, b = 1, lty = 2, col=4, lwd=3) # identical line
		if('h0' %in% show) abline(h=0, lty=2, col='red') # horizontal line at 0
		if('v0' %in% show) abline(v=0, lty=2, col='red') # vertical line at 0
		#browser()
		if('lm' %in% show){
			try_lm_line <- try(abline(lm(y~x), lty = 1, lwd=2), silent=T) # lm line
		} 
		if('anno' %in% show) {
			legend('bottomright', legend = c("Perfect concordance", "Observed linear trend"), lty = c(2,1), bty = "n", col=c(4, 1), lwd=c(3, 2), text.col=col_legend)
			labCCC <- sprintf("-->CCC: %.3f, (95%% CI %.2f~%.2f)", res['ccc'], res['ccc_lo'], res['ccc_hi'])
			lab_s_shift <- sprintf("-->Scale shift: %.2f", res['s_shift'])
			lab_l_shift <- sprintf("-->Location shift: %.2f", res['l_shift'])
			legend('topleft', legend = c(labCCC, lab_s_shift, lab_l_shift), bty = "n", text.col=col_legend)		
		}
	}
	if(plotOutlier){
		#gtitle <- str_c(tag, sprintf('%d disconcordant outliers at significance level of %.3f', sum(dat_outlier$isOutlier, na.rm=TRUE), alpha_outlier), collapse='')
		#p <- ggplot(dat_outlier, aes(x=avg, y=dif))+geom_point(shape=1, size=6)+
		#	geom_text(data=subset(dat_outlier, isOutlier==TRUE), aes(label=sampleName))+
		#	geom_hline(yintercept = extremes, linetype=2, color='red')+
		#	geom_hline(yintercept = diff_mean, linetype=2, color='black')+
		#	xlab('Average measurement')+ylab('Difference between the two measurements')+ggtitle(gtitle)
		#print(p)
		#browser()
		gtitle <- str_c(tag, sprintf('%d disconcordant outliers at significance level of %.3f', sum(dat_outlier$isOutlier, na.rm=TRUE), alpha_outlier), collapse='')
		ylim <- c(min(c(extremes, dat_outlier$dif), na.rm=TRUE), max(c(extremes, dat_outlier$dif), na.rm=TRUE))
		with(dat_outlier, plot(avg, dif, cex=2, pch=1, ylim=ylim, xlab='Average measurement', ylab='Difference between the two measurements', main=gtitle))
		abline(h=extremes, lty=2, col='red')
		abline(h=diff_mean, lty=2, col='black')
		tt <- (data=subset(dat_outlier, isOutlier==TRUE))
		if(nrow(tt)>0){
			with(tt, calibrate:::textxy(avg, dif, labs=sampleName, cex=1.1, offset=0))
		}
		res <- list(CCCinfo=CCCinfo, dat_outlier=dat_outlier, diff_mean=diff_mean, diff_sd=diff_sd, indexOutlier=which(dat_outlier$isOutlier), indexNonOutlier=which(dat_outlier$isOutlier==FALSE))
	}	
	res
}

# plotCCC <- function(x, y, xlab=NA, ylab=NA, tag='', col_legend='red', pch=1, main=NA, cex.main=1, cex=2, cex.legend=1, plot=TRUE, ylim=NA, xlim=NA, maskBeyondLim=FALSE, plotOutlier=FALSE, sampleName=NA, alpha_outlier=0.01){
# 	require(epiR)
# 	require(calibrate) # text around points in base graphics
# 	if(is.na(xlab)) xlab <- deparse(substitute(x)) # should happen at very beginning
# 	if(is.na(ylab)) ylab <- deparse(substitute(y))
# 	#browser()
# 	if(!is.vector(x)) x <- as.vector(x)
# 	if(!is.vector(y)) y <- as.vector(y)
# 	if(length(x)!=length(y)) stop('Length of x and y is unequal!')
# 	if(!is.na(ylim[1]) & maskBeyondLim) y <- maskByLimit(y, Lower=ylim[1], Upper=ylim[2])
# 	if(!is.na(xlim[1]) & maskBeyondLim) x <- maskByLimit(x, Lower=xlim[1], Upper=xlim[2])
# 	## inf value will lead to NaN in CCC and cor: for simplicity, mask as NA
# 	isInf_x <- is.infinite(x)
# 	isInf_y <- is.infinite(y)
# 	if(sum(isInf_x)>0) warning(sprintf('%d infinite values observed in x, masked as NA in order to compute CCC and Cor!', sum(isInf_x)))
# 	if(sum(isInf_y)>0) warning(sprintf('%d infinite values observed in y, masked as NA in order to compute CCC and Cor!', sum(isInf_y)))
# 	x[isInf_x] <- NA
# 	y[isInf_y] <- NA
# 	## remove NA: so that cor and CCC will not be NA
# 	tdf <- noNA(data.frame(x=x, y=y))
# 	x <- tdf$x
# 	y <- tdf$y
# 	## the default epi.ccc use sd(y)/sd(x) and y-x as the scale and location shift!
# 	ccc_fit <- epi.ccc(x, y, ci = "z-transform", conf.level = 0.95) 
# 	############### outlier calculation
# 	dat_outlier <- data.frame(avg=(x+y)*0.5, dif=y-x)
# 	### so inconsistent!: this is x-y!; use my own calc!
# 	#diff_mean <- mean(ccc_fit$blalt$delta)
# 	#diff_sd <- sqrt(var(ccc_fit$blalt$delta))
# 	diff_mean <- mean(dat_outlier$dif, na.rm=TRUE)
# 	diff_sd <- sd(dat_outlier$dif, na.rm=TRUE)
# 	### basic statistics that may be used for filtering
# 	mean_x <- mean(x, na.rm=TRUE)
# 	mean_y <- mean(y, na.rm=TRUE)
# 	sd_x <- sd(x, na.rm=TRUE)
# 	sd_y <- sd(y, na.rm=TRUE)
# 	#browser()
# 	extremes <- qnorm(c(alpha_outlier/2, 1-alpha_outlier/2), diff_mean, diff_sd) # lower and upper bound for outliers
# 	dat_outlier$isOutlier <- dat_outlier[, 'dif'] < extremes[1] | dat_outlier[, 'dif'] > extremes[2] 
# 	if(is.na(sampleName[1])) sampleName <- 1:nrow(dat_outlier)
# 	dat_outlier$sampleName <- sampleName
# 	#diff_mean-2*diff_sd
# 	#diff_mean+2*diff_sd
# 	############### CCC calculation
# 	tt <- unlist(ccc_fit$rho.c)
# 	s_shift <- ccc_fit$s.shift # scale shift , that is sigma1/sigma2
# 	l_shift <- ccc_fit$l.shift # location shift, that is (mu1-mu2)/sqrt(sigma1*sigma2)
# 	ccc <- tt['est'] # CCC estimate
# 	ccc_lo <- tt['lower'] # 95% CI lower
# 	ccc_hi <- tt['upper'] # 95% CI upper
# 	corr <- cor(x, y, use='complete.obs')
# 	Cb <- ccc_fit$C.b # bias correction factor; cor*Cb=CCC
# 	CCCinfo0 <- c(ccc, s_shift, l_shift, ccc_lo, ccc_hi, Cb, corr)
# 	names(CCCinfo0) <- c('ccc', 's_shift', 'l_shift', 'ccc_lo', 'ccc_hi', 'Cb', 'corr') 
# 	#res <- list(CCCinfo=res, mean_x=mean_x, mean_y=mean_y, sd_x=sd_x, sd_y=sd_y)
# 	aux <- c(mean_x=mean_x, mean_y=mean_y, sd_x=sd_x, sd_y=sd_y)
# 	CCCinfo <- c(CCCinfo0, aux)
# 	res <- CCCinfo
# 	#browser()
# 	if(plot){
# 		## when x is a mat, this xlab is the element of the mat and really long!
# 		if(length(xlab)>2) xlab <- 'x'
# 		if(length(ylab)>2) ylab <- 'y'
# 		#browser()
# 		if(is.na(main)) main <- sprintf('Concordance plot: CCC=%.3f, Corr=%.3f', ccc, corr)
# 		rr <- range(c(x, y), na.rm=TRUE) # using largest range possible so that we have a symmetric view
# 		if(!is.na(ylim[1])){
# 			ylim <- ylim
# 		} else {
# 			ylim <- rr
# 		}
# 		if(!is.na(xlim[1])){
# 			xlim <- xlim
# 		} else {
# 			xlim <- rr
# 		}
# 		plot(x, y, xlab=xlab, ylab=ylab, main=str_c(tag, '\n', main, collapse=''), pch=pch, cex=cex, xlim=xlim, ylim=ylim, cex.main=cex.main)
# 		#abline(0, 1, col=80, lty=2, lwd=2)
# 		abline(a = 0, b = 1, lty = 2, col=4, lwd=3) # identical line
# 		#browser()
# 		# may be a constant x: get an error if abline(lm(y~x), lty = 1, lwd=2) # lm line
# 		tp=try(abline(lm(y~x), lty = 1, lwd=2), silent=TRUE) # lm line
# 		legend('bottomright', legend = c("Perfect concordance", "Observed linear trend"), lty = c(2,1), bty = "n", col=c(4, 1), lwd=c(3, 2), cex=cex.legend, text.col=col_legend)
# 		labCCC <- sprintf("-->CCC: %.3f, (95%% CI %.2f~%.2f)", res['ccc'], res['ccc_lo'], res['ccc_hi'])
# 		lab_s_shift <- sprintf("-->Scale shift: %.2f", res['s_shift'])
# 		lab_l_shift <- sprintf("-->Location shift: %.2f", res['l_shift'])
# 		legend('topleft', legend = c(labCCC, lab_s_shift, lab_l_shift), bty = "n", cex=cex.legend, text.col=col_legend)
# 	}
# 	if(plotOutlier){
# 		#gtitle <- str_c(tag, sprintf('%d disconcordant outliers at significance level of %.3f', sum(dat_outlier$isOutlier, na.rm=TRUE), alpha_outlier), collapse='')
# 		#p <- ggplot(dat_outlier, aes(x=avg, y=dif))+geom_point(shape=1, size=6)+
# 		#	geom_text(data=subset(dat_outlier, isOutlier==TRUE), aes(label=sampleName))+
# 		#	geom_hline(yintercept = extremes, linetype=2, color='red')+
# 		#	geom_hline(yintercept = diff_mean, linetype=2, color='black')+
# 		#	xlab('Average measurement')+ylab('Difference between the two measurements')+ggtitle(gtitle)
# 		#print(p)
# 		#browser()
# 		gtitle <- str_c(tag, sprintf('%d disconcordant outliers at significance level of %.3f', sum(dat_outlier$isOutlier, na.rm=TRUE), alpha_outlier), collapse='')
# 		ylim <- c(min(c(extremes, dat_outlier$dif), na.rm=TRUE), max(c(extremes, dat_outlier$dif), na.rm=TRUE))
# 		with(dat_outlier, plot(avg, dif, cex=2, pch=1, ylim=ylim, xlab='Average measurement', ylab='Difference between the two measurements', main=gtitle))
# 		abline(h=extremes, lty=2, col='red')
# 		abline(h=diff_mean, lty=2, col='black')
# 		tt <- (data=subset(dat_outlier, isOutlier==TRUE))
# 		if(nrow(tt)>0){
# 			with(tt, textxy(avg, dif, labs=sampleName, cex=1.1, offset=0))
# 		}
# 		res <- list(CCCinfo=CCCinfo, dat_outlier=dat_outlier, diff_mean=diff_mean, diff_sd=diff_sd, indexOutlier=which(dat_outlier$isOutlier), indexNonOutlier=which(dat_outlier$isOutlier==FALSE))
# 	}	
# 	res
# }


#' plot a table with colour indicated
#'
#' Given a table object, this function plot the count and indicate the actual count
#'
#' modified on 2014/01/29: fisher test might fail; add try-catch (i.e. user chisq or NA); default xlab is related to tab dimnames now.
#' modified on 2014/01/30: add ordinal logistic regression info
#' modified on 2014/11/26: when chisq is successfully fitted, cell will be colored by standardized residual to understand which cell is enriched
#'
#' @param x a factor
#' @param y a factor
#' @param xlab xlab
#' @param ylab ylab
#' @param OrdReg logical indicating if an ordinal regression is to be fitted
#' @param plot logical indicating if a plot is to be generated
#' @param plotType type of plot. stdres for coloring by std residuals. count for coloring by count. dull for no coloring
#' @param verbose logical indicating whether print model fitting result on screen
#' @param tag the tag on the title; a fisher exact test p value will be appended
#' @param low.color the color for smallest value
#' @param high.color the color for largest value
#' @param text.size text size for the count
#' @examples
#' x <- sample(LETTERS[1:4], 100, replace=TRUE)
#' y <- sample(LETTERS[5:8], 100, replace=TRUE)
### plotHeatmapTable(x, y)	
#' @import plyr
#' @export
plotHeatmapTable <- function(x, y, xlab, ylab, tag=NULL, main=NULL, verbose=TRUE,
	OrdReg=TRUE, plot=TRUE, plotType='stdres', low.color='white', high.color='blue', text.size=8, stdres_range=c(-4, 4),
	sizex=12, sizey=12, sizetitlerel=0.8,sizetitleXrel=1, sizetitleYrel=1, # pars when type='jitter'
	anglex=NULL, angley=0){
	# retired: boundVofP=0.001, 
	#browser()
	tab <- table(y, x) # so as to make x on xaxis
	tnames <- names(dimnames(tab))
	if(missing(xlab)){
		xlab <- deparse(substitute(x))
	}
	if(missing(ylab)){
		ylab <- deparse(substitute(y))
	}
	mtab <- melt(tab)
	colnames(mtab) <- c("y", "x", "count") # for melt data, x is x from table(x, y)
	#mtab <- mutate(mtab, x=as.factor(x), y=as.factor(y)) # in case 0/1 will be interpretated as continous, giving 0.5 in the axis
	#browser()
	mtab$x <- orderFactor(mtab$x, new.order=levels(x), pmatch=FALSE)
	mtab$y <- orderFactor(mtab$y, new.order=levels(y), pmatch=FALSE)
	# in some special cases, Fisher test fails, i.e. Error in fisher.test(tab) : FEXACT error 6. LDKEY is too small for this problem. Try increasing the size of the workspace.
	fitP <- try(fisher.test(tab)$p.value, silent=TRUE)
	getChiOutMelt <- function(var='stdres'){
				#browser()
				tab_std <- csq[[var]]
				# trun value
				tab_std <- truncMat(tab_std, scale=F, Lower=stdres_range[1], Upper=stdres_range[2])
				#tab_std_m <- melt(tab_std[rownames(tab), colnames(tab)]) this is wrong
				tab_std_m <- melt(tab_std) # both table are table(x, y) order using melt, giving y, x, ...
				tab_std_m$value
	}
	if(class(fitP)=='try-error'){
		#csq <- chisq.test_boundV(tab, boundVofP=boundVofP)
		csq <- chisq.test(tab)
		fitP <- try(chisq.test(tab)$p.value, silent=TRUE)
		test <- 'Chi-squre'
		if(class(fitP)=='try-error'){
			fitP <- NA
		} else {
			# standardized residuals
			#tab_std <- chisq.test(tab)$stdres
			#tab_std_m <- melt(tab_std[rownames(tab), colnames(tab)])
			#browser()
			mtab$stdres <- getChiOutMelt(var='stdres')
			mtab$residuals <- getChiOutMelt(var='residuals')
			mtab$expected <- getChiOutMelt(var='expected')
		}
	} else {
		#browser()
		test <- 'Fisher exact'
		csq <- chisq.test(tab)
		# add stdres table
		#tab_std <- chisq.test(tab)$stdres
		#tab_std_m <- melt(tab_std[rownames(tab), colnames(tab)])
		#mtab$stdres <- tab_std_m$value
		mtab$stdres <- getChiOutMelt(var='stdres')
		mtab$residuals <- getChiOutMelt(var='residuals')
		mtab$expected <- getChiOutMelt(var='expected')
	}
	if(is.null(main)){
		main <- sprintf('%s test: P=%.3e', test, fitP)
	}
	if(!is.null(tag)) {
		main <- sprintf('%s\n%s test: P=%.3e', tag, test, fitP)
	}
	res <- list(test=test, Pval=fitP, mtab=mtab) # save test and Pval
	if(OrdReg){
		#browser()
		mydf <- data.frame(x=x, y=y)
		# deceive the variable name: to allow spaces etc
		xlab_ <- make.names(xlab); ylab_ <- make.names(ylab)
		colnames(mydf) <- c(xlab_, ylab_)
		fomula <- as.formula(sprintf('%s~%s', ylab_, xlab_)) # this is to preserve the names
		fit <- try(MASS::polr(fomula, data = mydf, Hess = TRUE), silent=TRUE)
		if(class(fit)!='try-error'){
			fitSmry <- summary(fit)
			ctable <- coef(summary(fit))
			p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
			## combined table
			coef <- cbind(ctable, "p value" = p)
			#res <- list(fit=fit, coef=coef)
			res$fit <- fit
			res$coef <- coef
		} else {
			#res <- list(fit=NA, coef=NA)		
			res$fit <- NA
			res$coef <- NA
		}
		## printing
		if(verbose){
			cat('Ordinal logistic regression\n')
			if(class(fit)!='try-error'){
			cat('########\n')
			cat('Warning\n')
			cat('########\n')
			cat('if the levels of the categorical variables are not set in a meaningful way, the result would be meaningless\n')
			cat('------------------------\n')
			cat('Model summary\n')
			print(fitSmry)
			cat('------------------------\n')
			cat('Coefficients and P values\n')
			print(coef)
			} else {
				cat('The model fails due to the nature of data. We catched the following information for further diagnosis\n')
				cat('------------------------\n')
				print(fit)
				cat('------------------------\n')
			}
		}
	}
	#browser()
	scale_axis_size <- theme(axis.text.x=element_text(size=sizex, colour = "black", angle=anglex, hjust=0.5, vjust=0.5), 
		axis.text.y=element_text(colour = "black", size=sizey, angle=angley),
		plot.title = element_text(size = rel(sizetitlerel)), 
		axis.title.y = element_text(size = rel(sizetitleYrel)),
		axis.title.x = element_text(size = rel(sizetitleXrel)))
	if(class(fitP)!='try-error' & plotType=='stdres'){
		# display residuals
		pal <- gplots:::bluered(31)
		p <- ggplot(mtab, aes(x=x, y=y, fill=stdres))+ 
			geom_tile()+
			geom_text(aes(label=count), size=text.size)+
			scale_fill_gradient2("Enrichment\n(Residual)", low=pal[1],mid=pal[17],high=pal[31],midpoint=0) + 
			xlab(xlab) + ylab(ylab) + ggtitle(main)+ 
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks=element_blank())
	} else if(plotType=='count'){
		# coloring by count
		p <- ggplot(mtab, aes(x=x, y=y, fill=count))+ 
			geom_tile()+
			geom_text(aes(label=count), size=text.size)+
			scale_fill_gradient("Count", low=low.color, high=high.color) + 
			xlab(xlab) + ylab(ylab) + ggtitle(main)+ 
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks=element_blank())
	} else {
		# dull no coloring
		p <- ggplot(mtab, aes(x=x, y=y, fill=count))+ 
			geom_tile()+
			geom_text(aes(label=count), size=text.size)+
			scale_fill_gradient("", low='white', high='white', guide=FALSE)+ 
			xlab(xlab) + ylab(ylab) + ggtitle(main)+ 
			theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.ticks=element_blank())
	}
	#theme(legend.position="top")
	p <- p+scale_axis_size
	res$p <- p
	#browser()
	if(plot==TRUE){
		return(p)
	} else {
		return(res)
	}	
}
#with(df_info, plotHeatmapTable(y=KLKPstatus, x=cl6_MUT_CN)) # cat-->cont

#with(df_info, plotHeatmapTable(cl6_MUT, cl6_MUT_CN))

#' plot cat-cat
#'
#' a wrappter to the plotHeatmapTable function
#' TODO: add options for type=c('barprop', 'barN', 'piefacet', 'residual')
#' @param x a factor
#' @param y a factor
#' @param xlab xlab
#' @param ylab ylab
#' @export
plotCatCat <- function(...){
	ll <- list(...)
	x <- ll[[1]]
	y <- ll[[2]]
	if(!is.factor(x) | !is.factor(y)) stop('both variable must be a factor!\n')
	#browser()
	plotHeatmapTable(...)
}
#with(df_info, plotCatCat(x=cl6_MUT, y=cl6_MUT_CN)) # cat-->cat
cat_cat_bar <- function(x, y, showProp=TRUE, xlab='', ylab=NULL, main='', legendTitle='', pallete="PRGn"){
	#browser()
	df <- data.frame(x=x, y=y)
	if(showProp){
		if(is.null(ylab)) ylab <- 'Percentage'
		dfs <- as.data.frame(with(df, prop.table(table(y, x), margin = 2)))
		p <- ggplot(dfs, aes(x, Freq, fill=y))+geom_bar(stat='identity', width=.5)+xlab(xlab)+ylab(ylab)+
		 ggtitle(main)+scale_fill_brewer(legendTitle, palette = pallete)+
		 theme_base(base_size = 14, base_family = "serif")	
	} else {
		if(is.null(ylab)) ylab <- '# of Samples'
		dfs <- as.data.frame(with(df, table(y, x), margin = 2))
		p <- ggplot(dfs, aes(x, Freq, fill=y))+geom_bar(stat='identity', width=.5)+xlab(xlab)+ylab(ylab)+
		 ggtitle(main)+scale_fill_brewer(legendTitle, palette = pallete)+
		 theme_base(base_size = 14, base_family = "serif")
	}
	p
}
cat_cat_pie <- function(x, y, showPropText=FALSE, xlab='', ylab=NULL, main='', legendTitle='', pallete="PRGn"){
	#browser()
	df <- data.frame(x=x, y=y)
	dfs <- as.data.frame(with(df, prop.table(table(y, x), margin = 2)))
	p0 <- ggplot(data = dfs, aes(x = factor(1), y = Freq, fill = y)) + 
		geom_bar(width = 1, stat = "identity") + scale_fill_brewer(palette =pallete) + 
		coord_polar(theta = "y") + facet_wrap( ~ x) + 
		xlab(xlab) + labs(fill = legendTitle)		
	if(showPropText){
		p <- p0+ scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) + 
		theme(panel.margin = unit(0.3, "lines"), legend.position = "right", 
  		axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), 
  		axis.ticks = element_blank(), panel.border = element_blank())
	} else {
		#browser()
		p <- p0 + blank_theme()+ 
			theme(strip.background=element_rect(colour="gray", fill="#ccd4e2"), axis.text.x=element_blank(), axis.text.y = element_blank())
	}
	p
}

#' guess angle for axis angle
#' @export
guess_angle <- function(catVec, catVecU=NULL, ncharMax=30, angle=45){
	#browser()
	if(!is.null(catVecU)) zz <- catVecU
	else zz <- as.character(unique(catVec))
	res <- ifelse(sum(nchar(zz))>=ncharMax, angle, 0)
	res
}

#' plot Cont-Cat variable 
#'
#' for cont~cat, we need a boxplot with jittered points. An ANOVA p value is needed. A summary of anova model can be available.
#' 
#' modified on 2014/01/29: x can be cont, y can be cat. This can be implemented by cord_flip()
#' modified on 2014/02/13: add the option ignoreNA. This happens when we want to remove the bar of NA values. Also, n_unique(x) should not count NA as a level when claiming a test name.
#' 	the reason is that when only 2 unique levels plus NA, we should say this is a t.test
#' 
#' @param x required to be a factor due to TukeyHSD(aov(y~x, data=mydf)) requires factor x. It is also necessary for the geom_boxplot otherwise a single box would be plotted when x is treated as continuous variable
#' @param y continuous vector
#' @param col a vector of color specifying the color of box and points. Color correspondance is through scale_color_manual(breaks=levels(x), values=col)
#' @param colI a vector of color of length(x) specifying the color of points. 
#' @param tag used to add a tag to the title
#' @param xlab xlab
#' @param ylab ylab
#' @param useRankTest whether to use wilcox rank test in place of t test; only relevant if only two categories is found. 
#' @param ignoreNA whether to exclude the bar of NAs in the data. Default is to leave it there. Notice this will not affect the P value
#' since NAs is not treated as a separate group but deleted
#' @param order whether order factor by median (so that to match with the median bar in boxplot) of y values
#' @param decreasing how the factor level are ordered, decreasing TRUE or FALSE
#' @param theme_bw logical if bw theme is imposed
#' @param barwidth barwidth when type='barplot'
#' @param barcol barcol when type='barplot'
#' @param cex.axis size of index text when type='barplot'
#' @param lx, ly legend position when type='barplot'
#' @param decreasing order direction of the bars
#' @param labels index text, e.g. cell line labels to show. default is integer index. 
#' @param near0 when type='barplot', bar at 0 will not be shown due to 0 length; we can use a bar hight of near0 to represent it. 
#' @param type plot type, either 'jitter' by geom_jitter (y is fixed) , 'dotplot' by geom_dotplot (y is binned), or 'barplot' by geom_barplot()
#' @param nybins number of bins on y axis if type=='dotplot'; this controls how y is grouped; dot size is determined through size; 
#' @param size points size
#' @param sizex xaxis text size
#' @param sizey yaxis text size
#' @param anglex xaxis text rotation
#' @param angley yaxis text rotation
#' @param shape points shape
#' @param width boxplot width; default is 1, which should be the ggplot2 default
#' @param plot if FALSE, no figure will be generated. If TRUE, figure and print will be added
#' @param seed jitter is random; add seed to make it reproducible
#' @param expand whether to expand the size and axis title (useful if figure is to be copied in slides with multiple panels)
#' @param base only relevant if categorical variable is binary which is used to calculate fold change. Default is Null, which will not enable FC calculation. 
#' @param FCbyRatio only active if two groups present, if specified as TRUE, need also to specify base=1 to be active
#' @param levels only relevant if categorical variable is binary which is used to calculate fold change and mean diff. Default is NULL, which will be using alphabetical order, handled by ttestP()
#' @param wrap_width_x use str_wrap to break long axis into multiple lines; default is not to break by setting this as Inf. 
#' @param saveTable whether to save anova result in a csv file
#' @param pathTable path to the saved csv file. default is NULL, which will be set as file.path(getwd(), 'OutputData/scatterTables')
#' @param print print text on console for test result
#' @export
plotContCat <- function(x, y, main, tag, col=NULL, colI=NULL, colVar=NULL, colVar_pal=NULL, 
	type='jitter', 
	xlab='', ylab='', useRankTest=FALSE, print=TRUE, 
	ignoreNA=TRUE, order=FALSE, decreasing_order=FALSE,  
	theme_bw=TRUE, ylim=NULL, jitter_w=0.2, nybins=100, 
	barwidth=0.2, barcol=NULL, cex.axis=NULL, lx=NULL, ly=NULL, decreasing=TRUE, labels=NULL, near0=0.01, showLegend=TRUE, # pars when type='barplot'
	width=0.8, size=4, sizeForce=NULL, sizex=12, wrap_width_x=Inf, sizey=12, sizetitlerel=0.8,sizetitleXrel=1, sizetitleYrel=1, # pars when type='jitter'
	anglex=NULL, angley=0, shape=20, 
	plot=TRUE, seed=1000, expand=TRUE,
	base=NULL, FCbyRatio=FALSE, levels=NULL, saveTable=FALSE, pathTable0=NULL, pathTable=NULL){
	set.seed(seed)
	#browser()
	if(missing(xlab))
		xlab <- deparse(substitute(x)) # deparse to make it a character so that it shows as it is ddf$clinical, not $(ddf, clinical)
	if(missing(ylab))
		ylab <- deparse(substitute(y))
	if(missing(tag)) {
		tag <- ''
	} else {
		tag <- sprintf('%s\n', tag) # add line break so that the user don't have to
	}
	#browser()
	# incase x is cont, y is cat, we use a coordinate flip to get same result
	flipped <- is.factor(y) & is.numeric(x)
	if(flipped){
		xold <- x
		yold <- y
		swap(x, y)
		swap(xlab, ylab)
		#nlevel <- n_unique(y, ignoreNA=TRUE)
	}
	nlevel <- n_unique(x, ignoreNA=TRUE) # unique levels should not count NA!
	if(is.null(col)) col <- rep('black', nlevel) # color is related to levels(x)
	mydf <- data.frame(x=x, y=y)
	# for label index in type='barplot'
	if(is.null(labels)) labels <- 1:nrow(mydf)
	mydf$labels <- labels
	#mydf <- mutate(mydf, cols=vec2colvec(vec=x, colpal=col))
	#browser()
	if(is.null(anglex)) {
		#browser()
		anglex <- ifelse(sum(nchar(as.character(unique(noNA(mydf$x)))))>16, 45, 0)
	}
	#browser()
	## when NA is to be deleted, we just extract complete cases
	if(ignoreNA){
		ind <- noNA(mydf, returnIndex=T)
		mydf <- mydf[ind, ]
	} else {
		# bug fix on 2014/02/20
		# need to conver NA to character 'NA'; otherwise, NA as missing will be just excluded
		NA2Char <- function(x) {
			x[is.na(x)] <- 'NA'
			x
		}
		mydf <- mutate(mydf, x=NA2Char(x))
	}
	aovP <- with(mydf, anovaP(expr=y, label=x, base=base))
	tukey <- with(mydf, TukeyHSD(aov(y~x, data=mydf)))
	#browser()
	if(plot & nlevel>2 & print){
		# only verbose the tukey result when plotting
		cat('Post test with Tukey HSD Test\n')
		print(tukey)
	}
	testName <- ifelse(nlevel>2, 'ANOVA', 't')
	#browser()
	if(missing(main)) {
			gtitle <- sprintf('%s%s test P value: %.3g\n', tag, testName, aovP['anovaP'])
			#browser()
			if(nlevel==2 & !is.null(base) & useRankTest==FALSE){ # base!=1 to have FC
				# FC is only available if nlevel=2 and base not null
				if(base!=1) {
					gtitle <- sprintf('%s\n%s test P value: %.3g, FC=%.3f\n', tag, testName, aovP['anovaP'], aovP['FC'])
				} else { 
					# use mean diff, base = 1
					ttestRes <- with(mydf, ttestP(data=y, class=x, levels=levels, base=base, more=T, FCbyRatio=FCbyRatio))
					#browser()
					if(FCbyRatio) {
						# show FC
						gtitle <- sprintf('%s\n%s test P value: %.3g, FC=%.3f\n', tag, testName, ttestRes['pval'], ttestRes['FC'])
					} else {
						# show mDiff
						gtitle <- sprintf('%s\n%s test P value: %.3g, meanDiff=%.3f, Effect size=%.3f\n', tag, testName, ttestRes['pval'], ttestRes['mDiff'], ttestRes['EffSize'])
					}	
				}	
			}	
			if(nlevel==2 & useRankTest==TRUE){
				res_rank <- wilcoxP(data=mydf$y, class=mydf$x)
				#browser()
				gtitle <- sprintf('%s%s test P value: %.3g\n', tag, 'Mann-Whitney rank', res_rank)
				# add FC
				if(!is.null(base)){
					gtitle <- sprintf('%s%s test P value: %.3g, FC=%.3f\n', tag, 'Mann-Whitney rank', res_rank, aovP['FC'])
				}
			}
	} else {
		gtitle <- main
	}	
	if(expand){
		#size <- 8
		sizex <- sizex*2
		sizey <- sizey*2
		# force dot size when EXPAND
		if(!is.null(sizeForce)) size <- sizeForce
		#if(nrow(mydf)>50) size <- 4
	}
	#browser()
	# scale_axis_size <- theme(axis.text.x=element_text(size=sizex, colour = "black", face ='bold', angle=anglex, hjust=0.5, vjust=0.5), 
	# 	axis.text.y=element_text(colour = "black", face ='bold', size=sizey, angle=angley),
	# 	plot.title = element_text(size = rel(sizetitlerel)), 
	# 	axis.title.y = element_text(size = rel(sizetitleYrel), face ='bold'),
	# 	axis.title.x = element_text(size = rel(sizetitleXrel), face ='bold'))
	scale_axis_size <- theme(axis.text.x=element_text(size=sizex, colour = "black", angle=anglex, hjust=0.5, vjust=0.5), 
		axis.text.y=element_text(colour = "black", size=sizey, angle=angley),
		plot.title = element_text(size = rel(sizetitlerel)), 
		axis.title.y = element_text(size = rel(sizetitleYrel)),
		axis.title.x = element_text(size = rel(sizetitleXrel)))
	#browser()
	#### p1 <- p0+scale_axis_size
	if(is.finite(wrap_width_x)){
		# wrap x text to form beautiful x-axis labels
		mydf0 <- mydf
		nlv <- str_wrap(levels(mydf0$x), width=wrap_width_x)
		mydf <- mutate(mydf, x=str_wrap(x, width=wrap_width_x))	
		mydf$x <- orderFactor(mydf$x, nlv)
		if(!is.null(names(col))) {
			# named col also fixed
			names(col) <- str_wrap(names(col), width=wrap_width_x)
		}
	}
	if(order==TRUE){
		mydf <- mutate(mydf, x=orderVecByVal(vec=x, val=y, func=function(x){median(x, na.rm=TRUE)}, decreasing=decreasing_order))
	}
	if(type=='jitter'){
		#browser()
		#geom_boxplot(outlier.size=0, width=width, aes(col=x), outlier.shape=NA)+
		if(!is.null(colI)){
		p0 <- ggplot(mydf, aes(x, y))+ylab(ylab)+xlab(xlab)+ggtitle(gtitle)+
			geom_boxplot(outlier.size=0, width=width, outlier.shape=NA)+
			geom_jitter(position=position_jitter(w=jitter_w, h=0), size=size, shape=shape, col=I(colI))
		} else if(!is.null(colVar)){
		p0 <- ggplot(mydf, aes(x, y))+ylab(ylab)+xlab(xlab)+ggtitle(gtitle)+
			geom_boxplot(outlier.size=0, width=width, outlier.shape=NA)+
			geom_jitter(aes(col=colVar), position=position_jitter(w=jitter_w, h=0), size=size, shape=shape)+
			scale_color_manual('', breaks=names(colVar_pal), values=colVar_pal)
		} else {
		p0 <- ggplot(mydf, aes(x, y))+ylab(ylab)+xlab(xlab)+ggtitle(gtitle)+
			geom_boxplot(outlier.size=0, width=width, aes(col=x), outlier.shape=NA)+
			geom_jitter(aes(col=x), position=position_jitter(w=jitter_w, h=0), size=size, shape=shape)+ 
			scale_color_manual(breaks=levels(x), values=col)
		}
		if(theme_bw) 
			p0 <- p0+theme_base()	
		if(is.null(lx)) # no legend
			p0 <- p0+theme(legend.position="none")	
		else
			p0 <- p0+theme(legend.position=c(lx, ly))	
		#browser()		
	} else if(type=='dotplot'){
		#browser()
		yrange <- range(mydf$y, na.rm=T)
		xrange <- c(1, n_unique(mydf$x))
		yxratio <- diff(yrange)/diff(xrange)
		#nybins <- 100
		# to ensure a consistent size of dot across figures, assume: ybinwidth * dotsize * yxratio *8000 = size; increasing size will increase actual size linearly 
		# however, this formula seems to miss something and does not work
		ybinwidth <- diff(yrange)/nybins
		#dotsize <- size/ybinwidth/3200/yxratio
		dotsize <- size/ybinwidth/44*yxratio
		#print(sprintf('%s; yxratio=%.2f; ybinwidth=%.3e; dotsize=%.3e\n', tag, yxratio, ybinwidth, dotsize))
		#browser()
		p0 <- ggplot(mydf, aes(x, y))+ylab(ylab)+xlab(xlab)+ggtitle(gtitle)+
			geom_boxplot(outlier.size=0, width=width, aes(col=x))+
			geom_dotplot(aes(col=x), binaxis = "y", stackdir = "center", dotsize=dotsize, binwidth=ybinwidth)+
			scale_color_manual(breaks=levels(x), values=col)+theme(legend.position="none")
		if(theme_bw) 
			p0 <- p0+theme_base()+theme(legend.position="none")		
		#geom_dotplot(aes(col=x), binaxis = "y", stackdir = "center", dotsize=size/ybinwidth/yxratio/8000, binwidth=ybinwidth)
	} else {
		# type == 'barplot'
		#browser()
		if(is.null(lx)){
			if(decreasing){
				lx <- 0.3; ly=0.2
			} else {
				lx <- 0.7; ly=0.2
			}
		}
		if(is.null(barcol)) barcol <- colorpalette2colvec('set1')[1:nlevel]
		if(is.null(cex.axis)) cex.axis <- guessCEX(length(mydf$x))
		mydf <- mutate(mydf, labels=orderFactor(labels, new.order=mydf$labels[order(mydf$y, decreasing=decreasing)]))
		#browser()
		ltitle <- xlab
		xlab <- ''
		if(!is.null(near0)){
			mydf[mydf$y==0, 'y'] <- near0
		}
		#browser()
		p0 <- ggplot(mydf, aes(labels, y, fill=x))+ylab(ylab)+xlab(xlab)+ggtitle(gtitle)+
			geom_bar(stat = "identity", position="dodge", width=barwidth)+
			scale_fill_manual(name=ltitle, breaks=levels(x), values=barcol)
		if(theme_bw) 
			p0 <- p0+theme_base()		
		p0 <- p0+theme(legend.position=c(lx, ly))
		if(!showLegend) p0 <- p0+theme(legend.position="none")
		## change scale
		scale_axis_size <- theme(axis.text.x=element_text(size=cex.axis*12, colour = "black", face ='bold', angle=90, hjust=1, vjust=0.5), 
		axis.text.y=element_text(colour = "black", face ='bold', size=sizey, angle=angley),
		plot.title = element_text(size = rel(sizetitlerel), face ='bold'), 
		axis.title.y = element_text(size = rel(sizetitleYrel)),
		axis.title.x = element_text(size = rel(sizetitleXrel)))
	}	
	p1 <- p0+scale_axis_size
	#browser()
	if(flipped)
		p1 <- p1+coord_flip()
	if(!is.null(ylim))
		p1 <- p1+ylim(ylim)
	### output table when necessary
	#browser()
	if(saveTable){
		#pathTable = pathTable0 + 'OutputData'
		if(is.null(pathTable0)) pathTable0 <- file.path(getwd()) # path
		if(is.null(pathTable)) pathTable <- file.path(pathTable0, 'scatterTables')
		createFolder(pathTable)
		#browser()
		write.csv(aovP, file.path(pathTable, make.names(sprintf('%s.csv', tag))))
	}
	# plot figure or return result
	if(plot){
		return(p1)
	} else {
		if(useRankTest){
			res <- list(anovaP=aovP, tukey=tukey, p=p1, testName=testName, gtitle=gtitle, rankTestP=res_rank)	
		} else {
			res <- list(anovaP=aovP, tukey=tukey, p=p1, testName=testName, gtitle=gtitle)
		}
		return(res)
	}	

}
#plotContCat(x=factor(pData$pathology), y=ic50_aug[1, ], ignoreNA=TRUE, plot=F, anglex=30)$p

#plotContCat(x=IC50group, y=mat_diff[1, ], ignoreNA=TRUE)
#plotContCat(x=ddf$clinical_stage_formated, y=myy, tag=sprintf("%s, %s", myyName, datName))

#with(df_info, plotContCat(y=as.factor(cl6_MUT_CN), x=NmutationInFM))
#plotContCat(x=Tstage, y=AXLmRNA)

#' Calculate statistics for continuous variable x and y
#' @param x variable 1
#' @param y variable 2
#' @return a list
#' @export
statContCont <- function(x, y){
	#browser()
	cor_ps <- corTest(x, y, method='pearson', minN=5)
	cor_sp <- corTest(x, y, method='spearman', minN=5)
	lmres <- lm_xy(x, y)
	list(corRes_pearson=cor_ps, corRes_spearman=cor_sp, lmRes=lmres)
}
#statContCont(x=fragSum_iso, y=rgInfo$total_mass)

#' plot Cont-Cont variable 
#'
#' Scatter plot for continuous variables and associated statistics
#'
#' for cont~cont, we need a scatter plot. An lm fit with p value and R^2. A trend fit with smoothing can be added. Spearman corr. and Pearson Corr.
#' 
#' @param x continuous vector
#' @param y continuous vector
#' @param tag used to add a tag to the title
#' @param xlab xlab
#' @param ylab ylab
#' @param alpha alpha transparancy for the points
#' @param size size passed to geom_point()
#' @param sizeExpand size when expanded
#' @param shape shape passed to geom_point()
#' @param theme_bw logical if bw theme is imposed
#' @param lm if to add a lm fit line
#' @param smooth if to add a smooth line
#' @param sizelab size for label under each data point
#' @param sizeaxis axis text size (tick text)
#' @param sizeaxistitlerel relative axis title size (tick text)
#' @param sizetitlerel relative size for main title
#' @param ylim ylim
#' @param xlim xlim
#' @param expand whether to expand the size and axis title (useful if figure is to be copied in slides with multiple panels)
#' @param labels a vector of IDs to be added to each data point
#' @param plot.margin.cm the right margin might be too small to show full axis text; this specifies margins in top, right, bottom, left
#' @export
plotContCont <- function(x, y, main=NULL, tag='',
	xlab='', ylab='',
	theme_bw=TRUE,
	alpha=1, size=5, sizeExpand=6, shape=20, sizeaxis=12, sizeaxistitlerel=1.2, sizetitlerel=0.8,
	lm=FALSE, lmcol='blue2', lmsize=2, lmlty=2, smooth=FALSE, 
	plot=TRUE, ylim=NULL, xlim=NULL, labels=NULL, sizelab=3, colab='purple', vjust=1.9, hjust=0.5, expand=NULL,
	plot.margin.cm = c(0.2,1.2,0.2,0.2)){
	if(is.null(expand)) {
		if(testObject(EXPAND)) {
			expand <- EXPAND # use global variable
		} else {
			expand <- TRUE
		}
	}
	if(missing(xlab))
		xlab <- deparse(substitute(x))
	if(missing(ylab))
		ylab <- deparse(substitute(y))
	if(!is.null(labels) & length(x)!=length(labels)) 
		stop('Length of labels not equal to length of x!')
	if(is.null(labels)) 
		text <- 1:length(x)
	else 
		text <- labels	
	addText <- ifelse(is.null(labels), FALSE, TRUE)
	#browser()
	# make sure the ylim focus on effective x and y
	td <- data.frame(x=x, y=y, label=text)
	td <- td[noNA(td, returnIndex=TRUE), ]
	if(is.null(ylim)) ylim <- range(pretty(td$y), na.rm=TRUE)
	if(is.null(xlim)) xlim <- range(pretty(td$x), na.rm=TRUE)
	#browser()
	#mydf <- data.frame(x=x, y=y) # before 2015/1015
	mydf <- td
	statL <- statContCont(x, y)
	if(plot){
		# only verbose the tukey result when plotting
		cat('Linear model fit:\n')
		smry <- try(summary(lm(y~x)), silent=TRUE)
		print(smry)
	}
	if(is.null(main)) {
			#gtitle <- sprintf('%s\n Linear model: P=%.3g, \n Coefficient= %.3g, Adjusted R Square=%.3g\n Pearson Correlation: corr=%.3g, P=%.3e\n Spearman Correlation: corr=%.3g, P=%.3e\n', 
			#	tag, statL$lmRes['pval'], statL$lmRes['coef'], statL$lmRes['adjRsquare'], 
			#	statL$corRes_pearson['cor'], statL$corRes_pearson['pval'],
			#	statL$corRes_spearman['cor'], statL$corRes_spearman['pval'])
			# gtitle <- sprintf('%s\n Pearson Correlation: corr=%.3g, P=%.3g\n Spearman Correlation: corr=%.3g, P=%.3g\nLinear model: P=%.3g, Adjusted R Square=%.3g', 
			# 	tag, 
			# 	statL$corRes_pearson['cor'], statL$corRes_pearson['pval'],
			# 	statL$corRes_spearman['cor'], statL$corRes_spearman['pval'],
			# 	statL$lmRes['pval'], statL$lmRes['adjRsquare'])
			gtitle <- sprintf('%s\n Pearson Correlation: corr=%.3g, P=%.4g\n Spearman Correlation: corr=%.3g, P=%.4g\nAdjusted R Square=%.3g', 
				tag, 
				statL$corRes_pearson['cor'], statL$corRes_pearson['pval'],
				statL$corRes_spearman['cor'], statL$corRes_spearman['pval'],
				statL$lmRes['adjRsquare'])
		} else {
		gtitle <- main
	}	
	if(expand) {
		size <- sizeExpand
		sizeaxis <- 24
	}
	#browser()
	scale_axis_size <- theme(axis.text=element_text(size=sizeaxis, colour = "black", face ='bold', hjust=0.5, vjust=0.5), 
		axis.title.x = element_text(size = rel(sizeaxistitlerel)),
		axis.title.y = element_text(size = rel(sizeaxistitlerel)),
		plot.title = element_text(size = rel(sizetitlerel)))
	p0 <- ggplot(mydf, aes(x, y))+ylab(ylab)+xlab(xlab)+ggtitle(sprintf('%s\n', gtitle))+ylim(ylim)+xlim(xlim)+
		geom_point(alpha=alpha, size=size, shape=shape)
	if(lm) {
	p0 <- p0+stat_smooth(method = "lm", se=FALSE, col=lmcol, size=lmsize, lty=lmlty)
	}
	if(smooth){
	p0 <- p0+stat_smooth(col='blue')
	}
	if(theme_bw) 
			p0 <- p0+theme_base()		
	#browser()
	if(addText){
		p0 <- p0+geom_text(aes(x, y, label=label), hjust=hjust, vjust=vjust, size=sizelab, colour=colab)
	}
	#browser()
	p0 <- p0+theme(plot.margin = unit(plot.margin.cm, "cm")) 
	p1 <- p0+scale_axis_size#+ theme(plot.title = element_text(family='mono', hjust = 0))
	if(plot){
		return(p1)
	} else {
		return(list(statContCont=statL, p=p1))
	}	
}
#plotContCont(x=fragSum_iso, y=rgInfo$total_mass)

#with(df_info, plotContCont(x=packYear, y=NmutationInFM))



#' flexible scatter plot (a wrapper)
#'
#' scatter for a cont-cont, cont-cat or cat-cat variable pair. 
#'
#' the program detects x, y as categorical if it is a factor. Then it wraps plotContCont plotContCat plotCatCat
#' @param x  either numeric or factor or numeric
#' @param y  either numeric or factor or numeric
#' @param ... other parameters passed to either plotContCont plotContCat plotCatCat
#' @export
plotScatter <- function(...){
	#browser()
	ll <- list(...) # a trick so that substitute(x) != 'x'; we need to avoid passing x!
	xnumeric <- is.numeric(ll[[1]]) # make sure x= some var is passed; otherwise, it is NULL!
	ynumeric <- is.numeric(ll[[2]])
	#browser()
	if(xnumeric){
		# x numeric
		if(ynumeric){
			plotContCont(...)
		} else {
			plotContCat(...)
		}
	} else {
		# x categorical
		if(ynumeric){
			plotContCat(...)
		} else {
			#browser()
			plotCatCat(...)
		}
	}
}
#plotScatter(x=fragSum_iso, y=rgInfo$total_mass)

#with(df_info, plotScatter(x=cl6_MUT, y=packYear)) # cont-->cont
#with(df_info, plotScatter(x=packYear, y=cl6_MUT_CN)) # cont-->cont

##' multiple scatter plot for a given data frame
##'
##' not implemented due to technical difficulty in passing arguments in a nested function wrapper!
##plotScatterMultiple <- function(dat, xvar, yvar, ...){
##


#
# another wrapper that unifies plotScatter and plotKM
#


#' smart library function
#'
#' load a library; when not present, install it
#'
#' @param ... a character string, no need to use c()
#' @param mute whether to mute
#' @export
mylibrary <- function(..., mute = FALSE) {
#browser()
	packChars <- deparse(substitute(...))
    temp <- sapply(packChars, function(x) {
     if(!require(x, quietly = mute, character.only = TRUE))
          install.packages(x, quiet = mute)
     library(x, character.only = TRUE, quietly = FALSE)
  })
}
#mylibrary(googleVis)





###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################
###############################################################################################################################





#' plot mean vs sd for a matrix for quick QC/filtering
#'
#' @param mat a matrix
#' @param type compute summary on rows or columns
#' @param useMAD whether to use MAD (more robust) rather than SD. Default is FALSE
#' @param alpha alpha blending for the points
#' @param qx quantile to be show for Mean (0~1)
#' @param qy quantile to be show for SD (0~1)
#' @param main main title
#' @param quantileCol line colors for quantiles
#' @param plot plot or return the summary
#' @export
plotMeanSD4mat <- function(mat, type='row', useMAD=FALSE, alpha=NA, qx=c(0.3, 0.5), qy=c(0.3, 0.5), main=NA, quantileCol='#B0AEA9', plot=TRUE){
	if(type!='row'){
		mat <- t(mat)
	}
	m <- apply(mat, 1, mean, na.rm=TRUE)
	if(useMAD){
	ylab <- 'MAD'
	sd <- apply(mat, 1, mad, na.rm=TRUE)
	} else {
	ylab <- 'SD'
	sd <- apply(mat, 1, sd, na.rm=TRUE)
	}
	if(is.na(alpha)) alpha <- guessAlpha(length(m))
	if(is.na(main)) {
		input <- deparse(substitute(mat))
		p <- plotScatter(x=m, y=sd, xlab='Mean', ylab=ylab, lm=FALSE, smooth=TRUE, plot=FALSE, alpha=alpha, tag=input)$p
	} else {
		p <- plotScatter(x=m, y=sd, xlab='Mean', ylab=ylab, lm=FALSE, smooth=TRUE, plot=FALSE, alpha=alpha, main=main)$p
	}
	if(!is.na(qx)[1]){
		xs <- quantile(m, prob=qx, na.rm=TRUE); names(xs) <- qx
		myy <- max(sd, na.rm=TRUE)*seq(0.2, 0.9, length.out=length(xs))
		p <- p+geom_vline(xintercept = xs, linetype=2, color=quantileCol)+
		 annotate("text", x = xs, y = jitter(myy), label = str_c('q=', qx), angle=-90, colour='blue')
	}
	if(!is.na(qy)[1]){
		ys <- quantile(sd, prob=qy, na.rm=TRUE); names(ys) <- qy
		myx <- max(m, na.rm=TRUE)*seq(0.2, 0.9, length.out=length(ys))
		p <- p+geom_hline(yintercept = ys, linetype=2, color=quantileCol)+
		 annotate("text", x = jitter(myx), y = ys, label = str_c('q=', qy), angle=0, colour='blue') 
	}
	if(plot) {
	res <- p
	} else {
		res <- list(p=p, Mean=m, SD=sd, qx=xs, qy=ys)
	}
	res
	#browser()
}
#' select genes based on mean~SD plot from plotMeanSD4mat
#' @param ll the list returned by plotMeanSD4mat()
#' @param m mean cutoff, usually provided in ll
#' @param s SD cutoff, usually provided in ll
#' @return a data frame with Index, Mean and SD. Index tells what rows are selected
#' @export
selectGeneFromMeanSD <- function(ll, m, s){
	#browser()
	mydf <- with(ll, data.frame(Index=1:length(Mean), Mean=Mean, SD=SD))
	if(!missing(m)) rdf <- subset(mydf, Mean>=m)
	if(!missing(s)) rdf <- subset(mydf, SD>=s)
	rdf
}
#sdf <- selectGeneFromMeanSD(pp1, m=pp1$qx['0.5'])