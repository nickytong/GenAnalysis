##
## this script is downloaded from https://github.com/kmiddleton/rexamples/blob/master/qplot_survival.R
## the purpose is to plot survival curve. It is better than the ususal plot(surfit(Surv(...))) in the 
## sense that we do not need to specify the color and legend ourselves!
##

### create survival/cumulative incidence plot
require(survival)
require(ggplot2)
#####
# this function is the core of adding N to the legends; 
# WARNING: this is not proved to be correct; actually in practive,there might be a difference of 1 or 2!!!
#			should find a more error prone way to do it
#####
# define custom function to create a survival data.frame
createSurvivalFrame0 <- function(f.survfit){
  # initialise frame variable
  f.frame <- NULL
  
  # check if more then one strata
  if(length(names(f.survfit$strata)) == 0){ # no strata?
    # create data.frame with data from survfit
    f.frame <- data.frame(time=f.survfit$time, 
                          n.risk=f.survfit$n.risk, 
                          n.event=f.survfit$n.event, 
                          n.censor = f.survfit$n.censor, 
                          surv=f.survfit$surv, 
                          upper=f.survfit$upper, 
                          lower=f.survfit$lower)
    # create first two rows (start at 1)
    f.start <- data.frame(time=c(0, f.frame$time[1]), 
                          n.risk=c(f.survfit$n, f.survfit$n), 
                          n.event=c(0,0), 
                          n.censor=c(0,0), 
                          surv=c(1,1), 
                          upper=c(1,1), 
                          lower=c(1,1)) 
    # add first row to dataset
    f.frame <- f.frame_raw <- rbind(f.start, f.frame)
    
    # remove temporary data
    rm(f.start)
  } 
  else {
    # create vector for strata identification
    f.strata <- NULL
    for(f.i in 1:length(f.survfit$strata)){
      # add vector for one strata according to number of rows of strata
      f.strata <- c(f.strata, rep(names(f.survfit$strata)[f.i], f.survfit$strata[f.i]))
    }
    # create data.frame with data from survfit (create column for strata)
    # f.frame_raw : nothhing added yet
    f.frame <- f.frame_raw <- data.frame(time=f.survfit$time, n.risk=f.survfit$n.risk, n.event=f.survfit$n.event, n.censor = f.survfit$n.censor, 
		surv=f.survfit$surv, upper=f.survfit$upper, lower=f.survfit$lower, strata=factor(f.strata))
    
    # remove temporary data
    #rm(f.strata)
    # create first two rows (start at 1) for each strata
    for(f.i in 1:length(f.survfit$strata)){
      
      # take only subset for this strata from data
      f.subset <- subset(f.frame, strata==names(f.survfit$strata)[f.i])
      
      # create first two rows (time: 0, time of first event)
      f.start <- data.frame(time=c(0, f.subset$time[1]), n.risk=rep(f.survfit[f.i]$n, 2), n.event=c(0,0), n.censor=c(0,0), surv=c(1,1), upper=c(1,1), lower=c(1,1), strata=rep(names(f.survfit$strata)[f.i],2))	
      
      # add first two rows to dataset
      f.frame <- rbind(f.start, f.frame)
      
      # remove temporary data
      #rm(f.start, f.subset)
      
    }
    
    # reorder data
    f.frame <- f.frame[order(f.frame$strata, f.frame$time), ]
    #browser()
    # rename row.names
    rownames(f.frame) <- NULL   
  }
  #browser()
  # return frame
  res <- list(f.frame=f.frame, f.frame_raw=f.frame_raw )
  return(res)
}
createSurvivalFrame <- function(f.survfit){
    createSurvivalFrame0(f.survfit)$f.frame
}
# define custom function to draw kaplan-meier curve with ggplot
qplot_survival <- function(f.frame, f.CI="default", f.shape=3){  
  # use different plotting commands dependig whether or not strata's are given
  if("strata" %in% names(f.frame) == FALSE){
    
    # confidence intervals are drawn if not specified otherwise
    if(f.CI=="default" | f.CI==TRUE ){
      #browser()
      # create plot with 4 layers (first 3 layers only events, last layer only censored)
      # hint: censoring data for multiple censoring events at timepoint are overplotted
      # (unlike in plot.survfit in survival package)
      ggplot(data=f.frame) + geom_step(aes(x=time, y=surv), direction="hv") + geom_step(aes(x=time, y=upper), direction="hv", linetype=2) + geom_step(aes(x=time,y=lower), direction="hv", linetype=2) + geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)
      
    }
    else {
      
      # create plot without confidence intervalls
      ggplot(data=f.frame) + geom_step(aes(x=time, y=surv), direction="hv") +  geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)
      
    }
    
  }
  else {
    
    if(f.CI=="default" | f.CI==FALSE){
      
      # without CI 
      ggplot(data=f.frame, aes(group=strata, colour=strata)) + geom_step(aes(x=time, y=surv), direction="hv") + geom_point(data=subset(f.frame, n.censor>0), aes(x=time, y=surv), shape=f.shape)
      
    }
    else {
      # with CI (hint: use alpha for CI)
      ggplot(data=f.frame, aes(colour=strata, group=strata)) + geom_step(aes(x=time, y=surv), direction="hv") + geom_step(aes(x=time, y=upper), directions="hv", linetype=2, alpha=0.5) + geom_step(aes(x=time,y=lower), direction="hv", linetype=2, alpha=0.5) + geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape) 
    } 
  } 
}

#' plot just KM
#'
#' no covariate KM plot
#'
#' @param OS Overall survival
#' @param OScensoring censoring status
#' @param f.CI whether to add C.I.
#' @param xlab x label
#' @param ylab y label 
#' @param main main title
#' @param f.shape shae
#' @param legendText text for legend
#' @param lx x position of legend
#' @param ly y position of legend
#' @param hcolor horizontal line color 
#' @param showAnno whether to show annotation
#' @export 
plotKM0 <- function(OS, OScensoring, f.CI=TRUE, xlab='Survival time', f.shape=3, 
        ylab='Survival Probability', main='Kaplan-Meier Curve', hcolor='#a2abb5',
        legendText='Death status', tag='', lx=0.2, ly=c(0.1, 0.05), showAnno=TRUE, cex.main=0.8){
    sfit <- survfit(Surv(OS, OScensoring)~1)
	  sframe <- createSurvivalFrame(sfit)
    main <- str_c(tag, '\n', main)
    #browser()
    p <- qplot_survival(sframe, f.CI=f.CI, f.shape=f.shape)+geom_hline(yintercept = 0.5, linetype=2, color=hcolor)+
        xlab(xlab)+ylab(ylab)+ggtitle(main)+ylim(c(0,1))+theme(plot.title = element_text(size = rel(cex.main), colour = "black"))
    if(showAnno){
        tt <- table(OScensoring)
        legendText <- str_c(legendText, ' ', names(tt), ': N=', tt)
        #browser()
        rr <- range(OS, na.rm=TRUE)
        xp <- rr[1]+diff(rr)*lx
        yp <- ly
        p <- p+ annotate("text", x = xp, y = ly, label = legendText)
    }
    p
}
#with(clin1u, plotKM0(OS=OS, OScensoring=OScensoring, lx=0.15))
# # create frame from survival class (survfit)
# t.survfit <- survfit(t.Surv~1, data=lung)
# t.survframe <- createSurvivalFrame(t.survfit)
# 
# # create kaplan-meier-plot with ggplot
# qplot_survival(t.survframe)
# 
# 
# # drawing survival curves with several strata
# t.Surv <- Surv(lung$time, lung$status)
# t.survfit <- survfit(t.Surv~sex, data=lung)
# plot(t.survfit)
# 
# # two strata
# t.survframe <- createSurvivalFrame(t.survfit)
# qplot_survival(t.survframe)
# # with CI
# qplot_survival(t.survframe, TRUE)
# # add ggplot options, use different shape
# qplot_survival(t.survframe, TRUE, 1) + theme_bw() + scale_colour_manual(value=c("green", "steelblue")) + opts(legend.position="none")
# 
# # multiple stratas
# t.survfit <- survfit(t.Surv~ph.karno, data=lung)
# t.survframe <- createSurvivalFrame(t.survfit)
# qplot_survival(t.survframe)
# # plot without confidence intervals and with different shape
# qplot_survival(t.survframe, FALSE, 20)

## add coxph for continuous and logrank test for binary category.

#' plot OS against a continuous variable by breaking the continuous variable
#'
#' for a covariate x and survival, break x into groups and plot KM curves
#'
#' updated on 2014/03/13: when x is missing, this means just plot KM for OS. 
#' This implementation is built upon plotKMcat. That is, we first breaks the continuous variable to categorical and feed it into plotKMcat().
#' How is the color decided; how does it corresponds to the broken intervals? This is actually decidied by plotKMcat() which use aes(group, colour) to
#' make sure color and group are correct. The labels in categorical variable is the label in legend to make sure we're safe. If one want to change this,
#' just manipulate the categorical variable. For the continuous variable, we can just specify the labels in the order of resulting increasing intervals from cut().
#'
#' updated on 2016/01/08: breaks and labels works with cut and good for >2 groups; However, sometimes dichotomized group for extreme values is needed (using dichotomizeExpr).
#' To fully use use_dichotomizeExpr, we add use_dichotomizeExpr to indicate if this function is activated and use par_dichotomizeExpr as a list to hold all parameters passed to it.
#' @param x a continuous variable
#' @param OS survival time
#' @param OScensoring censoring status, 1 is event
#' @param main title
#' @param tag a tag added at the start of title
#' @param xlab x label
#' @param ylab y label
#' @param use_dichotomizeExpr whether to use dichotomizeExpr to dichotomize the data; when this is activated, cut is inactivated
#' @param par_dichotomizeExpr a list of parameters passed to dichotomizeExpr
#' @param labels labels for each broken categories, corresponding to values (intervals) from small to large
#' @param legendTitle legend title name
#' @param nrowL legend row number
#' @param theme_bg background theme
#' @param cols either a named vector (names from x) to accurately specify the color or a color vecto which is to be input in the alphabetical order as the x values
#' @param legend_within whether to put legend within the figure
#' @param lx relative x position of legend (between 0~1)
#' @param ly relative y position of legend (between 0~1)
#' @param addEventCount addEventCount whether to add event count in the legend
#' @export 
ggKM <- function(x, OS, OScensoring, breaks, labels, use_dichotomizeExpr=FALSE, par_dichotomizeExpr=list(method='quantile', q1=NA, q2=NA, Lower=NULL, Upper=NULL, labels=c('Low', 'High'), level_labels=c('Low', 'High'), plot=FALSE), legendTitle, nrowL=NA, main, cex.main=0.8, tag,
	xlab='Survival time', ylab='Survival probability',
	cols, legend_within=TRUE, lx=0.8, ly=0.88, xlim=NULL, 
	theme_bg=theme_base(), sizeaxis=12, sizeaxistitlerel=1.2, sizetitlerel=0.8, 
	plot=TRUE, addEventCount=F, rmMissingEarly=TRUE){
	#browser()
  # dichotomization before removing missing value or after matters!
  # on 2018/02/22, we start removing NA value before everything starts. This Robert thinks is less confusing
  # since median cutoff will always be symmetrical. 
  if(missing(x)){
        p <- plotKM0(OS=OS, OScensoring=OScensoring, cex.main=cex.main)
        return(p)
  }
	if(rmMissingEarly){
      tp <- noNA(data.frame(x=x, OS=OS, OScensoring=OScensoring))
      x <- tp$x
      OS <- tp$OS
      OScensoring <- tp$OScensoring
  }
  #browser()
  if(missing(breaks)){
		breaks <- c(-Inf, median(x, na.rm=TRUE), Inf)
		#if(missing(labels)){
			#labels <- c('Low', 'High') # alphabetical order
		#	stop("Need to specify labels")
		#}
	}	
	if(missing(legendTitle)){
			legendTitle <- substitute(x)
	}
	#browser()
	if(missing(tag)){
			tag <- sprintf('KM Curve for %s', deparse(substitute(x)))
	}
	# for missing labels, use the intervals
	if(missing(labels)){
		labels <- levels(cut(x, breaks=breaks))
	}
	#browser() # Error in cut.default(x, breaks = breaks, labels = labels) : lengths of 'breaks' and 'labels' differ
  if(!use_dichotomizeExpr){
    xGrouped <- cut(x, breaks=breaks, labels=labels)
    nCurve <- length(breaks)-1
  } else {
    xGrouped0 <- with(par_dichotomizeExpr, dichotomizeExpr(x, method=method, q1=q1, q2=q2, Lower=Lower, Upper=Upper, labels=labels, plot=plot))
    xGrouped <- orderFactor(xGrouped0, new.order=par_dichotomizeExpr$level_labels)
    nCurve <- 2
  }
	if(missing(cols)){
		cols <- brewer.pal(8, "Set1")[1:nCurve]
		## override the color for aethestics when 2 curves
		if(nCurve==2) cols <- c("red1","blue3")
	}
	#browser()
	# p value without dichotomizing
	coxfit <- coxphP2(expr=x, OS=OS, OScensoring=OScensoring)
	contP <- coxfit$P
	paredText <- parse_coxphP2(coxfit, type='HR')
	## fit on grouped
	coxfitGrouped <- coxphP2(expr=xGrouped, OS=OS, OScensoring=OScensoring)
	groupedP <- coxfitGrouped$P
	paredTextGrouped <- parse_coxphP2(coxfitGrouped)
	#browser()
	#browser()
	if(missing(main)) {
		gtitle <- sprintf('%s\nLogrank test (overall) P: %.3g\n%s\n%s\ncox PH on continuous scale P: %.3g\n%s, %s', tag, groupedP, 
				paredTextGrouped$coefExp_text, paredTextGrouped$CI_text, contP, paredText$coefExp_text, paredText$CI_text)
	} else {
	gtitle <- main
	}
	# below gtitle has everything, so disable tag
	p1 <- ggKMcat(x=xGrouped, OS=OS, OScensoring=OScensoring, cols=cols, legend_within=legend_within, lx=lx, ly=ly, nrowL=nrowL, addEventCount=addEventCount, legendTitle=legendTitle, theme_bg=theme_bg, tag='', plot=FALSE, cex.main=cex.main, sizeaxis=sizeaxis, sizeaxistitlerel=sizeaxistitlerel, sizetitlerel=sizetitlerel)$p+ggtitle(gtitle)+xlab(xlab)+ylab(ylab)
	if(!is.null(xlim))
    p1 <- p1+xlim(xlim)
  #browser()
  if(plot){
		cat('------------------------\n')
		cat('Cox PH fit on continuous covariate:\n')
		print(coxfit$fit)
		#p0 <- qplot_survival(sframe, f.CI=F)+ ggtitle(gtitle)+ylim(c(0,1))+xlab(xlab)+ylab(ylab)
		#browser()
		return(p1)
	} else {
		#browser()
    return(list(groupedP=groupedP, contP=contP, p=p1, 
				coxfitGrouped=coxfitGrouped, coxfit=coxfit, 
				paredText=paredText, paredTextGrouped=paredTextGrouped))
	}
}

#ggKM(x=x, OS=OS, OScensoring=OScensoring, breaks=breaks, labels=labels, tag=tag, cols=cols)

#plotKM(x=myy, OS=ddf$OS, OScensoring=ddf$OScensoring, breaks=c(-Inf, median(myy, na.rm=TRUE), Inf), 
#	labels=c('low', 'high'), tag=sprintf("%s High vs Low, %s", myyName, datName), legendTitle=sprintf('%s Low/High', myyName))

#plotKM(x=myy, OS=ddf$OS, OScensoring=ddf$OScensoring, breaks=c(-Inf, median(myy, na.rm=TRUE), Inf), labels=c('low', 'high'), tag=sprintf("%s, %s", myyName, datName), legendTitle=sprintf('%s Low/High', myyName))

#plotKM(x=myy, OS=ddf$OS, OScensoring=ddf$OScensoring, breaks=c(-Inf, median(myy, na.rm=TRUE), Inf), labels=c('low', 'high'), tag=sprintf("%s, %s", myyName, datName), legendTitle=sprintf('%s Low/High', myyName))

#with(ddf, plotKM(x=AXL_rppa, OS=OS, OScensoring=OScensoring, breaks=c(-Inf, median(AXL_rppa, na.rm=TRUE), Inf), labels=c('low', 'high'), tag=sprintf("AXL RPPA, %s", datName), legendTitle='AXL Low/High'))

#ggKM(x=AXLmRNA, OS=os_MDACC_HN$OS, OScensoring=os_MDACC_HN$OScensoring, breaks=c(-Inf, median(AXLmRNA), Inf), labels=c('low', 'high'))

#' guess nrowL from a specified number of categories
#'
#' @param nCategory number of categories
#' @param useMax whether to use maximum value
guess_nrowL <- function(nCategory, useMax=FALSE){
    if(useMax) return (nCategory)
    if(nCategory==2) {
        return (1)
    } else if (nCategory==3) {
        return (2)
    } else {
        return (4)
    }
}
	
#' for categorical covariate x, plot KM curve
#'
#' This wraps qplot_survival() function for KM curve generation.
#' 
#' Note no sample size indicated in text: NA omitted in x, the strata var, however, OS can still be missing
#' so this is nominal N, effective N due to OS missing is not computed; this is good for tracking sample size and avoid confusion. We just
#' need to bear in mind that some samples will miss OS. 
#'
#' modified 2014/02/14: when x has only one unique value, get error: Error in `$<-.data.frame`(`*tmp*`, "strata", value = character(0)) : 
#'  replacement has 0 rows, data has 68
#' What color corresponds to what label? The labels by default are extracted from levels(x); The cols and computed labels (after
#' adding N) are 1-to-1 mapping. So it is a matter wheather the color is pointing to the same curve, which depends on qplot_survival() function.
#' I cannot prove qplot_survival() gets the color correct. But it seems work for categorical variable at this moment.
#' Updates: qplot_survival() computed sframe which is a df format for ggplot2 where a variable strata is formed which makes it possible
#' to use aes(group=strata, colour=strata). By default, the labels in strata is in alphabetical order;  
#' ------------------->
#' Now, ggKMcat does not specify the labels. Instead, if one wants to change the labels, just rename the values in x through mapNames. 
#' Since the color and curves are created by aes(group, colour), there is always correct correspondance between curve and colour. The cols can be specified through a named
#' vector so that one can accurately control the color for each curve. 
#' To do: add N to the legend
#'
#' To check with traditional surv curve, do:
#' 		mysurv <- with(clin1u, survfit(Surv(PFS, PFScensoring) ~ GENDER_PHENO)) 
#		plot(mysurv, lty = 2:3) 
#' @param x a factor
#' @param OS survival time
#' @param OScensoring censoring status, 1 is event
#' @param main title
#' @param tag a tag added at the start of title
#' @param xlab x label
#' @param ylab y label
#' @param labels this is deprecated since labels are just values from x
#' @param legendTitle legend title name
#' @param nrowL legend row number
#' @param theme_bg background theme
#' @param cols either a named vector (names from x) to accurately specify the color or a color vecto which is to be input in the alphabetical order as the x values 
#' @param legend_within whether to put legend within the figure
#' @param lx relative x position of legend (between 0~1)
#' @param ly relative y position of legend (between 0~1)
#' @param verbose whether to print warning due to survframe count incompatible with observed.
#' @param addEventCount addEventCount whether to add event count in the legend
#' @export 
ggKMcat <- function(x, OS, OScensoring, main, cex.main=0.8, tag,
	xlab='Survival time', ylab='Survival probability',
	labels, legendTitle, nrowL=NA, sizeaxis=12, sizeaxistitlerel=1.2, sizetitlerel=0.8, 
	theme_bg=theme_base(),
	cols, legend_within=TRUE, lx=0.8, ly=0.88, plot=TRUE, verbose=FALSE, addEventCount=F, countFromInput=TRUE){
	#browser()
	#if(missing(labels)){
	#		labels <- sort(levels(as.factor(x)), decreasing=FALSE)
	#}	
	if(missing(legendTitle)){
			legendTitle <- deparse(substitute(x))
	}
	if(missing(tag)){
			tag <- sprintf('KM Curve for %s', legendTitle) ## legendTitle might be specified and more informative
	}
	xGrouped <- x
	nCurve <- n_unique(xGrouped[!is.na(xGrouped)])
	if(missing(cols)){
		cols <- brewer.pal(8, "Set1")[1:nCurve]
		## override the color for aethestics when 2 curves
		if(nCurve==2) cols <- c("red1","blue3")
	}
	#browser()
	# p value without dichotomizing
	#contP <- coxphP(expr=x, OS=OS, OScensoring=OScensoring)
	# p value after categorizing
	#groupedP <- coxphP(expr=xGrouped, OS=OS, OScensoring=OScensoring)
	coxfitGrouped <- coxphP2(expr=xGrouped, OS=OS, OScensoring=OScensoring)
	paredTextGrouped <- parse_coxphP2(coxfitGrouped)
	#browser()
	groupedP <- coxfitGrouped$P
	#coefExp <- coxfit$coef_exp
	#browser()
	#coefExp_text <- coefExp2text(coefExp)
	#browser()
    # added noNA on 2014/03/14 so that when NA is present, the count N is more accurate
    # # the problem is not because of NA, it is because of createSurvivalFrame() which adds two rows starting at time=0 for each group
    # therefore, added createSurvivalFrame0 so that N is from sframe built before adding the rows for plotting
    tdat <- data.frame(OS=OS, OScensoring=OScensoring, xGrouped=xGrouped)
    tdat <- noNA(tdat)
	sfit <- with(tdat, survfit(Surv(OS, OScensoring)~xGrouped))
	#browser()
  temp <- createSurvivalFrame0(sfit)
	sframe_raw <- temp$f.frame_raw
	sframe <- temp$f.frame
    #sframe <- sframe0 <- createSurvivalFrame(sfit)
	## 
    # here is the magic: replace the automatically generated group labels with something we want (i.e. through user specification)
    # since xGrouped= is added everywhere, just replace it with ''
	##
    ## when xGrouped has one unique value, strata is null; add this to avoid error
    if(!is.null(sframe$strata)){
	sframe$strata <- str_replace(sframe$strata, sprintf("%s=", as.character(quote(xGrouped))), '') # back to the origianl xGrouped values
	sframe_raw$strata <- str_replace(sframe_raw$strata, sprintf("%s=", as.character(quote(xGrouped))), '') # back to the origianl xGrouped values
	}
	## when labels is specified, it is a named vector
	#browser()
	if(!missing(labels) & !is.null(sframe$strata)) {
		#browser()
		sframe$strata <- mapNames(sframe$strata, labels)
		sframe_raw$strata <- mapNames(sframe_raw$strata, labels)
	}
    # on 2014/03/06: futher add N=xx info
  ttab <- table(sframe_raw$strata) # this is already trying to be the correct number; but in practive difference of 1 or 2 is observed!!!
	# added on 2018/02/22: force the count to be from the input data
  if(countFromInput){
    ttab <- table(tdat$xGrouped)[names(ttab)] # this may be more accurate? 
  }
  # now force it to be the same as table(x)!
	# NA omitted in x, the strata var, however, OS can still be missing
  # so this is nominal N, effective N due to OS missing is not computed
  tab_x <- table(noNA(x)) 
  tab_x_censoring <- table(x, OScensoring) 
  tab_event <- tab_x_censoring[, '1']
	# now updating!
  # ttab1 not used on 02/22/2018: to make nominal N always shown in figure legend!
  # <--- can be omitted
  ttab1 <- ttab 
	if(all(names(ttab1) %in% names(tab_x))){
	for(e in names(ttab1)){
		if(verbose){
			if(ttab1[e] != tab_x[e]) 
				warning(sprintf('For level %s, count from survframe %d != %d observed in data', e, ttab1[e], tab_x[e]))
		}
		ttab1[e] <- tab_x[e]
	}
	}
  # can be omitted----->
  #browser()
  tab_event <- tab_event[names(ttab)]
  if(addEventCount) {
    strataAugMap <- data.frame(names(ttab), str_c(names(ttab), ' (N=', ttab, ', events=', tab_event, ')', sep='')) 
  } else {
    strataAugMap <- data.frame(names(ttab), str_c(names(ttab), ' (N=', ttab, ')', sep=''))  
  }
  #browser()
  sframe$strata <- mapNames(sframe$strata, strataAugMap)
  # trying to preserve the order of levels: first passed to ttab, than to strataAugMap
  # this will affect the legend ordering
  sframe <- mutate(sframe, strata=orderFactor(strata, new.order=strataAugMap[, 2])) 
	#browser()
  if(missing(main)) {
		gtitle <- sprintf('%s\nLogrank test (overall) P: %.3g\n%s\n%s', 
			tag, groupedP, paredTextGrouped$coefExp_text, paredTextGrouped$CI_text)
	} else {
		gtitle <- main
	}
    # more data adaptive nrowL
    if(is.na(nrowL)){
        nrowL <- guess_nrowL(length(ttab), useMax=legend_within)
    }
	#browser()
	p0 <- qplot_survival(sframe, f.CI=F)+ ggtitle(gtitle)+ylim(c(0,1))+theme_bg+
		xlab(xlab)+ylab(ylab)+theme(plot.title = element_text(size = rel(cex.main), colour = "black"))
	#if(theme_bw) p0 <- p0+theme_bw()	
    #browser()
	scale_axis_size <- theme(axis.text=element_text(size=sizeaxis, colour = "black", face ='bold', hjust=0.5, vjust=0.5), 
    axis.title.x = element_text(size = rel(sizeaxistitlerel)),
    axis.title.y = element_text(size = rel(sizeaxistitlerel)),
    plot.title = element_text(size = rel(sizetitlerel)),
    legend.title = element_text(size = rel(sizetitlerel)),
    legend.text = element_text(size = rel(sizetitlerel)))
  p1 <- p0+guides(color = guide_legend(nrow = nrowL))+theme(legend.position="bottom")+
		scale_colour_manual(name=legendTitle,values=cols)+scale_axis_size
	if(legend_within){
        # nrowL is overriden now
        p1 <- p1+theme(legend.position=c(lx,ly))+guides(color = guide_legend(nrow = nrowL))
    }
    if(plot){
		#coxphfit <- try(coxph(Surv(OS, OScensoring)~xGrouped), silent=TRUE)
		cat('------------------------\n')
		cat('Cox PH fit on categorical covariate:\n')
		print(coxfitGrouped$fit)
		#browser()
		return(p1)
	} else {
    return(list(groupedP=groupedP, p=p1, coxfitGrouped=coxfitGrouped, paredTextGrouped=paredTextGrouped, N1=unname(ttab[levels(x)[1]]), N2=unname(ttab[levels(x)[2]])))
	}
}
#with(dat_SETD2_moreClinics, plotKM(x=reorder(as.factor(stat_CN_MUT), new.order=lv_stat_CN_MUT), OS=OS, OScensoring=OScensoring, addEventCount=T, lx=0.03, ly=0.23, xlab=KM_xlab, legendTitle='', tag='stat_CN_MUT'))

#with(clin1u_mut, ggKMcat(x=myx, OS=PFS, OScensoring=PFScensoring, xlab='Survival in months', tag=tag))+theme_title


		
#with(tdf, ggKM(x=EMTscore, OS=OS, OScensoring=OScensoring, breaks=c(-Inf, 0, Inf), 
#		labels=c('E', 'M'), tag=sprintf("%s", myhist), legendTitle='EMT Status'))
#with(aml, ggKMcat(x=x2, OS=time, OScensoring=status, cols=c("cyan", "purple")))
#with(aml, ggKMcat(x=x2, OS=time, OScensoring=status, cols=c("a_Nonmaintained"="cyan", "b_Maintained"="purple")))

#plotKM(x=cut(myy, breaks=c(-Inf, median(myy, na.rm=TRUE), Inf)), OS=ddf$OS, OScensoring=ddf$OScensoring, labels=c('low', 'high'), tag=sprintf("%s, %s", myyName, datName), legendTitle=sprintf('%s Low/High', myyName))

#ggKMcat(x=x, OS=OS, OScensoring=OScensoring, tag=tag, legendTitle='T Stage')

#' flexible KM curve wrapper
#'
#' plot KM curve for a cont/cat x. 
#'
#' the program detects x as categorical if it is a factor or it is a character variable. 
#' @param x the stratification variable, either numeric or factor or character
#' @param ... other parameters passed to either ggKMcat or ggKM 
#' @export
plotKM <- function(x, legendTitle, tag, ...){
	#browser()
	if(missing(x)){
        p <- plotKM0(...)
        return(p)
    }
    if(missing(legendTitle)){ # otherwise, it always is x when not specified
			legendTitle <- substitute(x)
	}	
	if(missing(tag)){
			tag <- sprintf('KM Curve for %s', deparse(substitute(x)))
	}	
	if(is.character(x) | is.factor(x)){
   # browser()
		return(ggKMcat(x=x, legendTitle=legendTitle, tag=tag, ...))
	} else {
		return(ggKM(x=x, legendTitle=legendTitle, tag=tag, ...))
	}
}

#with(df_info, plotKM(x=NmutationInFM, OS=OS, OScensoring=OScensoring)) # cont-->OS
