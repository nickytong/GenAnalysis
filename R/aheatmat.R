if(FALSE){
# cd /home/ptong1/Backup/Package/; R --vanilla
library(roxygen2)
roxygenize("GenAnalysis")

library(devtools)
build('GenAnalysis')
build('GenAnalysis')
install('GenAnalysis')

load_all('GenAnalysis')
load_all('/home/ptong1/Backup/Package/GenAnalysis')

build_win('GenAnalysis')

##
detach("package:GenAnalysis", unload=TRUE)
library(GenAnalysis)

pdf(file=file.path(dirSource, 'mypal.pdf') )
showmypal()
dev.off()
    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/GenAnalysis/R/aheatmap.R'))



pdf(file=file.path(dirSource, 'mypal2.pdf'), height=7*length(get_pal())/15)
showmypal2()
dev.off()


    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/GenAnalysis/R/aheatmap.R'))

## get distinct colors using Kevin's Polychrome
library(Polychrome)
lastpal <- createPalette(15, colorpalette2colvec('set1'), M=10000)
lastpal <- createPalette(33, colorpalette2colvec('set1'), M=10000)
swatch(lastpal)


}
	#### aheatmap using aheatmat with better design
##
## personal color pallete
##
colpal_cat <- c('set1', 'dark2', 'accent', 'set2', 'set3')
colpal_cont_2sided <- c('bl2gr', 'cy2yl', 'gr2rd', 'ma2gr', 'bluered', 'blueorange', 'jet', 'bl2yl')
colpal_cont_1sided <- c('heat.colors', 'heat', 'redscale', 'greenscale')
colpal_cont <- c(colpal_cont_2sided, colpal_cont_1sided)

#' convert a character vector to a color vector by specified col. 
#' @param vec a vector of values. it can be numeric or categorical. 
#' @param colpal a color vector. If NA is present, it needs to to be taken into account. 
#' @return a color vecotor same length as vec
#' @export
vec2colvec <- function(vec, colpal=NULL){
	isnum <- is.numeric(vec)
	#browser()
	if(!isnum){
		ncat <- n_unique(vec, ignoreNA=FALSE)
		vec <- as.character(vec) # force to character, even logical, otherwrise, error come up
	}
	if(is.null(colpal)){
		if(isnum) {
			colpal <- colorpalette2colvec('bluered')
		} else {
			if(ncat<=9)	colpal <- colorpalette2colvec('set1')[1:ncat]
			else colpal <- colorpalette2colvec('bluered', nrcol=ncat)
		}
	}	
	#browser()
	pheno <- data.frame(x=vec)
	colpal <- list(x=colpal)
	colbar <- prepcolbar(pheno, colpal=colpal)[, 1]
	colbar
}


#' extracts color given a specified pallete (string) possibly with alpha blending and number of colors desired
#'
#' this function is originally from LSD but adapted a little (adding more colors). This function will return the color if nrcol
#' is not specified; otherwise, color interperlation through grDevices::colorRampPalette will be used. 
#' 
#' @param pal string specifying pallete (predefined)
#' @param nrcol number of colors desired
#' @param alpha alpha blending
#' @return a vector of colors
#' @import Polychrome
#' @export
colorpalette2colvec <- function (pal='brewer.comb', nrcol = NULL, alpha = NULL, seed=NULL) {
	# allow 'set1_dark2'
	isComposite <- str_detect(pal, '_')
	if(length(isComposite)==0) isComposite <- FALSE # pal is NULL
	#browser()
	if(!isComposite){
		palette <- colorpalette2colvec0(pal=pal, seed=seed, nrcol = nrcol, alpha = alpha)
	} else {
		#browser()
		pals <- unlist(str_split(pal, '_'))
		palette <- foreach(x=pals, .combine='c') %do% {
			# notice we disable nrcol and alpha at this stage; 
			colorpalette2colvec0(pal=x, seed=seed, nrcol = NULL, alpha = NULL)
		}	
	}
	if (!is.null(nrcol)) {
    	if(length(palette)>=nrcol){
    		palette <- palette[1:nrcol]
    	} else {
	    	# interpolation
	        colfct = grDevices::colorRampPalette(palette)
	        palette = colfct(nrcol)
		}
    }
	palette
}

palette_L <- list(standard = c("white", "light yellow", 
            "yellow", "gold", "brown"), heat = c("grey", "dark blue", 
            "red", "orange", "gold"), red = c("#E0C6C6", "red"), 
            crazyred = c("#940000", "#A50000", "#FF5C5C", "#FFB9B9"), 
            green = c("#C6E0C6", "green"), crazygreen = c("dark green", 
                "#009700", "green", "#C0F5D0"), blue = c("#C6C6E0", 
                "blue"), crazyblue = c("dark blue", "blue", "#7390EE", 
                "light blue"), black = c("#E5E5E5", "black"), 
            mountain = c("light green", "dark green", "black", 
                "dark grey", "#F0F0F0"), girly = c("violet", 
                "violetred", "violetred1", "purple", "purple3"), 
            jamaica = c("red", "yellow", "green"), wird = c("white", 
                "darkred"), wiorrd = c("white", "darkorange", 
                "darkred"), wiocrd = c("white", "darkorchid4", 
                "darkred"), wiorocgn = c("white", "orange", "darkorange", 
                "darkorchid1", "darkorchid4", "darkorchid1", 
                "green", "darkgreen", "darkgreen"), wicymgyl = c("white", 
                "cyan", "darkcyan", "magenta", "darkmagenta", 
                "magenta", "yellow", "yellow2", "lawngreen"), 
            wibugr = c("white", "skyblue1", "skyblue4", "slateblue1", 
                "slateblue4", "slateblue1", "slategray1", "slategray2", 
                "slategray4"), boxes = c("darkgreen", "cornflowerblue", 
                "yellow", "darkred", "darkorchid"), pies = c("red", 
                "yellow", "green", "purple", "blue"), colorblind = c("#000000", 
                "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                "#D55E00", "#CC79A7"), ylorrd = RColorBrewer::brewer.pal(9, 
                "YlOrRd"), ylorbr = RColorBrewer::brewer.pal(9, "YlOrBr"), 
            ylgnbu = RColorBrewer::brewer.pal(9, "YlGnBu"), ylgn = RColorBrewer::brewer.pal(9, 
                "YlGn"), reds = RColorBrewer::brewer.pal(9, "Reds"), rdpu = RColorBrewer::brewer.pal(9, 
                "RdPu"), purples = RColorBrewer::brewer.pal(9, "Purples"), 
            purd = RColorBrewer::brewer.pal(9, "PuRd"), pubugn = RColorBrewer::brewer.pal(9, 
                "PuBuGn"), pubu = RColorBrewer::brewer.pal(9, "PuBu"), orrd = RColorBrewer::brewer.pal(9, 
                "OrRd"), oranges = RColorBrewer::brewer.pal(9, "Oranges"), 
            greys = RColorBrewer::brewer.pal(9, "Greys"), greens = RColorBrewer::brewer.pal(9, 
                "Greens"), gnbu = RColorBrewer::brewer.pal(9, "GnBu"), bupu = RColorBrewer::brewer.pal(9, 
                "BuPu"), bugn = RColorBrewer::brewer.pal(9, "BuGn"), blues = RColorBrewer::brewer.pal(9, 
                "Blues"), spectral = RColorBrewer::brewer.pal(11, "Spectral"), 
            rdylgn = RColorBrewer::brewer.pal(11, "RdYlGn"), rdylbu = RColorBrewer::brewer.pal(11, 
                "RdYlBu"), rdgy = RColorBrewer::brewer.pal(11, "RdGy"), rdbu = RColorBrewer::brewer.pal(11, 
                "RdBu"), puor = RColorBrewer::brewer.pal(11, "PuOr"), prgn = RColorBrewer::brewer.pal(11, 
                "PRGn"), piyg = RColorBrewer::brewer.pal(11, "PiYG"), brbg = RColorBrewer::brewer.pal(11, 
                "BrBG"), set1 = RColorBrewer::brewer.pal(9, "Set1")[-6], set2 = RColorBrewer::brewer.pal(8, 
                "Set2"), set3 = RColorBrewer::brewer.pal(12, "Set3")[-2], pastel1 = RColorBrewer::brewer.pal(9, 
                "Pastel1"), pastel2 = RColorBrewer::brewer.pal(8, "Pastel2"), 
            paired = RColorBrewer::brewer.pal(12, "Paired"), dark2 = RColorBrewer::brewer.pal(8, 
                "Dark2"), accent = RColorBrewer::brewer.pal(8, "Accent"), standardterrain = grDevices::terrain.colors(9), 
            standardtopo = grDevices::topo.colors(9), standardheat = grDevices::heat.colors(9), 
            standardrainbow = grDevices::rainbow(9, start = 0.7, end = 0.1), 
            standardcm = grDevices::cm.colors(9), bl2gr = colorRamps::blue2green(11), 
            bl2gr2rd = colorRamps::blue2green2red(11), bl2rd = colorRamps::blue2red(11), 
            bl2yl = colorRamps::blue2yellow(11), cy2yl = colorRamps::cyan2yellow(11), 
            gr2rd = colorRamps::green2red(11), ma2gr = colorRamps::magenta2green(11), 
            matlablike = colorRamps::matlab.like(11), matlablike2 = colorRamps::matlab.like2(11), 
            primarycolors = colorRamps::primary.colors(11), ygob = colorRamps::ygobb(11), # starts with added color by myself
			blueyellow=oompaBase::blueyellow(11),
			greenred=oompaBase::redgreen(11),
			cyblyl=c("cyan3", 'black', "yellow"), # cyan to black to yellow
			jetcolors=oompaBase::jetColors(11),
			jet=squash::jet(11),
			redscale=oompaBase::redscale(11),
			greenscale=oompaBase::greenscale(11),
			bluescale=oompaBase::bluescale(11),
			bluered=gplots::bluered(11),
			heat.colors = rev(heat.colors(9)),
			blueorange=squash::blueorange(11),
			darkbluered=squash::darkbluered(11),
			heat=squash::heat(11),
			NMF=NMF:::ccRamp('-RdYlBu'), # NMF col
			mutblack=c('gray', 'black'), # mutation color 0/1
			mutred=c('#CCCCFF', '#FF0000'), # mutation color 0/1
			mutblue=c('#CCCCFF', '#0000FF'), # mutation color 0/1
			brewer.comb=c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Dark2"),
				RColorBrewer::brewer.pal(8, "Accent"), RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(9, "Set3"))
			)
#' obtain color pallete defined in this package
#' @param pal the name for color pallete; if NULL, all will be returned in a list
#' @export
get_pal <- function(pal=NULL){
	if(is.null(pal))
		return(palette_L)
	if(!pal %in% names(palette_L))
		stop(sprintf('%s is not in exisiting pallete: \n%s', pal, str_c(names(palette_L), sep='|\n')))
	return(palette_L[[pal]]	)
}	
colorpalette2colvec0 <- function (pal='brewer.comb', nrcol = NULL, alpha = NULL, seed=NULL) {
    #browser()
    if(!is.null(seed)){
    	# using polychrome
    	palette <- Polychrome::createPalette(nrcol, seed, M=10000)
    } else {
    	if (length(pal) > 1) {
    		# palette is the color vector; if pal length >1, means it is a color vector; no processing
    	    palette = pal
    	} else {
    		# in case pal not specified
    	    palette = get_pal(pal)
    	}
    }	
    if (is.null(palette)) {
        stop(paste(pal, "is not a valid palette name for colorpalette"))
    }
	#
	# if nrcol is NULL, the specified color will be returned without interpolating; 
	# if nrcol is specified, the color will always be interpolated. In most cases, the color would be changed and only overlap if nrcol is special (through command divider)
	#
    if (!is.null(nrcol)) {
    	if(length(palette)>=nrcol){
    		palette <- palette[1:nrcol]
    	} else {
	    	# interpolation
	        colfct = grDevices::colorRampPalette(palette)
	        palette = colfct(nrcol)
		}
    }
    if (!is.null(alpha)) {
        palette = LSD::convertcolor(palette, alpha)
    }
    return(palette)
}
#show_col(colorpalette2colvec('set1'))
#show_col(colorpalette2colvec('brewer.comb'))
# show_col(colorpalette2colvec(c('#FF0505', '#920FDA', '#FFE923', '#20E469'), nrcol=15))

#' check if a vector is categorical or continuous
#'
#' @param vec a vector
#' @return logical
#' @export
iscat <- function(vec){
	if(is.factor(vec) | is.character(vec)){
		res <- TRUE
	} else {
		res <- FALSE
	}
	res
}
#' check columnwise of a data frame whether it is categorical (character or factor) or continuous
#'
#' this is a colwise application of iscat()
#'
#' @param dat a data frame
#' @return logic vector
#' @export
iscatDF <- function(dat){
	#colwise(iscat)(dat)
	unlist(colwise(iscat)(dat))
	#browser()
}
# iscatDF(annCol4NMF)
# extract color for a given vec and pal
ispaloverflow <- function(pal){
	if(is.null(pal)){
		# NULL, not overflow as in usedCat=TRUE: economic usage of color
		return(FALSE)
	} else if (is.na(pal[1])){
		return(TRUE) # overflow
	} else {
		return(TRUE)
	}
}
# auto-generating color vector given specified color pallete name
autocolpal <- function(vec, pal=NULL, alpha = NULL, usedCat=NULL){
	#
	# warning: this will alway generate lots of colors for cont var! (n_unique is used); may need improvement, e.g. separating cont and cat
	#
	nu <- n_unique(vec)
	#browser()
	if(nu<=8){
		if(ispaloverflow(pal)) pal <- 'set1' # in case colpallete is not enough
		# special treatment for categorical vars
		res <- colorpalette2colvec(pal, alpha=alpha)[1:nu]
	} else {
		if(ispaloverflow(pal)) pal <- 'greenred' # in case colpallete is not enough
		res <- colorpalette2colvec(pal, nrcol=nu, alpha=alpha)
	}
	# overwride when ecomonic usage of brewer color is activated
	if(!is.null(usedCat)) {
		# for economic usa of brewer color: pick colors starting from last leftover in brewer.comb
		res <- colorpalette2colvec('brewer.comb', alpha=alpha)[(usedCat+1):(usedCat+nu)]
	} 
	res
}
#' get default colorpal (color vector) for a data frame (colwise)
#' 
#' this is auto color pallete; for better vis, the user need to specify it explicitely by modifying the returned list or start from scratch
#' 
#' @param dat a data frame
#' @param pal optional user specified palette to update the final colpal
#' @return a (named) list of color vectors that can be used to build colormap and colormat
#' @export
autocolpal4DF <- function(dat, alpha=NULL, economicBrewer=TRUE, pal=NULL){
	if(is.null(dat)){ # when no pheno, for convenience purpose
		return(NULL)
	}
	#browser()
	dat <- data.frame(dat, check.names=F)
	catstatus <- iscatDF(dat)
	index_cont <- 1
	index_cat <- 1
	usedCat <- 0 # economic usage of brewer col
	if(economicBrewer){
		res <- foreach(i=1:ncol(dat)) %do% {
			#browser()
			#if(i==13) browser()
			if(catstatus[i]==TRUE){
				tt <- autocolpal(vec=dat[, i], alpha=alpha, usedCat=usedCat)
				index_cat <- index_cat+1
				usedCat <- usedCat+n_unique(dat[, i])
			} else {
				tt <- autocolpal(vec=dat[, i], pal=colpal_cont[index_cont], alpha=alpha)
				index_cont <- index_cont+1	
			}
			tt
		} 
	} else {
		res <- foreach(i=1:ncol(dat)) %do% {
			#browser()
			if(catstatus[i]==TRUE){
				tt <- autocolpal(vec=dat[, i], pal=colpal_cat[index_cat], alpha=alpha)
				index_cat <- index_cat+1
			} else {
				tt <- autocolpal(vec=dat[, i], pal=colpal_cont[index_cont], alpha=alpha)
				index_cont <- index_cont+1	
			}
			tt
		}
	}	
	#browser()
	names(res) <- colnames(dat)
	# add meta info for col type
	attr(res, 'iscat') <- catstatus
	attr(res, 'pal_len') <- sapply(res, length)
	attr(res, 'class') <- c('colpal4DF', 'list')
	if(!is.null(pal)){
		#browser()
		res <- updatePal(oldpal=res, newpal=pal)
	}
	res
}
#colpal_rel1 <- autocolpal4DF(dat=pheno)

#autocolpal4DF(dat=annCol4NMF)
#' plot method for colpal4DF
#'
#' @method plot colpal4DF
#' @param obj a colpal4DF class usually returned by autocolpal4DF()
#' @param nr number of colors to show
#' @param mar mar in par() setting
#' @export
plot.colpal4DF <- function(obj, nr=12, mar=c(2, 9, 2, 1)){
	npal = length(obj)
	op <- par(mar=mar)
    plot(1, 1, xlim = c(0, nr), ylim = c(0, npal), type = "n", 
        axes = FALSE, bty = "n", xlab = "", ylab = "", main = "view of specified colpal4DF")
    for (i in 1:npal) {
    	nr_effective <- length(obj[[i]])
    	if(nr_effective<=nr){
    		nr_pl <- nr_effective
    		pal_pl <- obj[[i]]
    	} else {
    		nr_pl <- nr
    		pal_pl <- grDevices::colorRampPalette(obj[[i]])(nr) # down-interpolation for continuous pal when it is too many
    	}
        rect(xleft = 0:(nr_pl - 1), ybottom = i - 1, xright = 1:nr_pl, 
            ytop = i - 0.2, col = pal_pl, 
            border = "light grey")
    }
    #browser()
    text(rep(-0.1, npal), (1:npal) - 0.6, labels = names(obj), xpd = TRUE, adj = 1)
    par(op)
}
#plot.colpal4DF(colpal_rel1)
#plot(autocolpal4DF(dat=annCol4NMF))

#' get method for class colpal4DF
#' 
#' this function is named getter, not get to minimize conflict with the usual get function
#'
#' @param obj an colpal4DF object
#' @param what specify what to get from augMat, i.e. base, var.equal
#' @method getter colpal4DF
#' @return an object as requested
#' @export
getter.colpal4DF <- function(obj, what='avail'){
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
#getter(autocolpal4DF(dat=annCol4NMF))

# a vetor to a named list
vec2namedList <- function(vec, name){
	res <- list(vec)
	names(res) <- name
	#browser()
	res
}
#' build a map given a vector of values and color pallete
#'
#' for continuous variable, makecmap is used; for categorical variable, buildcolorfacs is used
#'
#' @param vec a vector of values to be mapped
#' @param pal color pallete
#' @name name of the vec as in the named list of colpal
#' @param iscat logical indicating if cont or cat
#' @param n only used to specify number of colors generated for continious vec
#' @return a color map
#' @export
vecpal2map <- function(vec, pal, name, iscat, n=59){
	if(iscat){
		vec_l <- vec2namedList(vec=vec, name=name)
		pal_l <- vec2namedList(vec=pal, name=name)
		#browser()
		map <- buildcolorfacs(barVecs=vec_l, barCols=pal_l)	
	} else {
		map <- makecmap(vec, colFn=grDevices::colorRampPalette(pal), n=n) 

	}
	#browser()
	map
}
# check if map is cont or categorical
ismapcat <- function(map){
	#browser()
	ifelse('include.lowest' %in% names(map), FALSE, TRUE)
}
#ismapcat(z)
#' build colmat/colmap given pheno (colwise of data frame) and colpal
#'
#' building colmap depends on iscat status of each column, catstat. catstat is first searched with attributes; if it is not specified, iscatDF is used to calculate it. 
#' therefore, the default assumes the user has prepared pheno in a way that categorical vars are either character or factor and cont vars are numeric. 
#'
#' @param pheno a data frame where columns represent different column bars (in a heatmap)
#' @param colpal color pallete, a named list of color vectors. When not specified, autocolpal4DF will be used to generate the default.
#' @return a list with elements mapL (a list of map for color legend) and colmat (color matrix for colbar). 
#'  input pheno and colpal are attached as attributes
#' @export
prepcolbar <- function(pheno, colpal=NULL){
	# force data frame
	if(is.null(pheno)) {
		return (NULL)
	}
	pheno <- data.frame(pheno, check.names=F)
	#browser()
	if(is.null(colpal)){
		colpal <- autocolpal4DF(dat=pheno)
	} else {
		if(!inherits(colpal, 'colpal4DF')){
			# in this case, colpal is just color palette name, not color vector
			colpal <- autocolpal4DF(dat=pheno, pal=colpal)
		}
	}
	#browser()
	zz <- colpal[match(colnames(pheno), names(colpal))]
	attributes(zz) <- attributes(colpal)
	colpal <- zz # so that all attributes are still preserved
	if(!all(colnames(pheno) %in% names(colpal))) stop(sprintf('Some colnames(pheno) are not present in names(colpal)!\n'))
	catstat <- getter(colpal, 'iscat')
	pal_len <- getter(colpal, 'pal_len')
	if(is.null(catstat)) catstat <- iscatDF(pheno)
	mapL <- vector('list')
	colmat <- NULL
	#browser()
	for(i in 1:ncol(pheno)) {
		coln <- colnames(pheno)[i]
		vec <- pheno[, i]
		#browser()
		# get col map
		map <- vecpal2map(vec=vec, pal=colpal[[coln]], name=coln, iscat=catstat[i], n=pal_len[i])
		#mapL[[i]] <- map 
		if(catstat[i]){
			mapL[[i]] <- map[[1]] # for cat, to preserve the structure: no added list nesting!
		} else {
			mapL[[i]] <- map
		}		
		#browser()
		# get col vector
		if(catstat[i]){
			# is cat
			cl <- buildcolorMat(map)
		} else {
			cl <- cmap(vec, map=map)
		}	
		colmat <- cbind(colmat, cl)
	}	
	names(mapL) <- colnames(colmat) <- names(colpal)
	rownames(colmat) <- rownames(pheno)
	#res <- list(mapL=mapL, colmat=colmat)
	res <- colmat
	attr(res, 'mapL') <- mapL
	attr(res, 'colpal') <- colpal
	attr(res, 'pheno') <- pheno
	attr(res, 'class') <- c('prepcolbar', 'list')
	res

}
#colbar <- prepcolbar(pheno, colpal=colpal)

#prepcolbar(annCol4NMF)

#' plot method for prepcolbar class
#'
#' @method plot prepcolbar
#' @param obj a prepcolbar class (a data frame)
#' @param cex.col cex.col
#' @param cex.row cex.row
#' @param mar margin setting
#' @export
plot.prepcolbar <- function(obj, cex.row=NULL, cex.col=NULL, mar=c(2, 9, 2, 1)){
	#browser()
	colmat <- obj
	cex.row <- guessCEX(ncol(colmat))
	cex.col <- guessCEX(nrow(colmat))
	op <- par(mar=mar)
    cimage(colmat, border = NA, add =F, xlab='', ylab='', axes=FALSE)
	axis(1, at = 1:nrow(colmat), labels = rownames(colmat), tick = FALSE, line = -0.5, cex.axis = cex.col, las=2)
	axis(2, at = 1:ncol(colmat), labels = colnames(colmat), tick = FALSE, line = -0.5, cex.axis = cex.row, las=2)
	par(op)
}
#plot.prepcolbar(prepcolbar(annCol4NMF))


#' subset main expression mat
#'
#' this can be used to subset genes, extract effective columns in augMat, and further subset columns. 
#'
#' Users can subset samples through two selectors: sSel0 and sSel1. However, sSel0 is primarily preserved to select effective samples after
#' data augmentation. sSel1 is for subsetting. colOrderIndex can be used to sort the samples (columns). The final samples presented
#' is the intersection between sSel0, sSel1 and colOrderIndex. Notice that when index is provided, 
#' it is assumed to be global index for easy tracking (both gSel and sSel0, sSel1 and colOrderIndex). Provoding sample name or gene name
#' in the form of rownames or colnames also works. 
#'
#' @param mat matrix for heatmap
#' @param gSel gene selection (global index or gene names). This can be index or gene names. Internally, both will be converted to integer based index (global or local)
#' @param sSel0 sample selection filter 0; by default, this deals with augMat to remove all NA columns
#'  the user can override this by supplying a vector of indeces or sample names.
#' @param cluster_rows whether to draw row clusters
#' @param cluster_cols whether to draw column clusters
#' @param sSel1 sample selection (index or sample names). This is useful to select a subset of samples based
#' on phenotype, i.e. some mutation, like ++ to be removed
#' @param colOrderIndex index/sample names to order column (global index or sample names); The actual number of samples surviving is an intersect between colOrderIndex and sSel (from sSel0 and sSel1).
#'  if not specified, column clustering will be instructed to construct; otherwise, no clumn clustering will be done later on
#' @param rowOrderIndex index (integer or gene names) to order rows; The actual genes surviving is an intersect between rowOrderIndex and gSel.
#'  if not specified, row clustering will be instructed to construct; otherwise, no row clustering will be done later on
#' @return expression matrix with attributes gSel, sSel and so on
#' @export
subsetMainMat <- function(mat, gSel=NULL, sSel0=NULL, cluster_rows, cluster_cols, sSel1=NULL, colOrderIndex=NULL, rowOrderIndex=NULL) {
	## 
	## unifying index: global integer index
	##
	# deal with augMat
	#browser()
	if(inherits(mat, 'augMat')){
		sSel0 <- getter(mat, 'filledColInd') # global integer indexing
	} else {
	    sSel0 <- 1:ncol(mat)
	}
	# override by user specification
	if(!is.null(sSel0)){
	  sSel0 <- sSel0
	  sSel0 <- ixunify2index(ix=sSel0, name_global=colnames(mat)) # global integer index
	}
	if(!is.null(sSel1)){
	  sSel1 <- ixunify2index(ix=sSel1, name_global=colnames(mat)) # global integer index
	} else {
	  sSel1 <- 1:ncol(mat)
	}
	sSel <- intersect(sSel0, sSel1)
	#browser()
	# column ordering
	if(is.null(colOrderIndex)){
		# no column ordering
		columnOrder <- sSel
		colcl <- TRUE
	} else {
		# color ordered and subsetted
		columnOrder <- sSel[noNA(match(colOrderIndex, sSel))]
		colcl <- FALSE # no column clustering
	}
	# gSel
	if(is.null(gSel)) gSel <- rownames(mat) # default uses all genes
	# global integer index for gene
	gSel <- ixunify2index(ix=gSel, name_global=rownames(mat)) 
	# row ordering
	if(is.null(rowOrderIndex)){
		# row ordering not specified: draw row clustering
		# no row ordering by default since cluster_rows=F
		rowOrder <- gSel
		#rowcl <- TRUE
	} else {
		# color ordered and subsetted
		rowOrderIndex <- ixunify2index(ix=rowOrderIndex, name_global=rownames(mat)) # ix unifying
		rowOrder <- gSel[noNA(match(rowOrderIndex, gSel))]
		#rowcl <- FALSE # no row clustering
	}
	# mat subsetted
	# all index output are integer based
	smat <- mat[rowOrder, columnOrder]
	attr(smat, 'gSel_g') <- rowOrder # index, based on global reference
	attr(smat, 'sSel_g') <- columnOrder # index, based on global reference
	attr(smat, 'gSel_l') <- ixglobal2local(rowOrder, rowOrder) # index, based on local reference
	attr(smat, 'sSel_l') <- ixglobal2local(columnOrder, columnOrder) # index, based on local reference
	attr(smat, 'rowcl') <- cluster_rows
	attr(smat, 'colcl') <- cluster_cols
	#browser()
	# lets see if it is really necessary
	#attr(smat, 'rawmat') <- mat # this is needed to reconstruct matrix or do plotting with pheno; this should be more convenient outside aheatmat
	smat
}  


#z = subsetMainMat(mat=rna, gSel=gsel_ttestRes_NF1vsKRAS, sSel0=NULL, sSel1=which(!is.na(pData$NF1vsKRAS)))
#
#' transform a matrix (center/scaling; truncation; clustering;) for heatmap; strongly recommends first apply subsetMainMat
#' to get input mat since some attributes are used
#'
#' @param mat a matrix, usually from subsetMainMat()
#' @param scale selection from c("none", "row", "column")
#' @param clusterWithScaledData logical indicating if use scaled data for clustering; default is FALSE
#' @param cluster_rows whether to cluster rows
#' @param cluster_cols whether to cluster columns 
#' @param clustering_distance_rows distance metric for rows
#' @param clustering_distance_cols distance metric for cols
#' @param clustering_method clustering method
#' @param truncate logical indicating if truncation is needed; default is NA will enable truncate if scale!='none'
#' @param Lower parameter Lower to truncByQuantile()
#' @param Upper parameter Upper to truncByQuantile()
#' @param q1 parameter q1 to truncByLimit
#' @param q2 parameter q2 to truncByLimit
#' @return a matrix with elements gSel, sSel, tree_row, tree_col and so on
#' @export 
transform_smat <- function(mat, scale = "row", clusterWithScaledData=FALSE, cluster_rows = TRUE, 
	cluster_cols = TRUE, 
	clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
	clustering_method = "ward.D",
	truncate=NA, q1=0, q2=1, Lower=NULL, Upper=NULL) {
	# tree_row = cluster_mat(matRaw, distance = clustering_distance_rows, method = clustering_method)
	# tree_col = cluster_mat(t(matRaw), distance = clustering_distance_cols, method = clustering_method)
	#matRaw = as.matrix(Getter(mat, 'rawmat')) # matRaw is the input matrix which is desired to do clustering and stuff 
	#attrL <- attributes(mat)
	matRaw = as.matrix(mat) # this is submatrix! Be sure to use loccal indexing!!!
	matMain <- matRaw # this is the matrix for heatmap, maybe scaled later
	rowcl <- Getter(mat, 'rowcl')
	colcl <- Getter(mat, 'colcl')
	# when no attributes was set (must not built from subsetMainMat), use user specified case
	# otherwise, attributes from mat will override
	if(is.null(rowcl)) rowcl <-cluster_rows
	if(is.null(colcl)) colcl <-cluster_cols
	if(is.na(truncate)){
		truncate <- ifelse(scale=='none', FALSE, TRUE)
	}
	if(scale != "none"){
		#browser()
		# this modifies mat; need to inherits its attributes
		#mostattributes(mat) <- attrL
		matMain = scale_mat(matRaw, scale)
	}	
	if(clusterWithScaledData){
		matRaw_backup <- matRaw # we still preserve a copy of the raw mat
		matRaw <- matMain
	}
	#browser()
	# clustering of rows
	if(rowcl){
		tree_row = cluster_mat(matRaw, distance = clustering_distance_rows, method = clustering_method)
		rowOrder <- tree_row$order
	}
	else{
		tree_row = NULL
		#rowOrder <- Getter(mat, 'gSel')
		rowOrder <- Getter(mat, 'gSel_l') # local indexing!
	}
	# clustering of cols
	if(colcl){
		#browser()
		tree_col = cluster_mat(t(matRaw), distance = clustering_distance_cols, method = clustering_method)
		columnOrder <- tree_col$order # index based on local reference
	}
	else{
		#tree_col = NA # NULL is better than missing and NA
		tree_col = NULL
		columnOrder <- Getter(mat, 'sSel_l') # index based on local reference
	}
	# truncate mat for vis
	if(truncate){
		if(!is.null(Lower) & !is.null(Upper)){
			temps <- truncByLimit(matMain, Lower=Lower, Upper=Upper)
		} else {
			#browser()
			temps <- truncByQuantile(matMain, q1=q1, q2=q2)
		}
	} else {
		temps <- matMain
	}
	#res <- list(mat4heat=temps, tree_row=tree_row, tree_col=tree_col)
	res <- temps # scaled matrix
	attr(res, 'submat') <- matRaw # raw matrix without scaling/ordering; submatrix
	attr(res, 'tree_row') <- tree_row # dendrogram
	attr(res, 'tree_col') <- tree_col # dendrogram
	attr(res, 'rowcl') <- rowcl # logical
	attr(res, 'colcl') <- colcl # logical
	attr(res, 'rowOrder') <- rowOrder # row clustering order or dummy: index, based on local reference; might be modified by clustering
	attr(res, 'colOrder') <- columnOrder # column clustering order or dummy: index, based on local reference
	attr(res, 'gSel_g') <- Getter(mat, 'gSel_g') # index, based on global reference
	attr(res, 'sSel_g') <- Getter(mat, 'sSel_g') # index, based on global reference
	# aux/meta
	attr(res, 'clustering_distance_rows') <- clustering_distance_rows 
	attr(res, 'clustering_distance_cols') <- clustering_distance_cols 
	attr(res, 'clustering_method') <- clustering_method 
	attr(res, 'scale') <- scale 
	#browser()
	res
}
# z2 <- transform_smat(z)
#tx_L <- transform_smat(smat, Lower=-2.4, Upper=2.4)

tree2dendro <- function(tree){
	if(is.null(tree)){
		res <- FALSE
	} else {
		res <- as.dendrogram(tree)
	}
	res
}
# ensure index of clustering and matrix
ensureClusterIndex <- function(rowcl, colcl, rowOrder, colOrder){
	#
	# to ensure the ordering among: temp, colMat_All, gc and sc
	# when gc or sc is specified, rows or columns will be ordered to construct dendrogram
	# when either one is disabled, need to mannually order rows or columns
	#
	if(rowcl){
		if(colcl){
			# no ordering: driven by dendrogram
			dendrogram <- 'both'
			rowOO <- 1:length(rowOrder) # gc, sc
			colOO <- 1:length(colOrder) # gc, sc
		} else {
			# ordered column: row is not ordered and should be driven by dendrogram 
			dendrogram <- 'row'
			rowOO <- 1:length(rowOrder) 
			colOO <- colOrder
		}
	} else {
		if(colcl){
			# ordered rows: column should not be reordered so as to be dirven by dendrogram
			dendrogram <- 'column'
			rowOO <- rowOrder
			colOO <- 1:length(colOrder)
		} else {
			# both ordered
			dendrogram <- 'none'
			rowOO <- rowOrder
			colOO <- colOrder
		}
	}
	#browser()
	# rowOO: 
	res <- list(dendrogram=dendrogram, rowOO=rowOO, colOO=colOO)
	res
}
#' present a heatmap
#'
#' Constructing a heatmap can be divided into 3 steps:
#'  (1) prepare column bars using prepcolbar() function: specify colorbar color schemes
#'  (2) subsetting main matrix for gene selection and sample selection using subsetMainMat()
#'  (3) transform main matrix including scaling, gene/sample clustering, and ordering. this step follows (2) with smat 
#'      and can be done with transform_smat()
#'  (4) present heatmap (specify main heatmap colors)
#'
#' @param colbar column bar object as returned by prepcolbar(): need to add option for colbar=NULL, no colbar case
#' @param rowbar row bar object as returned by prepcolbar(). Default is NULL, meaning no rowbar
#' @param color color vector for expression data
#' @param ncolor number of colors to be interpolated based on color parameter
#' @param colBarSel a vector to select a subset of column bars; this should be a subset of colnames of pheno data (also colnames of colbar)
#' @param oma oma passed to par()
#' @param plot whether to plot the heatmap
#' @param ... additional parameters to heatmat
#' @return mapL list with additional attributes to reconstruct heatmap
#'  to reconstruct, use:
#'  with(attributes(phm), ...)
#' @export
plotHeatmap <- function(smat_tx, colbar, rowbar=NULL, color=colorpalette2colvec('blueyellow'), 
	ncolor=60, colBarSel=NULL,
	oma=c(2,1,3,11)+0.1, plot=TRUE, cexRow=NULL, cexCol=NULL, labCol = NULL, labRow = NULL, 
	Kcol=NULL, Krow=NULL, labColpal='set1', labRowpal='set2', labRowcolor=NULL, 
	...) {
	#browser()
	mat0 <- smat_tx
	colMat_All <- colbar
	mapL <- Getter(colbar, 'mapL')
	rowcl <- Getter(smat_tx, 'rowcl')
	colcl <- Getter(smat_tx, 'colcl')
	rowOrder <- Getter(smat_tx, 'rowOrder') # this is local indexing
	colOrder <- Getter(smat_tx, 'colOrder') # this is local indexing!
	colFn <- grDevices::colorRampPalette(color)
	ensure <- ensureClusterIndex(rowcl, colcl, rowOrder, colOrder)
	rowOO <- ensure$rowOO #
	colOO <- ensure$colOO
	# trees
	tree_col <- Getter(smat_tx, 'tree_col')
	tree_row <- Getter(smat_tx, 'tree_row')
	#colOrder_g <- ixlocal2global(colOrder, Getter(smat_tx, 'sSel_g')) # this is global indexing
	colOrder_g <- ixlocal2global(colOO, Getter(smat_tx, 'sSel_g')) # this is global indexing; for subsetting samples
	#rowOrder_g <- ixlocal2global(rowOrder, Getter(smat_tx, 'gSel_g')) # this is global indexing
	if(is.null(colBarSel)) {
		colBarSel <- colnames(colMat_All)
	}
	if(!is.null(rowbar)){
		#rowMat <- rowbar
		# currently only deals with one row bar; need to modify heatmat for more rowbars
		#RowSideColors <- rowbar[rowOrder, 1] # notice we use local order here
		RowSideColors <- as.matrix(rowbar[rowOrder, ]) # notice we use local order here (because rowbar is local indexing)
		mapL_row <- Getter(rowbar, 'mapL')
	} else {
		RowSideColors <- NULL
		mapL_row <- NULL
	}
	temp <- mat0[rowOO, colOO]
	#browser()
	#add map for expr bar
	mapExpr <- makecmap(temp, colFn=colFn, n=ncolor)
	if(!is.null(mapL_row)){
		mapL[["geneGroup"]] <- mapL_row[[1]]
	}
	mapL[["Expression"]] <- mapExpr
	Rowv <- tree2dendro(Getter(smat_tx, 'tree_row'))
	Colv=tree2dendro(Getter(smat_tx, 'tree_col'))
	#browser()
	ColSideColors <- as.matrix(colMat_All)[colOrder_g, rev(colBarSel), drop=F] # drop=F very important for 1-column bar!!!
	if(is.null(cexRow))	cexRow <- guessCEX(nrow(temp))
	if(is.null(cexCol))	cexCol <- guessCEX(ncol(temp))
	if(is.null(labCol))	labCol=as.character(colnames(temp))
	if(is.null(labRow))	labRow=as.character(rownames(temp))
	dendrogram=ensure$dendrogram
	col=colFn(ncolor)
	if(!is.null(Kcol)){
		#browser()
		#identical(names(cl_col), labCol) # TRUE
		cl_col <- cutree(tree_col, Kcol)
		labCol <- str_c(labCol, cl_col, sep=' ||')
		dd <-data.frame(cluster=cl_col)
		rownames(dd) <- names(cl_col)
		#browser()
		labColcolor <- prepcolbar(dd, colpal=labColpal)[colnames(temp), 1][tree_col$order]	# this color is in absolute order, not to be ordered by tree		
	} else {
		labColcolor <- NULL
	}
	if(!is.null(Krow)){
		#browser()
		#identical(names(cl_col), labCol) # TRUE
		cl_row <- cutree(tree_row, Krow)
		labRow <- str_c(names(cl_row), cl_row, sep=' ||')
		dd <-data.frame(cluster=cl_row)
		rownames(dd) <- names(cl_row)
		#browser()
		labRowcolor <- prepcolbar(dd, colpal=labRowpal)[rownames(temp), 1][tree_row$order]	# this color is in absolute order, not to be ordered by tree	
	}
	#browser()
	if(plot){
		op <- par(oma=oma) 
		if(identical(colnames(ColSideColors), 'PC1')){
			# no colbar is specified: PC1 is arbitrary contrast
			heatmat(temp, Rowv = Rowv, 
				Colv=Colv, 
				dendrogram=dendrogram,
	        	col=col,
	        	labCol=labCol,
	        	labRow=labRow,
	        	labColcolor=labColcolor, 
	        	RowSideColors=RowSideColors, labRowcolor=labRowcolor,
	        	trace='none',cexRow=cexRow, scale='none', cexCol=cexCol, ...) 
		} else {
			heatmat(temp, Rowv = Rowv, 
				Colv=Colv, 
				dendrogram=dendrogram,
	       		col=col,
	       		labCol=labCol,
	       		labRow=labRow,
	       		ColSideColors=ColSideColors, labColcolor=labColcolor, 
	       		RowSideColors=RowSideColors, labRowcolor=labRowcolor,
	       		trace='none',cexRow=cexRow, scale='none', cexCol=cexCol, ...) 
		}		
		par(op)
	}
	res <- mapL
	pheno <- Getter(colbar, 'pheno')
	gpheno <- Getter(rowbar, 'pheno') # maybe NULL
	#browser()
	# BUG fixed 10/03/2017: why it is FALSE now?
	# if(!is.null(Colv)) { 
	if(class(Colv)=='dendrogram') {
		hc <- as.hclust(Colv)
		colo <- hc$order
	} else {
		colo <- 1:ncol(temp)
	}	
	#browser()
	#bug fix: if(!is.null(Rowv)){
	if(class(Rowv)=='dendrogram'){
		hr <- as.hclust(Rowv)
		rowo <- rev(hr$order)
	} else {
		rowo <- 1:nrow(temp)
	}
	hdat <- temp[rowo, colo]
	#browser()
	attr(res, 'pheno') <- pheno
	attr(res, 'gpheno') <- gpheno
	attr(res, 'tree_row') <- tree_row
	attr(res, 'tree_col') <- tree_col
	attr(res, 'oma') <- oma
	attr(res, 'temp') <- temp
	attr(res, 'hdat') <- hdat # data showing in heatmap
	attr(res, 'Rowv') <- Rowv
	attr(res, 'Colv') <- Colv
	attr(res, 'dendrogram') <- ensure$dendrogram
	attr(res, 'col') <- col
	attr(res, 'labCol') <- labCol
	attr(res, 'labRow') <- labRow
	attr(res, 'ColSideColors') <- ColSideColors
	attr(res, 'ColSideColorsIndex_g') <- colOrder_g # global index for colside colors matrix; this is useful if we want to testing association with colside bar in present of sample subsetting
	attr(res, 'RowSideColorsIndex_l') <- rowOrder # local indexing
	attr(res, 'RowSideColors') <- RowSideColors
	attr(res, 'cexRow') <- cexRow
	attr(res, 'cexCol') <- cexCol
	attr(res, 'scale') <- 'none'
	## add command for reconstruction
	attr(res, 'cmd') <- "
	op <- par(oma=oma) 
	heatmat(temp, Rowv = Rowv, 
			Colv=Colv, 
			dendrogram=dendrogram,
	        col=col,
	        labCol=labCol,
	        labRow=labRow,
	        ColSideColors=ColSideColors,
	        RowSideColors=RowSideColors,
	        trace='none',cexRow=cexRow, scale='none', cexCol=cexCol) 
	par(op)
	"
	#attr(res, 'class') <- 'mapL'
	attr(res, 'class') <- c('phm', 'list')
	res
}
#phm <- plotHeatmap(obj_tx, colbar=colbar_rel1)

#' testree method for class phm
#' 
#' this function test row tree or column tree with specified pheno/gpheno data
#'
#' @param obj an phm object
#' @param orientation specify what tree to test with, row or column
#' @param bar column name to specify the pheno data (column of pheno or gpheno)
#' @param k k as in cutree
#' @param h h as in cutree
#' @param plot whether to draw a figure or just test result
#' @method testree phm
#' @return test result or nothing
#' @export
testree.phm <- function(obj, orientation='column', bar=NULL, k=3, h=NULL, plot=TRUE){
	if(orientation=='row'){
		tree <- Getter(obj, 'tree_row')
		pheno <- Getter(obj, 'gpheno')
		sel <- Getter(obj, 'RowSideColorsIndex_l')
	} else {
		tree <- Getter(obj, 'tree_col')
		pheno <- Getter(obj, 'pheno')
		sel <- Getter(obj, 'ColSideColorsIndex_g')
	}
	if(is.null(bar)) bar <- colnames(pheno)[1]
	#browser()
	# bug! subsetting is not consistent: vec <- pheno[, bar]	
	vec <- pheno[sel, bar]
	test_tree_split(tree=tree, vec=vec, k=k, h=h, plot=plot, ylab=bar, tag='Association with cluster membership')
}
#' test association between tree cut and a vector
#'
#' @param tree a dendrogram 
#' @param k k as in cutree
#' @param h h as in cutree
#' @param plot whether to draw a figure or just test result
#' @param vec value to be tested
#' @param ... additional parameters passed to plotScatter
#' @export
test_tree_split <- function(tree, k=3, h=NULL, vec, plot=TRUE, ...){
	group <- as.factor(cutree(tree, k=k, h=h))
	cat(sectionString('Summary of clusters by cutree'))
	print(table(group))
	cat('\n')
	cat(sectionString('Hypothesis test result'))
	#browser()
	if(plot){
		print(plotScatter(group, vec, xlab='Clusters', ...))
	} else {
		z <- plotScatter(group, vec, xlab='Clusters', tag='', plot=F, ...)
		res <- z[! names(z) %in% 'p']
		return(res)
	}	
}
#test_tree_split(tree=Getter(phm, 'tree_row'), vec=Getter(phm, 'gpheno')$group)
#test_tree_split(tree=Getter(phm, 'tree_col'), vec=Getter(phm, 'pheno')$LKB1expr)
#
# hkey2 stretch is dependent on number of breaks; here we empirically determins it
# the rule is: n_breaks * stretch = 10

#' compute stretch from a continuous map
#'
#' @param map a map
#' @return a value
#' @export
get_stretch <- function(map, ref=10){
	ref/length(map$breaks)
}
#' print method for phm class
#'
#' @method print phm
#' @param obj a phm class (a data frame)
#' @param ... additional parameters to plot4rowanovaMatVec
#' @export
print.phm <- function(obj){
	cat('Attributes included in phm class:\n')
	cat(Getter(obj))
	cat('\n')
	#obj[1:length(obj)]	
	cat(sprintf('cexRow=%.4f; cexCol=%.4f\n', Getter(obj, 'cexRow'), Getter(obj, 'cexCol')))
	cat(sprintf('cmd=%s\n', Getter(obj, 'cmd')))
}

#' plot method for phm class
#'
#' @method plot phm
#' @param obj a phm class (a data frame)
#' @param ... additional parameters to plot4rowanovaMatVec
#' @export
plot.phm <- function(obj, y0=1.0, x0=0.8, yd=0.17) {
	#browser()
	mapL <- obj[1:length(obj)]
	# 2015/08/25: conditionally decide if to show NA in the legend; this should only affects the legend, not other parts
	ColSideColorsIndex_g <- Getter(obj, 'ColSideColorsIndex_g') # this index is good for map$fac
	nLegend <- length(mapL)
	lorder <- c(nLegend, 1:(nLegend-1))
	i <- 1
	for(lo in lorder){
		ln <- names(mapL)[lo]
		map <- mapL[[lo]]
		catmap <- ismapcat(map)
		y <- y0-yd*(i-1)
		i <- i+1
		#browser()
		#if(lo==4) browser()
		if(catmap){
			#browser()
			# show NA only if NA is present in the column bar
			showNA <- any(is.na(map$fac[ColSideColorsIndex_g]))
			# 2015/08/31do not want to show the artificial PC1 column bar
			if(ln!='PC1'){
				vlegend(map, ln, y=y, x=x0, side=4, cex=0.6, stretch=0.8, showNA=showNA)
				#browser()
				cat(sprintf('vlegend(phm[[\'%s\']], \'%s\', y=%.1f, x=%.1f, side=4, cex=0.6, stretch=0.8, showNA=%s) \n', ln, ln, y, x0, showNA))
			}
		} else {
			#browser()
			hkey2(map, ln, y=y, x=x0, stretch=get_stretch(map), side=1)  
			cat(sprintf('hkey2(phm[[\'%s\']], \'%s\', y=%.1f, x=%.1f, stretch=%.3f, side=1) \n', ln, ln, y, x0, get_stretch(map)))

		}
		#browser()
	}
}
#plot.mapL(mapL)
#plotHeatmap(obj_tx, colbar=colbar_rel1, color=colorpalette2colvec('cyblyl'))

#plotHeatmap(obj_tx, colbar=colbar_rel1, color=colorpalette2colvec('bluered'), colBarSel=colBarOrder_all_addPredCl)

#plotHeatmap(obj_tx, colbar=colbar_rel1)

#' update color pal (as input for prepcolbar)
#'
#' given an old pal and a list of newly specified pal, update when compatible; when no shared names, no updates will be made
#' @param oldpal a list of colors, usually returned by autocolpal4DF
#' @param newpal a list of colors, usually user specified; default is NULL, which means no updates at all.
#' @return a named list of colors 
#' @export
updatePal <- function(oldpal, newpal=NULL)	{
	#browser()
	spalnames <- intersect(names(newpal), names(oldpal))
	for(name in spalnames) {
		oldpal[[name]] <- newpal[[name]]
		#browser()
	}
	oldpal
}
#updatePal(oldpal=colpal_rel1, newpal=newpal)

#' annotated heatmap using heatmat engine
#'
#' @param mat matrix for heatmap
#' @param pheno a data frame where columns represent different column bars (in a heatmap)
#'  notice: always do not subset pheno! 
#' @param gSel gene selection. This can be index or gene names. Internally, both will be converted to integer based index (global or local)
#' @param sSel0 sample selection filter 0; by default, this deals with augMat to remove all NA columns
#'  the user can override this by supplying a vector of indeces or sample names.
#' @param sSel1 sample selection (index or sample names). This is useful to select a subset of samples based
#' on phenotype, i.e. some mutation, like ++ to be removed
#' @param colOrderIndex index/sample names to order column; The actual number of samples surviving is an intersect between colOrderIndex and sSel (from sSel0 and sSel1).
#'  if not specified, column clustering will be instructed to construct; otherwise, no clumn clustering will be done later on
#' @param rowOrderIndex index (integer or gene names) to order rows; The actual genes surviving is an intersect between rowOrderIndex and gSel.
#'  if not specified, row clustering will be instructed to construct; otherwise, no row clustering will be done later on. To disable row clustering, need to specify rowOrderIndex and cluster_rows=FALSE simultaneously
#' @param scale selection from c("none", "row", "column")
#' @param clusterWithScaledData logical indicating if use scaled data for clustering; default is FALSE
#' @param cluster_rows whether to cluster rows
#' @param cluster_cols whether to cluster columns 
#' @param clustering_distance_rows distance metric for rows
#' @param clustering_distance_cols distance metric for cols
#' @param clustering_method clustering method
#' @param truncate logical indicating if truncation is needed; default is NA will enable truncate if scale!='none'
#' @param Lower parameter Lower to truncByQuantile()
#' @param Upper parameter Upper to truncByQuantile()
#' @param q1 parameter q1 to truncByLimit
#' @param q2 parameter q2 to truncByLimit
#' @param colbar column bar object as returned by prepcolbar(): need to add option for colbar=NULL, no colbar case
#' @param rowbar row bar object as returned by prepcolbar(). Default is NULL, meaning no rowbar
#' @param color color vector for expression data
#' @param ncolor number of colors to be interpolated based on color parameter
#' @param colBarSel a vector to select a subset of column bars; this should be a subset of colnames of pheno data (also colnames of colbar)
#' @param oma oma passed to par()
#' @param plot whether to plot the heatmap
#' @param plotLegend whether to plot color legend
#' @param  colpal user specified column-wise color palette (a named list of color vector); only shared names will be used to update default colpal. Default is NULL
#'   thus no updates at all, only focusing on default
#' @param  rowpal user specified row-wise color palette (a named list of color vector); only shared names will be used to update default rowpal. Default is NULL
#'   thus no updates at all, only focusing on default
#' @param gpheno a one-column gene annotation data frame. (need more work to include multiple row bars and check indexing)
#' @param ... additional parameters to heatmat
#' @return mapL list with additional attributes to reconstruct heatmap (returned by plotHeatmap, a phm class)
#' @export
aheatmat = function(mat, pheno=NULL, # starting pars in subsetMainMat
	gSel=NULL, sSel0=NULL, sSel1=NULL, 
	colOrderIndex=NULL, rowOrderIndex=NULL, # starting pars in transform_smat, need to be specified, not null to disable row tree
	scale = "row", clusterWithScaledData=FALSE, cluster_rows = TRUE, 
	cluster_cols = TRUE, 
	clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
	clustering_method = "ward.D",
	truncate=NA, q1=0.01, q2=0.99, Lower=NULL, Upper=NULL, 
	cexRow=NULL, cexCol=NULL, labCol = NULL, labRow = NULL, # starting pars in plotHeatmap
	colbar=NULL, rowbar=NULL, color=colorpalette2colvec('bluered'), 
	ncolor=60, colBarSel=NULL, # starting plot.mapL pars
	y0=1.0, x0=0.8, yd=0.17, 
	oma=c(0,1,3,11)+0.1, plot=TRUE, # added pars
	colpal=NULL, rowpal=NULL, gpheno=NULL, plotLegend=TRUE, labRowcolor=NULL, 
	...
	){
	#if(is.null(rownames(mat))) rownames(mat) <- 1:nrow
	#if(is.null(colnames(mat))) colnames(mat) <- 1:ncol
	mat <- namedmat(mat)
	if(is.null(pheno)){
		#PC1 <- calcpca(mat[, scale = T, center = T)$projmat[, 1]
		pheno <- data.frame(PC1=rep('Class', ncol(mat)))
		rownames(pheno) <- colnames(mat)
	}
	if(!identical(colnames(mat), rownames(pheno))) {
		stop(sprintf('colnames(mat) and rownames(pheno) does not match!\ncolnames(mat): \n%s\n\nrownames(pheno): \n%s\n', 
			str_c(colnames(mat), collapse=', '), str_c(rownames(pheno), collapse=', ')))
	}
	colpal0 <- autocolpal4DF(dat=pheno) # default
	# update colpal0 from specified pal (named list)
	if(!is.null(colpal)){
		#browser()
		colpal0 <- updatePal(oldpal=colpal0, newpal=colpal)	
	}
	colbar0 <- prepcolbar(pheno, colpal=colpal0)
	# in case the user has specified colbar, this overrides all !!!
	if(!is.null(colbar)){
		colbar0 <- colbar 
	}
	# gene bar
	rowpal0 <- autocolpal4DF(dat=gpheno) # default
	# update colpal0 from specified pal (named list)
	rowpal0 <- updatePal(oldpal=rowpal0, newpal=rowpal)	
	#browser()
	rowbar0 <- prepcolbar(gpheno, colpal=rowpal0)
	# in case the user has specified rowbar, this overrides all !!!
	if(!is.null(rowbar)){
		rowbar0 <- rowbar 
	}
	#browser()
	# colbar0: subsetting
	## subset data
	smat <- subsetMainMat(mat=mat, gSel=gSel, sSel0=sSel0, sSel1=sSel1, cluster_cols=cluster_cols, cluster_rows=cluster_rows, colOrderIndex=colOrderIndex, rowOrderIndex=rowOrderIndex)
	## transform data
	#browser()
	obj_tx <- transform_smat(smat, scale = scale, clusterWithScaledData=clusterWithScaledData, cluster_rows = cluster_rows, 
		cluster_cols = cluster_cols, 
		clustering_distance_rows = clustering_distance_rows, clustering_distance_cols = clustering_distance_cols, 
		clustering_method = clustering_method,
		truncate=truncate, q1=q1, q2=q2, Lower=Lower, Upper=Upper)
	## plot heatmap
	phm <- plotHeatmap(obj_tx, colbar=colbar0, rowbar=rowbar0, color=color, ncolor=ncolor, colBarSel=colBarSel, oma=oma, plot=plot, 
		cexRow=cexRow, cexCol=cexCol, labCol = labCol, labRow = labRow, labRowcolor=labRowcolor, ...)
	#browser()
	if(plotLegend) {
		plot(phm, y0=y0, x0=x0, yd=yd)
	}
	#browser()
	phm
}	
#aheatmat(mat=rna, gSel=rownames(gpheno), pheno=pheno_rel1, q1=0.1, q2=0.9, clustering_distance_rows = "euclidean", 
#	clustering_distance_cols = "euclidean", clustering_method = "ward.D", colbar=colbar_rel1, rowbar=rowbar)

#' plot correlation heatmap for rows from a matrix
#'
#' @param mat a matrix where correlation will be computed from its rows
#' @param method correlation method
#' @param main main title
#' @param ykey y value for correlation key
#' @export
aheatmatcor <- function(mat, method='pearson', CexRow=NULL, CexCol=NULL, oma=c(2,1,4,3)+0.1, dendrogram='row', main='', ykey=1.05, xkey=0.05, stretch=0.08, colFn=bluered, ncol=128, return=FALSE){
	#browser()
	cormat <- cor(t(mat), method=method, use="pairwise.complete.obs")
	## a computation engine for heatmap
	phm <- aheatmat(cormat, scale='none', plotLegend=F, plot=F, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")	
	#colFn <- bluered; n <- 128
	n <- ncol
	### breaks around 0
	n1 <- floor(n/2); n2 <- n-n1
	temp <- Getter(phm, 'temp')
	breaks <- c(seq(min(temp, na.rm=T), 0, length=n1+1), seq(0, max(temp, na.rm=T), length=n2+1)[-1])
	#oma=c(2,1,6,3)+0.1
	op <- par(oma=oma) 
	with(attributes(phm),
	{
	if(is.null(CexRow)) 
		cexRow <- cexRow*0.9
	else cexRow <- CexRow
	if(is.null(CexCol)) 
		cexCol <- cexCol*0.9
	else cexCol <- CexCol
	heatmat(temp, Rowv = Rowv, 
					Colv=Colv, 
					dendrogram='row',
					col=colFn(n),
					labCol=labCol,
					labRow=labRow,
					labColcolor=NULL, 
					RowSideColors=NULL, labRowcolor=NULL, main=main,
					trace='none',cexRow=cexRow, scale='none', cexCol=cexCol, breaks=breaks) 
	}
	)
	hkey2(makecmap(temp, colFn=colFn, n=n, breaks= breaks), y=ykey, x=xkey, stretch=stretch, title="Correlation", side=1, digits=2) # make sure breaks
	par(op)	
	if(return)
		return(list(cormat=cormat, phm=phm))
}
#aheatmatcor(mat_helicase)