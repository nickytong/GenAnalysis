####
#### to build
####
if(FALSE){
# cd /home/ptong1/Backup/Package/; R --vanilla
library(roxygen2)
#library(roxygen) # not working
roxygenize("GenAnalysis")

library(devtools)
build('GenAnalysis')
build('GenAnalysis')
install('GenAnalysis')

load_all('GenAnalysis')

devtools::load_all('/home/ptong1/Backup/Package/GenAnalysis')
source(file.path(osbaseZdrive, 'ptong1/Backup/RSource/heatmat.R')) # source my utility function: i.e. heatmat, createFolder

build_win('GenAnalysis')

##
detach("package:GenAnalysis", unload=TRUE)
library(GenAnalysis)

    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/GenAnalysis/R/aheatmap.R'))

	
	
###
### compare to pheatmap
###
cd /workspace/ptong1/Projects/PackageSourceCode/
library(roxygen2)
#library(roxygen) # not working
roxygenize("pheatmap")

library(devtools)
load_all('pheatmap')
build('pheatmap')
install('pheatmap')
}
####
#### to devel
####
if(FALSE){
# (1) disable the library(GenAnalysis) command in heatmat.R file
# (2) it's a good practice to save the old version at: Z:\Backup\Package\GenAnalysisHistory
# (3) in the script calling this package, source the R file, located at /home/ptong1/Backup/Package/
    if(FALSE){
    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/GenAnalysis/R/aheatmap.R'))
    source(file.path(osbaseZdrive, 'ptong1/Backup/Package/GenAnalysis/R/myplots_p.R'))
	}
# (4) interactively modify the script

}


######################
lo = function(rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, shrink=0.95, treeheight_col, treeheight_row, legend, annotation, annotation_colors, 
	annotation_legend, main, fontsize, fontsize_row, fontsize_col, cexRow=cexRow, cexCol=cexCol, ...){
	# Get height of colnames and length of rownames
	if(!is.null(coln[1])){
		longest_coln = which.max(strwidth(coln, units = 'in'))
		gp = list(fontsize = fontsize_col, ...)
		coln_height = unit(1, "grobheight", textGrob(coln[longest_coln], rot = 90, gp = do.call(gpar, gp))) + unit(5, "bigpts")
	}
	else{
		coln_height = unit(5, "bigpts")
	}
	# when there is rowname, get the longest one as a reference
	if(!is.null(rown[1])){
		longest_rown = which.max(strwidth(rown, units = 'in'))
		gp = list(fontsize = fontsize_row, ...)
		#rown_width = unit(1, "grobwidth", textGrob(rown[longest_rown], gp = do.call(gpar, gp))) + unit(10, "bigpts")
		# now scales the width so no extra space will be wasted
		rown_width = unit(1, "grobwidth", textGrob(rown[longest_rown], gp = do.call(gpar, gp)))*cexRow + unit(10, "bigpts")
	}
	else{
		rown_width = unit(5, "bigpts") # when no rownames
	}
	
	gp = list(fontsize = fontsize, ...)
	# Legend position
	if(!is.na(legend[1])){
		longest_break = which.max(nchar(names(legend)))
		longest_break = unit(1.1, "grobwidth", textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))
		title_length = unit(1.1, "grobwidth", textGrob("Scale", gp = gpar(fontface = "bold", ...)))
		legend_width = unit(12, "bigpts") + longest_break * 1.2
		legend_width = max(title_length, legend_width)
	}
	else{
		legend_width = unit(0, "bigpts")
	}
	
	# Set main title height
	if(is.na(main)){
		main_height = unit(0, "npc")
	}
	else{
		main_height = unit(1.5, "grobheight", textGrob(main, gp = gpar(fontsize = 1.3 * fontsize, ...)))
	}
	
	# Column annotations
	if(!is.na(annotation[[1]][1])){
		# Column annotation height : this is for the col-bar height
		annot_height = unit(ncol(annotation) * (8 + 2) + 2, "bigpts")
		# Width of the correponding legend; annot_legend_width: this is for the legend width
		longest_ann = which.max(nchar(as.matrix(annotation))) # this can be increased by passing a parameter!!!
		annot_legend_width = unit(1.2, "grobwidth", textGrob(as.matrix(annotation)[longest_ann], gp = gpar(...))) + unit(12, "bigpts")
		if(!annotation_legend){
			annot_legend_width = unit(0, "npc")
		}
	}
	else{
		annot_height = unit(0, "bigpts")
		annot_legend_width = unit(0, "bigpts")
	}
	
	# Tree height
	treeheight_col = unit(treeheight_col, "bigpts") + unit(5, "bigpts")
	treeheight_row = unit(treeheight_row, "bigpts") + unit(5, "bigpts") 
	
	# Set cell sizes
	if(is.na(cellwidth)){
		matwidth = unit(1, "npc") - rown_width - legend_width - treeheight_row - annot_legend_width
	}
	else{
		matwidth = unit(cellwidth * ncol, "bigpts")
	}
	
	if(is.na(cellheight)){
		matheight = unit(1, "npc") - main_height - coln_height - treeheight_col - annot_height
	}
	else{
		matheight = unit(cellheight * nrow, "bigpts")
	}	
	matwidth <- matwidth*shrink
	matheight <- matheight*shrink
	
	# Produce layout()
	# this sets up the global layout that should obey; further function calls to vplayout is based on this!!!
	pushViewport(viewport(layout = grid.layout(nrow = 5, ncol = 5, widths = unit.c(treeheight_row, matwidth, rown_width, legend_width, annot_legend_width), heights = unit.c(main_height, treeheight_col, annot_height, matheight, coln_height)), gp = do.call(gpar, gp)))
	
	# Get cell dimensions: 4, 2 is the column legend
	pushViewport(vplayout(4, 2))
	cellwidth = convertWidth(unit(0:1, "npc"), "bigpts", valueOnly = T)[2] / ncol
	cellheight = convertHeight(unit(0:1, "npc"), "bigpts", valueOnly = T)[2] / nrow
	upViewport()
	
	# Return minimal cell dimension in bigpts to decide if borders are drawn
	mindim = min(cellwidth, cellheight)  # what does this do
	return(mindim)
}

### the dendrogram is draw with grid system from scratch!!!
### default is to draw column dendrogram
draw_dendrogram = function(hc, horizontal = T){
	h = hc$height / max(hc$height) / 1.05
	m = hc$merge # this indicates at each iteration what sample is merged: n-1 by 2 matrix
	o = hc$order
	n = length(o) # sample size

	m[m > 0] = n + m[m > 0]  # for positive values, add n to it since it merges with a conceptual node
	m[m < 0] = abs(m[m < 0]) # for negative values, use its abs value since it is singleton merging with actual nodes
	### stores x (positon) and y (height)
	dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y"))) 
	# for 1:n rows, dist[, 1] has min: 1/n/2, max: (2n-1)/(2n)   -------> this is the middle; thus the n leaves at the bottome
	# for these bottome leaves, y=0
	dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)
	### 2*n - 1 , but only updates rows (n+1):(2n-1)
	for(i in 1:nrow(m)){
		### computes the middle; for node 1:n, since x1, x2 are equal, there is no effect; y is also both 0
		### where there is a conceptual node, it gets averaged
		### it may seem a conceptual node have x=0; but this will not happen since n+1 will always be averaged from previous computation! Thus it is always an average
		### another way to say it, starting from the first few rows with hc$merge negative, the rows after row n will always be filled with nonzero!!!
		dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2 
		dist[n + i, 2] = h[i]
	}
	### this piece draws a |^^|, that is, three segments connected. IN total, there is n-1 such shapes.
	### we need to specify x1, x2, y1, y2 that tells us the two x position and two y positions
	### x1, x2 represents the sample index; y1, y2 relates to the height!
	draw_connection = function(x1, x2, y1, y2, y){
		grid.lines(x = c(x1, x1), y = c(y1, y))
		grid.lines(x = c(x2, x2), y = c(y2, y))
		grid.lines(x = c(x1, x2), y = c(y, y))
	}
	
	if(horizontal){
		for(i in 1:nrow(m)){
			draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
		}
	}	
	else{
		gr = rectGrob()
		# the trick: rotate 90 degree!
		# notice gr fills the current viewport fully; height and width computes based on this gr
		pushViewport(viewport(height = unit(1, "grobwidth", gr), width = unit(1, "grobheight", gr), angle = 90))
		dist[, 1] = 1 - dist[, 1] 
		for(i in 1:nrow(m)){
			draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
		}
		upViewport() # go back to parent vp where we starts
	}
}
### matrix is processed into a color matrix; since using grid.rect, it is easy to add a border!!!
### it is also easy to add text by directly overlay on it
### added rot: rotation
draw_matrix = function(matrix, border_color, fmat, fontsize_number, rot=0){
	n = nrow(matrix)
	m = ncol(matrix)
	x = (1:m)/m - 1/2/m
	y = 1 - ((1:n)/n - 1/2/n)
	for(i in 1:m){
		grid.rect(x = x[i], y = y[1:n], width = 1/m, height = 1/n, gp = gpar(fill = matrix[,i], col = border_color))
		if(attr(fmat, "draw")){
			grid.text(x = x[i], y = y[1:n], label = fmat[, i], gp = gpar(col = "grey30", fontsize = fontsize_number), rot=rot)
		}
	}
}

## draws colnames at the bottom
draw_colnames = function(coln, ...){
	m = length(coln)
	x = (1:m)/m - 1/2/m # min: 1/(2m); max: 1-1/(2m)
	### x is easy; y is always set at 0.96npc! we can specify rotation!
	grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = 0.5, hjust = 0, rot = 270, gp = gpar(...))
}
## draw rownames at the right
draw_rownames = function(rown, ...){
	n = length(rown)
	y = 1 - ((1:n)/n - 1/2/n)
	#browser()
	grid.text(rown, x = unit(0.04, "npc"), y = y, vjust = 0.5, hjust = 0, gp = gpar(...))	
}

draw_legend = function(color, breaks, legend, ...){
	height = min(unit(1, "npc"), unit(150, "bigpts"))
	pushViewport(viewport(x = 0, y = unit(1, "npc"), just = c(0, 1), height = height))
	legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
	breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
	h = breaks[-1] - breaks[-length(breaks)]
	grid.rect(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
	grid.text(names(legend), x = unit(12, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))
	upViewport()
}

#' convert annotation (data frame) and annotation_colors (a list) to color matrix
#'
#' @param annotation a data frame of continuous or categorical variables
#' @param annotation_colors a list of color pallete; for continuous variable, specify a vector of colors; for categorical or factor
#' variable, speacify a named vector; for NA, the named vector should contain 'NA'
#' @return a data frame of colors
#' @export
convert_annotations = function(annotation, annotation_colors){
	new = annotation
	for(i in 1:ncol(annotation)){
		a = annotation[, i]
		b = annotation_colors[[colnames(annotation)[i]]]
		if(is.character(a) | is.factor(a)){
			a = as.character(a)
			a[is.na(a)] <- 'NA' # to match with 'NA' col
			if(length(setdiff(a, names(b))) > 0){
				stop(sprintf("Factor levels on variable %s do not match with annotation_colors", colnames(annotation)[i]))
			}
			new[, i] = b[a]
		}
		else{
			a = cut(a, breaks = 100)
			new[, i] = colorRampPalette(b)(100)[a]
		}
	}
	return(as.matrix(new))
}

# draws columnbars, the annotation by passing a color matrix
draw_annotations = function(converted_annotations, border_color){
	n = ncol(converted_annotations) # number of col bars
	m = nrow(converted_annotations)
	# x is the center of each cell
	x = (1:m)/m - 1/2/m # min: 1/(2m); max: 1-1/(2m) where m is number of sample
	### don't fully understand this. So each column bar has 8 bigpts, but separated by 2bigpoints (note the whiteline inbetween bars!)
	### using cumsum since at each different y, we'll draw a rect
	y = cumsum(rep(8, n)) - 4 + cumsum(rep(2, n)) ## this is further transformed into big points. 
	#browser()
	## iterate each sample, draws all bars for that sample simultaneously! 
	## rect.width = 1/m, filling the entire region, while height is 8 bigpts
	for(i in 1:m){
		###
		## this is the order from last col to first column; but can be changed by y[n:1]
		###
		#grid.rect(x = x[i], unit(y[1:n], "bigpts"), width = 1/m, height = unit(8, "bigpts"), gp = gpar(fill = converted_annotations[i, ], col = border_color))
		grid.rect(x = x[i], unit(y[n:1], "bigpts"), width = 1/m, height = unit(8, "bigpts"), gp = gpar(fill = converted_annotations[i, ], col = border_color))
	}
}

draw_annotation_legend = function(annotation, annotation_colors, border_color, ...){
	y = unit(1, "npc") # from the top of this vp
	text_height = unit(1, "grobheight", textGrob("FGH", gp = gpar(...)))
	for(i in names(annotation_colors)){
		grid.text(i, x = 0, y = y, vjust = 1, hjust = 0, gp = gpar(fontface = "bold", ...))
		y = y - 1.5 * text_height # y decreases a bit
		if(is.character(annotation[, i]) | is.factor(annotation[, i])){
			for(j in 1:length(annotation_colors[[i]])){
				grid.rect(x = unit(0, "npc"), y = y, hjust = 0, vjust = 1, height = text_height, width = text_height, gp = gpar(col = border_color, fill = annotation_colors[[i]][j]))
				grid.text(names(annotation_colors[[i]])[j], x = text_height * 1.3, y = y, hjust = 0, vjust = 1, gp = gpar(...))
				y = y - 1.5 * text_height # y decreases for each color drawn
			}
		}
		else{
			# for continous col legend, decrease is differnt
			yy = y - 4 * text_height + seq(0, 1, 0.02) * 4 * text_height
			h = 4 * text_height * 0.02
			grid.rect(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = h, width = text_height, gp = gpar(col = "#FFFFFF00", fill = colorRampPalette(annotation_colors[[i]])(50)))
			txt = rev(range(grid.pretty(range(annotation[, i], na.rm = TRUE))))
			yy = y - c(0, 3) * text_height
			grid.text(txt, x = text_height * 1.3, y = yy, hjust = 0, vjust = 1, gp = gpar(...))
			y = y - 4.5 * text_height
		}
		y = y - 1.5 * text_height # this affects actually the gap between
	}
}

draw_main = function(text, ...){
	grid.text(text, gp = gpar(fontface = "bold", ...))
}

vplayout = function(x, y){
	return(viewport(layout.pos.row = x, layout.pos.col = y))
}

## 2 step procedure to get rowname or colname: first test if labRow/labCol is specified; then conditional choose rowname or colname
getName <- function(mat, name, type='row'){
	if(is.null(name)) return(name) # instantly return NULL if it is NULL which is due to show_rownames/show_colnames=FALSE
	if(!is.na(name)[1]) {
		res <- name
	} else {
		if(type=='row') {
			res <- rownames(mat)
		} else {
			res <- colnames(mat)
		}		
	}
	res
}
# filename=NA by default; what does this mean? We do see a figure window poped out! Does this mean nothing is generated?
# yes: when it is NA, there will be the figure drawed in the R window; if not, then the heatmap_motor() called itself by setting filename=NA, smart!!!!!!
# the input matrix is now a color matrix!!! (by scale_colours() function)
# major problem of pheatmap is that when there are many genes or samples, the smallest fontsize is 1, leading to text size to large or line widt too large.
# based on gpar() and get.gpar(), The size of text is fontsize*cex. The size of a line is fontsize*cex*lineheight. 
# added: cexRow, cexCol to enable shrinking rowname and colnames
# added labRow, labCol: this is passed with the ordering of image. Thus, if ordered by dendrogram this should also be ordered, not the original lab!
# added annoHeaderPosition: left, which is the old one that adds annotation bar names at the left of the bars (c(3, 1)); now we add right so that the barname shows to the right. This is needed
#       since when no row dendrogram, the header will not show completely
heatmap_motor = function(matrix, border_color, cellwidth, cellheight, shrink=0.95, tree_col, tree_row, treeheight_col, treeheight_row, filename, 
	width, height, breaks, color, legend, annotation, annotation_colors, annotation_legend, main, 
	fontsize, fontsize_row, fontsize_col, fmat, fontsize_number, cexRow, cexCol, labRow=NA, labCol=NA, annoHeaderPosition=NA, rot=0, ...){
	grid.newpage()	
	# Set layout
	rown <- getName(mat=matrix, name=labRow, type='row')
	coln <- getName(mat=matrix, name=labCol, type='column')
	# decide annoHearder position
	if(is.na(annoHeaderPosition)){
		if(treeheight_row==0) {
		annoHeaderPosition <- 'right' # norw row tree: thus add to the right
		} else {
		annoHeaderPosition <- 'left'
		}
	}
	if(!is.na(labCol)[1]) coln <- labCol
	#mindim = lo(coln = colnames(matrix), rown = rownames(matrix), nrow = nrow(matrix), ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, cexRow=cexRow, cexCol=cexCol,  ...)
	mindim = lo(coln = coln, rown = rown, nrow = nrow(matrix), ncol = ncol(matrix), cellwidth = cellwidth, shrink=shrink, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, cexRow=cexRow, cexCol=cexCol,  ...)
	if(!is.na(filename)){
		pushViewport(vplayout(1:5, 1:5))
		
		if(is.na(height)){
			height = convertHeight(unit(0:1, "npc"), "inches", valueOnly = T)[2]
		}
		if(is.na(width)){
			width = convertWidth(unit(0:1, "npc"), "inches", valueOnly = T)[2]
		}		
		# Get file type
		r = regexpr("\\.[a-zA-Z]*$", filename)
		if(r == -1) stop("Improper filename")
		ending = substr(filename, r + 1, r + attr(r, "match.length"))		
		f = switch(ending,
			pdf = function(x, ...) pdf(x, ...),
			png = function(x, ...) png(x, units = "in", res = 300, ...),
			jpeg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
			jpg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
			tiff = function(x, ...) tiff(x, units = "in", res = 300, compression = "lzw", ...),
			bmp = function(x, ...) bmp(x, units = "in", res = 300, ...),
			stop("File type should be: pdf, png, bmp, jpg, tiff")
		)
		## smart: calls filename=NA!!!
		# print(sprintf("height:%f width:%f", height, width))
		f(filename, height = height, width = width)
		heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, cexRow=cexRow, cexCol=cexCol, ...)
		dev.off()
		upViewport()
		return()
	}	
	# Omit border color if cell size is too small 
	if(mindim < 3) border_color = NA	
	# Draw title
	if(!is.na(main)){
		pushViewport(vplayout(1, 2))
		draw_main(main, fontsize = 1.3 * fontsize, ...)
		upViewport()
	}	
	#browser()
	# Draw tree for the columns
	if(!is.na(tree_col[[1]][1]) & treeheight_col != 0){
		pushViewport(vplayout(2, 2))
		draw_dendrogram(tree_col, horizontal = T)
		upViewport()
	}	
	# Draw tree for the rows
	if(!is.na(tree_row[[1]][1]) & treeheight_row != 0){
		pushViewport(vplayout(4, 1))
		draw_dendrogram(tree_row, horizontal = F)
		upViewport()
	}	
	# Draw matrix
	#browser()
	pushViewport(vplayout(4, 2))
	draw_matrix(matrix, border_color, fmat, fontsize_number, rot=rot)
	upViewport()	
	#browser()
	# Draw colnames: note we add cexCol for enhanced power
	#if(length(colnames(matrix)) != 0){
	if(length(coln) != 0){
		# already taken care in aheatmap
		#if(is.na(cexCol))
		#	cexCol <- guessCEX(ncol(matrix))
		pushViewport(vplayout(5, 2))
		pars = list(coln, fontsize = fontsize_col, cex=cexCol, ...)
		if (fontsize_col!=0)
			do.call(draw_colnames, pars)
		upViewport()
	}
	# Draw rownames
	#if(length(rownames(matrix)) != 0){
	if(length(rown) != 0){
		#if(missing(cexRow))
		#	cexRow <- guessCEX(nrow(matrix))
		pushViewport(vplayout(4, 3))
		pars = list(rown, fontsize = fontsize_row, cex=cexRow, ...)
		if(fontsize_row!=0)
			do.call(draw_rownames, pars)
		upViewport()
	}	
	# Draw annotation tracks: this is the column-bars
	# converted_annotation is a column matrix, N*ncolBars, still using the grid.rect
	if(!is.na(annotation[[1]][1])){
		pushViewport(vplayout(3, 2))
		converted_annotation = convert_annotations(annotation, annotation_colors)
		draw_annotations(converted_annotation, border_color)
		upViewport()
	}	
	# Draw annotation labels at (3, 1) 
	# needs to compute position based on its length if were to align to the right
	# quick and dirty: center
	##
	# to: (1) calculate grid.text fontsize; (2) align to the right
	##
	#browser()
	###
	# there is nothing happening in vplayout(3, 1) from pheatmap: so this section is added by me
	if(annoHeaderPosition=='left'){
		if(!is.na(annotation[[1]][1])){
		pushViewport(vplayout(3, 1))
		nTT = ncol(annotation) # number of col bars
		yTT = cumsum(rep(8, nTT)) - 4 + cumsum(rep(2, nTT)) ## this is further transformed into big points. 
		barNames <- colnames(annotation)
		# may need bigpts to fontsize: 8 bigpts per row
		grid.text(barNames, x=0.5, y=unit(yTT[nTT:1], "bigpts"), just='centre') ## this is the order from last col to first column
		upViewport()
		}	
	} else {
		# add to the right: c(3, 3)
		if(!is.na(annotation[[1]][1])){
		pushViewport(vplayout(3, 3))
		nTT = ncol(annotation) # number of col bars
		yTT = cumsum(rep(8, nTT)) - 4 + cumsum(rep(2, nTT)) ## this is further transformed into big points. 
		barNames <- colnames(annotation)
		# may need bigpts to fontsize: 8 bigpts per row
		grid.text(barNames, x=0.5, y=unit(yTT[nTT:1], "bigpts"), just='centre') ## this is the order from last col to first column
		upViewport()
		}	
	}	
	# Draw annotation legend
	if(!is.na(annotation[[1]][1]) & annotation_legend){
		if(length(rown) != 0){
			#pushViewport(vplayout(4:5, 5))
			pushViewport(vplayout(2:5, 5)) # starts at the very top to save space
		}
		else{
			#pushViewport(vplayout(3:5, 5))
			pushViewport(vplayout(2:5, 5))
		}
		draw_annotation_legend(annotation, annotation_colors, border_color, fontsize = fontsize, ...)
		upViewport()
	}	
	# Draw legend
	if(!is.na(legend[1])){
		### need this??? potentially mistake: length(colnames(matrix))
		if(length(rown) != 0){
			pushViewport(vplayout(4:5, 4))
		}
		else{
			pushViewport(vplayout(3:5, 4))
		}
		draw_legend(color, breaks, legend, fontsize = fontsize, ...)
		upViewport()
	}	
}

##' when center=TRUE, the broken bars at 0 has equal length of value and color
##' This is the old generate_breaks function which does not allow for midpoint
generate_breaks0 = function(x, n, center = F){
	if(center){
		m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
		res = seq(-m, m, length.out = n + 1)
	}
	else{
		res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
	}
	
	return(res)
}

#' generate breaks for specified start, end position and n, posibly with a midpoint
#'
#' when midpoint is not explicitely specified, use the middle point of st, end (the mean)
#'
#' @param st start position
#' @param ed end position
#' @param midpoint midpoint
#' @param n number of intervals to form. This actually needs n+1 breaks
#' @return a vector of length n+1 for breaks
#' @export
genbreaks <- function(st, ed, midpoint=NA, n){
	n1 <- floor(n/2); n2 <- n-n1
	# when midpoint is not explicitely specified, use the middle point of st, end (the mean)
	if(is.na(midpoint)) midpoint <- mean(c(st, ed))
	if(midpoint<st | midpoint>ed) 
		stop(sprintf('midpoint is not proper: start=%f, end=%f, midpoint=%f', st, ed, midpoint))
	breaks <- c(seq(st, midpoint, length=n1+1), seq(midpoint, ed, length=n2+1)[-1])
}
#' generate breaks for an input matrix or vector 
#'
#' this is an wrapper for genbreaks which is more straightforward
#'
#' enhenced with midpoint
#' @param x a matrix or vector
#' @param n desired number of intervals, which corresponds to n+1 breaks
#' @param midpoint midpoint. default is NA, which will be recomputed as the middle point of the value ladder (mean of c(min, max)). 
#'   This is equivalent to no midpoint at all, and hence the default way
#' @param center whether to generate breaks centering 0 by expanding the values and hence, waste some colors (not recommended if the data is skewed)
#' @return a vector of break points
#' @export
generate_breaks = function(x, n, midpoint=NA, center = F){
	if(center){
		m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
		#res = seq(-m, m, length.out = n + 1)
		res = genbreaks(st=-m, ed=m, midpoint=midpoint, n=n)
	}
	else{
		#res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
		res = genbreaks(st=min(x, na.rm = T), ed=max(x, na.rm = T), midpoint=midpoint, n=n)
	}
	
	return(res)
}

scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
	return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

scale_colours = function(mat, col = rainbow(10), breaks = NA){
	mat = as.matrix(mat)
	return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

#' rcorr to compute correlation coefficient among columns.  
#'
#' Missing values are deleted in pairs rather than deleting all rows of x having any missing variables.
#'
#' rcorr requires at least 5 rows. So when nrow<5, use cor function.
#' @param mat a matrix
#' @param method "pearson","spearman"
#' @return correlation matrix
#' @export
rcor <- function(mat, method='pearson'){
	if(nrow(mat)>=5){
		res <- rcorr(x=mat, type=method)$r
	} else {
		res <- cor(mat, method=method)
	}
	res
}
# try to modify based on Kevin Coombes distanceMatrix() function can deal with NA 
# it seems this function has been updated so that by default through using the use parameter, missing will be solved by pairwise complete deletion (less conservative than complete.obs across all cases)
# copied here to be independent from oompa
distanceMatrix2 <- function(dataset, metric, use='pairwise.complete.obs', ...) {
  if(inherits(dataset, "ExpressionSet")) {
    dataset <- exprs(dataset)
  }
  METRICS <- c('pearson', 'sqrt pearson', 'spearman', 'weird',
               'absolute pearson', 'uncentered correlation')
  m <- pmatch(metric, METRICS, nomatch=NA)
  if (is.na(m)) {
    distance <- dist(t(dataset), method=metric, ...)
  } else {
    metric <- METRICS[m]
    uncent <- function(dataset) {
      temp <- sqrt(apply(dataset^2, 2, mean))
      dataset <- sweep(dataset, 2, temp, '/')
      temp <- t(dataset) %*% dataset/nrow(dataset)
      as.dist(1 - temp)
    }
#browser()
    distance <- switch(metric,
                       pearson = as.dist((1-cor(dataset, use=use))/2),
                       'sqrt pearson' = as.dist(sqrt(1-cor(dataset, use=use))),
                       weird = dist(cor(dataset, use=use)),
                       'absolute pearson' = as.dist((1-abs(cor(dataset, use=use)))),
                       spearman = as.dist((1-cor(apply(dataset, 2, rank), use=use))/2),
                       'uncentered correlation' = uncent(dataset))
  }
  distance
}
#' hierarchical clustering with rows using given distance and linkage
#'
#' @param mat input matrix
#' @param distance distance metric
#' @param method linkage method
#' @return hclust object
#' @export
cluster_mat = function(mat, distance, method){
	distance_coombes <- c('pearson', 'sqrt pearson', 'spearman', 'weird',
               'absolute pearson', 'uncentered correlation')
	distance_pheatmap <- c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
	distance_all <- c(distance_coombes, distance_pheatmap)
	if(!(method %in% c("ward", "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))){
		stop("clustering method has to one form the list: 'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
	}
	if(!(distance[1] %in% distance_all) & class(distance) != "dist"){
		print(!(distance[1] %in% c('pearson', 'sqrt pearson', 'spearman', 'weird',
               'absolute pearson', 'uncentered correlation', "correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) | class(distance) != "dist")
		stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'pearson', 'sqrt pearson', 'spearman', 'weird',
               'absolute pearson', 'uncentered correlation', 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
	}
	if(method=='ward') {
		method <- 'ward.D' # this is the old ward linkage, not exact implementation of ward; doing this will remove the warning
	}
	if(distance[1] == "correlation"){
		#d = as.dist(1 - cor(t(mat)))
		#cor0 <- (1 - cor(t(mat)))/2
		#cor1 <- (1 - rcor(t(mat)))/2
		#d = as.dist((1 - cor(t(mat)))/2) # oompa compliant
		d = as.dist((1 - rcor(t(mat)))/2) # to deal with NA: pair wise deletion; added on 2014/02/26
		#browser()
	} else if(distance[1] %in% distance_coombes){
		# implements Coombes pearson approach dealing with NA
		#d = distanceMatrix2(mat, distance[1]) 
		d = distanceMatrix2(t(mat), distance[1]) # notice this is clustering on rows!!!
	} else{
		if(class(distance) == "dist"){
			d = distance
		}
		else{
			# "Missing values are allowed, and are excluded from all computations involving the rows within which they occur."
			# this means for rows with NA, the distance is computed from other columns without NA---> seems not the case
			d = dist(mat, method = distance)
		}
	}
	#browser()
	#return(stats::hclust(d, method = method))
	return(fastcluster::hclust(d, method = method)) # faster
}

scale_rows = function(x){
	m = apply(x, 1, mean, na.rm = T)
	s = apply(x, 1, sd, na.rm = T)
	return((x - m) / s)
}
#' scale a matrix by row, by column or none
#'
#' notice both center and scale are simultaneously activated
#' @param mat a matrix
#' @param scale selection from c("none", "row", "column")
#' @return scaled matrix
#' @export 
scale_mat = function(mat, scale){
	if(!(scale %in% c("none", "row", "column"))){
		stop("scale argument shoud take values: 'none', 'row' or 'column'")
	}
	mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
	return(mat)
}

generate_annotation_colours = function(annotation, annotation_colors, drop){
	if(is.na(annotation_colors)[[1]][1]){
		annotation_colors = list()
	}
	count = 0
	for(i in 1:ncol(annotation)){
		if(is.character(annotation[, i]) | is.factor(annotation[, i])){
			if (is.factor(annotation[, i]) & !drop){
				count = count + length(levels(annotation[, i]))
			}
			else{
				count = count + length(unique(annotation[, i]))
			}
		}
	}
	
	factor_colors = hsv((seq(0, 1, length.out = count + 1)[-1] + 
      0.2)%%1, 0.7, 0.95)
	
	set.seed(3453)
	
	for(i in 1:ncol(annotation)){
		if(!(colnames(annotation)[i] %in% names(annotation_colors))){
			if(is.character(annotation[, i]) | is.factor(annotation[, i])){
				n = length(unique(annotation[, i]))
				if (is.factor(annotation[, i]) & !drop){
					n = length(levels(annotation[, i]))
				}
				ind = sample(1:length(factor_colors), n)
				annotation_colors[[colnames(annotation)[i]]] = factor_colors[ind]
				l = levels(as.factor(annotation[, i]))
				l = l[l %in% unique(annotation[, i])]
				if (is.factor(annotation[, i]) & !drop){
					l = levels(annotation[, i])
				}
				names(annotation_colors[[colnames(annotation)[i]]]) = l
				factor_colors = factor_colors[-ind]
			}
			else{
				r = runif(1)
				annotation_colors[[colnames(annotation)[i]]] = hsv(r, c(0.1, 1), 1)
			}
		}
	}
	return(annotation_colors)
}

kmeans_aheatmap = function(mat, k = min(nrow(mat), 150), sd_limit = NA, ...){
	# Filter data
	if(!is.na(sd_limit)){
		s = apply(mat, 1, sd)
		mat = mat[s > sd_limit, ]	
	}
	
	# Cluster data
	set.seed(1245678)
	km = kmeans(mat, k, iter.max = 100)
	mat2 = km$centers
	
	# Compose rownames
	t = table(km$cluster)
	rownames(mat2) = sprintf("cl%s_size_%d", names(t), t)
	
	# Draw heatmap
	aheatmap(mat2, ...)
}

 
#' A function to draw clustered heatmaps.
#' 
#' A function to draw clustered heatmaps where one has better control over some graphical 
#' parameters such as cell size, etc. 
#' 
#' Default distance is correlation for both row and column; linkage is set to ward by default. Scaling is applied to rows by default.
#' Modification has been made so that scaling only applies how the matrix is colored; the old matrix would be used to do the clustering.
#' For the annotation data frame to show as a column bar, categorical variables need to be a factor; the function can calculate default colors
#' But do not forget to have rownames of annotation data frame set to colnames of the matrix!
#' 
#' Added features: 
#' (1) truncate=NA, q1=0.01, q2=0.99, Lower, Upper: this is to truncate the valus to save color space. if truncate=FALSE, no truncation; if NA, it will be internally
#' set to TRUE if scale!='none'
#' (2) ...
#' 
#' The function also allows to aggregate the rows using kmeans clustering. This is 
#' advisable if number of rows is so big that R cannot handle their hierarchical 
#' clustering anymore, roughly more than 1000. Instead of showing all the rows 
#' separately one can cluster the rows in advance and show only the cluster centers. 
#' The number of clusters can be tuned with parameter kmeans_k.
#'
#' @param mat numeric matrix of the values to be plotted.
#' @param color vector of colors used in heatmap.
#' @param kmeans_k the number of kmeans clusters to make, if we want to agggregate the 
#' rows before drawing heatmap. If NA then the rows are not aggregated.
#' @param breaks a sequence of numbers that covers the range of values in mat and is one 
#' element longer than color vector. Used for mapping values to colors. Useful, if needed 
#' to map certain values to certain colors, to certain values. If value is NA then the 
#' breaks are calculated automatically.
#' @param midpoint whether explicitely to match a mid point of the values to a mid-colour (usually white color) for the main image. 
#' This is useful when one what to hide result around i,e. 0 which
#' is not an interesting result.  This value is based on the matrix inputted for color matrix. Therefore, if scaled, it is on the scaled scale. 
#' @param border_color color of cell borders on heatmap, use NA if no border should be 
#' drawn.
#' @param cellwidth individual cell width in points. If left as NA, then the values 
#' depend on the size of plotting window.
#' @param cellheight individual cell height in points. If left as NA, 
#' then the values depend on the size of plotting window.
#' @param shrink sometimes no space is left for margin; set shrink (percentage of current figure) will the image
#' @param scale character indicating if the values should be centered and scaled (first subtracts mean and then scaled by sd) in 
#' either the row direction or the column direction, or none. Corresponding values are 
#' \code{"row"}, \code{"column"} and \code{"none"}; By default, scale by row is applied which is mostly done in microarray data (scale each gene)
#' @param clusterWithScaledData this parameter controls whether the scaling of input matrix would affect clustering. By default this is set to be FALSE, which means
#' clustering only affects the color, not affecting the dendrogram and the numbers showing; if set to be TRUE, the scaled data would be used for clustering and coloring (as well
#' as the numbers shown) if necessary. It is recommended to set this as FALSE since scaling would attenuate the difference among samples and hence hide important pattern.
#' @param cluster_rows boolean values determining if rows should be clustered,
#' @param cluster_cols boolean values determining if columns should be clustered.
#' @param clustering_distance_rows distance measure used in clustering rows. Possible 
#' values are \code{"correlation"} for Pearson correlation and all the distances 
#' supported by \code{\link{dist}}, such as \code{"euclidean"}, etc. If the value is none 
#' of the above it is assumed that a distance matrix is provided.
#' @param clustering_distance_cols distance measure used in clustering columns. Possible 
#' values the same as for clustering_distance_rows.
#' @param clustering_method clustering method used. Accepts the same values as 
#' \code{\link{hclust}}.
#' @param treeheight_row the height of a tree for rows, if these are clustered. 
#' Default value 50 points.
#' @param treeheight_col the height of a tree for columns, if these are clustered. 
#' Default value 50 points.
#' @param legend logical to determine if legend should be drawn or not.
#' @param legend_breaks vector of breakpoints for the legend.
#' @param legend_labels vector of labels for the \code{legend_breaks}.
#' @param annotation data frame that specifies the annotations shown on top of the 
#' columns. Each row defines the features for a specific column. The columns in the data 
#' and rows in the annotation are matched using corresponding row and column names. Note 
#' that color schemes takes into account if variable is continuous or discrete.
#' @param annotation_colors list for specifying annotation track colors manually. It is 
#' possible to define the colors for only some of the features. Check examples for 
#' details. For categorical column bar, define a named vector so that color -category
#' will be correctly mapped; for continuous column bar, the specified color vector will be used
#' to interpolate colors for values inbetween.
#' @param annoHeaderPosition if left, which is the old one that adds annotation bar names at the left of the bars (c(3, 1)); now we add right so that the barname shows to the right. This is needed
#'       when no row dendrogram where he header will not show completely.
#' @param annotation_legend boolean value showing if the legend for annotation tracks 
#' should be drawn. 
#' @param drop_levels logical to determine if unused levels are also shown in the legend
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @param main the title of the plot. set to "" so that there will be one line blank on top for good visual
#' @param fontsize base fontsize for the plot 
#' @param fontsize_row fontsize for rownames (Default: fontsize) 
#' @param fontsize_col fontsize for colnames (Default: fontsize) 
#' @param display_numbers logical determining if the numeric values are also printed to 
#' the cells. 
#' @param number_format format strings (C printf style) of the numbers shown in cells. 
#' @param number_rotation rotation for number text display. default is 0
#' For example "\code{\%.2f}" shows 2 decimal places and "\code{\%.1e}" shows exponential 
#' notation (see more in \code{\link{sprintf}}).    
#' @param fontsize_number fontsize of the numbers displayed in cells
#' @param filename file path where to save the picture. Filetype is decided by 
#' the extension in the path. Currently following formats are supported: png, pdf, tiff,
#'  bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is 
#' calculated so that the plot would fit there, unless specified otherwise.
#' @param width manual option for determining the output file width in inches.
#' @param cexRow scale rownames by this factor.
#' @param cexCol scale colnames by this factor.
#' @param labRow labRow as in heatmap.2 Notice this is specified as the input matrix; if ordered, i.e. by clustering, the names should be the original order since it will be internally updated.
#' @param labCol labCol as in heatmap.2 Notice this is specified as the input matrix; if ordered, i.e. by clustering, the names should be the original order since it will be internally updated.
#' @param height manual option for determining the output file height in inches.
#' @param returnTree whether to return the clustering tree back. Default is FALSE
#' @param truncate 
#' @param Lower parameter Lower to truncByQuantile()
#' @param Upper parameter Upper to truncByQuantile()
#' @param q1 parameter q1 to truncByLimit
#' @param q2 parameter q2 to truncByLimit
#' @param \dots graphical parameters for the text used in plot. Parameters passed to 
#' \code{\link{grid.text}}, see \code{\link{gpar}}. 
#' 
#' @return return the following list if returnTree is TRUE
#' \itemize{
#' 	\item \code{tree_row} the clustering of rows as \code{\link{hclust}} object 
#' 	\item \code{tree_col} the clustering of columns as \code{\link{hclust}} object
#' 	\item \code{kmeans} the kmeans clustering of rows if parameter \code{kmeans_k} was 
#' specified 
#' }
#' 
#' @author  Pan Tong \email{nickytong@@gmail.com} 
#' @export
aheatmap = function(mat, color = bluered(100), kmeans_k = NA, breaks = NA, midpoint=NA, 
	border_color = "grey60", cellwidth = NA, cellheight = NA, shrink=0.92, scale = "row", clusterWithScaledData=FALSE, cluster_rows = TRUE, cluster_cols = TRUE, 
	clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", clustering_method = "ward",  
	treeheight_row = ifelse(cluster_rows, 50, 0), treeheight_col = ifelse(cluster_cols, 50, 0), legend = TRUE, legend_breaks = NA, 
	legend_labels = NA, annotation = NA, annotation_colors = NA, annoHeaderPosition=NA, annotation_legend = TRUE, drop_levels = TRUE, show_rownames = T, 
	show_colnames = T, main = "", fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = F, number_format = "%.2f", number_rotation=0, 
	fontsize_number = 0.8 * fontsize, filename = NA, width = NA, height = NA, cexRow, cexCol, labRow=NA, labCol=NA, 
	truncate=NA, q1=0.01, q2=0.99, Lower, Upper, 
	returnTree=FALSE, ...){
	rown <- getName(mat=mat, name=labRow, type='row')
	coln <- getName(mat=mat, name=labCol, type='column')
	if(is.na(truncate)){
		truncate <- ifelse(scale=='none', FALSE, TRUE)
	}
	breaksSpecified <- ifelse(is.na(breaks)[1], FALSE, TRUE) # add this indicator to track if user has specified breaks
	### breaks: (1) generated by default for scaled data
	# Preprocess matrix
	matRaw = as.matrix(mat) # matRaw is the input matrix which is desired to do clustering and stuff like kmeans aggregation
	# note that the matrix is scaled at the very begining; this means the clustering would also be affected!
	if(scale != "none"){
		mat = scale_mat(matRaw, scale) # scaled mat
		if(!breaksSpecified){
			# it is a good specification for center=TRUE since this is scaled data. 
			breaks = generate_breaks(mat, length(color), midpoint=midpoint, center = T) # breaks based on the scaled matrix, ok
		}
	}
	####
	# we use mat to track the matrix for display (used to calculate breaks the color matrix); we use matRaw to track the matrix for clustering
	####
	# Kmeans: input is the raw matrix; the output mat which would be used to do clustering later needs to be scaled when necessary
	if(!is.na(kmeans_k)){
		# Cluster data
		km = kmeans(matRaw, kmeans_k, iter.max = 100)
		matCenters = km$centers
		matRaw <- matCenters # this will be used for clustering the aggregated data
		# needs to scale the data after aggregation for display
		if(scale != "none"){
			mat = scale_mat(matCenters, scale) # scaled mat
			if(is.na(breaks[1])){ # since breaks is now a vector computed before; we need to update it
				breaks = generate_breaks(mat, length(color), midpoint=midpoint, center = T) # breaks based on the scaled matrix for aggregated data, ok
			}
		}
		# Compose rownames
		t = table(km$cluster)
		rownames(mat) = sprintf("cl%s_size_%d", names(t), t)
	}
	else{
		km = NA
	}
	### deal with clustering: either use scaled matrix for clustering or use raw matrix, depending on clusterWithScaledData
	# now mat is the scaled mat; matRaw is the raw mat (might after aggregation)
	# but if the user request to cluster with scaled data, we are still not too late to do it! just set the matRaw as mat
	if(clusterWithScaledData){
		matRaw_backup <- matRaw # we still preserve a copy of the raw mat
		matRaw <- mat
	}
	# Do clustering: tree_row is returned by  return(hclust(d, method = method))
	if(cluster_rows){
		tree_row = cluster_mat(matRaw, distance = clustering_distance_rows, method = clustering_method)
		mat = mat[tree_row$order, , drop = FALSE] # reorder matrix based on clustering; no drop!!! ---> in pheatmap, this means the color matrix is also ordered
		matRaw = matRaw[tree_row$order, , drop = FALSE] # this is needed so that when show raw values, they are ordered accordingly (2014/02/18)
		rown <- rown[tree_row$order] # reorder the rown correspondingly
	}
	else{
		tree_row = NA
		treeheight_row = 0
	}
	
	if(cluster_cols){
		tree_col = cluster_mat(t(matRaw), distance = clustering_distance_cols, method = clustering_method)
		mat = mat[, tree_col$order, drop = FALSE]
		matRaw = matRaw[, tree_col$order, drop = FALSE] ## fix a bug so that digit numbers are ordered accordingly (2014/02/18)
		coln <- coln[tree_col$order]
	}
	else{
		tree_col = NA
		treeheight_col = 0
	}
	###### 
	# Format numbers to be displayed in cells: the number is the raw matrix; scaled matrix is meaningless!
	if(display_numbers){
		#fmat = matrix(sprintf(number_format, matRaw), nrow = nrow(mat), ncol = ncol(mat)) # there is a bug here: when row/column ordered, the digits are not ordered!
		fmat = matrix(sprintf(number_format, matRaw), nrow = nrow(mat), ncol = ncol(mat)) # we have also ordered matRaw accordingly
		attr(fmat, "draw") = TRUE
	}
	else{
		fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
		attr(fmat, "draw") = FALSE
	}
	
	
	# Colors and scales
	if(!is.na(legend_breaks[1]) & !is.na(legend_labels[1])){
		if(length(legend_breaks) != length(legend_labels)){
			stop("Lengths of legend_breaks and legend_labels must be the same")
		}
	}
	#### in some cases, the mat ought to be truncated. Notice the truncation does not affect clustering
	#### but we need to truncate the data before generating breaks and legend labels and color mat
	if(truncate){
		if(!missing(Lower) & !missing(Upper)){
			mat <- truncByLimit(mat, Lower=Lower, Upper=Upper)
		} else {
			mat <- truncByQuantile(mat, q1=q1, q2=q2)
		}
		# also update the breaks if data is truncated
		if(!breaksSpecified){
		breaks = generate_breaks(mat, length(color), midpoint=midpoint, center = T) # breaks based on the scaled matrix, ok
		}
	}
	### breaks: (2) generated by default for scale='none'; notice if scale!='none', breaks must have been generated and not NA, hence this code is not executed for scaled mat
	if(is.na(breaks[1])){
			breaks = generate_breaks(as.vector(mat), length(color), midpoint=midpoint, center=F) # do not center since data is not scaled so as not to waste color		
	}
  if (legend & is.na(legend_breaks[1])) {
      legend = grid.pretty(range(as.vector(breaks)))
			names(legend) = legend
  }
	else if(legend & !is.na(legend_breaks[1])){
		legend = legend_breaks[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
		
		if(!is.na(legend_labels[1])){
			legend_labels = legend_labels[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
			names(legend) = legend_labels
		}
		else{
			names(legend) = legend
		}
	}
  else {
      legend = NA
  }
#browser()
	################# scales mat to a color matrix: from hereafter, mat is a colormat!!!! ###############
	mat = scale_colours(mat, col = color, breaks = breaks)
	
	# Preparing annotation colors
	if(!is.na(annotation[[1]][1])){
		annotation = annotation[colnames(mat), , drop = F]
		annotation_colors = generate_annotation_colours(annotation, annotation_colors, drop = drop_levels)
	}
	
	if(!show_rownames){
		#rown <- NULL # this is not working; a bug
		# now I'm using cheating: overwrite fontsize so that it is not shown
		fontsize_row <- 0
	}
	
	if(!show_colnames){
		#coln <- NULL
		fontsize_col <- 0
	}
	if(missing(cexCol))
			cexCol <- guessCEX(ncol(mat))
	if(missing(cexRow))
			cexRow <- guessCEX(nrow(mat))
	# Draw heatmap
	# so from hereafter, cexCol and cexRow will never be missing or NA
	# passing a color matrix mat for imaging
	heatmap_motor(mat, labRow=rown, labCol=coln, annoHeaderPosition=annoHeaderPosition, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight, shrink=shrink, treeheight_col = treeheight_col, treeheight_row = treeheight_row, tree_col = tree_col, tree_row = tree_row, filename = filename, width = width, height = height, breaks = breaks, color = color, legend = legend, annotation = annotation, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, cexRow=cexRow, cexCol=cexCol, rot=number_rotation, ...)
	res <- list(tree_row = tree_row, tree_col = tree_col, kmeans = km)
	if(returnTree){
		return(res)
	}
}
#set.seed(100)
#test = matrix(rnorm(200), 20, 10)
# test[1:3, 1:3] <- -10
#aheatmap(test, scale='none', clustering_distance_rows='euclidean', clustering_distance_cols='euclidean', clustering_method='average', display_numbers=TRUE)
#aheatmap(test, scale='none', clustering_distance_rows='euclidean', midpoint=0, clustering_distance_cols='euclidean', clustering_method='average', display_numbers=TRUE)


#' Prepare annotation data frame so that aheatmap can show it
#'
#' aheatmap requires the annotation data frame to have rownames identical to colnames of the heatmap matrix; otherwise, nothing is drawn
#'
#' the main goal for this function is to (1) add rownames (2) convert categorical variables to factor so that aheatmap can draw it;
#' rownames can be passed to set rownames; of no rownames, aheatmap will not draw the column bar
#' the conversion from categorical to factor is achieved using formatVec() as the engine.
#' 
#' @param dat annotation data frame, one row per sample
#' @param sampleNames the colnames of the matrix in heatmap
#' @param minUV minimum unique value to not be a factor
#' @return a formated data frame
#' @export
formatAnnotation4aheatmap <- function(dat, sampleNames, minUV=12){	
	res <- colwise(formatVec)(dat)
	if(!missing(sampleNames)) rownames(res) <- sampleNames
	res
}


#' format a vector by forcing into a factor if not numeric or unique value less than some cutoff
#' @param x a vector
#' @param minUV minimum unique value to not be a factor
#' @return a vector formated
#' @export
formatVec <- function(x, minUV=12) {
	ux <- sort(unique(x), decreasing=FALSE) # alphabet for character
	#browser()
	if(class(x)[1]!="numeric" & length(ux)<minUV){
		res <- factor(x, levels=ux)
	} else {
		res <- x
	}
	res	
}
#' format a data frame using formatVec
#' @param dat a data frame
#' @param minUV minimum unique value to not be a factor
#' @return a data frame
#' @export
formatDF <- function(dat, minUV=12){
	#apply(2, dat, function(x) formatVec(x, minUV))
	tt <- colwise(formatVec, minUV=minUV)(dat)
	rownames(tt) <- rownames(dat)
	tt
}
#' index of columns of a data frame is continuous (numeric)
#' @param dat a data frame or a vector
#' @param minUV passed to formatVec
#' @return logical
#' @export
isCont <- function(dat, minUV=12){
	if(!is.vector(dat)){
	## is a data frame
	if(length(minUV)==1){
	res <- apply(dat, 2, function(x) is.numeric(formatVec(x, minUV)))
	} else {
	res <- apply(dat, 2, function(x, y) is.numeric(formatVec(x, y)), minUV)
	}
	} else {
	res <- is.numeric(formatVec(dat, minUV))
	}
	res
}
#' which columns of a data frame is continuous (numeric)
#' @param ... parameters passed to isCont
#' @return a vector of index
#' @export
whichCont <- function(...){
	which(isCont(...))
}
#' indicator if columns are categorical
#' @export
isCat <- function(...) !isCont(...)
#' which columns are categorical
#' @export
whichCat <- function(...) which(isCat(...))
#whichCont(phenoGanglia)


#' plot matrix values
#'
#' We may need to plot correlation, expression matrix and so on which is not large. We thus can show the actual value. No clustering is needed.
#'
#' @param mat a matrix to be plotted
#' @param color color pallete ( a vector of colors)
#' @param display_numbers whether to display numbers
#' @param scale ways to scale the values
#' @param ... additional parameters to be passed to aheatmap (pheatmap)
#' @export 
aMatrix <- function(mat, color=bluered(50), display_numbers=TRUE, scale='none', isint=FALSE, ...){
	#browser()
	if(isint) {
		number_format = "%d"
	}  else {
		number_format = "%.2f"
	}
	aheatmap(mat, scale=scale, display_numbers=display_numbers, color = color, cluster_rows=FALSE, cluster_cols=FALSE, number_format=number_format, ...) # cluster_rows=FALSE, cluster_cols=FALSE, 
}
  #aMatrix(cor_imTs, fontsize_number=12)

#aMatrix(cor(dfRBM10), color=redscale(30), fontsize_number=12)

#' 2-step procedure to decide ordering
#'
#' when order is specified, use this order; otherwise, try to use clustering; if no clustering is explicitely specified, do not do any ordering
#'
#' @param ord order vector
#' @param mat a matrix
#' @param distance clustering distance
#' @param direction when specified as row, do clustering on rows and calculate row ordering; when specified as column, do the same thing but with trnasformed mat
#' @return a permutation
decideOrder <- function(ord, mat, distance, method, direction='row'){
	if(direction=='column') {
		mat <- t(mat) # when performed on column, transform it
	}
	#browser()
	if(is.null(ord)){
		# when Rowv not specified, try to use clustering
		if(!missing(distance)){	
		tree <- cluster_mat(mat, distance = distance, method = method)
		oo <- tree$order
		} else {
		# when even no clustering is known, no ordering performed
		oo <- 1:nrow(mat)	
		}
	} else {
	# when Rowv is specified, this is the primary order
	if(length(ord)!=nrow(mat)) stop('Length of order integers is not compatible!\n')
	oo <- ord
	}
	oo	
}
	
#' plot an ordered (sorted) matrix
#'
#' this is an enhancement of aMatrix where one can specify Rowv and Colv as in heatmap.2
#' typical application is that one can need to order according to a phenotype, i.e. EMT score
#'
#' The working horce is aMatrix (more precisely its code). We only need to computer row order and column order and pass the reordered matrix to aMatrix.
#' Ordering is determined first by looking at Rowv and Colv. If either one is not specified, order is computed based on clustering with a given distance and linkage rule. 
#' If even clustering is not specified, then there is no ordering.
#'
#' @param mat a matrix to be plotted
#' @param color color pallete ( a vector of colors)
#' @param display_numbers whether to display numbers
#' @param scale ways to scale the values
#' @param Rowv either NULL or a vector of integers for row ordering
#' @param Colv either NULL or a vector of integers for column ordering
#' @param clustering_distance_rows distance measure used in clustering rows. Possible 
#' values are \code{"correlation"} for Pearson correlation and all the distances 
#' supported by \code{\link{dist}}, such as \code{"euclidean"}, etc. If the value is none 
#' of the above it is assumed that a distance matrix is provided.
#'
#' @param clustering_distance_cols distance measure used in clustering columns. Possible 
#' values the same as for clustering_distance_rows.
#'
#' @param clustering_method clustering method used. Accepts the same values as 
#' \code{\link{hclust}}.
#'
#' @param annotation data frame that specifies the annotations shown on top of the 
#' columns. Each row defines the features for a specific column. The columns in the data 
#' and rows in the annotation are matched using corresponding row and column names. Note 
#' that color schemes takes into account if variable is continuous or discrete.
#'
#' @param return whether to return a list of orders
#' @param ... additional parameters to be passed to aheatmap (pheatmap)
#' @export 
aoMatrix <- function(mat, color=bluered(50), display_numbers=TRUE, scale='none', Rowv=NULL, Colv=NULL, clustering_distance_rows, 
		clustering_distance_cols, clustering_method='ward.D2', return=FALSE, 
		annotation = NA, annotation_colors = NA, labRow=NA, labCol=NA, ...){
	rown <- getName(mat=mat, name=labRow, type='row')
	coln <- getName(mat=mat, name=labCol, type='column')
	#browser()
	# below we use a 2-step procedure to decide row and column order
	rowO <- decideOrder(Rowv, mat, clustering_distance_rows, clustering_method, direction='row')
	colO <- decideOrder(Colv, mat, clustering_distance_cols, clustering_method, direction='column')
	# matO <- mat; annotationO <- annotation
	#browser()
	if(!is.na(annotation)[1]){ # this requires annotation to be directly observable!!! hence need to have this parameter copied!
		annotation <- annotation[colO, , drop=FALSE] # annotation needs to be ordered with colO
	}
	mat <- mat[rowO, colO]
	rown <- rown[rowO]
	coln <- coln[colO]
	#browser()
	aheatmap(mat, scale=scale, color=color, display_numbers=display_numbers, cluster_rows=FALSE, cluster_cols=FALSE, 
		annotation=annotation, annotation_colors=annotation_colors, labRow=rown, labCol=coln, ...)
	if(return){
		return(list(rowO=rowO, colO=colO))
	}
}
if(FALSE){
aoMatrix(mat, annotation=pheno_f, annotation_colors = annoColors, display_numbers=FALSE, scale='row', Rowv=OL_1$rowO, Colv=order(pheno$EMTscore, decreasing=FALSE), 
	color = exprColFn(31), border_color=NA, 
		cexCol=0.8, cexRow=0.8, main=tumor, labRow=rownames(cor_imTs))
}
#aoMatrix(mat=test, clustering_distance_rows='correlation', clustering_distance_cols='correlation', clustering_method='ward')

########################################################################################################################################
# BUM
########################################################################################################################################

#' Create a bum 'object'
#'
#' Given a vector of p values, fit BUM and calculate FDR table 
#' 
#' Where the pval is highly significant, FDR tends to be small and the proportion of H0 genes is also small (pi0). It turns out the estimate of pi0
#' is still consistant (very similar to the qvalue package) and thus the Bum model fits really well even in this case. The problem arises when we want to achieve
#' a FDR larger than that is possible. i.e. maximum FDR=0.2 but specified FDR=0.3, in which case, the  cutoffSignificant() function will return a value larger than 1.
#' This is not a problem as long as we truncate this cutoff to 1 (notice that even using qvalue package will not achieve FDR=0.3). Finally, Bum model is very consistent with
#' the qvalue package. 
#'
#' @param pval pvalue vector
#' @param alphas alpha values to create FDR table
#' @param FDR FDR value to be used for gene selection; if NA, guess it from FDR table
#' @param pcutoff cutoff of p value to select indeces based on p values; default is 0.05
#' @param maskPthr when the P value contains a lot of 1's or close to 1, BUM model might fails. Thus we can mask out these P values
#' as NA so that we get a good estimate and selection for significant events. Default is 1, which means no masking is done. 
#' @return a list
#' @export
create_bum <- function(pval, alphas=c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5), FDR=NULL, maskPthr=1, pcutoff=0.05) {
	if(maskPthr<1) warning(sprintf('P<%.3f are masked as NA: %d P values masked!', maskPthr, sum(pval>maskPthr, na.rm=TRUE)))
	pval[which(pval > maskPthr)] <- NA
	#browser()
	bum <- Bum(pval)
	nSig <- sapply(alphas, function(x) countSignificant(bum, alpha=x, by = "FDR"))
	pCut <- sapply(alphas, function(x) cutoffSignificant(bum, alpha=x, by = "FDR"))
	if(any(pCut>1)) {
		warning('cutoffSignificant returns p value >1; truncated at 1!')
		pCut[pCut>1] <- 1
	}
	tabFDR <- data.frame(FDR=alphas, Selected=nSig, Pvalue=pCut)
	# FDR based on qvalue
	require(qvalue)
	# if p value rang is not [0, 1], give error. so nned to specify lambda as the range of p value
	# https://support.bioconductor.org/p/74637/
	#browser()
	lambda_default <- seq(0.05, 0.95, 0.05)
	lambda_max <- max(noNA(pval))
	qobj <- try(qvalue(p = noNA(pval), lambda = lambda_default[lambda_default<=lambda_max]), silent=TRUE)
	#browser()
	while(class(qobj)=='try-error'){
		lambda_max <- 0.9 * lambda_max
		cat(sprintf('--->lambda_max iteration: %.2f\n', lambda_max))
		qobj <- try(qvalue(p = noNA(pval), lambda = lambda_default[lambda_default<=lambda_max]), silent=TRUE)
		if(lambda_max<0.05) {
			warning('Qvalue estimation failed!\n')
			break
		}
	}
	if(class(qobj)!='try-error'){
		qvalues <- qobj$qvalues
		dfq <- data.frame(p=qobj$pvalues, q=qobj$qvalues)
		dfq <- dfq[order(dfq$p, decreasing=FALSE), ]
		# none-smooth estimation
		tabQ <- foreach(x=alphas, .combine='rbind') %do% {
			tp <- subset(dfq, q<=x)
			c(Qvalue=x, Selected=nrow(tp), Pvalue=max(c(tp$p, 0)))
		}
		tabQ <- data.frame(tabQ)
		rownames(tabQ) <- NULL
	# table of selection by P
	tabP <- foreach(x=c(0.001, 0.01, 0.05, 0.1, 0.2), .combine='rbind') %do% {
		tp <- subset(dfq, p<=x)
		c(Pvalue=x, Selected=nrow(tp), FDR=max(c(tp$q, 0)))
	}
	tabP <- data.frame(tabP)
	rownames(tabP) <- NULL
	}
	if(is.null(FDR)) {
		FDR <- guessFDRfromFDRtable(tabFDR)
		warning(sprintf('FDR not specified. Using guessed FDR=%.3f\n', FDR))
	}
	isFlatBUM <- ifelse(tabFDR[tabFDR[, 'FDR']==0.3, 'Selected']<2, TRUE, FALSE) # guess if this is a flat FDR
	if(length(isFlatBUM)==0) isFlatBUM <- FALSE # in case 0.3 is not specified, it is likely very significant as smaller cutoff is used
	gSelectInd <- which(selectSignificant(bum, alpha=FDR, by = "FDR"))
	# order by decreasing P value
	gSelectInd <- gSelectInd[order(pval[gSelectInd], decreasing=FALSE)]
	gSelectInd_all <- order(pval, decreasing=FALSE) # select all genes but returns in order
	# select all genes with P<0.05
	#browser()
	gSelectInd_P05 <- logic2sortedIndex(pval<0.05, pval, decreasing=FALSE)
	gSelectInd_byPcutoff <- logic2sortedIndex(pval<pcutoff, pval, decreasing=FALSE)
	if(class(qobj)!='try-error'){
		res <- list(pval=pval, FDR=FDR, dfq=dfq, 
			 bum=bum, tabFDR=tabFDR, tabQ=tabQ, tabP=tabP, gSelectInd=gSelectInd, gSelectInd_all=gSelectInd_all, isFlatBUM=isFlatBUM, gSelectInd_P05=gSelectInd_P05, 
			 gSelectInd_byPcutoff=gSelectInd_byPcutoff, pcutoff=pcutoff, fitQ=TRUE)
	} else {
		res <- list(pval=pval, FDR=FDR, 
			 bum=bum, tabFDR=tabFDR, gSelectInd=gSelectInd, gSelectInd_all=gSelectInd_all, isFlatBUM=isFlatBUM, gSelectInd_P05=gSelectInd_P05, 
			 gSelectInd_byPcutoff=gSelectInd_byPcutoff, pcutoff=pcutoff, fitQ=FALSE) 	
	}
	res
	
}
#bumFDR <- create_bum(ll_ttestRes[[i]][, 'pval'])
#gIndexSel <- bumFDR$gSelectInd_P05


# not working: getParName <- function(x){deparse(substitute(x))}

#' plot bum with FDR table on it
#'
#' this usually couples with  create_bum()
#' @param bumFDR a list with at least a bum model and tabFDR
#' @param main main title
#' @param ylab ylab
#' @param tableType not implemented
#' @param shown number of tables to be shown
#' @param breaks number of breaks for the histogram
#' @param digits digits of P value cutoff
#' @param tag used to add a tag to the title
#' @export
plotBumFDR <- function(bumFDR, xpos=0.5, yadj=0.7, main='', xlab='p Values', ylab='Density', tableType=NA, 
	digits=4, tag='', breaks=100, shown=2){
	if(is.null(main)){
		main <- str_c(tag, sprintf('Bum plot for %s', deparse(substitute(bumFDR))), sep='\n')
	}
	#if(is.na(tableType)) tableType <- ifelse(bumFDR$isFlatBUM, 'P', 'FDR') # guess table type
	# nclass=100 added based on oompa: this is to make sure it is actually the ymax
	#browser()
	bum <- bumFDR$bum	
	if(!is.null(bumFDR$fitQ)){
		dfq <- bumFDR$dfq
		reso <- breaks
		#browser()
		rhist <- hist(dfq$p, breaks=reso, plot=F)
		xvals <- (0:reso)/reso
		fit <- bum@lhat + (1-bum@lhat)*dbeta(xvals, bum@ahat, 1)
		#lines(xvals, fit, col='#025627', lwd=2)
		#abline(h=bum@pihat, col='#1109ed', lwd=2)
		dff <- data.frame(x=xvals, y=fit)  
		dff <- subset(dff, y<Inf) # no overshoot of the figure
		p1 <- ggplot(dfq, aes(x=p)) + 
				geom_histogram(aes(y=..density..), breaks=rhist$breaks, fill = "white", colour='black') + 
				xlim(c(0, 1))+xlab(xlab)+ylab(ylab)+ggtitle(main)+
				geom_line(data=dff, aes(x, y), colour='#025627', size=0.9)+
				geom_segment(aes(x = 0, y = bum@pihat, xend = 1, yend = bum@pihat),  colour='#1109ed', size=0.4)+
				theme(panel.border = element_blank(),
				axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
				axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black"))
		
		tabFDR <- bumFDR$tabFDR
		tabQ <- bumFDR$tabQ
		tabQ1 <- tabQ
		colnames(tabQ1)[1] <- 'FDR'
		tabP <- bumFDR$tabP
		formatTab <- function(x, digits=digits){
			x[, 3] <- signif(x[, 3], digits=digits)
			x
		}
		tblTheme <- gridExtra::ttheme_default()
		tbl_FDR <- gridExtra::tableGrob(formatTab(tabFDR), rows=NULL, theme=tblTheme)
		tbl_Q <- gridExtra::tableGrob(formatTab(tabQ), rows=NULL, theme=tblTheme)
		tabQ1 <- gridExtra::tableGrob(formatTab(tabQ1), rows=NULL, theme=tblTheme)
		tbl_P <- gridExtra::tableGrob(formatTab(tabP), rows=NULL, theme=tblTheme)
		
		
		tbl_FDR$widths <- unit(rep(1/3.1, 3), "npc")
		tbl_Q$widths <- unit(rep(1/3.1, 3), "npc")
		tabQ1$widths <- unit(rep(1/3.1, 3), "npc")
		tbl_P$widths <- unit(rep(1/3.1, 3), "npc")
		#blank <- grid.rect(gp=gpar(col="white"))
		if(shown==3){
			gridExtra::grid.arrange(p1, tbl_FDR, tbl_Q, tbl_P,
					nrow=2, layout_matrix = rbind(c(1,1,1), c(2,3,4)), 
					as.table=TRUE,
					heights=c(3,1.5))				 
		}
		if(shown==2){
			# this is not consistent!
			# grid.arrange(p1, tbl_FDR, tbl_P,
			# 		nrow=2, layout_matrix = rbind(c(1,1), c(2,3)), 
			# 		as.table=TRUE,
			# 		heights=c(3,1.5))	
			gridExtra::grid.arrange(p1, tabQ1, tbl_P,
					nrow=2, layout_matrix = rbind(c(1,1), c(2,3)), 
					as.table=TRUE,
					heights=c(3,1.5))				 
		}	
	} else {
		hist(bum, main=main)
		#browser()
	}	
}

### old version
plotBumFDR0 <- function(bumFDR, xpos=0.5, yadj=0.7, main, xlab='p Values', ylab='Density', cex.axis=1, cex.lab=1, 
	showFDRtable=TRUE, showPcutoff=TRUE, tableType=NA, P_thr=c(1e-4, 1e-3, 1e-2, 0.05), digits=4, tag=NULL, breaks=NULL){
	if(missing(main)){
		main <- str_c(tag, sprintf('Bum plot for %s', deparse(substitute(bumFDR))), sep='\n')
	}
	if(is.na(tableType)) tableType <- ifelse(bumFDR$isFlatBUM, 'P', 'FDR') # guess table type
	# nclass=100 added based on oompa: this is to make sure it is actually the ymax
	yPos <- max(hist(bumFDR$bum@pvals, nclass=100, plot=FALSE)$density)*yadj
	#browser()
	if(is.null(breaks)){
		hist(bumFDR$bum, main=main, ylab=ylab, xlab=xlab, cex.lab=cex.lab, cex.axis=cex.axis, xlim=c(0, 1))	
	} else {
		hist(bumFDR$bum, main=main, ylab=ylab, xlab=xlab, cex.lab=cex.lab, cex.axis=cex.axis, breaks=breaks, xlim=c(0, 1))
	}
	
	#browser()
	if(showFDRtable){
		if(tableType=='FDR'){
			tpTabFDR <- bumFDR$tabFDR
			tpTabFDR[, 3] <- signif(tpTabFDR[, 3], digits=digits)
			if(showPcutoff){
			addtable2plot(xpos, yPos, tpTabFDR[, 1:3],bty="o",display.rownames=FALSE,hlines=TRUE,
                         vlines=TRUE,title=NULL, text.col='red')
			} else {
			addtable2plot(xpos, yPos, tpTabFDR[, 1:2],bty="o",display.rownames=FALSE,hlines=TRUE,
                         vlines=TRUE,title=NULL, text.col='red')
			}
		} else {
			tpTabByP <- getPvecSmry(P=bumFDR$bum@pvals, thr=P_thr)$smryTab
			colnames(tpTabByP) <- c('P cutoff', 'Selected')
			addtable2plot(xpos, yPos, tpTabByP[, 1:2],bty="o",display.rownames=FALSE,hlines=TRUE,
                         vlines=TRUE,title=NULL, text.col='red')
		}
	}					 
}
#plotBumFDR(t_bum, xpos=0.6, main='', showPcutoff=F, showFDRtable=F, ylab='Probability density', xlab='p Values')

########################################################################################################################################
# evaluation/validation of hierarchical clustering stability
########################################################################################################################################


#' calculate cophenetic coefficient for a hierarchical clustering result
#' 
#' To compute the Cophenetic Correlation Coefficient of hierarchical clustering, we need two informations: Distance matrix  and Cophenetic Matrix 
#' Cophenetic Correlation Coefficient is simply correlation coefficient between distance matrix and Cophenetic matrix
#' As the value of the Cophenetic Correlation Coefficient is quite close to 1, we can say that the clustering is quite fit. A random data would have Cophenetic Correlation Coefficient
#' near 0.
#'
#' @param dist a dist object which can be constructed by dist, ClassDiscovery::distanceMatrix, ade4::dist.binary or my binaryDist function
#' @param hclust a hclust object returned by hclust function (either from stats::hclust or fastcluster::hclust
#' @return a correlation coeficient scalar
#' @export
copheneticCoef <- function(dist, hclust){
	d1 <- dist
	d2 <- cophenetic(hclust)
	res <- cor(d1, d2)
	res
}
# copheneticCoef(dist=distO, hclust=hclustO)




#' extract silhouette  information for clustering
#'
#' the input is the same as copheneticCoef
#' @param dist a dist object which can be constructed by dist, ClassDiscovery::distanceMatrix, ade4::dist.binary or my binaryDist function
#' @param hclust a hclust object returned by hclust function (either from stats::hclust or fastcluster::hclust
#' @param ... additional parameters to be passed to cutree, i.e. k and h
#' @return an object sil
#' @export
buildSilhouette <- function(dist, hclust, ...){
	# in case the user specified a membership vector
	if(is.vector(dist)) {
		membership <- dist
	} else {
	# in this case, use cutree to generate membership
		membership <- cutree(hclust, ...)
	}
	sk <- silhouette(membership, dist)
	sk
}
#sil <- buildSilhouette(dist=distO, hclust=hclustO, k=6)

#####
# seems this is not needed due to the powerful ConsensusClusterPlus package!
#####
#' Consensus clustering based on BootstrapClusterTest from OOMPA
#'
#' Bootstrap genes to calculate a consensus matrix
#'
#' There are generally two doable approach to obtain a consensus matrix for clustering samples: 
#' (1) subsampling samples. (2) bootstrap genes. Bootstrap samples is not recommended since when a sample
#' is sampled multiple times, it will form a cluster by itself.
#'  This function applies bootstrap on the genes (rows of a matrix) so as to perturb the data when generating a consensus matrix 
#' 
#' @param data A data matrix, numerical data frame, or ExpressionSet object.
#' @param FUN A function that, given a data matrix, returns a vector of cluster assignments based on the columns
#' Examples of functions with this behavior are cutHclust, cutKmeans, cutPam, and cutRepeatedKmeans.
#' @param subsetSize An optional integer argument used to select a subset of genes. If present, each iteration of the bootstrap selects
#'  subsetSize rows from the original data matrix. If missing, each bootstrap contains the same number of rows as the original data matrix.
#' @param nTimes The number of bootstrap samples to collect.
#' @param verbose A logical flag
#' @param ... Additional arguments passed to the classifying function, FUN.
BootstrapGenes <- function (data, FUN, subsetSize, nTimes = 100, verbose = TRUE, 
    ...) 
{
    call <- match.call()
    if (inherits(data, "ExpressionSet")) {
        data <- exprs(data)
    }
    N <- ncol(data)
    if (missing(subsetSize)) {
        subsetSize <- nrow(data)
    }
    subsetSize <- as.integer(subsetSize)
    bootMatch <- matrix(0, nrow = N, ncol = N)
    for (i1 in 1:nTimes) {
        tempData <- data[sample(nrow(data), subsetSize, replace = TRUE), 
            ]
        if (verbose) {
            cat(paste("[", i1, "] ", nrow(tempData), " ", sep = ""))
            if (i1%%10 == 0) 
                cat("\n")
        }
        tempCut <- FUN(tempData, ...)
        K <- max(tempCut)
        tempMatch <- matrix(0, nrow = N, ncol = N)
        for (i2 in 1:K) {
            tempMatch[tempCut == i2, tempCut == i2] <- 1
        }
        bootMatch <- bootMatch + tempMatch
    }
    dimnames(bootMatch) <- list(colnames(data), colnames(data))
    if (verbose) 
        cat("\n")
    testResult <- new("ClusterTest", call = call, result = bootMatch/nTimes)
    new("BootstrapClusterTest", testResult, f = FUN, nTimes = nTimes, 
        subsetSize = subsetSize)
}

#' calculate dispersion from a consensus matrix
#' extracted from NMF package
#' @param object a consensus matrix between 0 and 1
#' @return a scaler between 0 and 1
#' @export
dispersionFromConsMat <- function(object){
	stopifnot( nrow(object) == ncol(object) )
	sum( 4 * (object-1/2)^2 ) / nrow(object)^2
}


#' calculate cophenetic coefficient for a consensus matrix
#' extracted from NMF package
#' @param object a consensus matrix between 0 and 1
#' @param linkage linkage rule to calculate cophenetic distance
#' @return a scaler between 0 and 1
#' @export
cophCorFromConsMat <- function(object, linkage='average'){
         # check for empty matrix
          if( nrow(object)==0  || ncol(object)==0 )
          {
               warning("NA produced [input matrix is of dimension ", nrow(object), "x", ncol(object), "]"
                         , call.=FALSE)
               return(NA)
          }         
          # safe-guard for diagonal matrix: to prevent error in 'cor'
          if( all(object[upper.tri(object)]==0) && all(diag(object)==object[1,1]) )
               return(1)
          # convert consensus matrix into dissimilarities
          d.consensus <- as.dist(1 - object)
          # compute cophenetic distance based on these dissimilarities
          hc <- hclust(d.consensus, method=linkage)
          d.coph <- cophenetic(hc)
          # return correlation between the two distances
          res <- cor(d.consensus, d.coph, method='pearson')
          return(res)
}

#' plot the fit from ConsensusClusterPlus()
#'
#' @param fit result returned by ConsensusClusterPlus()
#' @param type which plot type, either "cophenetic" or "dispersion"
#' @param plot logical if plot is generated
#' @return a data frame when plot=FALSE
#' @export
plotSummaryFromConsensusClusterPlus <- function(fit, type='cophenetic', plot=TRUE, ...){
	type <- match.arg(type, choices=c("cophenetic", "dispersion"))
	#browser()
	Ks <- 2:length(fit)
	coph <- sapply(fit[Ks], function(x) cophCorFromConsMat(x$consensusMatrix))
	dispersion <- sapply(fit[Ks], function(x) dispersionFromConsMat(x$consensusMatrix))
	if(plot){
		switch(type,
			   cophenetic=plot(Ks, coph, type='b', ...),
			   dispersion=plot(Ks, dispersion, type='b', ...))
	} else {
		res <- data.frame(K=Ks, coph=coph, dispersion=dispersion)
		return(res)
	}
}











########################################################################################################################################
# color mapping utility
########################################################################################################################################


#' Generate a named color vector
#'
#' ggplot2 use named color vector for scale_colour_manual. This is an efficient way to generate colors. the namedColVec function facilitates
#' this by requring a color vec and a factor. The number of levels in factor should equal the number of colors. Note this is different from the 
#' buildcolorfacs function
#'
#' when NA is present, we deliberately tell the function what color and what name the user want it to be. To extract just the unique values without NA, we combine
#' unique+na.omit or levels+na.omit
#'
#' @param col color vector
#' @param fac the factor as a reference
#' @param by colors to be named in the order of levels in the factor or the alphabetical order (in this case, fac is not required to be a factor)
#' @param na.col color for NA. it is named by na.name
#' @param na.name name for NA, na.col
#' @param return a named color vector. color for NA is always at the end
#' @export  
#' @examples
#' vv <- letters[c(5, 1, 3, 2, 26)]
#' vvNA <- c(vv, NA)
#' vvfac <- factor(vv)
#' vvNAfac <- factor(vvNA)
#' cols <- c('black', 'red', 'green', 'blue', 'yellow')
#' namedColVec(cols[1:length(unique(na.omit(vv)))], vv, by='alphabetical')
#' namedColVec(cols[1:length(unique(na.omit(vvNA)))], vvNA, by='alphabetical')
#' namedColVec(cols[1:length(unique(na.omit(vvfac)))], vvfac, by='level')
#' namedColVec(cols[1:length(unique(na.omit(vvNAfac)))], vvNAfac, by='level')
namedColVec <- function(col, fac, by='level', na.col='gray', na.name='NA'){
	#browser()
	by <- match.arg(by, choices=c("level", "alphabetical"))
	if(by=="level"){
		## make sure is a factor
		if(!is.factor(fac)) stop('fac needs to be a factor!\n')
		lv <- levels(na.omit(fac)) # when NA is present, it is excluded from lv
		## make sure length is compatible
		if(length(col)!=length(lv)) stop(sprintf('%d cols specified but factor has %d levels\n', length(col), length(lv)))
		colName <- lv
	}
	if(by=="alphabetical"){
		# extracting unique value. first converting to a char vec when it is a factor
		fac <- as.character(fac)
		uvSort <- sort(unique(na.omit(fac)), decreasing=FALSE) # when NA is present, NA is excluded from uvSort
		# make sure length is compatible
		if(length(col)!=length(uvSort)) stop(sprintf('%d cols specified but only observe %d unique values\n', length(col), length(uvSort)))
		colName <- uvSort
	}
	res <- col
	names(res) <- colName
	# handling NA
	if(any(is.na(fac)))
		res[na.name] <- na.col
	res
}




#' prepare a matrix for plotting by converting it to a color matrix
#'
#' this is an application of mapNames4mat() function that maps values to a color space.
#'
#' this is a safe way to map value to color 1-to-1 since we can specify the old value and the new value
#' deliberately.
#' NA can be handled easily
#' @param mat the original matrix
#' @param colmap either a named vector or a 2 column matrix specifying old/new mapping. This might be constructed from namedColVec which handles NA
#' @param na.name when there is NA, colmap must have specified how to map na.name to its color. So to index mat where there is missing val, we replace missing val with na.name
#' @return the color matrix as well as colorfacs for legend returned as a list
#' @export
makeColMat <- function(mat, colmap, na.name='NA'){
	# when colmap is a named vector, build the colmap
	if(is.vector(colmap)) {
		#colvec <- colmap # stores the old colors
		colmap <- cbind(names(colmap), colmap)
	}
	if(any(is.na(mat))) {
		# now we have NA observed. first check na.name is in colmap
		if(!na.name %in% colmap[, 1]) {
			stop(sprintf("na.name %s is not in colmap!\n", na.name))
		}
		# convert NA to na.name
		mat[is.na(mat)] <- na.name
	}	
	cmat <- mapNames4mat(mat, colmap) # colormap
	# prepare a colorfacs for legend
	bc <- colmap[order(colmap[, 1]), 2] # this is due to the artifact of buildcolorfacs that specify colors in the order of alphabetical of the barVecs
	colorfacs <- buildcolorfacs(barVecs=list(values=colmap[, 1]),
                  barCols=list(bc))
	#browser()
	list(cmat=cmat, colorfacs=colorfacs)
}
#makeColMat(mat=mat_str_sf[gOrdered, ], colmap=colPalMutCN)
#makeColMat(mat, colmap=colmap)

#' for a matrix, i.e. 0, 1, we may want to format it as WT, MUT
#'
#' @param mat a matrix to be formated
#' @param map 2-column data frame where first column is input names; 2nd column for output names
#' @param fac2string whether to treat factor as string
#' @return a formated data frame
#' @export
mapNames4mat <- function(mat, map, fac2string=TRUE){
	res <- foreach(i=1:ncol(mat), .combine='cbind') %do% {
		mapNames(mat[, i], map, fac2string=fac2string)
	}
	dimnames(res) <- dimnames(mat)
	res
}

#' mapping a vec to another name space
#'
#' given a vector and a 2-column map, map names for vector based on the map.
#' This is useful if we want to update labels according to a table. 
#'
#' This function returns a new label vector with instructions in the map if a hit is found in the map
#' @param vec the input vector
#' @param map 2-column data frame where first column is input names; 2nd column for output names
#' @param fac2string whether to treat factor as string
#' @return a vector same length as vec
#' @export
#### modified on 2013/12/03: if colnames of map does not include 'old' and 'new', use first and second columns instead
mapNames <- function(vec, map, fac2string=TRUE){
	#browser()
	if(fac2string){
		vec <- makestring(vec)
		map <- colwise(makestring)(data.frame(map))
	}
	res <- vec
	names(res) <- vec
	if(!all(c('old', 'new') %in% colnames(map))) {
		index <- 1:2
	} else {
		index <- match(c('old', 'new'), colnames(map))
	}
	for(i in 1:nrow(map)){
		old <- map[i, index[1]]
		new <- map[i, index[2]]
		if(!is.na(new)){
			ind <- which(vec==old)
			#browser()
			res[ind] <- new
		}
	}
	res
}
#mydf$Condition <- mapNames(mydf$X2, map=data.frame(rownames(sampleAnno), sampleAnno$Condition))

#' impose as.character
#' @param vec a vector, or a factor
#' @return a character vector
#' @export 
makestring <- function(vec){
	if(is.factor(vec)) 
		vec <- as.character(vec)
	vec
}


#' given a mat and phenoDat and barColList, calculates sampleAnno and annoCols that can
#'  be directly passed to aheatmap
#' @param mat the main mat for heatmap
#' @param phenoDat pheno data frame
#' @param barColList column bar color list. This should be a list of color vectors (the color pallete). By default, this is set as NULL,
#' which means getDefaultBarColList() will be called to set the default values.
#' @return a list for aheatmap visualization
#' @export
prepare4aheatmap <- function(mat, phenoDat, barColList=NULL){
	#browser()
	if(ncol(mat)!=nrow(phenoDat)) stop(sprintf('mat has %d columns; phenoDat has %d rows', ncol(mat), nrow(phenoDat)))
	#if(length(barColList)!=ncol(phenoDat)) stop(sprintf('%d columns in phenoDat; %d barColList elements', length(barColList), ncol(phenoDat)))
	#if(!identical(sort(colnames(phenoDat)), sort(names(barColList)))) stop('colnames of phenoDat does not match with names in barColList')
	sampleAnno <- formatAnnotation4aheatmap(phenoDat, sampleNames=colnames(mat))
	annoCols <- prepare4ColList(barColList=barColList, phenoDat=phenoDat)
	res <- list(mat=mat, sampleAnno=sampleAnno, annoCols=annoCols)
	res
}
#' a wrapper to aheatmap which enables passing phenoData. This serves as a quick aheatmat version. 
#' @param phenoDat a data frame for pheno data
#' @export
aheatmap2 <- function(mat, phenoDat, barColList=NULL, ...){
	#browser()
	L <- prepare4aheatmap(mat=mat, phenoDat=phenoDat, barColList=barColList)
	temp <- L$mat
	aheatmap(temp, annotation=L$sampleAnno, 
	annotation_colors = L$annoCols, ...) 
}
#aheatmap2(mat=matGanglia[sig2$sym, ], phenoDat=data.frame(REST=phenoGanglia[, 1]), col=exprColors)
#aheatmap2(mat=matGanglia[sig2$sym, ], phenoDat=phenoGanglia, col=exprColors, barColList=list(REST=bluered(30), RESTscore=coolheat(30)))

#' sometimes we just want a default barColList rather than specify it everytime
#' We first guess which columns of phenoDat is categorical and numeric and than assign colors correspondingly
#' @param phenoDat phenoDat
#' @return a barColList
#' @export
getDefaultBarColList <- function(phenoDat){
	nCol <- 30
	iCatCol <- iContCol <- 1
	res2 <- foreach(i=1:ncol(phenoDat)) %do% {
		tt <- phenoDat[, i]
		if(isCat(tt)){
			res <- myCatColList[[iCatCol]][1:n_unique(tt)]
			iCatCol <- iCatCol+1
		} else {
			res <- myContColList[[iContCol]]
			iContCol <- iContCol+1
		}
		res
	}
	names(res2) <- colnames(phenoDat)
	res2
}
#getDefaultBarColList(phenoGanglia)

#' given a list of barColors and phenoData, prepare a list of colors that can be passed to aheatmap
#' @param barColList a named list of colors
#' @param phenoDat a data frame for pheno data
#' @return a list of vectors each with length nrow(phenoDat)
#' @export
prepare4ColList <- function(barColList=NULL, phenoDat){
	if(is.null(barColList)){
		barColList <- getDefaultBarColList(phenoDat)
	}
	if(length(barColList)!=ncol(phenoDat)) stop(sprintf('%d columns in phenoDat; %d barColList elements', length(barColList), ncol(phenoDat)))
	if(!identical(sort(colnames(phenoDat)), sort(names(barColList)))) stop('colnames of phenoDat does not match with names in barColList')
	Lnames <- names(barColList)
	res <- vector('list')
	#browser()
	res <- foreach(nn=Lnames) %do% {
		if(isCat(phenoDat[, nn])){
			res <- namedColVec(col=barColList[[nn]], fac=phenoDat[, nn]) 
		} else {
			res <- barColList[[nn]]
		}
		res
	}
	names(res) <- Lnames
	res
}
#prepare4ColList(phenoDat=phenoGanglia)
#prepare4ColList(barColList=list(REST=bluered(30), RESTscore=coolheat(30)), phenoDat=phenoGanglia)
#tt <- prepare4aheatmap(mat=matGanglia[sig2$sym, ], phenoDat=phenoGanglia, mainCol=exprColors, barColList=list(REST=bluered(30), RESTscore=coolheat(30)))
