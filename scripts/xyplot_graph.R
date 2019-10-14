#! /usr/bin/env Rscript

library(ggplot2)
library(ggExtra)
library(optparse)

################################################################
## FUNCTIONS
################################################################
ggplotRegression <- function (data, x_axis, y_axis, category) {
	fit <-lm(data[[opt$y_column]] ~ data[[opt$x_column]], data = data)
	# Code from https://susanejohnston.wordpress.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/
	g <- ggplot(data, aes_string(x = x_axis, y = y_axis, colour = category)) +
	  geom_point(alpha=.2) +
	   scale_color_manual(values=c("#FF0000", "#00FF00", "#0000FF")) +
	  #stat_smooth(method = "lm", col = "red") +
	  # labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
	  #                    "Intercept =",signif(fit$coef[[1]],5 ),
	  #                    " Slope =",signif(fit$coef[[2]], 5),
	  #                    " P =",signif(summary(fit)$coef[2,4], 5))) +
	  xlab(x_axis) +
	  ylab(y_axis) + 
	  theme(legend.position='bottom')
	ggMarginal(g, type= "boxplot", groupColour = TRUE, groupFill = TRUE)
}


################################################################
## OPTPARSE
################################################################
option_list <- list(
	make_option(c("-d", "--data_file"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-o", "--output"), type="character", default="results",
		help="Output figure file"),
	make_option(c("-x", "--x_column"), type="character", 
		help="Name of column to be used for X dimension"),
	make_option(c("-y", "--y_column"), type="character", 
		help="Name of column to be used for Y dimension"),
	make_option(c("-g", "--group_column"), type="character", 
		help="Name of the group column")

)
opt <- parse_args(OptionParser(option_list=option_list))


################################################################
## MAIN
################################################################

data <- read.table(opt$data_file, sep="\t", header=TRUE)
pdf(paste(opt$output, '.pdf', sep=""))
	# ggplot(data, aes(x=data[[opt$x_column]], y=data[[opt$y_column]])) +
 	#  	geom_point(shape=1) +
	# 	xlab(opt$x_column) +
	# 	ylab(opt$y_column) +
	# 	geom_smooth(method=lm)
	g <- ggplotRegression(data,
		opt$x_column,
		opt$y_column,
		opt$group_column
	 )
	g
dev.off()
