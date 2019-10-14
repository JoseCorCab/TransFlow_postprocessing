#! /usr/bin/env Rscript

library(ggplot2)
library(optparse)

################################################################
## OPTPARSE
################################################################

option_list <- list(
	make_option(c("-d", "--data_file"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-o", "--output"), type="character", default="results",
		help="Output figure file"),
	make_option(c("-x", "--species"), type="character", 
		help="Name of column with the measurement categories"),
	make_option(c("-t", "--title"), type="character", default="Boxplot",
		help="Figure title"),
	make_option(c("-y", "--values"), type="character", 
		help="Name of column to be used for Y dimension")


)
opt <- parse_args(OptionParser(option_list=option_list))

################################################################
## MAIN
################################################################

data <- read.table(opt$data_file, header = TRUE , sep="\t")
data[[opt$species]] <- factor(data[[opt$species]], levels=c("Unmapped TTs","Atypical TTs", "Low Quality TTs", "Reliable TTs"))
#if(!is.na(opt$desviation){
#data <- summarySE(data, measurevar=opt$values, groupvars=c(opt$conditions))
pdf(paste(opt$output, '.pdf', sep=""))
		obj <- ggplot(data, aes(y=data[[opt$values]], x=data[[opt$species]])) +
			geom_boxplot() +
 			labs(title=opt$table_name, y="Length (bp)") +
			ggtitle(opt$title) + 
 			theme(axis.text.y = element_text(size=15, face="bold"),
			axis.text.x = element_text(size=11, face="bold"), 
			axis.title=element_text(size=17 ,face="bold"),
			axis.title.x=element_blank(),
			plot.title=element_text(size=20, face="bold")) +
			ylim(0, 25000)
			 

	

		ggplot_build(obj)
    	print(ggplot_build(obj))


dev.off()
