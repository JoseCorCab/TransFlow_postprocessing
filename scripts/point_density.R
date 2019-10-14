#! /usr/bin/env Rscript

library(ggplot2)
library(optparse)
library(RColorBrewer)

################################################################
## OPTPARSE
################################################################

option_list <- list(
	make_option(c("-d", "--data_file"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-n", "--table_name"), type="character",
                help="Table name"),
	make_option(c("-o", "--output"), type="character", default="results",
		help="Output figure file"),
	make_option(c("-c", "--conditions"), type="character", 
		help="Name of column with the conditions names"),
	make_option(c("-l", "--legend"), type="logical", default = FALSE, action="store_true",   
		help="Use this flag for skip legend from graph %default"),
	make_option(c("-Y", "--ymax"), type="double", 
		help="Max limit in plot"),
	make_option(c("-x", "--species"), type="character", 
		help="Name of column with the measurement categories"),
	make_option(c("-y", "--values"), type="character", 
		help="Name of column to be used for Y dimension")
	

)
opt <- parse_args(OptionParser(option_list=option_list))

################################################################
## MAIN
################################################################

data <- read.table(opt$data_file, header = TRUE , sep="\t")

sorted <- subset(data, data[[opt$species]] == "Raw")
sorted2 <- sorted[order(-sorted[[opt$values]]) , c(2)] 
data[[opt$conditions]] <- factor(data[[opt$conditions]], levels=c("Artifacts","Unknown","Complete", "C-terminal","Internal","ncRNAs","N-terminal","New coding"))
print(data)
data[[opt$species]] <- factor(data[[opt$species]], levels=c("Raw TTs","Typical TTs", "Reliable TTs"))
#data[[opt$species]] <- factor(data[[opt$species]], levels=c("Raw","Typical", "Reliable"))
#if(!is.na(opt$desviation){
#data <- summarySE(data, measurevar=opt$values, groupvars=c(opt$conditions))
print("test")
pdf(paste(opt$output, '.pdf', sep=""))

		obj <- ggplot(data, aes(x=data[[opt$species]], y=data[[opt$values]])) +
			#geom_bar( position=position_dodge(), stat="identity") +
 			#geom_point(aes(color = data[[opt$values]], size=data[[opt$conditions]])) + 
			geom_area(aes(group = data[[opt$conditions]], fill = data[[opt$conditions]])) +
			labs(x=opt$table_name, y=opt$values, fill=opt$conditions) 
 			cbp2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
 			          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

 			obj <- obj + 
 			scale_fill_manual(values = cbp2) +
 			theme_classic() +	
 			theme(
			plot.title = element_text(size = 20, face = "bold"),
			axis.text.y = element_text(size=15),
 			axis.text.x = element_text(size=15, angle=20, hjust=0.8),
 			axis.title=element_text(size=17 ,face="bold"), 
 			legend.text=element_text(size=15), 
 			legend.title=element_text(size=15, face="bold"),
 			legend.position = "bottom") +
 			guides(fill=guide_legend(title.position="top", title.hjust =0.5)) +
 			scale_x_discrete(expand = c(0.1,0.0)) 
 			#  scale_y_discrete(expand = c(0,0))
 			if(opt$legend){
 				obj <- obj + theme(legend.position = "none") 
 			}
 			if(!is.null(opt$ymax)){
 				obj <- obj + ylim(0,opt$ymax)
 			}
 			#ylim(0, 90000)
 	#	obj <- obj + scale_fill_brewer(palette="Paired")
	

		obj

dev.off()
