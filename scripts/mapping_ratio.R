#! /usr/bin/env Rscript

library(ggplot2)
library(optparse)
library(gridExtra)
library(dplyr)
library(ggthemes)



################################################################
## OPTPARSE
################################################################

option_list <- list(
	make_option(c("-b", "--barplot"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-s", "--scatter"), type="character",
		help="Tabulated file with information about each sample"),
	make_option(c("-o", "--output"), type="character", default="results",
		help="Output figure file")


)
opt <- parse_args(OptionParser(option_list=option_list))

################################################################
## MAIN
################################################################

barplot_data <- read.table(opt$barplot, header = TRUE , sep="\t")
scatter_data <- read.table(opt$scatter, header = TRUE , sep="\t")

#data[[opt$species]] <- factor(data[[opt$species]], levels=c("Raw","Filtered"))


# data <- transform(data, rate = Aligned_reads / Size)
# cdata <- ddply(data, c("Transcriptome", "Stage"), summarise,
#                N    = length(rate),
#                mean = mean(rate),
#                sd   = sd(rate),
#                se   = sd / sqrt(N)["mapped_transcripts"],["mapped_transcripts_v4"]
# )
#["(%) Mapped reads"],["Reads per TTs"],["(%) Mapped TTs"],["(%) Mapped TTs [v4.0]"]]
barplot_data$metric <- factor(barplot_data$metric, levels=c("Mapped reads", "Mapped TTs", "Mapped TTs [v4.0]"))
barplot_data$transcriptome <- factor(barplot_data$transcriptome, levels=c("v4.0", "v4.0*", "Oases_k45*", "Min2_Oases_Cap3*"))


scatter_data$transcriptome <- factor(scatter_data$transcriptome, levels=c("v4.0", "v4.0*", "Oases_k45*", "Min2_Oases_Cap3*"))


print(barplot_data)
print(scatter_data)
#quit(status = 1)

 pdf(paste(opt$output, '.pdf', sep=""))
	slc <- 100
   	cbp2 <- c("#F0E442", "#D55E00", "#56B4E9","#0072B2", "#009E73", "#000000", "#E69F00", "#CC79A7")
	obj <- ggplot(data=barplot_data 
			) + 
	geom_bar(aes(y = mean,x = transcriptome, fill = metric),
			position=position_dodge(),
			stat="identity") + 
   	scale_fill_manual(values = cbp2) +
	geom_errorbar(data=barplot_data,aes( x = transcriptome,ymin=mean-se, ymax=mean+se, group = metric), width=.2, position=position_dodge(0.9)) +

   	 theme_classic()+ 
   	# scale_colour_hc()+
   	geom_errorbar(data=scatter_data, aes(x=transcriptome, ymin=(mean-se)/slc, ymax=(mean+se)/slc), width=.2,position = position_nudge(x = +0.1)) +
	geom_line(data = scatter_data, aes(x=transcriptome, y = mean/slc, group = 1, linetype = metric), colour = "black", position = position_nudge(x = +0.1)) +
   	geom_point(data=scatter_data, aes(x=transcriptome, y = mean/slc, shape = metric), colour = "black", size = 1, position = position_nudge(x = +0.1),colour = "red") +	
   	scale_y_continuous(breaks = seq(0,1,0.2),sec.axis = sec_axis(~.*slc, name = "Reads per TTs", breaks = seq(0, max(scatter_data$mean) +20, 20))) + 
   	theme(legend.position = "top", 
   		legend.title = element_blank(),
   		axis.text.y = element_text(size = 12),
   		axis.text.x = element_text(size = 12),
   		legend.text=element_text(size=12),
   		panel.grid.major.y = element_line(colour = "grey", linetype = "dashed"),
   		axis.title = element_text(size= 14, face = "bold"))+
   	labs(x="Transcriptome", y = "Ratio") 
   	# 			geom_bar( position=position_dodge(), stat="identity") +
#  			theme(axis.text.x = element_text(face="bold", size=10, angle=20, hjust=0.8)) +
#  			ylim(0, 40000)
 		
# 		if(!is.na(opt$desviation){
# 			data <- summarySE(data, measurevar=opt$values, groupvars=c(opt$conditions))
# 		obj <- ggplot(cdata, aes(x=Transcriptome, y=mean, fill=Stage)) + 
# 		 		scale_fill_grey(start=0.7, end=0.3) + theme_classic() +

# 		theme(axis.text.x = element_text(face="bold", size=10, angle=5, hjust=0.8))
# 		geom_boxplot() +
# ###################### histogram
# 		scale_color_grey(start=0.6, end=0.2) + theme_classic() +
#  	 	geom_line(aes(color=data[[opt$conditions]], linetype=data[[opt$conditions]])) +
#   		geom_point(aes(color=data[[opt$conditions]])) +
#  		labs(title='', x=opt$species, y='Ratio', fill=opt$conditions) +
# 		theme(axis.text.y = element_text(size=15), axis.text.x = element_text(size=15), axis.title=element_text(size=17 ,face="bold"), legend.text=element_text(size=15), legend.title=element_text(size=15, face="bold"), legend.position="top")

#   		ylim(0, 12000)
#  			labs(title='', x='Stages', y='', group="Transcriptome")+
#  			  theme(legend.position="bottom") 
 			  

# 		}
	  

		obj

dev.off()
