########################################################################################################
#' Detection and prevalence of SARS-CoV-2 co-infections during the Omicron variant circulation, France, December 2021 - February 2022

#' Antonin Bal, Bruno Simon, Gregory Destras, Richard Chalvignac, Quentin Semanas, Antoine Oblette, Gregory Queromes, Remi Fanget, 
#' Hadrien Regue, Florence Morfin, Martine Valette, Bruno Lina, Laurence Josset 
#' 
#' (c) 2022 Laurence Josset, <laurence.josset(at)univ-lyon1.fr>
#' 
#' MIT LICENSE

########################################################################################################

options(stringsAsFactors = FALSE) 
library(dplyr)
require(data.table)
library(ggplot2)
library(reshape)
library(cowplot)

########################################################################################################
## load source data
########################################################################################################

load("SOURCE_Variant_Table_MIXES.Rdata")
load("SOURCE_seq_metrics_MIXES.Rdata")


########################################################################################################
### EXCLUDE SAMPLES OF RUN 220530 (PERFORMED TO TEST THE DETECTION LIMIT) : To include all samples, comment the 2 following lines
########################################################################################################
selected_samples = seq_metrics$sample[seq_metrics$RUN != "220530_13146"]
seq_metrics = seq_metrics[seq_metrics$RUN != "220530_13146",]

vcf = vcf[is.element(vcf$sample,selected_samples),]


########################################################################################################
## Annotate the sequencing metrics file for the plots
########################################################################################################
seq_metrics$Coinf = seq_metrics$Delta > 0 & seq_metrics$Omicron > 0

seq_metrics$Delta = as.numeric(seq_metrics$Delta) ## % of Delta strain will be used for the color gradient

seq_metrics$coinf_min_match[seq_metrics$coinf_min_match==""] <- "None"

########################################################################################################
## Either AY.102 or B.1.617.2 identified as they have the same list of defining polymorphism, except for the S:G142D found in AY.102 only
########################################################################################################
seq_metrics$coinf_maj_match = gsub("AY.102","B.1.617.2",seq_metrics$coinf_maj_match)
seq_metrics$coinf_min_match = gsub("AY.102","B.1.617.2",seq_metrics$coinf_min_match)

########################################################################################################
## Fig 2.A : Identification of MAIN LINEAGE
########################################################################################################
p1 <- ggplot(data=seq_metrics, aes(x=factor(coinf_maj_match), y=coinf_maj_ratio ,colour=Delta)) 

p1 <-	p1 +
		geom_jitter(aes(shape=Coinf),alpha=0.9,width = 0.25, height = 0) +
		facet_grid(.~proto) + 
		ylab("Main Lineage Mutation Ratio") + xlab("Main Lineage Identification") +
		ylim(c(0.5,1)) +
		scale_color_gradient2(low = "steelblue1", mid = "cyan4", high = "tomato",
			breaks=c(0,10,20,30,40,50,60,70,80,90,100), midpoint = 50,
			name = "Delta:Omicron", labels = c("0:100","10:90","20:80", "30:70","40:60","50:50","60:40","70:30", "80:20","90:10","100:0")) +
		scale_shape(name = "Experimental\nCoinfection") +
		theme_bw()  + theme( axis.text.x  = element_text(vjust=0.5, size=10,angle=45))


########################################################################################################
##  Fig 2.B : Identification of SECONDARY LINEAGE
########################################################################################################
p2 <- ggplot(data=seq_metrics, aes(x=factor(coinf_min_match), y=coinf_min_ratio ,colour=Delta)) 

p2 <-	p2 +
		geom_jitter(aes(shape=Coinf),alpha=0.9,width = 0.25, height = 0) +
		facet_grid(.~proto) + 
		ylab("Secondary Lineage Mutation Ratio") + xlab("Secondary Lineage Identification") +
		ylim(c(0,1)) +
		scale_color_gradient2(low = "steelblue1", mid = "cyan4", high = "tomato",
			breaks=c(0,10,20,30,40,50,60,70,80,90,100), midpoint = 50,
			name = "Delta:Omicron", labels = c("0:100","10:90","20:80", "30:70","40:60","50:50","60:40","70:30", "80:20","90:10","100:0")) +
		scale_shape(name = "Experimental\nCoinfection") +
		theme_bw()  + theme( axis.text.x  = element_text(vjust=0.5, size=10,angle=45))


########################################################################################################
## Fig 2
########################################################################################################
plot_grid(p1, p2, ncol = 1,labels = c("A","B"), align = 'v',axis = "rl")
ggsave("Fig2.pdf",width = 9 ,height=5.3 )
ggsave("Fig2.png",width = 9 ,height=5.3 )