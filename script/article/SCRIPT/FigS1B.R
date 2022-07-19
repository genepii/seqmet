############################################################################################################################
#' Detection and prevalence of SARS-CoV-2 co-infections during the Omicron variant circulation, France, December 2021 - February 2022

#' Antonin Bal, Bruno Simon, Gregory Destras, Richard Chalvignac, Quentin Semanas, Antoine Oblette, Gregory Queromes, Remi Fanget, 
#' Hadrien Regue, Florence Morfin, Martine Valette, Bruno Lina, Laurence Josset 
#' 
#' (c) 2022 Laurence Josset, <laurence.josset(at)univ-lyon1.fr>
#' 
#' MIT LICENSE

############################################################################################################################

############################################################################################################################
#### CODE USED TO GENERATE FIG S1
############################################################################################################################

options(stringsAsFactors = FALSE) 

library(ggplot2)
library(plyr)
library(dplyr)

library(cowplot)

############################################################################################################################
## DATA SOURCE : SEQUENCING METRICS (seq_metrics) for n=21448 sequencing results (21387 samples + 61 duplicate sequencing) 
############################################################################################################################
load("FORMATED_seq_metrics.Rdata")
 
#############################################################################################################################
## Fig 1B Calculate Distribution of the number of secondary lineage specific mutations (coinf_min_common)
#############################################################################################################################

# for samples sequenced in duplicate, determine which is the first and second passage based on Year of sequencing, and Plate number (Pl)
first_replicate = ddply(seq_metrics, .(ID),function(x) x[order(x$Year,x$pl),][1,])

seq_metrics$Replicate = NA
seq_metrics$Replicate[is.element(seq_metrics$sample,first_replicate$sample)] = "First"
seq_metrics$Replicate[!is.element(seq_metrics$sample,first_replicate$sample)] = "Second"

## Distrib VARCOUNT

p <- ggplot(seq_metrics[seq_metrics$Replicate=="First",], aes(x=coinf_min_common)) 

p <- p + geom_histogram(,binwidth=1)+
	ggtitle("All samples") +
	xlab("Number of specific-variants from secondary lineage") +
	 theme_bw() 

p_inset <- ggplot(seq_metrics[seq_metrics$Replicate=="First"& seq_metrics$coinf_min_common>1,], aes(x=coinf_min_common)) 

p_inset <- p_inset + geom_histogram(,binwidth=1)+
	ggtitle("Samples with more than 1 specific-variants\nfrom secondary lineage") +
	xlab("Number of specific-variants from secondary lineage") +
	 theme_bw() 

ggdraw(p + theme_half_open(12)) +
  draw_plot(p_inset, .25, .25, .75, .75) 

ggsave("FigS1B.pdf",width=7,height=5)
