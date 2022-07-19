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
#### CODE USED TO GENERATE FIG S7
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
 
# for samples sequenced in duplicate, determine which is the first and second passage based on Year of sequencing, and Plate number (Pl)
first_replicate = ddply(seq_metrics, .(ID),function(x) x[order(x$Year,x$pl),][1,])

seq_metrics$Replicate = NA
seq_metrics$Replicate[is.element(seq_metrics$sample,first_replicate$sample)] = "First"
seq_metrics$Replicate[!is.element(seq_metrics$sample,first_replicate$sample)] = "Second"

#############################################################################################################################
## IMPACT OF CHANGING THE TRESHOLD ON PREVALENCE
#############################################################################################################################
#### Compute theoritical prevalence for each treshold
ntot = length(first_replicate$sample)

prev = first_replicate %>%
	  group_by(coinf_min_common) %>%
	  dplyr::summarise(
		n = n())


prev = prev[prev$coinf_min_common>0,]	
cs = prev[order(prev$coinf_min_common,decreasing = TRUE),] %>% mutate(cs = cumsum(n))

cs$prev_estimate = sapply(cs$cs , function(x) binom.test(x,ntot)$estimate)	
cs$CI1 = sapply(cs$cs , function(x) binom.test(x,ntot)$conf.int[1])
cs$CI2 = sapply(cs$cs , function(x) binom.test(x,ntot)$conf.int[2])



p <- ggplot(cs, aes(x=prev_estimate,y=coinf_min_common)) 

p <- p + 
	geom_vline(xintercept = as.numeric(cs[cs$coinf_min_common==6,c("CI1","prev_estimate","CI2")]),linetype =  "dotted", color = "red") +
	geom_point() + 
	geom_errorbarh(aes(xmin = CI1, xmax = CI2), height = 0.25) +
	xlab("Prevalence estimate (95% CI)") +
	ylab("Number of defining-variants of secondary lineage") +
	 theme_bw() 

### Stats : compare prevalence to prevalence found with current treshold (coinf_min_common=6)
n_6 = cs$cs[cs$coinf_min_common==6]

cs$Fisher = sapply(cs$cs , function(x) fisher.test(matrix(c(n_6,ntot-n_6,x,ntot-x),ncol=2,byrow=TRUE))$p.value)	

## inset for more than 1 specific mutation

p_inset <- ggplot(cs[cs$coinf_min_common>1,], aes(x=prev_estimate,y=coinf_min_common)) 

p_inset <- p_inset + 
	geom_vline(xintercept = as.numeric(cs[cs$coinf_min_common==6,c("CI1","prev_estimate","CI2")]),linetype =  "dotted", color = "red") +
	geom_point() + 
	geom_text(data=cs[cs$coinf_min_common>1 & cs$Fisher<0.05,],aes(x=CI2,y=coinf_min_common),label="*",nudge_x = 0.0001)+
	geom_errorbarh(aes(xmin = CI1, xmax = CI2), height = 0.25) +
	xlab("Prevalence estimate (95% CI)") +
	ylab("Number of defining-variants of secondary lineage") +
	 theme_bw() 

### save

ggdraw(p + theme_half_open(12)) +
  draw_plot(p_inset, .25, .25, .75, .75) 

ggsave("FigS7.pdf",width=7,height=7)

## which are the cut-off with no impact on prevalence? between 3-17
cs[cs$Fisher>0.05,]
