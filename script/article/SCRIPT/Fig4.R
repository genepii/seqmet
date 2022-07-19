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
#### CODE USED TO GENERATE FIG 4 AND CALCULATE PREVALENCE
############################################################################################################################

options(stringsAsFactors = FALSE) 

library(reshape)
library(plyr)
library(ggplot2)
library(ggsci)
library(cowplot)
library(reshape2)

############################################################################################################################
## DATA SOURCE : SOURCE_ANONYMISED_METADATA.Rdata with sampling dates and weeks, final results of WGS (in CLADE column) for n=21387 samples
############################################################################################################################
load("SOURCE_ANONYMISED_METADATA.Rdata")


############################################################################################################################
## CALCULATE AND WRITE PREVALENCE
############################################################################################################################

n_CLADE = data.frame(table(metadata[c("WEEK","CLADE")]))
n_CLADE_cast = dcast(n_CLADE, WEEK ~ CLADE , value.var = "Freq")

Perc_CLADE_cast = t(apply(n_CLADE_cast[,-1],1, function(x) x/sum(x)))
colnames(Perc_CLADE_cast) = paste("Ratio",colnames(Perc_CLADE_cast),sep="_")

CLADE_cast = cbind(n_CLADE_cast,Perc_CLADE_cast)
rownames(CLADE_cast) = CLADE_cast$WEEK
CLADE_cast = CLADE_cast[,colnames(CLADE_cast) != "WEEK"]
write.table(CLADE_cast,file="Prevalence.txt", sep="\t", quote=F,col.names=NA)

## Co-circulation Delta and BA.1 : weeks with >1% circulation BA.1 and Delta 
CLADE_cast_DO = CLADE_cast[CLADE_cast$Ratio_Delta>0.01 & CLADE_cast[,"Ratio_Omicron(BA.1)"]>0.01,]
## estimate Delta/BA.1 prevalence
ntot = sum(CLADE_cast_DO[,grep("Ratio",colnames(CLADE_cast_DO), invert=TRUE)])
binom.test(sum(CLADE_cast_DO[,"Delta/Omicron(BA.1)"]),ntot) 


## Co-circulation BA.1 and BA.2 : weeks with >1% circulation BA.1 and BA.2 
CLADE_cast_OO = CLADE_cast[CLADE_cast[,"Ratio_Omicron(BA.1)"]>0.01 & CLADE_cast[,"Ratio_Omicron(BA.2)"]>0.01,]
## estimate Delta/BA.1 prevalence
ntot = sum(CLADE_cast_OO[,grep("Ratio",colnames(CLADE_cast_OO), invert=TRUE)])
binom.test(sum(CLADE_cast_OO[,"Omicron(BA.1)/Omicron(BA.2)"]),ntot) 


############################################################################################################################
### PLOT FIG 4
############################################################################################################################

## PANEL A. NUMBER of WHOLE GENOME SEQUENCES WITH COINFECTION (n WGS)
metadata_coinf = metadata[metadata$CLADE=="Omicron(BA.1)/Omicron(BA.2)" | metadata$CLADE=="Delta/Omicron(BA.2)"| metadata$CLADE=="Delta/Omicron(BA.1)",]
metadata_coinf$CLADE = factor(metadata_coinf$CLADE, levels = c("Delta","Omicron(BA.1)","Delta/Omicron(BA.1)","Delta/Omicron(BA.2)","Omicron(BA.1)/Omicron(BA.2)","Omicron(BA.2)","Omicron(BA.3)","Other lineages"))  ## to match colors with Fig 3
																
p1 <- ggplot(metadata_coinf, aes(x=WEEK, group=CLADE))

p1 <- p1 + geom_bar(stat="count", aes(fill = CLADE)) +
    scale_x_discrete(breaks = factor(metadata$WEEK),limits = levels(factor(metadata$WEEK)) ) +## include all weeks
	ylab("n WGS") + xlab("Week") +
	scale_y_continuous(breaks=seq(0,10,2))+
	theme_bw() + theme( axis.text.x  = element_text(vjust=0.5, size=8,angle=45))

p1_npg = p1 + scale_fill_npg(drop=FALSE)
	
## PANEL B. PROPORTIONS OF WHOLE GENOME SEQUENCES (% WGS)
metadata$CLADE = factor(metadata$CLADE, levels = c("Delta","Omicron(BA.1)","Delta/Omicron(BA.1)","Delta/Omicron(BA.2)","Omicron(BA.1)/Omicron(BA.2)","Omicron(BA.2)","Omicron(BA.3)","Other lineages"))  ## to match colors with Fig 3
														
p2 <- ggplot(metadata, aes(x=WEEK, group=factor(CLADE)))

p2 <- p2 + geom_bar(stat="count", aes(fill = CLADE),position="fill") +
	# scale_fill_manual(values=cbPalette) +
	scale_y_continuous(limits=c(0,1),labels = scales::percent) +
	ylab("% WGS") + xlab("Week") +
	theme_bw() + theme( axis.text.x  = element_text(vjust=0.5, size=8,angle=45))
	
p2_npg = p2 + scale_fill_npg(drop=FALSE)

plot_grid(p1_npg, p2_npg, ncol = 1, align = 'v',axis = "rl",rel_heights=c(1,2))
ggsave("Fig4.pdf",width = 7 ,height=7 )
