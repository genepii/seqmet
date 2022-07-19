############################################################################################################################
#' Detection and prevalence of SARS-CoV-2 co-infections during the Omicron variant circulation, France, December 2021 - February 2022

#' Antonin Bal, Bruno Simon, Gregory Destras, Richard Chalvignac, Quentin Semanas, Antoine Oblette, Gregory Queromes, Remi Fanget, 
#' Hadrien Regue, Florence Morfin, Martine Valette, Bruno Lina, Laurence Josset 
#' 
#' (c) 2022 Laurence Josset, <laurence.josset(at)univ-lyon1.fr>
#' 
#' MIT LICENSE

############################################################################################################################

options(stringsAsFactors = FALSE) 
library(dplyr)
library(ggplot2)
library(reshape)
library("gridExtra")

library(data.table)

############################################################################################################################
## load source data
############################################################################################################################
load("SOURCE_seq_metrics_MIXES.Rdata")

load("SOURCE_Variant_Table_MIXES.Rdata")
vcf$nt_mut = paste(vcf$nt_pos,vcf$nt,sep=" ")


############################################################################################################################
##  SEQMET-db : list of defining polymorphism for each lineage
############################################################################################################################
annot <-as.data.frame(fread("profile_lineages.tsv"))

annot$nt_mut = paste(annot$pos.nt,annot$nt.mut,sep=" ")

### Select only BA.1 and B.1.617.2 SPECIFIC mutations (excluding shared positions)
annot = annot[annot$var == "BA_1"|annot$var == "B_1_617_2",]
annot = annot[!is.element(annot$nt_mut,annot$nt_mut[duplicated(annot$nt_mut)]),]


############################################################################################################################
### INCLUDE ONLY SAMPLES OF RUN 220530 (PERFORMED TO TEST THE DETECTION LIMIT) 
############################################################################################################################
selected_samples = seq_metrics$sample[seq_metrics$RUN == "220530_13146"]

vcf = vcf[is.element(vcf$sample,selected_samples),]

############################################################################################################################
### Fig S6 : number of secondary lineage-specific polymorphism in consensus sequence
############################################################################################################################

## For each variant, determine whether it is a Delta- or OMICRON-specific polymorphism :
vcf$VOC <- sapply(vcf$nt_mut, function(x) ifelse(is.element(x,annot$nt_mut),annot$var[x==annot$nt_mut],NA)) 

## Among MAJOR variants , how many samples have both Delta- and Omicron-specific mutations (excluding shared positions) :
vcf_MAJ_DO = vcf[is.element(vcf$VOC,c("BA_1","B_1_617_2")) & vcf$type == "MAJOR",]

nb_Omicron_Delta_MAJ = vcf_MAJ_DO %>%
  group_by(sample,VOC) %>%
  dplyr::summarise(n = n())
  
nb_MAJ = cast(nb_Omicron_Delta_MAJ,  sample ~ VOC, value.var = "n")

nb_MAJ[is.na(nb_MAJ)] <- 0

## identify whether consensus sequence is Delta or Omicron (based on presence of majority variants)
nb_MAJ$n_min_consensus = apply(nb_MAJ,1,function(x) min(x[c("B_1_617_2","BA_1")]))

nb_MAJ$min_consensus = apply(data.frame(nb_MAJ),1,function(x) ifelse(x["B_1_617_2"]<x["BA_1"],"B.1.617.2","BA.1"))
nb_MAJ$min_consensus[nb_MAJ$n_min_consensus==0] <- "None"

## Add coinf_min_ratio from the seq_metrics file
nb = merge(nb_MAJ,seq_metrics,by.x="sample",by.y="sample")

## Coinfection detected using seqmet?
nb$Positive_coinf_min_ratio <- nb$coinf_min_ratio >0 

## order mixes for plot
nb$Delta = factor(nb$Delta, levels = unique(nb$Delta)[order(as.numeric(unique(nb$Delta)))])

### Plot Fig S2
p <- ggplot(data=nb, aes(x=Delta, y=n_min_consensus)) 

p <-	p +
		geom_jitter(aes(colour=min_consensus,shape=Positive_coinf_min_ratio),alpha=0.9, width = 0.25, height = 0) +
		facet_grid(.~proto) + 
		facet_grid(.~proto) + 
		ylab("Number of contaminant lineage\nspecific mutations in the consensus") + xlab("Delta:Omicron mixes") +
		scale_x_discrete(breaks=c(1,5,10),
                      labels = c("1:99","5:95","10:90")) +
		scale_color_manual(values = c( "tomato","steelblue1","grey"),name = "Contaminant lineage", labels = c("B.1.617.2","BA.1","None")) +
		scale_shape_manual(values = c(8, 2),name = "Secondary lineage\nmutation ratio", labels = c("Null","Positive")) +
		theme_bw() + theme( axis.text.x  = element_text(vjust=1, angle=90, hjust=1)) 

p
ggsave(file="FigS6.pdf",width = 9 ,height= 4 )

