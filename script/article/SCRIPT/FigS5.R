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
library(gridExtra)


############################################################################################################################
## load source data
############################################################################################################################
load("SOURCE_seq_metrics_MIXES.Rdata")

## vcf from SEQMET PIPELINE
load("SOURCE_Variant_Table_MIXES.Rdata")
vcf$nt_mut = paste(vcf$nt_pos,vcf$nt,sep=" ")
vcf = vcf[,colnames(vcf)!="depth"]
vcf$pipeline = "SEQMET"

## vcf from DRAGEN PIPELINE
load("SOURCE_Variant_Table_MIXES_DRAGEN.Rdata")
vcf_DRAGEN$nt_mut = paste(vcf_DRAGEN$nt_pos,vcf_DRAGEN$nt,sep=" ")
vcf_DRAGEN$pipeline = "DRAGEN"

## vcf from IVAR PIPELINE
load("SOURCE_Variant_Table_MIXES_IVAR.Rdata")
vcf_IVAR$nt_mut = paste(vcf_IVAR$nt_pos,vcf_IVAR$nt,sep=" ")
vcf_IVAR$pipeline = "IVAR"


## combine all vcf
vcf = rbind(vcf,vcf_DRAGEN,vcf_IVAR)

## List of Delta and Omicron specific and shared variants based on covariants.org
annot <- read.delim("annot_vcf_DO.txt")
annot$nt_mut = paste(annot$pos.nt,annot$nt.mut,sep=" ")


############################################################################################################################
### EXCLUDE SAMPLES OF RUN 220530 (PERFORMED TO TEST THE DETECTION LIMIT) : To include all samples, comment the 2 following lines
############################################################################################################################
selected_samples = seq_metrics$sample[seq_metrics$RUN != "220530_13146"]

vcf = vcf[is.element(vcf$sample,selected_samples),]

############################################################################################################################
### Fig S5 : number of secondary lineage-specific polymorphism in consensus sequence
############################################################################################################################

## For each variant, determine whether it is a Delta- or OMICRON-specific polymorphism :
vcf$VOC <- sapply(vcf$nt_mut, function(x) ifelse(is.element(x,annot$nt_mut),annot$var[x==annot$nt_mut],NA)) 

## Among MAJOR variants , how many samples have both Delta- and Omicron-specific mutations (excluding shared positions) :
vcf_MAJ_DO = vcf[is.element(vcf$VOC,c("O","D")) & vcf$type == "MAJOR",]

nb_Omicron_Delta_MAJ = vcf_MAJ_DO %>%
  group_by(sample,pipeline,VOC) %>%
  dplyr::summarise(n = n())
  
nb_MAJ = cast(nb_Omicron_Delta_MAJ,  sample + pipeline ~ VOC, value.var = "n")

nb_MAJ[is.na(nb_MAJ)] <- 0

## identify the lineage with lowest number of specific mutations within consensus
nb_MAJ$n_min_consensus = apply(nb_MAJ,1,function(x) min(x[c("D","O")]))

nb_MAJ$min_consensus = apply(data.frame(nb_MAJ),1,function(x) ifelse(x["D"]<x["O"],"B.1.617.2","BA.1"))
nb_MAJ$min_consensus[nb_MAJ$n_min_consensus==0] <- "None"

nb_MAJ$is_chimeric <- nb_MAJ$n_min_consensus >0

### Annotate the file (Primers, Extraction, % Delta, % Omicron of mixes)
nb_MAJ$ID = sapply(nb_MAJ$sample, function(x)  strsplit(x,"_")[[1]][1])

nb_MAJ$proto = sapply(nb_MAJ$ID,function(x) strsplit(strsplit(x,"-")[[1]][2],"x")[[1]][1])
nb_MAJ$proto = gsub("MID","MIDNIGHT ",nb_MAJ$proto)
nb_MAJ$proto = gsub("V4","ARTIC V4",nb_MAJ$proto)
nb_MAJ$proto = gsub("V41","V4.1",nb_MAJ$proto)
nb_MAJ$proto = gsub("MIDNIGHT $","MIDNIGHT V1",nb_MAJ$proto)

nb_MAJ$proto = factor(nb_MAJ$proto,levels=c("MIDNIGHT V1","MIDNIGHT V2","ARTIC V4","ARTIC V4.1"))

nb_MAJ$ext = sapply(nb_MAJ$ID,function(x) strsplit(strsplit(x,"-")[[1]][2],"x")[[1]][2])
nb_MAJ$Delta = sapply(nb_MAJ$ID,function(x) strsplit(strsplit(x,"-")[[1]][2],"x")[[1]][3])
nb_MAJ$Omicron = sapply(nb_MAJ$ID,function(x) strsplit(strsplit(strsplit(x,"-")[[1]][2],"x")[[1]][4],"_")[[1]][1])
nb_MAJ$Omicron = sapply(nb_MAJ$Omicron,function(x) strsplit(x,"X")[[1]][1])

## order mixes for plot
nb_MAJ$Delta = factor(nb_MAJ$Delta, levels = unique(nb_MAJ$Delta)[order(as.numeric(unique(nb_MAJ$Delta)))])

### Plot Fig S5
p <- ggplot(data=nb_MAJ, aes(x=Delta, y=n_min_consensus)) 

p <-	p +
		geom_point(aes(colour=min_consensus,shape=is_chimeric),alpha=0.9) +
		facet_grid(pipeline~proto) + 
		ylab("Number of contaminant lineage specific mutations in the consensus") + xlab("Delta:Omicron mixes") +
		scale_x_discrete(breaks=c(0,1,5,10,20,30,40,50,60,70,80,90,100),
                      labels = c("0:100","1:99","5:95","10:90","20:80", "30:70","40:60","50:50","60:40","70:30", "80:20","90:10","100:0")) +
		scale_color_manual(values = c( "tomato","steelblue1","grey"),name = "Contaminant lineage", labels = c("B.1.617.2","BA.1","None")) +
		scale_shape_manual(values = c(17,16), name = "Chimeric consensus sequence") +
		theme_bw() + theme( axis.text.x  = element_text(vjust=1, angle=90, hjust=1)) 

p
ggsave(file="FigS5.pdf",width = 10 ,height=6 )

