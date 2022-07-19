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
library(tidyr)
library(ggplot2)
library(cowplot)
library(data.table)

############################################################################################################################
## load source data for experimental mixes
############################################################################################################################
load("SOURCE_seq_metrics_MIXES.Rdata")

load("SOURCE_Variant_Table_MIXES.Rdata") 
vcf$nt_mut = paste(vcf$nt_pos,vcf$nt,sep=" ")

load("SOURCE_depth_MIXES.Rdata")


############################################################################################################################
## List of Delta and Omicron specific and shared variants based on covariants.org
############################################################################################################################

annot_R <- read.delim("annot_vcf_DO.txt")
annot_R$nt_mut = paste(annot_R$pos.nt,annot_R$nt.mut,sep=" ")

### Select only Delta and Omicron SPECIFIC mutations (excluding shared positions)
annot_R = annot_R[annot_R$var!="DO",]

############################################################################################################################
##  SEQMET-db : list of defining polymorphism for each lineage
############################################################################################################################
annot_S <-as.data.frame(fread("profile_lineages.tsv"))

annot_S$nt_mut = paste(annot_S$pos.nt,annot_S$nt.mut,sep=" ")

### Select only BA.1 and B.1.617.2 SPECIFIC mutations (excluding shared positions)
annot_S = annot_S[annot_S$var == "BA_1"|annot_S$var == "B_1_617_2",]
annot_S = annot_S[!is.element(annot_S$nt_mut,annot_S$nt_mut[duplicated(annot_S$nt_mut)]),]

############################################################################################################################
## Compare annotations
############################################################################################################################
annot_R$source = "covariant"
annot_S$source = "seqmetdb"

annot = rbind(annot_R,annot_S)

### annotations in seqmet that were also in covariants
annot$source[is.element(annot$nt_mut,annot$nt_mut[duplicated(annot$nt_mut)])] <- "covariant&seqmetdb"

annot = annot[!duplicated(annot$nt_mut),]

############################################################################################################################
### EXCLUDE SAMPLES OF RUN 220530 (PERFORMED TO TEST THE DETECTION LIMIT) : To include all samples, comment the 2 following lines
############################################################################################################################
selected_samples = seq_metrics$sample[seq_metrics$RUN != "220530_13146"]

vcf = vcf[is.element(vcf$sample,selected_samples),]


############################################################################################################################
### Which Omicron-specific mutations are missed in each sample?
############################################################################################################################

## Subset vcf file to only Delta and/or Omicron-defining mutations => vcfDO

vcf$VOC <- sapply(vcf$nt_mut, function(x) ifelse(is.element(x,annot$nt_mut),annot$var[x==annot$nt_mut],NA)) ## For each variant, determine whether it is Delta- or OMICRON-specific
vcfDO = vcf[!is.na(vcf$VOC),]

##STOP if all defining mutations not found at least once in vcf
stopifnot(length(setdiff(annot$nt_mut,vcfDO$nt_mut))==0)

vcfDO$VOC_nt_mut = paste(vcfDO$VOC,vcfDO$nt_mut,sep="_")

table_S3 = spread(vcfDO[,c("sample","VOC_nt_mut","af")], key = VOC_nt_mut, value = af)
rownames(table_S3) = table_S3$sample
table_S3 = table_S3[,colnames(table_S3)!="sample"]

table_S3[is.na(table_S3)] <- 0
colnames(table_S3) = gsub("D|B_1_617_2","B.1.617.2",colnames(table_S3))
colnames(table_S3) = gsub("O|BA_1","BA.1",colnames(table_S3))

############################################################################################################################
### Differentiate mutations which are absent in the sample (=> af should be set to 0) from not covered positions (af should be set to NA)
### Indeed the vcf file list all variants and their allele frequency found for each sample. An absent variant may be : absent in the sample (ie wt position) OR not covered with enough reads to be called.
############################################################################################################################


nt_pos = as.numeric(sapply(colnames(table_S3),function(x) strsplit(strsplit(x,"_")[[1]][2]," ")[[1]][1]))

## determine not covered position to put NA in table_S3 :  
not_cov = depth[depth$depth<10,]

for (s in rownames(table_S3)) {
				print(s)
				position_to_mask = not_cov[not_cov$sample==s,"nt_pos"]
				# print(nt_pos[is.element(nt_pos,position_to_mask)])
				if(length(nt_pos[is.element(nt_pos,position_to_mask)])>0) {
					table_S3[s,is.element(nt_pos,position_to_mask)] <- NA
					print(nt_pos[is.element(nt_pos,position_to_mask)])
					}
				rm(position_to_mask)
				}


############################################################################################################################
### Plot heatmap (HM)
############################################################################################################################
library(gplots)
matHM = t(table_S3)
matHM[is.na(colnames(matHM))] <- 0
# colors
hmcols <- colorRampPalette(c("white","black"))(100)

# Order column based on  expected Delta frequency in the mix + colors : function expected Omicron frequency in the mix
colnames(matHM) = sapply(colnames(matHM),function(x) strsplit(x,"-")[[1]][2])

Delta = as.numeric(sapply(colnames(matHM),function(x) strsplit(strsplit(x,"x")[[1]][3],"_")[[1]][1]))
v =  as.numeric(as.factor(Delta))
matHM=matHM[,order(v)]
Delta = as.numeric(sapply(colnames(matHM),function(x) strsplit(strsplit(x,"x")[[1]][3],"_")[[1]][1]))
v =  as.numeric(as.factor(Delta))
col_column = sapply(v, function(x) colorRampPalette(c("steelblue1","cyan4","tomato"))(max(v)) [x])

# row class colors : Delta, and Omicron 
r = as.numeric(as.factor(sapply(rownames(matHM),function(x) strsplit(x,"_")[[1]][1])))
col_row = sapply(r, function(x) c("tomato","steelblue1")[x])

h=heatmap.2(as.matrix(matHM), Colv = F, dendrogram="row", scale="none" ,breaks=seq(0,100,1),col=hmcols,trace="none",cexCol=0.3, cexRow=0.3,margin=c(15, 15),
density.info="none",keysize=1,ColSideColors = col_column, RowSideColors = col_row, na.rm=FALSE, na.color = "yellow")

# dev.copy2pdf(file="FigS4_HM_yellow.pdf")


############################################################################################################################
### SAVE TABLE S3 with mutations ordered as in the HM
############################################################################################################################
t3 = as.data.frame(matHM[rev(h$rowInd),])

## add a column indicating if the mutations is found in covariants, seqmet-db or in both lists
t3$source = as.vector(sapply(rownames(t3), function(x){
	lab = strsplit(x,"_")[[1]][2]
	annot$source[annot$nt_mut==lab]
	}))
write.table(t3,file="table_S3.txt", sep="\t", quote=F,col.names=NA)

