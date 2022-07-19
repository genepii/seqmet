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
#### CODE USED TO GENERATE FIG S2
############################################################################################################################

options(stringsAsFactors = FALSE) 

library(reshape)
library(plyr)
library(ggplot2)
library(data.table)
library(gridExtra)
library(reshape2)

############################################################################################################################
## DATA SOURCE : SEQUENCING METRICS (seq_metrics) and VARIANT TABLE for n=21448 sequencing results (21387 samples + 61 duplicate sequencing) 
############################################################################################################################
load("SOURCE_Variant_Table.Rdata")

load("FORMATED_seq_metrics.Rdata")


############################################################################################################################
## DATA SOURCE BIS: SEQUENCING METRICS (seq_metrics_bis) for duplicates of n=7 samples with coinf_min_common = 4 or 5 at first replicate
############################################################################################################################
load("SOURCE_seq_metrics_bis.Rdata")
load("SOURCE_Variant_Table_bis.Rdata")


############################################################################################################################
##  SEQMET-db : list of defining polymorphism for each lineage
############################################################################################################################
annot <-as.data.frame(fread("profile_lineages.tsv"))

annot$nt_mut = paste(annot$pos.nt,annot$nt.mut,sep=" ")

############################################################################################################################
## PLOT ALLELE FREQUENCY FONCTION 
############################################################################################################################

plot_AF <- function(id, vcf_file,annot_file, seq_metrics,colors_bar = c("tomato","steelblue1")) {
	library(ggplot2)

	annot=annot_file
	
	### identify main and secondary lineages for the sample : 
	### for each sample choose main and secondary lineages based on the most coinfected duplicate (coinf_min_ratio maximum) : this is to have same mutations between duplicates to compare discordant samples
	seq_metrics_id = seq_metrics[seq_metrics$ID == id,]
	
	virus = as.vector(unlist(seq_metrics_id[which.max(seq_metrics_id$coinf_min_ratio),c("coinf_maj_match","coinf_min_match")]))
	
	### subset annotation file to only secondary and main lineages
	annot=annot[is.element(annot$var,virus),]
	annot$nt_mut = paste(annot$pos.nt,annot$nt.mut,sep="")
	
	## filter-out shared mutations
	annot = annot[!is.element(annot$nt_mut,annot$nt_mut[duplicated(annot$nt_mut)]),]

	##########################
	### annotate vcf
	##########################
	vcf=vcf_file[is.element(vcf_file$ID,id),]
	
	vcf$nt_mut = paste(vcf$nt_pos,vcf$nt,sep="")
	
	vcf$VOC <- sapply(vcf$nt_mut, function(x) ifelse(is.element(x,annot$nt_mut),annot$var[x==annot$nt_mut],NA)) ## For each variant, determine whether it is Delta- or OMICRON-specific

	vcf$AA_mut = paste(vcf$AA_pos,vcf$AA,sep="")
	vcf$AA_mut = apply(vcf[,c("prot","AA_mut")],1, function(x) paste(x,collapse=":"))
	
	vcf$nt_AA_mut = paste(vcf$nt_mut,vcf$AA_mut,sep=";")

	########## Select only specific mutations
	vcfDO = vcf[!is.na(vcf$VOC),]

	# Plot
	p <- ggplot(data=vcfDO, aes(x=nt_pos, y=af, group = sample )) 

	p <- p +
		geom_bar(aes(fill = VOC), stat="identity") + 
		geom_point(data=vcfDO, aes(x=nt_pos, y=af, group = sample,color = VOC),size=3,alpha=0.8)+
		geom_text(data=vcfDO, aes(x=nt_pos, y=af, group = sample ,label=nt_AA_mut, color = VOC), vjust=0, hjust=0,size=0.5,angle=45)+
		facet_grid( sample ~ID) + 
		ylab("Mutation frequency (%)") + xlab("Nucleotide position") +
		scale_color_manual(values=colors_bar,guide = "none") +
		scale_fill_manual(values=colors_bar,name = "Specific Mutations") +
		geom_hline(yintercept = c(5,50),linetype =  "dotted", alpha = 0.5, inherit.aes = FALSE) +
		theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	return(p)
}

############################################################################################################################
### SELECT SAMPLES 
############################################################################################################################
seq_metrics_cutoff4_5 = seq_metrics[seq_metrics$coinf_min_common==4 | seq_metrics$coinf_min_common==5,]
borderline_samples_ID = seq_metrics_cutoff4_5$ID

seq_metrics_cutoff4_5 = rbind(seq_metrics_cutoff4_5[,colnames(seq_metrics_bis)],seq_metrics_bis)

vcf = rbind(vcf,vcf_bis)



## For each sample : plot allele frequencies for the main and secondary lineages
list_plots <- lapply(as.vector(borderline_samples_ID), function(x) plot_AF(id = x, vcf_file=vcf,annot_file = annot, seq_metrics,colors_bar = c("tomato","steelblue1")))

#############################################################################################################################
### SAVE FIG S2
############################################################################################################################

pdf(file="FigS2.pdf",width = 15 ,height=12 )
grid.arrange(grobs =list_plots, ncol = 2)
dev.off()
