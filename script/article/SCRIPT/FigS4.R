############################################################################################################################
#' Detection and prevalence of SARS-CoV-2 co-infections during the Omicron variant circulation, France, December 2021 - February 2022

#' Antonin Bal, Bruno Simon, Gregory Destras, Richard Chalvignac, Quentin Semanas, Antoine Oblette, Gregory Queromes, Remi Fanget, 
#' Hadrien Regue, Florence Morfin, Martine Valette, Bruno Lina, Laurence Josset 
#' 
#' (c) 2022 Laurence Josset, <laurence.josset(at)univ-lyon1.fr>
#' 
#' MIT LICENSE

############################################################################################################################


library(dplyr)
require(data.table)
library(ggplot2)
library(reshape)

############################################################################################################################
## load source data
############################################################################################################################

load("SOURCE_Variant_Table_MIXES.Rdata")
load("SOURCE_seq_metrics_MIXES.Rdata")

vcf$nt_mut = paste(vcf$nt_pos,vcf$nt,sep=" ")

## List of Delta and Omicron specific and shared variants based on covariants.org
annot <- read.delim("annot_vcf_DO.txt")
annot$nt_mut = paste(annot$pos.nt,annot$nt.mut,sep=" ")

############################################################################################################################
### EXCLUDE SAMPLES OF RUN 220530 (PERFORMED TO TEST THE DETECTION LIMIT) : To include all samples, comment the 2 following lines
############################################################################################################################
selected_samples = seq_metrics$sample[seq_metrics$RUN != "220530_13146"]

vcf = vcf[is.element(vcf$sample,selected_samples),]

############################################################################################################################
### Fig S1 : ALLELE FREQUENCY PLOTS for ALL MIXES
############################################################################################################################

plot_vcf_all_mixes <- function(titre, vcf_file,annot_file, ncol=11,w=FALSE,h=FALSE) {
	require(data.table)
	library(ggplot2)

	annot=annot_file
	annot$nt_mut = paste(annot$pos.nt,annot$nt.mut,sep=" ")

	##########################
	### annotate vcf
	##########################
	vcf=vcf_file
	
	vcf$VOC <- sapply(vcf$nt_mut, function(x) ifelse(is.element(x,annot$nt_mut),annot$var[x==annot$nt_mut],NA)) ## For each variant, determine whether it is Delta- or OMICRON-specific

	vcf$nt_mut = paste(vcf$nt_pos,vcf$nt,sep="")
	vcf$AA_mut = paste(vcf$AA_pos,vcf$AA,sep="")
	vcf$AA_mut = apply(vcf[,c("prot","AA_mut")],1, function(x) paste(x,collapse=":"))
	
	vcf$nt_AA.mut = paste(vcf$nt_mut,vcf$AA_mut,sep=";")

	vcf$proto = sapply(vcf$ID , function(x) strsplit(x,"x")[[1]][1])
	vcf$proto = gsub("MID","MIDNIGHT ",vcf$proto)
	vcf$proto = gsub("V4","ARTIC V4",vcf$proto)
	vcf$proto = gsub("V41","V4.1",vcf$proto)
	vcf$proto = gsub("MIDNIGHT $","MIDNIGHT V1",vcf$proto)
	vcf$proto = factor(vcf$proto,levels=c("MIDNIGHT V1","MIDNIGHT V2","ARTIC V4","ARTIC V4.1"))
	
	#####Correct ordering of samples of the plot
	levels_samples = unique(vcf$sample)
	levels_samples_Delta = as.numeric(sapply(levels_samples,function(x) strsplit(strsplit(x,"x")[[1]][3],"x")[[1]][1]))	
	levels_samples_ext = sapply(levels_samples,function(x) strsplit(strsplit(x,"x")[[1]][2],"x")[[1]][1])
	vcf$sample = factor(vcf$sample,levels=levels_samples[order(levels_samples_ext,levels_samples_Delta)])

	########## plot AF bars of DELTA-specific and/or OMICRON-specific mutations
	vcfDO = vcf[!is.na(vcf$VOC),]

	###### SIZE of the plot	
	l = length(unique(vcf$ID))
	width = ifelse(!w,l * 1.2,w)
	height = ifelse(!h,l ,h)
	
	# Plot
	p <- ggplot(data=vcfDO, aes(x=nt_pos, y=af, group = sample )) 

	p +
		geom_bar(aes(fill = VOC), stat="identity") + 
		geom_point(data=vcfDO, aes(x=nt_pos, y=af, group = sample,color = VOC),size=5,alpha=0.8)+
		geom_text(data=vcfDO, aes(x=nt_pos, y=af, group = sample ,label=nt_AA.mut), color = "black", vjust=0, hjust=0,size=0.5,angle=45)+
		facet_wrap(proto ~ sample, ncol=ncol) + 
		ylab("Mutation frequency (%)") + xlab("Nucleotide position") +
		scale_color_manual(values=c("tomato","darkgrey","steelblue1"),name = "Mutations", labels = c("Delta-specific", "shared", "Omicron-specific")) +
		scale_fill_manual(values=c("tomato","darkgrey","steelblue1"),name = "Mutations", labels = c("Delta-specific", "shared", "Omicron-specific")) +
		geom_hline(yintercept = 50,linetype =  "dotted", alpha = 0.5, inherit.aes = FALSE) +
		theme_bw() + theme( axis.text.x  = element_text( size=20),axis.text.y  = element_text(size=20),axis.title.x = element_text(size = 50),axis.title.y = element_text(size = 50)) +
		theme(legend.title = element_text(size=50),legend.text = element_text(size=50),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

	ggsave(file = paste(titre,".pdf",sep=""),width = width, height = height,limitsize = FALSE)
	
}


plot_vcf_all_mixes(titre="FigS4", vcf_file=vcf,annot_file=annot, ncol=11)
