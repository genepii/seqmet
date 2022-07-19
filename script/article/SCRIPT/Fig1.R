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

## List of Delta and Omicron specific and shared variants based on covariants.org
annot <- read.delim("annot_vcf_DO.txt")
annot$nt_mut = paste(annot$pos.nt,annot$nt.mut,sep=" ")


########################################################################################################
### EXCLUDE SAMPLES OF RUN 220530 (PERFORMED TO TEST THE DETECTION LIMIT) : To include all samples, comment the 2 following lines
########################################################################################################

selected_samples = seq_metrics$sample[seq_metrics$RUN != "220530_13146"]

vcf = vcf[is.element(vcf$sample,selected_samples),]

########################################################################################################
### Fig 1.A : Ratio of Delta- and Omicron-specific mutations
########################################################################################################

vcf$nt_mut = paste(vcf$nt_pos,vcf$nt,sep=" ")

##### Select Delta- and Omicron-specific variants (i.e. Exclude variants shared between Delta and Omicron (labeled "DO"))
Delta_pos = annot$nt_mut[annot$var == "D"]
Omicron_pos = annot$nt_mut[annot$var == "O"]

### Select only specific-mutations present in each sample 
vcf_Delta = vcf[is.element(vcf$nt_mut,Delta_pos) ,]
vcf_Omicron = vcf[is.element(vcf$nt_mut,Omicron_pos) ,]

##### Calculate the ratio of Delta- and Omicron-specific variants found in each sample
nb_Delta = dplyr::summarise(group_by(vcf_Delta,sample),
		n = n())

nb_Omicron = dplyr::summarise(group_by(vcf_Omicron,sample),
		n = n())

nb_Omicron_Delta = merge(nb_Delta,nb_Omicron,by="sample",all=TRUE) ##TRUE ici on veut aussi garder les non coinf pour le plot!!!
nb_Omicron_Delta[is.na(nb_Omicron_Delta)] <- 0
colnames(nb_Omicron_Delta) = c("sample","nb_Delta","nb_Omicron")
	
nb_Omicron_Delta$ratio_Delta = nb_Omicron_Delta$nb_Delta / length(Delta_pos)
nb_Omicron_Delta$ratio_Omicron = nb_Omicron_Delta$nb_Omicron / length(Omicron_pos)

### Annotate the file (Primers, Extraction, % Delta, % Omicron of mixes)
nb_Omicron_Delta$ID = sapply(nb_Omicron_Delta$sample, function(x)  strsplit(x,"_")[[1]][1])

nb_Omicron_Delta$proto = sapply(nb_Omicron_Delta$ID,function(x) strsplit(strsplit(x,"-")[[1]][2],"x")[[1]][1])
nb_Omicron_Delta$proto = gsub("MID","MIDNIGHT ",nb_Omicron_Delta$proto)
nb_Omicron_Delta$proto = gsub("V4","ARTIC V4",nb_Omicron_Delta$proto)
nb_Omicron_Delta$proto = gsub("V41","V4.1",nb_Omicron_Delta$proto)
nb_Omicron_Delta$proto = gsub("MIDNIGHT $","MIDNIGHT V1",nb_Omicron_Delta$proto)

nb_Omicron_Delta$proto = factor(nb_Omicron_Delta$proto,levels=c("MIDNIGHT V1","MIDNIGHT V2","ARTIC V4","ARTIC V4.1"))

nb_Omicron_Delta$ext = sapply(nb_Omicron_Delta$ID,function(x) strsplit(strsplit(x,"-")[[1]][2],"x")[[1]][2])
nb_Omicron_Delta$Delta = sapply(nb_Omicron_Delta$ID,function(x) strsplit(strsplit(x,"-")[[1]][2],"x")[[1]][3])

nb_Omicron_Delta$Coinf = nb_Omicron_Delta$ratio_Delta > 0 & nb_Omicron_Delta$ratio_Omicron > 0

## PLOT FIG 1.A
nb_Omicron_Delta$Delta = as.numeric(nb_Omicron_Delta$Delta)
p1 <- ggplot(data=nb_Omicron_Delta, aes(x=ratio_Omicron, y=ratio_Delta ,colour=Delta)) 

p1 <-	p1 +
		geom_point(aes(colour=Delta,shape=Coinf),alpha=0.9) +
		facet_grid(.~proto) + 
		ylab("Delta-specific mutations detection rate") + xlab("Omicron-specific mutations detection rate") +
		scale_color_gradient2(low = "steelblue1", mid = "cyan4", high = "tomato",
			breaks=c(0,10,20,30,40,50,60,70,80,90,100), midpoint = 50,
			name = "Delta:Omicron", labels = c("0:100","10:90","20:80", "30:70","40:60","50:50","60:40","70:30", "80:20","90:10","100:0")) +
		scale_shape(name = "Experimental\nCoinfection") +
		geom_vline(xintercept = 0.25,linetype =  "dotted", alpha = 0.5, inherit.aes = FALSE) +
		geom_hline(yintercept = 0.90,linetype =  "dotted", alpha = 0.5, inherit.aes = FALSE) +
		theme_bw() 


########################################################################################################
### Fig 1.B : Measured vs expected frequency of Delta and Omicron in mixes
########################################################################################################

## CALCULATE MEDIAN AF + IQR of Delta- and Omicron-specific variants
med_Delta = dplyr::summarise(group_by(vcf_Delta,sample),
		med = median(af),
		iqr=IQR(af)
		)
boxplot(med_Delta$med)

med_Omicron = dplyr::summarise(group_by(vcf_Omicron,sample),
		med = median(af),
		iqr=IQR(af)
		)
boxplot(med_Omicron$med)


med_Omicron_Delta = merge(med_Delta,med_Omicron,by="sample",all=TRUE) ##TRUE : keep non co-infected mixes

colnames(med_Omicron_Delta) = c("sample","median_Delta","IQR_Delta","median_Omicron","IQR_Omicron")
med_Omicron_Delta$sample <- sapply(med_Omicron_Delta$sample, function(x) strsplit(x,"\tM")[[1]][1])

table_MAF = merge(nb_Omicron_Delta,med_Omicron_Delta,by="sample",all=TRUE)
table_MAF[is.na(table_MAF)] <- 0
save(table_MAF,file="table_MAF_MIXES.Rdata") ## saved for Table S2

### Plot Fig 1.B
p2 <- ggplot(data=table_MAF, aes(x=Delta, y=median_Delta)) 

p2 <-	p2 +
		geom_smooth(color="grey", fill="grey", linetype="blank") +
		geom_point(aes(colour=Delta,shape=Coinf)) +
		facet_grid(.~proto) + 
		ylab("Measured Frequency (%)") + xlab("Expected Frequency (%)") +
		geom_abline(intercept = 0, slope = 1,linetype =  "dotted", alpha = 0.5, inherit.aes = FALSE) +
		scale_color_gradient2(low = "steelblue1", mid = "cyan4", high = "tomato",
			breaks=c(0,10,20,30,40,50,60,70,80,90,100), midpoint = 50,
			name = "Delta:Omicron", labels = c("0:100","10:90","20:80", "30:70","40:60","50:50","60:40","70:30", "80:20","90:10","100:0")) +
		scale_shape(name = "Experimental\nCoinfection") +
		theme_bw() 

########################################################################################################
### Fig 1.C : ALLELE FREQUENCY PLOT for 20:80 MIX
########################################################################################################

plot_AF_fig1C <- function(vcf_file,annot_file, ncol=11) {
	require(data.table)
	library(ggplot2)

	annot=annot_file

	##########################
	### annotate vcf
	##########################
	vcf=vcf_file
	
	vcf$VOC <- sapply(vcf$nt_mut, function(x) ifelse(is.element(x,annot$nt_mut),annot$var[x==annot$nt_mut],NA)) ## For each variant, determine whether it is Delta- or OMICRON-specific

	vcf$proto = sapply(vcf$ID , function(x) strsplit(x,"x")[[1]][1])
	vcf$proto = gsub("MID","MIDNIGHT ",vcf$proto)
	vcf$proto = gsub("V4","ARTIC V4",vcf$proto)
	vcf$proto = gsub("V41","V4.1",vcf$proto)
	vcf$proto = gsub("MIDNIGHT $","MIDNIGHT V1",vcf$proto)
	vcf$proto = factor(vcf$proto,levels=c("MIDNIGHT V1","MIDNIGHT V2","ARTIC V4","ARTIC V4.1"))

	########## plot AF bars of DELTA-specific and/or OMICRON-specific mutations
	vcfDO = vcf[!is.na(vcf$VOC),]

	# Plot
	p <- ggplot(data=vcfDO, aes(x=nt_pos, y=af, group = sample )) 

	p <- p +
		geom_bar(aes(fill = VOC), stat="identity") + 
		geom_point(data=vcfDO[vcfDO$af>50,], aes(x=nt_pos, y=af, group = sample,color = VOC),size=1,alpha=0.8)+
		facet_wrap(.~proto, ncol=ncol) + 
		ylab("Mutation frequency (%)") + xlab("Nucleotide position") +
		scale_color_manual(values=c("tomato","darkgrey","steelblue1"),name = "Mutations in\nconsensus", labels = c("Delta-specific", "shared", "Omicron-specific")) +
		scale_fill_manual(values=c("tomato","darkgrey","steelblue1"),name = "Mutations", labels = c("Delta-specific", "shared", "Omicron-specific")) +
		geom_hline(yintercept = 50,linetype =  "dotted", alpha = 0.5, inherit.aes = FALSE) +
		theme_bw() + theme( axis.text.x  = element_text(angle=45, size=8,hjust =1, vjust=1)) +
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

	return(p)
	
}


p3 <- plot_AF_fig1C( vcf_file=vcf[grep("EMAGx20x80",vcf$sample),],annot_file=annot, ncol=4)


########################################################################################################
### SAVE FIG 1AB
########################################################################################################

plot_grid(p1, p2,p3, ncol = 1,labels = c("A","B","C"), align = 'v',axis = "rl",rel_heights=c(1,1,1.1))
ggsave("Fig1.pdf",width = 9 ,height=8 )
ggsave("Fig1.png",width = 9 ,height=8 )
