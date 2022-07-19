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
#### CODE USED TO GENERATE FIG 3 AND DEFINE CO-INFECTIONS FROM SEQUENCING RESULTS
############################################################################################################################

options(stringsAsFactors = FALSE) 

# library(reshape)
library(plyr)
library(ggplot2)
library(ggsci)
library(cowplot)

library(reshape2)

library(dplyr)
library(data.table)

############################################################################################################################
## DATA SOURCE : SEQUENCING METRICS (seq_metrics) and VARIANT TABLE for n=21448 sequencing results (21387 samples + 61 duplicate sequencing) 
############################################################################################################################
load("SOURCE_Variant_Table.Rdata")

load("FORMATED_seq_metrics.Rdata")

############################################################################################################################
##  SEQMET-db : list of defining polymorphism for each lineage
############################################################################################################################
annot <-as.data.frame(fread("profile_lineages.tsv"))

############################################################################################################################
#### FIG S8 : A. Reproducibility of coinf_min_ratio for 61 duplicates 
############################################################################################################################
## Select duplicated samples
seq_metrics_dup = seq_metrics[is.element(seq_metrics$ID,seq_metrics$ID[duplicated(seq_metrics$ID)]) ,]

# for samples sequenced in duplicate, determine which is the first and second Replicate based on Year of sequencing, and Plate number (Pl)
first_Replicate = ddply(seq_metrics_dup, .(ID),function(x) x[order(x$Year,x$pl),][1,])

seq_metrics_dup$Replicate = NA
seq_metrics_dup$Replicate[is.element(seq_metrics_dup$sample,first_Replicate$sample)] = "First"
seq_metrics_dup$Replicate[!is.element(seq_metrics_dup$sample,first_Replicate$sample)] = "Second"


# prepare data for ggplot2
figSA <- reshape2::dcast(seq_metrics_dup[,c("ID","coinf_min_ratio","COINF","Replicate")],ID + COINF  ~ Replicate ,value.var = "coinf_min_ratio")
 
figSA$matches = sapply(figSA$ID, function(x) seq_metrics_dup$matches[seq_metrics_dup$ID==x & seq_metrics_dup$Replicate=="Second"])## keep only 2nd Replicate as it is what was validated for discordant Replicates
figSA$matches = gsub(" $","",figSA$matches)
figSA$matches = factor(figSA$matches,levels = c("Delta","Omicron(BA.1)","Delta/Omicron(BA.1)","Delta/Omicron(BA.2)","Omicron(BA.1)/Omicron(BA.2)","Omicron(BA.2)","Omicron(BA.3)","Other lineages"))  ## to match colors with fig 3 and 4

pSA <- ggplot(figSA, aes(x=First, y=Second))
	
pSA <- pSA + 
	geom_point(aes(shape=COINF,color=matches),alpha=0.9) +
	xlim(c(0,1)) + ylim(c(0,1))+
	xlab("Secondary lineage mutation ratio\nin First Replicate") + ylab("Secondary lineage mutation ratio\nin Second Replicate") +
	geom_abline(intercept = 0, slope = 1,linetype =  "dotted", alpha = 0.5, inherit.aes = FALSE) +
	scale_shape(name = "Natural Coinfection") +
	theme_bw() + scale_color_npg(drop=FALSE)


#### Comparison of secondary mutation ratio in first replicate among confirmed and refuted coinfection
Ratio_not_coinf = figSA$First[!figSA$COINF]
Ratio_coinf = figSA$First[figSA$COINF]
median (Ratio_not_coinf)
median (Ratio_coinf)
kruskal.test(list(Ratio_not_coinf,Ratio_coinf))

###########################################################################################################
####  FIG S8 : B Distribution of relative abundance for the ‘minor virus’ in confirmed vs refuted coinfections
############################################################################################################

compute_rel_abund <- function(s, vcf_file,annot_file, seq_metrics) {
	require(dplyr)

	annot=annot_file
	
	### identify main and secondary lineages for the sample : 
	virus = as.vector(unlist(seq_metrics[seq_metrics$sample == s,c("coinf_maj_match","coinf_min_match")]))
	
	### subset annotation file to only secondary and main lineages
	annot=annot[is.element(annot$var,virus),]
	annot$nt_mut = paste(annot$pos.nt,annot$nt.mut,sep="")
	
	## filter-out shared mutations
	annot = annot[!is.element(annot$nt_mut,annot$nt_mut[duplicated(annot$nt_mut)]),]

	##########################
	### select vcf
	##########################
	vcf=vcf_file[is.element(vcf$sample,s),]
	
	vcf$nt_mut = paste(vcf$nt_pos,vcf$nt,sep="")
	
	vcf$VOC <- sapply(vcf$nt_mut, function(x) ifelse(is.element(x,annot$nt_mut),annot$var[x==annot$nt_mut],NA)) 
	
	########## Select only mutations specific of main and secondary lineage
	vcf_m_sec = vcf[is.element(vcf$VOC,virus),]

	#### return Median AF + IQR and number of specific mutations within the consensus for main and secondary lineages
	RelAb = vcf_m_sec %>%
	  group_by(sample,VOC) %>%
	  dplyr::summarise(
		med = median(af),
		iqr=IQR(af),
		coeff_variation = IQR(af)/median(af)
		)
		
	return(list(RelAb))

}

## select only samples sequenced in duplicates
dup = seq_metrics$sample[is.element(seq_metrics$ID,seq_metrics$ID[duplicated(seq_metrics$ID)])]

RelAb <- sapply(dup, function(s) compute_rel_abund(s, vcf,annot, seq_metrics))

RelAb_table = do.call(rbind,RelAb)


## identify lineage with lowest relative abundance = MINOR LINEAGE (/!\ this is not always the secondary lineage which is identified by the lowest ratio of present mutations)
RelAb_min = RelAb_table %>%
		group_by(sample) %>%
		slice(which.min(med))




## determine whether the least abundant lineage is also the secondary lineage (based on ratio of present mutation)
seq_metrics_dup$is_lowest_secondary <- seq_metrics_dup$lowest_lineage == seq_metrics_dup$coinf_min_match

seq_metrics_dup$lowest_lineage = sapply(seq_metrics_dup$sample, function(x) ifelse(length(RelAb_min$VOC[RelAb_min$sample==x])==1,RelAb_min$VOC[RelAb_min$sample==x],NA))
seq_metrics_dup$lowest_lineage_RelAb = sapply(seq_metrics_dup$sample, function(x) ifelse(length(RelAb_min$med[RelAb_min$sample==x])==1,RelAb_min$med[RelAb_min$sample==x],NA))


# prepare data for ggplot2
seq_metrics_dup$lowest_lineage_RelAb[seq_metrics_dup$coinf_min_ratio==0] = 0 ## because when no secondary lineage there is only one RelAb which is the RelAb of the main lineage

cdata <- reshape2::dcast(seq_metrics_dup[,c("ID","lowest_lineage_RelAb","COINF","Replicate")],ID + COINF  ~ Replicate ,value.var = "lowest_lineage_RelAb")
 
cdata$matches = sapply(cdata$ID, function(x) seq_metrics_dup$matches[seq_metrics_dup$ID==x & seq_metrics_dup$Replicate=="Second"])## keep only 2nd Replicate as it is what was validated for discordant Replicates
cdata$matches = gsub(" $","",cdata$matches)
cdata$matches = factor(cdata$matches,levels = c("Delta","Omicron(BA.1)","Delta/Omicron(BA.1)","Delta/Omicron(BA.2)","Omicron(BA.1)/Omicron(BA.2)","Omicron(BA.2)","Omicron(BA.3)","Other lineages"))  ## to match colors with fig 3 and 4

####### Relative abundance of contaminant lineage in the 8 samples => FigSxB

pSB <- ggplot(cdata, aes(x=COINF, y=First))
	
pSB <- pSB +
	geom_violin() +
	geom_jitter(aes(shape=COINF, color=matches),width = 0.25, height = 0) +
	xlab("Confirmed Co-infection in duplicates") + ylab("Relative abundance of minor lineage\nin First Replicate (%)") +
	geom_hline(yintercept = 5,linetype =  "dotted", alpha = 0.5) +
	scale_shape(name = "Natural Coinfection") +
	ylim(c(0,50))+
	theme_bw() + scale_color_npg(drop=FALSE)
	

#### SATS : relab of minor lineage (in first replicate) between refuted and confirmed coinf (at 2nd replicate) 

RelAb_Refuted = cdata$First[!cdata$COINF]
			
RelAb_Confirmed = cdata$First[cdata$COINF]

median (RelAb_Refuted)
median (RelAb_Confirmed)
kruskal.test(list(RelAb_Refuted,RelAb_Confirmed))


###########################################################################################################
####  SAVE FIG S8 
############################################################################################################
plot_grid(pSA, pSB,labels = c("A","B"), ncol=1,align = 'v',axis = "rl")
ggsave("FigS8.pdf",width = 6 ,height=7 )

