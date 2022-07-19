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

library(plyr)
library(ggplot2)
library(ggsci)
library(cowplot)
library(ggExtra)
library(reshape2)

library(dplyr)
library(data.table)

############################################################################################################################
## DATA SOURCE : SEQUENCING METRICS (seq_metrics) and VARIANT TABLE for n=21448 sequencing results (21387 samples + 61 duplicate sequencing) 
############################################################################################################################
load("SOURCE_Variant_Table.Rdata")

load("SOURCE_seq_metrics.Rdata")

############################################################################################################################
##  SEQMET-db : list of defining polymorphism for each lineage
############################################################################################################################
annot <-as.data.frame(fread("profile_lineages.tsv"))

############################################################################################################################
#### FORMAT SEQUENCING METRICS : this part of code generate the file "FORMATED_seq_metrics.txt" also provided as data source
############################################################################################################################
## Coinfection is defined as a positive coinfection minority ratio (coinf_min_ratio) in duplicate sequencing
Likely_COINF = seq_metrics[seq_metrics$coinf_min_ratio>0 ,]

Confirmed_COINF = Likely_COINF$ID[duplicated(Likely_COINF$ID)] ## coinfection is confirmed only if the score is positive in duplicate sequencing

seq_metrics$COINF = is.element(seq_metrics$ID,Confirmed_COINF) ## add a column to seq_metrics to indicate which samples are coinfected

## Format MAIN lineage match (coinf_maj_match)
seq_metrics$Maj_match = seq_metrics$coinf_maj_match
seq_metrics$Maj_match[grep("AY",seq_metrics$Maj_match)] = "Delta"
seq_metrics$Maj_match[grep("B_1_617_2",seq_metrics$Maj_match)] = "Delta"
seq_metrics$Maj_match[grep("BA_1",seq_metrics$Maj_match)] = "Omicron(BA.1)"
seq_metrics$Maj_match[grep("BA_2",seq_metrics$Maj_match)] = "Omicron(BA.2)"
seq_metrics$Maj_match[grep("BA_3",seq_metrics$Maj_match)] = "Omicron(BA.3)"

## Format SECONDARY lineage match (coinf_min_match)
seq_metrics$Min_match = seq_metrics$coinf_min_match
seq_metrics$Min_match[seq_metrics$coinf_min_ratio==0]="" ##for samples with negative coinfection minority ratio, set match to ""
seq_metrics$Min_match[grep("AY",seq_metrics$Min_match)] = "Delta"
seq_metrics$Min_match[grep("B_1_617_2",seq_metrics$Min_match)] = "Delta"
seq_metrics$Min_match[grep("BA_1",seq_metrics$Min_match)] = "Omicron(BA.1)"
seq_metrics$Min_match[grep("BA_2",seq_metrics$Min_match)] = "Omicron(BA.2)"
seq_metrics$Min_match[grep("BA_3",seq_metrics$Min_match)] = "Omicron(BA.3)"

## Concatenate main and secondary lineages
seq_metrics$matches = paste(seq_metrics$Maj_match,seq_metrics$Min_match)

seq_metrics$matches = gsub("Delta Omicron(BA.1)","Delta/Omicron(BA.1)",seq_metrics$matches,fixed=TRUE)
seq_metrics$matches = gsub("Delta Omicron(BA.2)","Delta/Omicron(BA.2)",seq_metrics$matches,fixed=TRUE)
seq_metrics$matches = gsub("Omicron(BA.1) Delta","Delta/Omicron(BA.1)",seq_metrics$matches,fixed=TRUE)
seq_metrics$matches = gsub("Omicron(BA.1) Omicron(BA.2)","Omicron(BA.1)/Omicron(BA.2)",seq_metrics$matches,fixed=TRUE)
seq_metrics$matches = gsub("Omicron(BA.2) Omicron(BA.1)","Omicron(BA.1)/Omicron(BA.2)",seq_metrics$matches,fixed=TRUE)
seq_metrics$matches = gsub("B_1_640_1","Other lineages",seq_metrics$matches,fixed=TRUE)
seq_metrics$matches = gsub(" $","",seq_metrics$matches)

table(as.factor(seq_metrics$matches[seq_metrics$COINF]))/2
      
# Delta/Omicron(BA.1) Delta/Omicron(BA.2)     Omicron/Omicron 
                 # 28                   1                  24 


## Annotate seq_metrics to identify Delta/Omicron coinfections and Omicron(BA.1)/Omicron(BA.2) coinfections
seq_metrics$COINF_Delta_Omicron_BA1 = is.element(seq_metrics$ID,Confirmed_COINF) & seq_metrics$matches=="Delta/Omicron(BA.1)"
seq_metrics$COINF_Delta_Omicron_BA2 = is.element(seq_metrics$ID,Confirmed_COINF) & seq_metrics$matches=="Delta/Omicron(BA.2)"
seq_metrics$COINF_Delta_Omicron = seq_metrics$COINF_Delta_Omicron_BA1 | seq_metrics$COINF_Delta_Omicron_BA2
seq_metrics$COINF_Omicron_Omicron = is.element(seq_metrics$ID,Confirmed_COINF) & seq_metrics$matches=="Omicron(BA.1)/Omicron(BA.2)"

## save
# save(seq_metrics, file="FORMATED_seq_metrics.Rdata")
# write.table(seq_metrics,file="FORMATED_seq_metrics.txt", sep="\t", quote=F,col.names=NA)


###########################################################################################################
#### FIG 3A VIOLIN PLOTS FOR MAIN AND SECONDARY VIRUS MUTATION RATIOS
############################################################################################################
mdata <- melt(seq_metrics[,c("sample","coinf_maj_ratio","coinf_min_ratio","COINF","matches")], id=c("sample","COINF","matches") )

mdata$matches = factor(mdata$matches,levels = c("Delta","Omicron(BA.1)","Delta/Omicron(BA.1)","Delta/Omicron(BA.2)","Omicron(BA.1)/Omicron(BA.2)","Omicron(BA.2)","Omicron(BA.3)","Other lineages"))  ## to match colors with fig 4

pA <- ggplot(mdata, aes(x=variable, y=value))
	
pA <- pA + 
	xlab("Score") + ylab("Mutation ratio") +
	geom_violin() +
	geom_jitter(aes(alpha=COINF,color=matches),width = 0.25, height = 0) +
	scale_alpha_discrete(range = c(0, 0.9)) +  ## alpha set to 0 for samples defined as non co-infected : only coinfected samples are plotted with points
	scale_x_discrete(labels=c("Main lineage","Secondary lineage"),name=NULL)+
	theme_bw() + scale_color_npg()


###########################################################################################################
#### FIG 3B : Distribution of relative abundance for the ‘minor virus’ in natural co-infections
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

## Select duplicated samples
seq_metrics_dup = seq_metrics[is.element(seq_metrics$ID,seq_metrics$ID[duplicated(seq_metrics$ID)]) ,]

# for samples sequenced in duplicate, determine which is the first and second Replicate based on Year of sequencing, and Plate number (Pl)
first_Replicate = ddply(seq_metrics_dup, .(ID),function(x) x[order(x$Year,x$pl),][1,])

seq_metrics_dup$Replicate = NA
seq_metrics_dup$Replicate[is.element(seq_metrics_dup$sample,first_Replicate$sample)] = "First"
seq_metrics_dup$Replicate[!is.element(seq_metrics_dup$sample,first_Replicate$sample)] = "Second"

seq_metrics_dup$lowest_lineage = sapply(seq_metrics_dup$sample, function(x) ifelse(length(RelAb_min$VOC[RelAb_min$sample==x])==1,RelAb_min$VOC[RelAb_min$sample==x],NA))
seq_metrics_dup$lowest_lineage_RelAb = sapply(seq_metrics_dup$sample, function(x) ifelse(length(RelAb_min$med[RelAb_min$sample==x])==1,RelAb_min$med[RelAb_min$sample==x],NA))

## determine whether the least abundant lineage is also the secondary lineage (based on ratio of present mutation)
seq_metrics_dup$is_lowest_secondary <- seq_metrics_dup$lowest_lineage == seq_metrics_dup$coinf_min_match

## Correlation
# prepare data for ggplot2
seq_metrics_dup$lowest_lineage_RelAb[seq_metrics_dup$coinf_min_ratio==0] = 0 ## because when no secondary lineage there is only one RelAb which is the RelAb of the main lineage

cdata <- reshape2::dcast(seq_metrics_dup[,c("ID","lowest_lineage_RelAb","COINF","Replicate")],ID + COINF  ~ Replicate ,value.var = "lowest_lineage_RelAb")
 
cdata$matches = sapply(cdata$ID, function(x) seq_metrics_dup$matches[seq_metrics_dup$ID==x & seq_metrics_dup$Replicate=="Second"])## keep only 2nd Replicate as it is what was validated for discordant Replicates
cdata$matches = factor(cdata$matches,levels = c("Delta","Omicron(BA.1)","Delta/Omicron(BA.1)","Delta/Omicron(BA.2)","Omicron(BA.1)/Omicron(BA.2)","Omicron(BA.2)","Omicron(BA.3)","Other lineages"))  ## to match colors with fig 3 and 4

cdata$maj_or_min = sapply(cdata$ID, function(x) 
	ifelse(all(seq_metrics_dup$is_lowest_secondary[seq_metrics_dup$ID==x ]),"Secundary",
	ifelse(all(!seq_metrics_dup$is_lowest_secondary[seq_metrics_dup$ID==x ]),"Main","Discordance"))
	)

#### Fig3B  :only confirmed coinfections in both duplicates 
p <- ggplot(cdata[cdata$COINF,], aes(x=First, y=Second))
	
p <- p +
	geom_point(aes(shape=maj_or_min,color=matches),alpha=0.9) +
	xlab("Relative abundance of minor lineage\nin First Replicate (%)") + ylab("Relative abundance of minor lineage\nin Second Replicate (%)") +
	geom_abline(intercept = 0, slope = 1,linetype =  "dotted", alpha = 0.5, inherit.aes = FALSE) +
	geom_hline(yintercept = 5,linetype =  "dotted", alpha = 0.5) +
	geom_vline(xintercept = 5,linetype =  "dotted", alpha = 0.5) +
	scale_shape_manual(values=c(17,16,18),name="Minor lineage identified\nin duplicates as",labels=c("discordant","secondary lineage","main lineage"))  +
	xlim(c(0,50)) + ylim(c(0,50))+
	theme_bw() + scale_color_npg(drop=FALSE) 
	
	
pB <- ggMarginal(p, type="boxplot",size=10) 

##### Distribution of relative abundance of minor lineage for CO-INFECTIONS
median(c(cdata$Second [cdata$COINF],cdata$First [cdata$COINF]))
# [1] 20

IQR(c(cdata$Second [cdata$COINF],cdata$First [cdata$COINF]))
# [1] 26.25

table(as.factor(cdata$maj_or_min[cdata$COINF]))
# Discordance   Secundary 
          # 8          45 

#### Comparison of relative abundance between coinfected samples where minor lineage is the secondary lineage in both duplicates and samples with discordance 
RelAb_discordant = c(cdata$Second[cdata$maj_or_min=="Discordance" & cdata$COINF],
			cdata$First[cdata$maj_or_min=="Discordance" & cdata$COINF])
RelAb_secondary = c(cdata$Second[cdata$maj_or_min=="Secundary" & cdata$COINF],
			cdata$First[cdata$maj_or_min=="Secundary" & cdata$COINF])
median (RelAb_discordant)
median (RelAb_secondary)
kruskal.test(list(RelAb_discordant,RelAb_secondary))


##########################################################################################################
#### Fig 3C : what was the rate of chimeric consensus sequences in natural co-infections? 
#############################################################################################################

compute_n_specific_mut_in_consensus <- function(s, vcf_file,annot_file, seq_metrics) {
	require(dplyr)
	
	annot=annot_file
	
	### identify main and secondary lineages for the sample : 
	virus = as.vector(unlist(seq_metrics[seq_metrics$sample == s,c("coinf_maj_match","coinf_min_match")]))
	
	secondary = seq_metrics[seq_metrics$sample == s,"coinf_min_match"]
	
	### subset annotation file to only secondary and main lineages
	annot=annot[is.element(annot$var,virus),]
	annot$nt_mut = paste(annot$pos.nt,annot$nt.mut,sep="")
	
	## filter-out shared mutations
	annot = annot[!is.element(annot$nt_mut,annot$nt_mut[duplicated(annot$nt_mut)]),]

	##########################
	### select vcf ONLY MAJOR MUTATIONS WHICH ARE IN THE CONSENSUS SEQUENCE
	##########################
	vcf=vcf_file[is.element(vcf$sample,s) & vcf$type == "MAJOR",]
	
	vcf$nt_mut = paste(vcf$nt_pos,vcf$nt,sep="")
	
	vcf$VOC <- sapply(vcf$nt_mut, function(x) ifelse(is.element(x,annot$nt_mut),annot$var[x==annot$nt_mut],NA)) 

	vcf = vcf[!is.na(vcf$VOC),]
	
	########## For each lineage, return number of specific mutations within the consensus
	
	n_mut_MAJ = vcf %>%
	  group_by(sample,VOC) %>%
	  dplyr::summarise(n = n())
  
	return(list(n_mut_MAJ))

}

n_mut_MAJ <- sapply(dup, function(s) compute_n_specific_mut_in_consensus(s, vcf,annot, seq_metrics))

n_mut_table = do.call(rbind,n_mut_MAJ)

## identify lineage with lowest number of SPECIFIC polymorphisms in the consensus sequence
n_min = n_mut_table %>%
		group_by(sample) %>%
		slice(if(length(n)==2) which.min(n)
			else (NA))

n_max = n_mut_table %>%
		group_by(sample) %>%
		slice(which.max(n))

seq_metrics_dup$n_min = sapply(seq_metrics_dup$sample, function(x) ifelse(length(n_min$n[RelAb_min$sample==x])==1,n_min$n[n_min$sample==x],NA))
seq_metrics_dup$n_min[is.na(seq_metrics_dup$n_min)] = 0

seq_metrics_dup$n_max = sapply(seq_metrics_dup$sample, function(x) ifelse(length(n_max$n[RelAb_min$sample==x])==1,n_max$n[n_max$sample==x],NA))

## categorize seq_metrics_dup$n_min
seq_metrics_dup$n_min_cat = cut(seq_metrics_dup$n_min, 
                   breaks=c(-Inf, 0, 2, Inf), 
                   labels=c("0","1-2",">2"))
seq_metrics_dup$n_min_cat = factor(seq_metrics_dup$n_min_cat,levels=c(">2","1-2","0"))

## Plot count of chimeric sequences
pC <- ggplot(seq_metrics_dup[seq_metrics_dup$COINF,], aes(x=lowest_lineage_RelAb)) 

pC <- pC + 
	geom_vline(xintercept = 5,linetype =  "dotted", alpha = 0.5) +
	geom_histogram(aes(fill=n_min_cat),binwidth=1)+
	scale_fill_manual(values=c("grey10","grey50","grey80"),name="Number of mutations\nspecific to minor lineage\nin consensus sequence")+
	xlab("Relative abundance of minor lineage") + ylab("Number of sequences") +
	xlim(c(0,50))+
	scale_y_continuous(breaks=seq(0,10,2))+
	theme_bw() 

### Rate of chimeric sequences
seq_metrics_dup$is_chim <- seq_metrics_dup$n_min>0
# FALSE  TRUE 
   # 37    69 

### Chimeric sequences in all sequences with minor lineage above 38%
median(seq_metrics_dup[seq_metrics_dup$COINF & seq_metrics_dup$lowest_lineage_RelAb>38,"n_min"])
IQR(seq_metrics_dup[seq_metrics_dup$COINF & seq_metrics_dup$lowest_lineage_RelAb>38,"n_min"])

########################################################################################################
### SAVE FIG 3
########################################################################################################

plot_grid(pA, pB,pC,labels = c("A","B","C"), ncol=1,align = 'v',axis = "rl")
ggsave("Fig3.pdf",width = 5.6 ,height=8 )
ggsave("Fig3.png",width = 5.6 ,height=8 )

