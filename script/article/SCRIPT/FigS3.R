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
#### CODE USED TO GENERATE FIG S3
############################################################################################################################

options(stringsAsFactors = FALSE) 

library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)

library(cowplot)
library(ggsci)

library(tidyr)

############################################################################################################################
##  SEQMET-db : list of defining polymorphism for each lineage
############################################################################################################################
annot <-as.data.frame(fread("profile_lineages.tsv"))
annot$var = gsub("B_1_617_2","AY_",annot$var)

############################################################################################################################
### FigS3A. number of defining-variants for each lineage
############################################################################################################################
n_defining_mut = annot %>%
	  group_by(var) %>%
	  dplyr::summarise(
		n = n())

n_defining_mut$matches = gsub("BA_","BA",n_defining_mut$var)
n_defining_mut$matches = sapply(n_defining_mut$matches, function(x) strsplit(x,"_")[[1]][1])

n_defining_mut$matches = gsub("AY","Delta",n_defining_mut$matches)
n_defining_mut$matches = gsub("BA1","Omicron(BA.1)",n_defining_mut$matches)
n_defining_mut$matches = gsub("BA2","Omicron(BA.2)",n_defining_mut$matches)
n_defining_mut$matches = gsub("BA3","Omicron(BA.3)",n_defining_mut$matches)
n_defining_mut$matches = gsub("^B","Other lineages",n_defining_mut$matches)
n_defining_mut$matches = factor(n_defining_mut$matches,levels = c("Delta","Omicron(BA.1)","Delta/Omicron(BA.1)","Delta/Omicron(BA.2)","Omicron(BA.1)/Omicron(BA.2)","Omicron(BA.2)","Omicron(BA.3)","Other lineages"))  ## to match colors with fig 3 and 4

p1 <- ggplot(n_defining_mut, aes(x=n,fill=matches)) 

p1 <- p1 +
	geom_vline(xintercept = 6,linetype =  "dotted", color = "red") +
	geom_histogram(,binwidth=1)+
	xlim(c(0,70))+
	xlab("Number of defining-mutations per lineage in seqmet-db") +
	 theme_bw() + scale_fill_npg(drop=FALSE) 

median(n_defining_mut$n)
n_defining_mut[which.min(n_defining_mut$n),]


### number of defining-variants for each lineage
lineage = unique(annot$var)

spe = matrix(NA, nrow= length(lineage),ncol=length(lineage),dimnames=list(lineage,lineage))


for (v1 in lineage){ ## row

	spe[v1,] = sapply(lineage, function(v2) 
		{ ## col
		### subset annotation file to only secondary and main lineages
		annot=annot[is.element(annot$var,c(v1,v2)),]
		annot$nt_mut = paste(annot$pos.nt,annot$nt.mut,sep="")
		
		## filter-out shared mutations
		annot = annot[!is.element(annot$nt_mut,annot$nt_mut[duplicated(annot$nt_mut)]),]
		return(length(annot$nt_mut[annot$var==v1]))
		})
	}

diag(spe) <- 0

write.table(spe,file="n_specific_mutations_per_lineage.txt", sep="\t", quote=F,col.names=NA)



############################################################################################################################
##  FigS3B.HISTO n_specific_mutations_per_lineage
############################################################################################################################

specific_mut = data.frame(spe)

specific_mut$secondary = rownames(specific_mut)

spe_pairwise = pivot_longer(specific_mut,!secondary, names_to = "main", values_to = "count")

## exclude potential coinf with same lineage
spe_pairwise = spe_pairwise[spe_pairwise$secondary != spe_pairwise$main,]

#### FORMAT FOR COLORS
## Format MAIN lineage match (coinf_maj_match)
spe_pairwise$Maj_match = spe_pairwise$main
spe_pairwise$Maj_match[grep("AY",spe_pairwise$Maj_match)] = "Delta"
spe_pairwise$Maj_match[grep("B_1_617_2",spe_pairwise$Maj_match)] = "Delta"
spe_pairwise$Maj_match[grep("BA_1",spe_pairwise$Maj_match)] = "Omicron(BA.1)"
spe_pairwise$Maj_match[grep("BA_2",spe_pairwise$Maj_match)] = "Omicron(BA.2)"
spe_pairwise$Maj_match[grep("BA_3",spe_pairwise$Maj_match)] = "Omicron(BA.3)"

## Format SECONDARY lineage match (coinf_min_match)
spe_pairwise$Min_match = spe_pairwise$secondary
spe_pairwise$Min_match[grep("AY",spe_pairwise$Min_match)] = "Delta"
spe_pairwise$Min_match[grep("B_1_617_2",spe_pairwise$Min_match)] = "Delta"
spe_pairwise$Min_match[grep("BA_1",spe_pairwise$Min_match)] = "Omicron(BA.1)"
spe_pairwise$Min_match[grep("BA_2",spe_pairwise$Min_match)] = "Omicron(BA.2)"
spe_pairwise$Min_match[grep("BA_3",spe_pairwise$Min_match)] = "Omicron(BA.3)"

## Concatenate main and secondary lineages
spe_pairwise$matches = apply(spe_pairwise[,c("Maj_match","Min_match")],1, function(x) paste(x[order(x)],sep="/",collapse="/"))

p2 <- ggplot(spe_pairwise, aes(x=count,fill=matches)) 

p2 <- p2 + 
	geom_vline(xintercept = 6,linetype =  "dotted", color = "red") +
	geom_histogram(,binwidth=1)+
	xlim(c(0,70))+
	xlab("Number of specific-mutations per lineage in all pairwise comparison of seqmet-db") +
	 theme_bw() 

dim(spe_pairwise[spe_pairwise$matches=="Delta/Delta" & spe_pairwise$count<6,])
dim(spe_pairwise[spe_pairwise$matches=="Delta/Delta" ,])

############################################################################################################################
##  SAVE FIG S3
############################################################################################################################
plot_grid(p1, p2, ncol = 1,labels = c("A","B"), align = 'v',axis = "rl")
ggsave("FigS3_seqmetdb.pdf",width = 7 ,height=7)
ggsave("FigS3_seqmetdb.png",width = 7 ,height=7)




