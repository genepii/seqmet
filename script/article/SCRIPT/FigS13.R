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
#### CODE USED TO GENERATE FIG S13 
############################################################################################################################
options(stringsAsFactors = FALSE) 

library(ggplot2)
library(ggsci)
library(gridExtra)
library(reshape2)

library(data.table)

library(stringr)
library(reshape2)

############################################################################################################################
##  SEQMET-db : list of defining polymorphism for each lineage
############################################################################################################################
lineages <-as.data.frame(fread("profile_lineages.tsv"))


############################################################################################################################
## FOR MIXES
############################################################################################################################
variant_table="SOURCE_Variant_Table_MIXES.Rdata"
seq_metrics="SOURCE_seq_metrics_MIXES.Rdata"

############################################################################################################################
## SELECT SPECIFIC MUTATIONS FOR EACH LINEAGE PRESENT WITHIN A SAMPLE/MIX FUNCTION
############################################################################################################################

load(paste0("~\\",variant_table))
vcf->study
load(paste0("~\\",seq_metrics))


merge(study,seq_metrics,by="sample")->gamme
sapply(1:nrow(gamme),function(x)paste(gamme$nt_pos[x],gamme$nt[x]))->gamme$mut

aftable<-gamme

aftable$V<-sapply(aftable$sample,function(x)str_split(str_split(x,"-")[[1]][2],"x")[[1]][1])
aftable$extract<-sapply(aftable$sample,function(x)str_split(str_split(x,"x")[[1]][2],"[0-9]|[0-9][0-9]|[0-9][0-9][0-9]x")[[1]][1])


colnames(aftable)[which(colnames(aftable)=="coinf_maj_match")]<-"lgmax"
colnames(aftable)[which(colnames(aftable)=="coinf_min_match")]<-"lgmin"
gsub("No match","",aftable$lgmin)->aftable$lgmin
melt(aftable,id=c(colnames(aftable))[-which(str_detect(colnames(aftable),c("lgmax|lgmin")))])->aftable2

aftable2$value<-sapply(aftable2$value,function(x)ifelse(x=="","no",x))

aftable2<-aftable2[aftable2$value!="no",]
aftable2[aftable2$sample%in%unique(aftable2$sample[aftable2$variable=="lgmin"&aftable2$coinf_min_ratio>0]),]->aftable2
gsub("\\.","_",aftable2$value)->aftable2$value
sapply(1:nrow(aftable2),function(x)ifelse(!is.na(unique(aftable2$value[aftable2$sample==aftable2$sample[x]&aftable2$variable=="lgmin"])),ifelse(str_detect(aftable2$mut[x],paste(setdiff(paste(lineages[lineages$var==aftable2$value[x],]$pos.nt,lineages[lineages$var==aftable2$value[x],]$nt.mut),intersect(paste(lineages[lineages$var==aftable2$value[x],]$pos.nt,lineages[lineages$var==aftable2$value[x],]$nt.mut),paste(lineages[lineages$var==unique(aftable2$value[aftable2$sample==aftable2$sample[x]])[-which(str_detect(unique(aftable2$value[aftable2$sample==aftable2$sample[x]]),aftable2$value[x]))],]$pos.nt,lineages[lineages$var==unique(aftable2$value[aftable2$sample==aftable2$sample[x]])[-which(str_detect(unique(aftable2$value[aftable2$sample==aftable2$sample[x]]),aftable2$value[x]))],]$nt.mut))),collapse="|")),ifelse(str_detect(paste(lineages[lineages$var==aftable2$value[x],]$pos.nt,lineages[lineages$var==aftable2$value[x],]$nt.mut,collapse="|"),aftable2$mut[x]),1,0),0),0))->aftable2$spe

aftable2[aftable2$spe==1,]->uniftable

############################################################################################################################
## REMOVE BIASES POSITIONS
############################################################################################################################

cleantable<-uniftable[(uniftable$nt_pos<22648|uniftable$nt_pos>23057)&(uniftable$nt_pos<26255|uniftable$nt_pos>27283)&uniftable$nt_pos!=28270&uniftable$nt_pos!=28271&uniftable$nt_pos!=21987&uniftable$nt_pos!=21986&uniftable$nt_pos!=28881&uniftable$nt_pos,]


############################################################################################################################
## COMPUTE ADJUSTED-R-SQUARED AND P.VALUE TESTING 0 INTERCEPT
############################################################################################################################

for(x in 1:nrow(cleantable)){
  
  data.frame(mut=cleantable$mut[cleantable$sample==cleantable$sample[x]&cleantable$variable==cleantable$variable[x]],cumsum=cumsum(cleantable$af[cleantable$sample==cleantable$sample[x]&cleantable$variable==cleantable$variable[x]]),NTPOS=c(1:length(cumsum(cleantable$af[cleantable$sample==cleantable$sample[x]&cleantable$variable==cleantable$variable[x]]))))->df

  if(is.na(tryCatch((lm(cumsum~NTPOS,data=df)),error=function(e){NA}))[1])
  {cleantable$adjR[x]<-NA}
  else {
    cleantable$adjR[x]<-summary(lm(cumsum~NTPOS,data=df))$adj.r.squared
    cleantable$intercept[x]<-summary(lm(cumsum~NTPOS,data=df))$coefficients[1,4]
  
    }
}

gammetable<-cleantable

############################################################################################################################
## FOR SAMPLES
############################################################################################################################
variant_table="SOURCE_Variant_Table.Rdata"
seq_metrics="FORMATED_seq_metrics.Rdata"
############################################################################################################################
## SELECT SPECIFIC MUTATIONS FOR EACH LINEAGE PRESENT WITHIN A SAMPLE/MIX FUNCTION
############################################################################################################################

load(paste0("~\\",variant_table))
vcf->study
load(paste0("~\\",seq_metrics))


merge(study,seq_metrics,by="sample")->gamme
sapply(1:nrow(gamme),function(x)paste(gamme$nt_pos[x],gamme$nt[x]))->gamme$mut

aftable<-gamme

aftable$V<-sapply(aftable$sample,function(x)str_split(str_split(x,"-")[[1]][2],"x")[[1]][1])
aftable$extract<-sapply(aftable$sample,function(x)str_split(str_split(x,"x")[[1]][2],"[0-9]|[0-9][0-9]|[0-9][0-9][0-9]x")[[1]][1])


colnames(aftable)[which(colnames(aftable)=="coinf_maj_match")]<-"lgmax"
colnames(aftable)[which(colnames(aftable)=="coinf_min_match")]<-"lgmin"
gsub("No match","",aftable$lgmin)->aftable$lgmin
melt(aftable,id=c(colnames(aftable))[-which(str_detect(colnames(aftable),c("lgmax|lgmin")))])->aftable2

aftable2$value<-sapply(aftable2$value,function(x)ifelse(x=="","no",x))

aftable2<-aftable2[aftable2$value!="no",]
aftable2[aftable2$sample%in%unique(aftable2$sample[aftable2$variable=="lgmin"&aftable2$coinf_min_ratio>0]),]->aftable2
gsub("\\.","_",aftable2$value)->aftable2$value
sapply(1:nrow(aftable2),function(x)ifelse(!is.na(unique(aftable2$value[aftable2$sample==aftable2$sample[x]&aftable2$variable=="lgmin"])),ifelse(str_detect(aftable2$mut[x],paste(setdiff(paste(lineages[lineages$var==aftable2$value[x],]$pos.nt,lineages[lineages$var==aftable2$value[x],]$nt.mut),intersect(paste(lineages[lineages$var==aftable2$value[x],]$pos.nt,lineages[lineages$var==aftable2$value[x],]$nt.mut),paste(lineages[lineages$var==unique(aftable2$value[aftable2$sample==aftable2$sample[x]])[-which(str_detect(unique(aftable2$value[aftable2$sample==aftable2$sample[x]]),aftable2$value[x]))],]$pos.nt,lineages[lineages$var==unique(aftable2$value[aftable2$sample==aftable2$sample[x]])[-which(str_detect(unique(aftable2$value[aftable2$sample==aftable2$sample[x]]),aftable2$value[x]))],]$nt.mut))),collapse="|")),ifelse(str_detect(paste(lineages[lineages$var==aftable2$value[x],]$pos.nt,lineages[lineages$var==aftable2$value[x],]$nt.mut,collapse="|"),aftable2$mut[x]),1,0),0),0))->aftable2$spe

aftable2[aftable2$spe==1,]->uniftable

############################################################################################################################
## REMOVE BIASES POSITIONS
############################################################################################################################

cleantable<-uniftable[(uniftable$nt_pos<22648|uniftable$nt_pos>23057)&(uniftable$nt_pos<26255|uniftable$nt_pos>27283)&uniftable$nt_pos!=28270&uniftable$nt_pos!=28271&uniftable$nt_pos!=21987&uniftable$nt_pos!=21986&uniftable$nt_pos!=28881&uniftable$nt_pos,]


############################################################################################################################
## COMPUTE ADJUSTED-R-SQUARED AND P.VALUE TESTING 0 INTERCEPT
############################################################################################################################

for(x in 1:nrow(cleantable)){
  
  data.frame(mut=cleantable$mut[cleantable$sample==cleantable$sample[x]&cleantable$variable==cleantable$variable[x]],cumsum=cumsum(cleantable$af[cleantable$sample==cleantable$sample[x]&cleantable$variable==cleantable$variable[x]]),NTPOS=c(1:length(cumsum(cleantable$af[cleantable$sample==cleantable$sample[x]&cleantable$variable==cleantable$variable[x]]))))->df
  
  if(is.na(tryCatch((lm(cumsum~NTPOS,data=df)),error=function(e){NA}))[1])
  {cleantable$adjR[x]<-NA}
  else {
    cleantable$adjR[x]<-summary(lm(cumsum~NTPOS,data=df))$adj.r.squared
    cleantable$intercept[x]<-summary(lm(cumsum~NTPOS,data=df))$coefficients[1,4]
    
  }
}

patientstable<-cleantable

gammetable[gammetable$V!="MID"&gammetable$V!="MIDV2",]->gammetableV
patientstablef<-data.frame(sample=patientstable$sample,value=patientstable$value,variable=patientstable$variable,af=patientstable$af,mut=patientstable$mut,adjR=patientstable$adjR,intercept=patientstable$intercept)
gammetableVf<-data.frame(sample=gammetableV$sample,value=gammetableV$value,variable=gammetableV$variable,af=gammetableV$af,mut=gammetableV$mut,adjR=gammetableV$adjR,intercept=gammetableV$intercept)

tablefini<-rbind(gammetableVf,patientstablef)
tablef<-unique(tablefini[,-(4:5)])
tablef$type<-sapply(tablef$sample,function(x)ifelse(str_detect(x,"22PlCo"),"gamme",ifelse(str_detect(x,"722000801801|021228537801"),"recomb","samples")))

FigS12<-ggplot(tablef,aes(x=adjR,y=log10(intercept)))+geom_point(aes(colour=type))+scale_colour_discrete(labels=c("mixes","recombinant","patients"))+geom_text(aes(label=ifelse(intercept<0.05&adjR<0.98,sample,"")),vjust=1.5,size=2.5)+geom_hline(yintercept=-1.3,linetype="dashed")+geom_vline(xintercept=0.98,linetype="dashed")+theme_bw()

pdf(file="Fig_S12.pdf",width = 10 ,height=15 )
FigS12
dev.off()





