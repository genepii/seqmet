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
#### CODE USED TO GENERATE FIG S14 S15
############################################################################################################################

options(stringsAsFactors = FALSE) 

library(ggplot2)
library(ggsci)
library(gridExtra)
library(reshape2)
library(data.table)

library(dplyr)

library(cowplot)

############################################################################################################################
## DATA SOURCE : sequencing metrics 
############################################################################################################################

load("FORMATED_seq_metrics.Rdata")

############################################################################################################################
## DATA SOURCE : VARIANT TABLE for recombinant sample 
############################################################################################################################
load("SOURCE_Variant_Table_RECOMB.Rdata")

load("SOURCE_Variant_Table.Rdata")
vcf$tech<-"SEQMET"
vcf$culture<-"P0"

vcf<-rbind(vcf,vcfrecomb)


############################################################################################################################
##  SEQMET-db : list of defining polymorphism for each lineage
############################################################################################################################
annot <-as.data.frame(fread("profile_lineages.tsv"))


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
  
  vcf$AA.mut = paste(vcf$AA_pos,vcf$AA,sep="")
  vcf$AA.mut = apply(vcf[,c("prot","AA.mut")],1, function(x) paste(x,collapse=":"))
  
  vcf$nt_AA.mut = paste(vcf$nt_mut,vcf$AA.mut,sep=";")
  
  ########## Select only specific mutations
  vcfDO = vcf[!is.na(vcf$VOC),]
  
  # Plot
  p <- ggplot(data=vcfDO, aes(x=nt_pos, y=af, group = sample )) 
  
  p <- p +
    geom_bar(aes(fill = VOC), stat="identity") + 
    geom_point(data=vcfDO, aes(x=nt_pos, y=af, group = sample,color = VOC),size=3,alpha=0.8)+
    geom_text(data=vcfDO, aes(x=nt_pos, y=af, group = sample ,label=nt_AA.mut, color = VOC), vjust=0, hjust=0,size=0.5,angle=45)+
    facet_grid( sample ~ID) + 
    ylab("Mutation frequency (%)") + xlab("Nucleotide position") +
    scale_color_manual(values=colors_bar,guide = "none") +
    scale_fill_manual(values=colors_bar,name = "Specific Mutations") +
    geom_hline(yintercept = c(5,50),linetype =  "dotted", alpha = 0.5, inherit.aes = FALSE) +
	xlim(c(0,30000))+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  return(p)
}

############################################################################################################################
## PLOT RECOMB FREQUENCY FONCTION 
############################################################################################################################

plot_recomb <- function(id,bp, chg,vcf_file,annot_file, virus,colors_bar = c("tomato","steelblue1")) {
  library(ggplot2)
  
  vcf_file->vcf
  annot->annot_file
  
  ### identify main and secondary lineages for the sample : 
  ### for each sample choose main and secondary lineages based on the most coinfected duplicate (coinf_min_ratio maximum) : this is to have same mutations between duplicates to compare discordant samples
  # seq_metrics_id = seq_metrics[seq_metrics$ID == id,]
  
  seq_metrics_id<-seq_metrics[seq_metrics$sample == id,]
  
  virus<-as.vector(unlist(seq_metrics_id[which.max(seq_metrics_id$coinf_min_ratio),c("coinf_maj_match","coinf_min_match")]))
  
  ### subset annotation file to only secondary and main lineages
  annot=annot[is.element(annot$var,virus),]
  annot$nt_mut = paste(annot$pos.nt,annot$nt.mut,sep="")
  
  ## filter-out shared mutations
  annot = annot[!is.element(annot$nt_mut,annot$nt_mut[duplicated(annot$nt_mut)]),]
  
  ##########################
  ### annotate vcf
  ##########################
  vcf<-vcf[is.element(vcf$sample,id),]
  
  vcf$nt_mut = paste(vcf$nt_pos,vcf$nt,sep="")
  
  vcf$VOC <- sapply(vcf$nt_mut, function(x) ifelse(is.element(x,annot$nt_mut),annot$var[x==annot$nt_mut],NA)) ## For each variant, determine whether it is Delta- or OMICRON-specific
  
  vcf$AA_mut = paste(vcf$AA_pos,vcf$AA,sep="")
  vcf$AA_mut = apply(vcf[,c("prot","AA_mut")],1, function(x) paste(x,collapse=":"))
  
  vcf$nt_AA_mut = paste(vcf$nt_mut,vcf$AA_mut,sep=";")
  
  ########## Select only specific mutations
  vcfDO = vcf[!is.na(vcf$VOC),]
  
  vcfDO$region<-cut(vcfDO$nt_pos, c(-Inf,c(bp), Inf),labels=c(1:(length(bp)+1)))
  
  mediandf<-as.data.frame(vcfDO%>%
                            group_by(region,VOC)%>%
                            summarise_at(vars(af),list(median=median)))
  
  for(n in 1:(length(bp)+1)){
    if(length(mediandf$region[mediandf$region==n])!=2){
      mediandf<-rbind(mediandf,data.frame(region=n,VOC=setdiff(unique(mediandf$VOC),mediandf$VOC[mediandf$region==n]),median=0))
    }
  }
  
  recomb=ifelse(length(mediandf$VOC)%%(length(bp)+1)==0,mean(unique(sapply(1:nrow(mediandf),function(x)max(mediandf$median[mediandf$VOC==mediandf$VOC[x]])-min(mediandf$median[mediandf$VOC==mediandf$VOC[x]])))),
                max(mediandf$median[mediandf$VOC==names(which.max(table(mediandf$VOC)))])-min(mediandf$median[mediandf$VOC==names(which.max(table(mediandf$VOC)))]))
  if(length(bp)>1) {
    recomb1=recomb-min(mediandf$median[mediandf$median!=0])
    recomb2=recomb-recomb1
  }else {recomb1=NA
  recomb2=NA}
  
  if(length(mediandf$VOC)%%2!=0){
    mediandf<-rbind(mediandf,data.frame(region=setdiff(mediandf$region[mediandf$VOC==names(which.max(table(mediandf$VOC)))],mediandf$region[mediandf$VOC!=names(which.max(table(mediandf$VOC)))]),VOC=unique(mediandf$VOC[mediandf$VOC!=names(which.max(table(mediandf$VOC)))]),median=100))}
  mediandf$virus<-sapply(1:nrow(mediandf),function(x)ifelse(mediandf$VOC[x]==mediandf$VOC[mediandf$median==min(mediandf$median[mediandf$median!=0])],"min","max"))
  
  if(is.na(recomb1)){
    
    dfplot<-rbind(data.frame(region=c(1:(length(bp)+1)),VOC="recomb",median=recomb,virus="recomb"),mediandf)
    dfplot$perc<-sapply(1:nrow(dfplot),function(x)ifelse(dfplot$virus[x]=="recomb",round(recomb),ifelse(dfplot$virus[x]=="min",min(dfplot$median[dfplot$virus=="min"]),100-round(recomb)-min(dfplot$median[dfplot$virus=="min"]))))
    dfplot$bp<-sapply(dfplot$region,function(x)ifelse(as.numeric(x)>length(bp),30000-as.numeric(bp[as.numeric(x)-1]),bp[as.numeric(x)]))
    dfplot$fill<-sapply(1:nrow(dfplot),function(x)ifelse(dfplot$VOC[x]=="recomb",ifelse(dfplot$region[x]==1,dfplot$VOC[dfplot$median==max(dfplot$median[dfplot$region==1&dfplot$VOC!="recomb"])&dfplot$VOC!="recomb"],dfplot$VOC[dfplot$median==min(dfplot$median[dfplot$region==1&dfplot$VOC!="recomb"])&dfplot$VOC!="recomb"]),dfplot$VOC[x]))
    
    p<-ggplot(dfplot,aes(x=bp,y=VOC,fill=fill))+geom_bar(stat="identity",position=position_stack(reverse=chg))+ theme(legend.background = element_rect(fill = "transparent"),legend.box.background = element_rect(fill = "transparent"),panel.background = element_rect(fill = "transparent"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "transparent",color = NA))+geom_text(aes(label=ifelse(bp>10000,paste(perc,"%"),"")))   
  }else{
    dfplot<-rbind(data.frame(region=rep(c(1:(length(bp)+1)),2),VOC=c(rep("recomb1",3),rep("recomb2",3)),median=c(rep(recomb1,3),rep(recomb2,3)),virus=c(rep("recomb1",3),rep("recomb2",3))),mediandf)
    dfplot$perc<-sapply(1:nrow(dfplot),function(x)ifelse(dfplot$virus[x]=="recomb1",round(recomb1),ifelse(dfplot$virus[x]=="recomb2",round(recomb2),ifelse(dfplot$virus[x]=="min",min(dfplot$median[dfplot$virus=="min"]),100-round(recomb)-min(dfplot$median[dfplot$virus=="min"])))))
    dfplot$bp<-sapply(dfplot$region,function(x)ifelse(x==2,as.numeric(bp[2]-bp[1]),ifelse(x==3,30000-as.numeric(bp[2]),bp[as.numeric(x)])))
    dfplot$fill<-sapply(1:nrow(dfplot),function(x)ifelse(dfplot$VOC[x]=="recomb2",ifelse(dfplot$region[x]==1,dfplot$VOC[dfplot$median==max(dfplot$median[dfplot$region==1&dfplot$VOC!="recomb1"&dfplot$VOC!="recomb2"])&dfplot$VOC!="recomb1"&dfplot$VOC!="recomb2"],ifelse(dfplot$region[x]==2,dfplot$VOC[dfplot$median==min(dfplot$median[dfplot$region==2&dfplot$VOC!="recomb1"&dfplot$VOC!="recomb2"])&dfplot$VOC!="recomb1"&dfplot$VOC!="recomb2"],dfplot$VOC[dfplot$median==max(dfplot$median[dfplot$region==3&dfplot$VOC!="recomb1"&dfplot$VOC!="recomb2"])&dfplot$VOC!="recomb1"&dfplot$VOC!="recomb2"])),ifelse(dfplot$VOC[x]=="recomb1",ifelse(dfplot$region[x]==1,dfplot$VOC[dfplot$median==max(dfplot$median[dfplot$region==1])],ifelse(dfplot$region[x]==2,dfplot$VOC[dfplot$median==max(dfplot$median[dfplot$region==2&dfplot$VOC!="recomb1"&dfplot$VOC!="recomb2"])&dfplot$VOC!="recomb1"&dfplot$VOC!="recomb2"],dfplot$VOC[dfplot$median==max(dfplot$median[dfplot$region==3&dfplot$VOC!="recomb1"&dfplot$VOC!="recomb2"])&dfplot$VOC!="recomb1"&dfplot$VOC!="recomb2"])),dfplot$VOC[x])))
    p <- ggplot(dfplot,aes(x=bp,y=VOC,fill=fill))+
		geom_bar(stat="identity",position=position_stack(reverse=chg))+ 
		xlim(c(0,30000))+
		theme(legend.background = element_rect(fill = "transparent"),legend.box.background = element_rect(fill = "transparent"),panel.background = element_rect(fill = "transparent"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.background = element_rect(fill = "transparent",color = NA))+geom_text(aes(label=ifelse(bp>10000,paste(perc,"%"),"")))
  }
  return(p)
}

#plot_recomb(id = sample,bp=c(breakpoint1,...), vcf_file=vcf_file,chg=TRUE,annot_file = annot, seq_metrics,colors_bar = c("tomato","steelblue1"))

############################################################################################################################
### FigS14
############################################################################################################################


##for each sample with recombinant suspicion, plot allele frequencies for the main and secondary lineages and the estimation of viral sub-populations frequencies

vcf_file<-vcf[vcf$ID=="021228537801"&vcf$culture=="P0"&vcf$sample=="Pl882-021228537801_S1293",]
P0<-plot_AF(id = "021228537801", vcf_file=vcf_file,annot_file = annot, seq_metrics,colors_bar = c("tomato","steelblue1"))

vcf_file<-vcf[vcf$ID=="021228537801"&vcf$culture=="P1"&vcf$sample=="ARTIC V4.1 SEQMET",]
vcf_file$sample<-"P1"
P1<-plot_AF(id = "021228537801", vcf_file=vcf_file,annot_file = annot, seq_metrics,colors_bar = c("tomato","steelblue1"))   

vcf_file<-vcf[vcf$ID=="021228537801"&vcf$culture=="P2"&vcf$sample=="ARTIC V4.1 SEQMET",]
vcf_file$sample<-"P2"
P2<-plot_AF(id = "021228537801", vcf_file=vcf_file,annot_file = annot, seq_metrics,colors_bar = c("tomato","steelblue1"))  


vcf_file<-vcf[vcf$ID=="021228537801"&vcf$culture=="P0",]
recP0<-plot_recomb(id = "Pl882-021228537801_S1293",bp=c(15000,20000), vcf_file=vcf_file,chg=TRUE,annot_file = annot, seq_metrics,colors_bar = c("tomato","steelblue1"))

vcf_file<-vcf[vcf$ID=="021228537801"&vcf$culture=="P1",]
vcf_file$sample<-"Pl882-021228537801_S1293"
recP1<-plot_recomb(id = "Pl882-021228537801_S1293",bp=c(15000,20000), vcf_file=vcf_file,chg=TRUE,annot_file = annot, seq_metrics,colors_bar = c("tomato","steelblue1"))

vcf_file<-vcf[vcf$ID=="021228537801"&vcf$culture=="P2"&vcf$sample=="ARTIC V4.1 SEQMET",]
vcf_file$sample<-"Pl882-021228537801_S1293"
recP2<-plot_recomb(id = "Pl882-021228537801_S1293",bp=c(15000,20000), vcf_file=vcf_file,chg=TRUE,annot_file = annot, seq_metrics,colors_bar = c("tomato","steelblue1"))


plotP0<-plot_grid(P0,recP0,align="v", axis="tblr",ncol=1,rel_heights = c(1,0.5))
plotP1<-plot_grid(P1,recP1,align="v", axis="tblr",ncol=1,rel_heights = c(1,0.5))
plotP2<-plot_grid(P2,recP2,align="v", axis="tblr",ncol=1,rel_heights = c(1,0.5))

vcf_file<-vcf[vcf$sample=="22Pl155a-722000801801_S551",]
p801<-plot_AF(id = "722000801801", vcf_file=vcf_file,annot_file = annot, seq_metrics,colors_bar = c("steelblue1","lightsteelblue3"))
recp801<-plot_recomb(id = "22Pl155a-722000801801_S551",bp=c(8000), vcf_file=vcf,chg=TRUE,annot_file = annot, seq_metrics,colors_bar = c("steelblue1","lightsteelblue3"))
plot801<-plot_grid(p801,recp801,align="v", axis="tblr",ncol=1,rel_heights = c(1,0.5))

vcf_file<-vcf[vcf$sample=="Pl910-021229656701_S869",]
p965<-plot_AF(id = "021229656701", vcf_file=vcf_file,annot_file = annot, seq_metrics,colors_bar = c("tomato","steelblue1"))
recp965<-plot_recomb(id = "Pl910-021229656701_S869",bp=c(23800), vcf_file=vcf,chg=TRUE,annot_file = annot, seq_metrics,colors_bar = c("tomato","steelblue1"))
plot965<-plot_grid(p965,recp965,align="v", axis="tblr", ncol=1,rel_heights = c(1,0.5))

FigS14<-plot_grid(plot_grid(plotP0,plotP1,plotP2,labels=c("A"),nrow=1),plot_grid(plot801,plot965,nrow=1,labels=c("B","C")),nrow=2)
FigS14
ggsave(file="Fig_S14.pdf",width = 15 ,height=7 )

plot_grid(plotP0,plotP1,plotP2,labels=c("A"),nrow=1)
ggsave(file="Fig_S14A.pdf",width = 17 ,height=7/2 )


############################################################################################################################
### FiGS15
############################################################################################################################

vcf_file<-vcf[vcf$ID=="021228537801"&vcf$culture=="P2",]
vcf_file$sample<-factor(vcf_file$sample,levels=c("ARTIC V4.1 SEQMET","ARTIC V4.1 COGUK","mNGS SEQMET"))
FigS15<-plot_AF(id = "021228537801", vcf_file=vcf_file,annot_file = annot, seq_metrics,colors_bar = c("tomato","steelblue1")) 

FigS15
ggsave(file="FigS15.pdf",width = 7 ,height=7 )