#!/usr/bin/env Rscript
#v0.0.1
library(optparse)
library(tidyr)
library(seqinr)
library(stringr)


option_list = list(
  make_option(c("--output"), type="character", default=NULL, 
              help="output file files.", metavar="character"),
  make_option(c("--coverage"), type="character", default=NULL, 
              help="coverage data from all_cov.cov file.", metavar="character"),
  make_option(c("--count-table"), type="character", default=NULL, 
              help="count data from all_counts.tsv file.", metavar="character"),
  make_option(c("--coinf-table"), type="character", default=NULL, 
              help="coinf data from coinf tsv file.", metavar="character"),
  make_option(c("--cons"), type="character", default=NULL, 
              help="fasta file containing all generated consensus.", metavar="character"),
  make_option(c("--posc"), type="character", default=NULL, 
              help="covidseq positive controls data from posc.tsv file.", metavar="character"),
  make_option(c("--conta"), type="character", default=NULL, 
              help="contamination results from contamination_common_poolt.tsv file.", metavar="character"),
  make_option(c("--mutscan"), type="character", default=NULL, 
              help="mutation screening from mutscan.tsv file.", metavar="character"),
  make_option(c("--runname"), type="character", default=NULL, 
              help="run name", metavar="character")  
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


checkPolySeg<-function(reference){
  split<-str_split(reference,coll("|"))
  if(length(split[[1]])>1) result<-as.character(split[[1]][[2]]) else result<-"FAIL"
  return(result)
}

getReference<-function(reference){
  split<-str_split(reference,coll("|"))
  result<-as.character(split[[1]][[1]])
  return(result)
}

## FASTA
if(!(is.null(opt$cons))){
  alnfa <- read.fasta(file = as.character(opt$cons))
  seqnames<-names(alnfa)
  lengthseq = unlist(lapply(alnfa, function(l) length(l)))
  length_N_tot = unlist(lapply(alnfa, function(l) length(grep("n",l))))
  #Sequence trimming for the internal N check
  alnfa_trim <- alnfa
  for (sequence in names(alnfa_trim)){
    consensus<-as.vector(unlist(alnfa_trim[[sequence]]))
    trimcons<-consensus[150:(length(consensus)-150)]
    alnfa_trim[[sequence]]<-as.list(trimcons)
  }
  Nseqtrim = unlist(lapply(alnfa_trim, function(l) length(grep("n",l))))
  N_check_results<-as.data.frame(cbind(seqnames,length_N_tot,lengthseq,Nseqtrim))
  rownames(N_check_results)<-NULL
  N_check_results$sample<-sapply(N_check_results$seqnames,getReference)
  N_check_results$seqtype<-sapply(N_check_results$seqnames,checkPolySeg)
  N_check_results$percCOV<-100-(as.numeric(as.character(N_check_results$length_N_tot))/as.numeric(as.character(N_check_results$lengthseq))*100)
  nb_seqtype<-length(unique(as.character(N_check_results[which(N_check_results$seqtype!="FAIL"),]$seqtype)))
  #check if polyseg
  N_check_results<-N_check_results[,-c(which(colnames(N_check_results)=="length_N_tot"),
                                       which(colnames(N_check_results)=="Nseqtrim"),
                                       which(colnames(N_check_results)=="lengthseq"))]
  if(nb_seqtype>1){
    N_check_results<-pivot_wider(N_check_results,id_cols = "sample",
                                 names_from = "seqtype",
                                 values_from = c("percCOV"))
    
    if(length(grep("FAIL",colnames(N_check_results)))>0) N_check_results<-N_check_results[,-grep("FAIL",colnames(N_check_results))]
    colnames(N_check_results)[2:length(colnames(N_check_results))]<-paste0("PercCOV_",colnames(N_check_results)[2:length(colnames(N_check_results))])
  }else{
    N_check_results<-N_check_results[,c("sample","percCOV")]
  }
}

## Coverage
if(!(is.null(opt$coverage))){
  cov_table<-read.table(paste0(as.character(opt$repres),as.character(opt$coverage)),header=F)
  colnames(cov_table)<-c("reference","seqtype","mean_depth","sample")
  cov_table$polyseg<-sapply(cov_table$seqtype,checkPolySeg)
  nb_seqtype<-length(unique(as.character(cov_table[which(cov_table$seqtype!="FAIL"),]$seqtype)))
  #check if polyseg
  if(nb_seqtype>1){
    mean_cov<-pivot_wider(cov_table,id_cols = c("sample","reference"),names_from = "polyseg",values_from = "mean_depth")
    if(length(grep("FAIL",colnames(mean_cov)))>0) mean_cov<-mean_cov[,-grep("FAIL",colnames(mean_cov))]
    colnames(mean_cov)[3:length(colnames(mean_cov))]<-paste0("mean_depth_",colnames(mean_cov)[3:length(colnames(mean_cov))])
  }else{
    mean_cov<-cov_table[,c("sample","reference","mean_depth")]
  }
}

## Count
if(!(is.null(opt$count))){
  count_table<-read.table(as.character(opt$count),header=F, sep="\t")
  colnames(count_table)<-c("reference","type","count","sample")
  count_table<-pivot_wider(count_table,id_cols = c("reference","sample"),names_from = "type",values_from = "count")
}


## Coinf
if(!(is.null(opt$coinf))){
  coinf_table<-read.table(as.character(opt$coinf),header=F, sep="\t")
  #count_table<-read.table("Z:/HadrienR/SUMMARY/000000_seqmet_varcall_fluabv/varcall/CY114381_varcount.tsv",header=F, sep="\t")
  colnames(coinf_table)<-c("reference","type","value","sample")
  coinf_table<-pivot_wider(coinf_table,id_cols = c("reference","sample"),names_from = "type",values_from = "value")
}

## posC
if(!(is.null(opt$posc))){
  posc_header<-read.table(as.character(opt$posc), nrows=1, header=F, stringsAsFactors=FALSE)
  posc_table<-read.table(as.character(opt$posc), skip=1, header=F)
  #posc_header<-read.table("Z:/HadrienR/SUMMARY/211018_seqmet_varcall_summaryNCOV/posc.tsv", nrows=1, header=F, stringsAsFactors=FALSE)
  #posc_table<-read.table("Z:/HadrienR/SUMMARY/211018_seqmet_varcall_summaryNCOV/posc.tsv", skip=1, header=F)
  #posc_header<-read.table("Z:/HadrienR/SUMMARY/000000_seqmet_varcall_fluabv/posc.tsv", nrows=1, header=F, stringsAsFactors=FALSE)
  #posc_table<-read.table("Z:/HadrienR/SUMMARY/000000_seqmet_varcall_fluabv/posc.tsv", skip=1, header=F)
  posc <- data.frame(sample = character(0), hasposc = character(0), stringsAsFactors=FALSE)
  for (i in 2:length(posc_header)) {posc[nrow(posc)+1,]<-c(posc_header[i], any(posc_table[,i] >= 0.80))}
  posc$hasposc<-gsub("FALSE", "FAILED", posc$hasposc)
  posc$hasposc<-gsub("TRUE", "OK", posc$hasposc)
}

## Conta
if(!(is.null(opt$conta))){
  conta_header<-read.table(as.character(opt$conta), nrows=1, header=F, stringsAsFactors=FALSE)
  conta_table<-read.table(as.character(opt$conta), skip=1, header=F)
  conta <- data.frame(sample = character(0), hasdp = character(0), stringsAsFactors=FALSE)
  for (i in 2:length(conta_header)) {conta[nrow(conta)+1,]<-c(conta_header[i], as.character(lapply(conta_table[,i][which.max(lapply(conta_table[,i], function(x) str_split(x, "/", n=Inf, simplify = TRUE)[,1]))], function(x) str_split(x, "/", n=Inf, simplify = TRUE)[,1])))}
}


## Merge
alldata<-merge(mean_cov,N_check_results,by.x = "sample",by.y = "sample", all.x=TRUE)
#optional
if(!(is.null(opt$conta))) alldata<-merge(alldata,conta,by = "sample", all.x=TRUE)
if(!(is.null(opt$posc))) alldata<-merge(alldata,posc,by = "sample",all.x=TRUE)
if(!(is.null(opt$mutscan))) alldata<-merge(alldata,t_mutscan,by = "sample" , all.x=TRUE)
if(!(is.null(opt$count))) alldata<-merge(alldata,count_table,by = c("sample","reference"), all.x=TRUE)
if(!(is.null(opt$coinf))) alldata<-merge(alldata,coinf_table,by = c("sample","reference"), all.x=TRUE)
if(!(is.null(opt$runname))) alldata$RUN<-as.character(opt$runname)

write.table(alldata,as.character(opt$output),quote=F,row.names = F,sep="\t")
