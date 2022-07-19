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


############################################################################################################################
## load source data
############################################################################################################################

load("SOURCE_seq_metrics_MIXES.Rdata")

load("table_MAF_MIXES.Rdata") ##From Fig1.r

VL = read.delim(file="SOURCE_MIXES_Viral_loads.txt")

############################################################################################################################
### EXCLUDE SAMPLES OF RUN 220530 (PERFORMED TO TEST THE DETECTION LIMIT) : To include all samples, comment the 2 following lines
############################################################################################################################
selected_samples = seq_metrics$sample[seq_metrics$RUN != "220530_13146"]

seq_metrics = seq_metrics[is.element(seq_metrics$sample,selected_samples),]

############################################################################################################################
## merge
############################################################################################################################

table_S2 = merge(table_MAF,seq_metrics,by="sample",all=TRUE)
table_S2 = merge(table_S2,VL,by="sample",all.x=TRUE,all.y=FALSE)


min(table_S2$percCOV)
# [1] 98.20574 

median(table_S2$percCOV)
# [1] 99.49816

median(table_S2$VL)

median(table_S2$mean_depth)
IQR(table_S2$mean_depth)

write.table(table_S2,file="table_S2.txt", sep="\t", quote=F,row.names=FALSE)

