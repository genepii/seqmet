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

VL = read.delim(file="SOURCE_MIXES_Viral_loads.txt")

############################################################################################################################
### SELECT ONLY SAMPLES OF RUN 220530 (PERFORMED TO TEST THE DETECTION LIMIT) 
############################################################################################################################
selected_samples = seq_metrics$sample[seq_metrics$RUN == "220530_13146"]

seq_metrics = seq_metrics[is.element(seq_metrics$sample,selected_samples),]

############################################################################################################################
## merge
############################################################################################################################
table_S4 = merge(seq_metrics,VL,by="sample",all.x=TRUE,all.y=FALSE)


min(table_S4$percCOV)
# [1] 97.88772


median(table_S4$percCOV)
# [1] 99.37071

median(table_S4$VL)

median(table_S4$mean_depth)
IQR(table_S4$mean_depth)

write.table(table_S4,file="table_S4.txt", sep="\t", quote=F,row.names=FALSE)

