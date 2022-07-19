Analysis code for:

Detection and prevalence of SARS-CoV-2 co-infections during the Omicron variant circulation, France, December 2021 - February 2022

Antonin Bal, Bruno Simon, Gregory Destras, Richard Chalvignac, Quentin Semanas, Antoine Oblette, Gregory Queromes, Remi Fanget, 
Hadrien Regue, Florence Morfin, Martine Valette, Bruno Lina, Laurence Josset 

published in medRxiv

Requirements

R software (tested on Windows 10 with v3.6.1 and v4.0.3).

Packages

	data.table
	plyr
	ggplot2
	ggsci
	cowplot
	ggExtra
	reshape2
	dplyr
	tidyr
	gridExtra
	
Usage

Scripts shoud be run individually in independent orders.
Table1.Rmd should be run with RStudio and knit.

Source data include :

- Tables for annotation of the mutations :

* annot_vcf_DO : List of Delta- and Omicron-specific and shared variants based on covariants.org
* profile_lineages : list of defining polymorphism for each lineage obtained with seqmet-db

- Tables related to the experimental Delta:Omicron mixes

* SOURCE_seq_metrics_MIXES : sequencing metrics obtained after running seqmet on the fastq files of Delta:Omicron mixes. fastq of the Delta:Omicron mixes were deposited on the SRA database under accession PRJNA817870 and PRJNA853723
* SOURCE_Variant_Table_MIXES : variant table describing the allele frequency (af) of each mutation found in each of the mixes by using seqmet.
* SOURCE_depth_MIXES : variant table describing the coverage depth (depth) for each position in each of the mixes by using seqmet.
* SOURCE_Variant_Table_MIXES_DRAGEN : variant table describing the allele frequency (af) of each mutation found in each of the mixes by using DRAGEN.
* SOURCE_Variant_Table_MIXES_IVAR : variant table describing the allele frequency (af) of each mutation found in each of the mixes by using https://github.com/connor-lab/ncov2019-artic-nf.
* SOURCE_MIXES_Viral_loads : viral loads in log10 copies/mL

- Tables related to the 21,387 samples sequenced between December 6th 2021 and February 27th 2022

* SOURCE_seq_metrics : sequencing metrics obtained after running seqmet on the fastq files of all samples. fastq of dehosted sequencing data of NPS with co-infections were deposited under accession PRJNA817806
* FORMATED_seq_metrics : formated sequencing metrics with mixes annotation (script for formating is described in Fig3.r)
* SOURCE_Variant_Table : variant table describing the allele frequency (af) of each mutation found in each of the mixes by using seqmet.
* ANONYMISED_METADATA : anonymised sample metadata including week of sampling, age (in year) and sex of patient, hospitalised (0 = outpatient, 1 = hospitalised, 2 = Heathcare worker), ICU = hospitalised in ICU, FIRSTSAMPLE = whether this is the first sample for each patient; together with final result of WGS

- Tables related to duplicate sequencing of samples with a number of secondary lineage-specific mutations between 4 and 5

* SOURCE_Variant_Table_bis : variant table describing the allele frequency (af) of each mutation found in each of the mixes by using seqmet.
* SOURCE_seq_metrics_bis: sequencing metrics obtained after running seqmet on the fastq files of all samples

- Tables related to sequencing of a co-infected samples with a recombination, and its isolation in cell culture (P1 : first passage; P2 : second passage)

* SOURCE_Variant_Table_RECOMB: variant table describing the allele frequency (af) of each mutation found in each of the mixes by using seqmet.


