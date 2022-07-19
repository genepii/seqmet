# seqmet

<p align="center">
  <img alt="HCL logo" src="https://github.com/genepii/seqmet/blob/main/doc/hcl_logo_full.png" width="30%">
&nbsp; &nbsp; &nbsp; &nbsp;
  <img alt="genEPII logo" src="https://github.com/genepii/seqmet/blob/main/doc/genepii_logo_full.png" width="63%">
</p>

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg)](https://sylabs.io/docs/)

## Introduction

**genepii/seqmet** is a bioinformatics pipeline designed for the analysis of next generation sequencing data (NGS) from viral pathogens.

The pipeline is built using [Nextflow](https://www.nextflow.io), is implemented with DSL2, rely on Singularity containers for its processes and on a json parameters file for flexibility and traceability. This framework permits to easily maintain the pipeline while keeping its portability and the reproducibility of its results. The pipeline is expected to be easy to install and be compliant with quality control demands of the French accreditation comity.

Currently it only permits to generate files for SARS-CoV-2, however flu and other viruses will be added progressively. In this first pre-release, it only manages paired-end reads obtained from Illumina sequencing, ONT sequencing and single-end reads should be available soon. De novo assembly and metagenomic approaches will be integrated in a not-so-distant future.

It generates for different purposes:
1. fastq files generated by multiple consecutive optional processes. Fastq files may be filtered, downsampled and trimmed during the pipeline with optional steps. Fastq files with humans reads filtered are compliant for submission on SRA.
2. bam files generated by alignment of the reads on the chosen reference. Multiple consecutive bam files are generated during the pipeline with optional steps, the last generated bam should be kept for visualizing the results in any complex cases.
3. vcf files reporting SNPs and indels as reported by minimap2/freebayes for each sample. Currently, minor SNPs corresponding to the reference base are also reported thanks to a second freebayes calling. Alleles frequences are correctly reported in the AF field of freebayes, usually missing in haploid mode. These files are used by seqmet to generate consensus sequences and detect potential co-infections or contamination of the sample (using a database generated with genepii/seqmet-db).
4. fasta files containing consensus sequences of each sample, only using non-ambiguous iupac codes. Low coverage, primer-related, and error-proned zones are masked in the sequences. These files are compliant for direct submission to GISAID/NCBI.
5. summary file reporting the mean sequencing depth, genome coverage percentage, results obtained for sequencing positive controls (mainly used in covidseqv3) and each main and secondary lineage suspected by the co-infection seeking script. This file provide an overview of the results and quality controls checks obtained for each sample.
6. nextclade file reporting the clade assignment and additional informations for each sample. This file is strictly identical to the one obtained with the official nextclade tool, and will be maintained to match the latest release without delays.
7. pangolin file reporting the pangolin assignment and additional informations for each sample. This file is strictly identical to the one obtained with the official pangolin tool, and will be maintained to match the latest release without delays.

## Pipeline summary

While the workflow and core processes of this pipeline are fixed, many variables can be tweaked in the parameters json file which permits a high flexibility but could also leads to highly variable results. A template of this json is provided, in the piperun folder, and should be used for testing before any issue is submitted. Our team will support any technical or methodological issues but won't be able to deal with any issue encountered when this json default parameters are changed.

1. Filter raw fastq file to exclude human reads for SRA submission ([`sra-human-scrubber`](https://github.com/ncbi/sra-human-scrubber)) ([`bbmap`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools)) - OPTIONAL
2. Downsample fastq pairs randomly, reproductible with a given seed ([`seqtk`](https://github.com/lh3/seqtk)) - OPTIONAL
3. Trim adapters from each read, drop any read pair if at least one read is shorter than a specified length, multiple adapters can be provided ([`cutadapt`](https://github.com/marcelm/cutadapt)) - OPTIONAL
4. Map reads on a set of reference, permit to choose an appropriate reference for each sample (in flu case for example), and verify sequencing success of positive controls ([`Minimap2`](https://github.com/lh3/minimap2))
5. Map reads on each reference passing a given threshold, multiple reference can be used ([`Minimap2`](https://github.com/lh3/minimap2))
6. Filter any duplicate reads as determined by Picard ([`Picard`](https://github.com/broadinstitute/picard)) - OPTIONAL
7. Realign soft-clipped reads, trying to align all reads of a given region in the same way, permit to improve detection of indels ([`abra2`](https://github.com/mozack/abra2)) - OPTIONAL
8. Soft-clip all reads in regions specified in a given bed file, then hard-clip all soft-clipped reads ([`samtools`](https://github.com/samtools/samtools)) ([`bcftools`](https://github.com/samtools/bcftools)) ([`jvarkit`](https://github.com/lindenb/jvarkit)) - OPTIONAL
9. Call, filter and format variants based on a given set of quality criteria, done separately for major and minor variants ([`freebayes`](https://github.com/freebayes/freebayes)) ([`bcftools`](https://github.com/samtools/bcftools)) ([`vt`](https://github.com/atks/vt)) ([`in-house script`])
10. Generate coverage and depth files ([`bedtools`](https://github.com/arq5x/bedtools2)) ([`samtools`](https://github.com/samtools/samtools))
11. Generate consensus sequences only using non-ambiguous iupac codes. Low coverage, primer-related, and error-proned zones are masked ([`bcftools`](https://github.com/samtools/bcftools)) ([`bedtools`](https://github.com/arq5x/bedtools2)) ([`in-house script`])
12. Compare vcf files to count minor/major variants in a sample matching the profile expected for each lineage, the best match will be considered "major lineage" of the sample. Then count minor variants matching the profile expected for each lineage, excluding variants matching the "major lineage", the best match will be considered "minor lineage" of the sample and indicate a potential co-infection. ([`in-house script`])
13. Merge various outputs to obtain comprehensive readable files ([`in-house script`])
14. Generate a comprehensive summary compiling samples quality data ([`in-house script`])
15. Generate a nextclade report ([`Nextclade`](https://github.com/nextstrain/nextclade))
16. Generate a pangolin report ([`pangolin`](https://github.com/cov-lineages/pangolin))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (`>=3.0`)

3. Clone this repository and prepare the nextflow config file and json parameters file:

   > - Build the singularity containers `nextclade`, `pangolin`, `samtools`, `varcall` using the definition files found in singularity/def.
   > - Edit nextflow.config to indicate the absolute path of the cloned repository, the absolute path of the singularity containes built previously, the maximum memory and threads to use, and the folder to mount in the containers (same as the cloned repository if working locally).
   > - Edit the json parameters template found in piperun/000000_seqmet_varcall_ncov to replace "/path/to/seqmet" by the absolute path of the cloned repository.

4. Prepare the database necessary to the worklow and indicate their path in the json parameters file:

   > - Prepare a vcf folder containing all vcf (lineage profiles) needed as reference when searching for co-infection/contamination. These vcfs can be obtained thanks to ([`genepii/seqmet-db`](https://github.com/genepii/seqmet-db)), by building it yourserlf or fetching a release.
   > - Prepare the human_filter.db provided by and for ([`sra-human-scrubber`](https://github.com/ncbi/sra-human-scrubber)). - OPTIONAL

5. Test the pipeline on the provided dataset with a single command using the json template:

   ```console
   cd seqmet/piperun
   ./launch_piperun.sh 000000_seqmet_varcall_ncov
   ```

6. Try on your on dataset, with different parameters:

   ```console
   cd seqmet/piperun
   mkdir date_seqmet_profiledb_ncov
   cp 000000_seqmet_profiledb_ncov/params_01.json date_seqmet_profiledb_ncov/
   # Edit date_seqmet_profiledb_ncov/params_01.json with your raw fastq absolute path in `fastq`, edit absolute result folder path in `result` and tweak any process parameters as suited
   ./launch_piperun.sh date_seqmet_varcall_ncov
   ```

## Documentation

Documentation is still work-in-progress and will be available at [wiki](https://github.com/genepii/seqmet/wiki).

## Credits

The pipeline is primarily designed, written, tested and maintained by Bruno Simon ([@sib0](https://github.com/sib0)) from [GenEPII Sequencing Platform, Institut des Agents Infectieux, Hospices Civils de Lyon, Lyon, France](https://genepii.univ-lyon1.fr/).

In-house scripts were written conjointly with Hadrien Regue ([@HadrienRegue](https://github.com/HadrienRegue)) and Theophile Boyer ([@BoyerTheo](https://github.com/BoyerTheo)) from [GenEPII Sequencing Platform, Institut des Agents Infectieux, Hospices Civils de Lyon, Lyon, France](https://genepii.univ-lyon1.fr/).

## Citations

If you use genepii/seqmet for your analysis, please cite it using the following doi: [medrxiv.org/content/10.1101/2022.03.24.22272871v1](https://doi.org/10.1101/2022.03.24.22272871)

You can cite the `seqmet` publication as follows:

> **Detection and prevalence of SARS-CoV-2 co-infections during the Omicron variant circulation, France, December 2021 - February 2022**
>
> Antonin Bal, Bruno Simon, Gregory Destras, Richard Chalvignac, Quentin Semanas, Antoine Oblette, Gregory Queromes, Remi Fanget, Hadrien Regue, Florence Morfin, Martine Valette, Bruno Lina, Laurence Josset.
>
> medRxiv 2022.03.24.22272871 doi: [10.1101/2022.03.24.22272871](https://doi.org/10.1101/2022.03.24.22272871).
