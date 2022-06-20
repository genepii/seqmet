def print_help() {
  log.info"""
  Usage:
    cd seqmet/piperun
    ./launch_piperun.sh [piperun_folder_name]

  Description:
    genepii/seqmet is a bioinformatics pipeline designed for the analysis of next generation sequencing data (NGS) from viral pathogens.

  Workflow options:
    All parameters can be changed in the params json file of the piperun.

  """.stripIndent()
}
