#!/bin/bash

PIPERUN_DIR="${PWD}"
NF_PATH="${PIPERUN_DIR%%/piperun}"

for PIPERUN in "$@"
do
echo "${PIPERUN}"
cd "${NF_PATH}/piperun/${PIPERUN}"
RUN_ID="${PIPERUN%%_*}_${PIPERUN##*_}"
DATA_TYPE="${PIPERUN##*_}"
    for JSON in *.json
    do
    rm -rf work .nextflow* trace.txt* report.html*
    LOG=""
    "${NF_PATH}/nextflow/nextflow" -C "${NF_PATH}/nextflow/nextflow.config" run "${NF_PATH}/nextflow/main.nf" -params-file "${JSON}" -with-trace -with-report --prefix "${RUN_ID}" || LOG="error"
    chmod -fR 777 "${NF_PATH}/piperun/${PIPERUN}"
    if [[ "${LOG}" != "" ]]; then echo "${RUN_ID} failed"; echo "${RUN_ID} failed - $(date)" >> "${NF_PATH}/piperun/launch_piperun.log"; continue; fi
    cp -r "${JSON}" trace.txt "${NF_PATH}/result/${PIPERUN}/"
        if [[ "${PIPERUN%%_*}" == "000000" ]]
        then
        echo "${DATA_TYPE}"
        #for f in ${NF_PATH}/result/${PIPERUN}/trimmed/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        #for f in ${NF_PATH}/result/${PIPERUN}/downsampled/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        #for f in ${NF_PATH}/result/${PIPERUN}/varcall/minimap2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        #for f in ${NF_PATH}/result/${PIPERUN}/varcall/picard/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        #for f in ${NF_PATH}/result/${PIPERUN}/varcall/abra2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        #cp -r -n "${NF_PATH}/result/${PIPERUN}/" "/srv/net/ngs-stockage.chu-lyon.fr/NGS_Virologie/seqmet/${DATA_TYPE}/" && rm -rf "${NF_PATH}/result/${PIPERUN}/"
        #chmod -fR 777 "/srv/net/ngs-stockage.chu-lyon.fr/NGS_Virologie/seqmet/${DATA_TYPE}/${PIPERUN}/"
        elif [[ ${DATA_TYPE} == "ncov" || ${DATA_TYPE} == "fluabv" || ${DATA_TYPE} == "hsv12" ]]
        then
        echo "${DATA_TYPE}"
        #for f in ${NF_PATH}/result/${PIPERUN}/trimmed/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        #for f in ${NF_PATH}/result/${PIPERUN}/downsampled/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        #for f in ${NF_PATH}/result/${PIPERUN}/varcall/minimap2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        #for f in ${NF_PATH}/result/${PIPERUN}/varcall/picard/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        #for f in ${NF_PATH}/result/${PIPERUN}/varcall/abra2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        #cp -r -n "${NF_PATH}/result/${PIPERUN}/" "/srv/net/ngs-stockage.chu-lyon.fr/NGS_Virologie/seqmet/${DATA_TYPE}/"
        #chmod -fR 777 "/srv/net/ngs-stockage.chu-lyon.fr/NGS_Virologie/seqmet/${DATA_TYPE}/${PIPERUN}/"
        else
        echo "${DATA_TYPE}"
        for f in ${NF_PATH}/result/${PIPERUN}/trimmed/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/downsampled/*.fastq.gz; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/minimap2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/picard/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        for f in ${NF_PATH}/result/${PIPERUN}/varcall/abra2/*/*.bam; do rm -rf "${f}"; touch "${f}"; done
        cp -r -n "${NF_PATH}/result/${PIPERUN}/" "/srv/net/ngs-stockage.chu-lyon.fr/NGS_Virologie/redmine/${DATA_TYPE}/"
        chmod -fR 777 "/srv/net/ngs-stockage.chu-lyon.fr/NGS_Virologie/redmine/${DATA_TYPE}/${PIPERUN}/"
        fi
    done
done
