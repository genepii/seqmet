process filter_sra {

    // Filter fastq pairs removing all reads identified as potentially of human origin

    label 'srahumanscrubber'
    storeDir params.result
    echo false
    tag { sampleId }

    when:
    params.filter_sra.todo == 1

    input:
    tuple(val(sampleId), path(R1), path(R2))
    
    output:
    tuple val(sampleId), path("filtered/${sampleId}_R1.fastq.gz"), path("filtered/${sampleId}_R2.fastq.gz")
    
    script:
    sampleId = (R1 =~ /(.+)_R([^\._]+)(.+)/)[0][1]
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    mkdir -p work filtered
    export TMPDIR=${params.tmp}

    if [[ ! -s ${R1} ]]; then
        touch "filtered/${sampleId}_R1.fastq.gz" "filtered/${sampleId}_R2.fastq.gz"
        
        exit 0
    fi

    zcat $R1 | OMP_NUM_THREADS=$task.cpus scrub.sh -o "work/${sampleId}_R1.fastq" -d ${params.filter_sra["db"]}
    zcat $R2 | OMP_NUM_THREADS=$task.cpus scrub.sh -o "work/${sampleId}_R2.fastq" -d ${params.filter_sra["db"]}

    repair.sh -Xmx${memory}G threads=$task.cpus in="work/${sampleId}_R1.fastq" in2="work/${sampleId}_R2.fastq" out="filtered/${sampleId}_R1.fastq" out2="filtered/${sampleId}_R2.fastq"

    pigz -p $task.cpus -0 "filtered/${sampleId}_R1.fastq"
    pigz -p $task.cpus -0 "filtered/${sampleId}_R2.fastq"
    #rm -rf work
    """
}

process downsample_seqtk {

    // Downsample fastq pairs randomly, reproductible with a given seed

    label 'varcall'
    storeDir params.result
    echo false
    tag { sampleId }

    when:
    params.downsample_seqtk.todo == 1

    input:
    tuple(val(sampleId), path(R1), path(R2))
    
    output:
    tuple val(sampleId), path("downsampled/${sampleId}_R1.fastq.gz"), path("downsampled/${sampleId}_R2.fastq.gz")
    
    script:
    sampleId = (R1 =~ /(.+)_R([^\._]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work downsampled

    if [[ ! -s ${R1} ]]; then
        touch "downsampled/${sampleId}_R1.fastq.gz" "downsampled/${sampleId}_R2.fastq.gz"
        
        exit 0
    fi

    seqtk sample -s 100 ${R1} ${params.downsample_seqtk["reads"]} | pigz -p $task.cpus -0 > "downsampled/${sampleId}_R1.fastq.gz"
    seqtk sample -s 100 ${R2} ${params.downsample_seqtk["reads"]} | pigz -p $task.cpus -0 > "downsampled/${sampleId}_R2.fastq.gz"
    rm -rf work
    """
}

process trim_cutadapt {

    // Trim adapters from each read, drop any read pair if at least one read is shorter than a specified length, multiple adapters can be provided

    label 'varcall'
    storeDir params.result
    echo false
    tag { sampleId }

    when:
    params.trim_cutadapt.todo == 1

    input:
    tuple(val(sampleId), path(R1), path(R2))
    
    output:
    tuple val(sampleId), path("trimmed/${sampleId}_R1.fastq.gz"), path("trimmed/${sampleId}_R2.fastq.gz")
    
    script:
    sampleId = (R1 =~ /(.+)_R([^\._]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work trimmed

    if [[ ! -s ${R1} ]]; then
        touch "trimmed/${sampleId}_R1.fastq.gz" "trimmed/${sampleId}_R2.fastq.gz"
        
        exit 0
    fi

    adapter_f=\$(for item in \$(echo ${params.trim_cutadapt["adapter_f"]} | tr -d '[\\[\\],\\"]'); do echo -n "-a \${item} "; done)
    adapter_r=\$(for item in \$(echo ${params.trim_cutadapt["adapter_r"]} | tr -d '[\\[\\],\\"]'); do echo -n "-A \${item} "; done)
    front_f=\$(for item in \$(echo ${params.trim_cutadapt["front_f"]} | tr -d '[\\[\\],\\"]'); do echo -n "-g \${item} "; done)
    front_r=\$(for item in \$(echo ${params.trim_cutadapt["front_r"]} | tr -d '[\\[\\],\\"]'); do echo -n "-G \${item} "; done)
    anywhere_f=\$(for item in \$(echo ${params.trim_cutadapt["anywhere_f"]} | tr -d '[\\[\\],\\"]'); do echo -n "-b \${item} "; done)
    anywhere_r=\$(for item in \$(echo ${params.trim_cutadapt["anywhere_r"]} | tr -d '[\\[\\],\\"]'); do echo -n "-B \${item} "; done)
    cutadapt --cores $task.cpus --pair-filter=any --minimum-length ${params.trim_cutadapt["minimum_length"]} --error-rate ${params.trim_cutadapt["error_rate"]} --nextseq-trim ${params.trim_cutadapt["quality"]} --trim-n --overlap ${params.trim_cutadapt["overlap"]} \${adapter_f}\${adapter_r}\${front_f}\${front_r}\${anywhere_f}\${anywhere_r} -o "trimmed/${sampleId}_R1.fastq" -p "trimmed/${sampleId}_R2.fastq" $R1 $R2
    pigz -p $task.cpus -0 "trimmed/${sampleId}_R1.fastq" "trimmed/${sampleId}_R2.fastq"
    rm -rf work
    """
}

process premap_minimap2 {

    // Map reads on a set of reference, permit to choose an appropriate reference for each sample, and verify sequencing success of positive controls
    // TODO Reformat depth tsv to keep positions as column while permitting merging

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.premap_minimap2.todo == 1

    input:
    tuple(val(sampleId), path(R1), path(R2))
    file 'ref/*'
    file 'posc/*'
    
    output:
    tuple val(sampleId), path("varcall/premap/depth/${sampleId}.tsv"), emit: depthtable
    tuple val(sampleId), path("varcall/premap/coverage/${sampleId}.tsv"), emit: coveragetable
    path("varcall/premap/posc/${sampleId}.tsv"), emit: poscsummary
    
    script:
    sampleId = (R1 =~ /(.+)_R([^\._]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p varcall/premap/depth varcall/premap/coverage varcall/premap/posc work/

    if [[ ! -s ${R1} ]]; then
        touch "varcall/premap/coverage/${sampleId}.cov" "varcall/premap/depth/${sampleId}.tsv" "varcall/premap/coverage/${sampleId}.tsv" "varcall/premap/posc/${sampleId}.tsv"
        
        exit 0
    fi

    for ref in ref/*; do bn=\$(basename "\$ref" | cut -d. -f1); grep -v "^>" "\${ref}" | tr -d '\\015' | tr -d '\\n' | \
    awk -v ref=\${bn%%.*} 'BEGIN { print ">"ref } { print }' >> ref/refpm.fna ; done;
    for ref in posc/*; do bn=\$(basename "\$ref" | cut -d. -f1); grep -v "^>" "\${ref}" | tr -d '\\015' | tr -d '\\n' | \
    awk -v ref=\${bn%%.*} 'BEGIN { print ">"ref } { print }' >> ref/refpm.fna ; done;
    for ref in posc/*; do bn=\$(basename "\$ref" | cut -d. -f1); echo \$bn >> posc.txt; done;
    minimap2 -t $task.cpus -a --sr --frag=${params.premap_minimap2["frag"]} -F${params.premap_minimap2["F"]} -k${params.premap_minimap2["k"]} -w${params.premap_minimap2["w"]} -A${params.premap_minimap2["A"]} -B${params.premap_minimap2["B"]} -O${params.premap_minimap2["O"]} -E${params.premap_minimap2["E"]} -r${params.premap_minimap2["r"]} -p${params.premap_minimap2["p"]} -N${params.premap_minimap2["N"]} -f${params.premap_minimap2["f"]} -n${params.premap_minimap2["n"]} -m${params.premap_minimap2["m"]} -s${params.premap_minimap2["s"]} -g${params.premap_minimap2["g"]} --heap-sort=${params.premap_minimap2["heap-sort"]} --secondary=${params.premap_minimap2["secondary"]} ref/refpm.fna $R1 $R2 | \
    samtools view --threads $task.cpus -b -F 4 - | \
    samtools sort --threads $task.cpus -o "ref/${sampleId}_sorted.bam"
    samtools index "ref/${sampleId}_sorted.bam"
    bedtools genomecov -ibam "ref/${sampleId}_sorted.bam" -d > "${sampleId}.tsv"
    posc=\$(cat posc.txt | perl -pe 'chomp if eof' | tr '\\n' '\\|'); egrep "\$posc" "${sampleId}.tsv" > "work/${sampleId}_posc_depth.tsv"; egrep -v "\$posc" "${sampleId}.tsv" > "work/${sampleId}_ref_depth.tsv"
    for ref in \$(cut -f 1 "work/${sampleId}_ref_depth.tsv" | sort -u); do cov=\$(bc <<< "scale=4; (\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref && \$3>=${params.premap_minimap2["min_depth"]} {bp+=1} {print bp}' "work/${sampleId}_ref_depth.tsv" | tail -1)/\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref {bp+=1} {print bp}' "work/${sampleId}_ref_depth.tsv" | tail -1))"); echo -e "\${ref}\\t\${cov}" >> "varcall/premap/coverage/${sampleId}.tsv"; done;
    for ref in \$(cut -f 1 "work/${sampleId}_posc_depth.tsv" | sort -u); do cov=\$(bc <<< "scale=4; (\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref && \$3>=${params.premap_minimap2["min_depth"]} {bp+=1} {print bp}' "work/${sampleId}_posc_depth.tsv" | tail -1)/\$(awk -v ref=\$ref 'BEGIN {bp=0} \$1==ref {bp+=1} {print bp}' "work/${sampleId}_posc_depth.tsv" | tail -1))"); echo -e "\${ref}\\t\${cov}" >> "varcall/premap/posc/${sampleId}.tsv"; done;
    awk -F'[\\t]' '{print \$1"_"\$2"\\t"\$3}' "${sampleId}.tsv" > "varcall/premap/depth/${sampleId}.tsv"
    #rm -rf ref
    """
}

process map_minimap2 {

    //  Map reads on each reference passing a given threshold, multiple reference can be used if max_match > 1

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.map_minimap2.todo == 1

    input:
    tuple(val(sampleId), path(R1), path(R2), path("coverage/${sampleId}.tsv"), path(ref))
    
    output:
    tuple val(sampleId), path("varcall/minimap2/${ref.simpleName}/${sampleId}.bam"), path("varcall/minimap2/${ref.simpleName}/${sampleId}.bam.bai"), path("varcall/ref/${ref.simpleName}.fna")
    
    script:
    sampleId = (R1 =~ /(.+)_R([^\._]+)(.+)/)[0][1]
    refName = (ref =~ /(.+)\.(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work varcall/ref/ varcall/minimap2/${ref.simpleName}/
    
    bn=\$(basename "${ref}" | cut -d. -f1)
    awk -F '[\\\t_]' '{print \$1"\\t"\$NF}' "coverage/${sampleId}.tsv" > "coverage/premap.tsv"

    if [[ ! -s ${R1} || ! \$bn =~ ^(\$(cat "coverage/premap.tsv" | sort -t\$'\\t' -k2 -nr | head -n ${params.map_minimap2["max_match"]} | awk -F'[\\t]' '{ORS="|"; if (\$2>${params.map_minimap2["min_cov"]}) print \$1 }'))\$ ]]; then
        touch "varcall/minimap2/${ref.simpleName}/${sampleId}.bam" "varcall/minimap2/${ref.simpleName}/${sampleId}.bam.bai"
        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    minimap2 -t $task.cpus -a --sr --frag=${params.map_minimap2["frag"]} -F${params.map_minimap2["F"]} -k${params.map_minimap2["k"]} -w${params.map_minimap2["w"]} -A${params.map_minimap2["A"]} -B${params.map_minimap2["B"]} -O${params.map_minimap2["O"]} -E${params.map_minimap2["E"]} -r${params.map_minimap2["r"]} -p${params.map_minimap2["p"]} -N${params.map_minimap2["N"]} -f${params.map_minimap2["f"]} -n${params.map_minimap2["n"]} -m${params.map_minimap2["m"]} -s${params.map_minimap2["s"]} -g${params.map_minimap2["g"]} --heap-sort=${params.map_minimap2["heap-sort"]} --secondary=${params.map_minimap2["secondary"]} --end-bonus=${params.map_minimap2["end-bonus"]} ${ref} $R1 $R2 | \
    samtools view --threads $task.cpus -b -F 4 - | \
    samtools sort --threads $task.cpus -o "varcall/minimap2/${ref.simpleName}/${sampleId}.bam"
    samtools index "varcall/minimap2/${ref.simpleName}/${sampleId}.bam"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    rm -rf work
    """
}

process postmap_picard {

    // Filter any duplicate reads as determined by Picard

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.postmap_picard.todo == 1

    input:
    tuple(val(sampleId), path(bam), path(bai), path(ref))
    
    output:
    tuple val(sampleId), path("varcall/picard/${ref.simpleName}/${sampleId}.bam"), path("varcall/picard/${ref.simpleName}/${sampleId}.bam.bai"), path("varcall/ref/${ref.simpleName}.fna")
    
    script:
    sampleId = (bam =~ /(.+)\.(.+)/)[0][1]
    refName = (ref =~ /(.+)\.(.+)/)[0][1]
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work/${ref.simpleName}/ varcall/ref/ varcall/picard/${ref.simpleName}/

    if [[ ! -s ${bam} ]]; then
        touch "varcall/picard/${ref.simpleName}/${sampleId}.bam" "varcall/picard/${ref.simpleName}/${sampleId}.bam.bai"
        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    java -XX:ParallelGCThreads=$task.cpus -Xmx${memory}G -jar /opt/conda/envs/varcall/share/picard-2.25.0-1/picard.jar MarkDuplicates I=${bam} O=varcall/picard/${ref.simpleName}/${sampleId}.bam M=marked_dup_metrics.txt REMOVE_DUPLICATES=true TMP_DIR=${params.tmp}
    samtools index "varcall/picard/${ref.simpleName}/${sampleId}.bam"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    rm -rf work
    """
}

process postmap_abra2 {

    // Realign soft-clipped reads, trying to align all reads of a given region in the same way, permit to improve detection of indels 

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.postmap_abra2.todo == 1

    input:
    tuple(val(sampleId), path(bam), path(bai), path(ref))
    
    output:
    tuple val(sampleId), path("varcall/abra2/${ref.simpleName}/${sampleId}.bam"), path("varcall/abra2/${ref.simpleName}/${sampleId}.bam.bai"), path("varcall/ref/${ref.simpleName}.fna")
    
    script:
    sampleId = (bam =~ /(.+)\.(.+)/)[0][1]
    refName = (ref =~ /(.+)\.(.+)/)[0][1]
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work varcall/ref/ varcall/abra2/${ref.simpleName}/

    if [[ ! -s ${bam} ]]; then
        touch "varcall/abra2/${ref.simpleName}/${sampleId}.bam" "varcall/abra2/${ref.simpleName}/${sampleId}.bam.bai"
        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    java -XX:ParallelGCThreads=$task.cpus -Xmx${memory}G -jar /opt/conda/envs/varcall/share/abra2-2.24-1/abra2.jar ${params.postmap_abra2["args"]} --amq ${params.postmap_abra2["amq"]} --ca ${params.postmap_abra2["ca"]} --mac ${params.postmap_abra2["mac"]} --mad ${params.postmap_abra2["mad"]} --mapq ${params.postmap_abra2["mapq"]} --maxn ${params.postmap_abra2["maxn"]} --mbq ${params.postmap_abra2["mbq"]} --mcl ${params.postmap_abra2["mcl"]} --mcr ${params.postmap_abra2["mcr"]} --mer ${params.postmap_abra2["mer"]} --mmr ${params.postmap_abra2["mmr"]} --mnf ${params.postmap_abra2["mnf"]} --mrn ${params.postmap_abra2["mrn"]} --mrr ${params.postmap_abra2["mrr"]} --msr ${params.postmap_abra2["msr"]} --rcf ${params.postmap_abra2["rcf"]} --in ${bam} --out "varcall/abra2/${ref.simpleName}/${sampleId}.bam" --ref ${ref} --threads $task.cpus --tmpdir ${params.tmp}
    samtools index "varcall/abra2/${ref.simpleName}/${sampleId}.bam"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    rm -rf work
    """
}

process postmap_ampliconclip {

    // Soft-clip all reads in regions specified in a given bed file, then hard-clip all soft-clipped reads

    label 'samtools'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.postmap_ampliconclip.todo == 1

    input:
    tuple(val(sampleId), path(bam), path(bai), path(ref))
    file 'ampliconclip/*'
    
    output:
    tuple val(sampleId), path("varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam"), path("varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam.bai"), path("varcall/ref/${ref.simpleName}.fna")
    
    script:
    sampleId = (bam =~ /(.+)\.(.+)/)[0][1]
    refName = (ref =~ /(.+)\.(.+)/)[0][1]
    memory = (task.memory =~ /([^\ ]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    mkdir -p work/${ref.simpleName}/ varcall/ref/ varcall/ampliconclip/${ref.simpleName}/

    if [[ ! -s ${bam} ]]; then
        touch "varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam" "varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam.bai"
        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi
    
    cat ampliconclip/${refName}* > amplicon.bed
    samtools ampliconclip -o "work/${ref.simpleName}/${sampleId}.bam" --both-ends --filter-len ${params.postmap_ampliconclip["filter_len"]} -b amplicon.bed ${bam}
    samtools sort --threads $task.cpus -o "work/${ref.simpleName}/${sampleId}_sorted.bam" "work/${ref.simpleName}/${sampleId}.bam"
    java -XX:ParallelGCThreads=$task.cpus -Xmx${memory}G -jar /opt/conda/envs/samtools/bin/biostar84452.jar --samoutputformat bam "work/${ref.simpleName}/${sampleId}_sorted.bam" > "varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam"
    samtools index "varcall/ampliconclip/${ref.simpleName}/${sampleId}.bam"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    rm -rf work
    """
}

process callvar_freebayes {

    // Call, filter and format variants based on a given set of quality criteria, done separately for major and minor variants
    // TODO When filtering minor variants, evaluate strand bias appropriately

    label 'freebayes'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.callvar_freebayes.todo == 1

    input:
    tuple(val(sampleId), path(bam), path(bai), path(ref))
    
    output:
    path "varcall/ref/${ref.simpleName}.fna"
    path "varcall/raw/${ref.simpleName}/${sampleId}.vcf"
    tuple val(refName), val(sampleId), path("varcall/ref/${ref.simpleName}.fna"), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: gentsv
    tuple val(refName), val(sampleId), path("varcall/ref/${ref.simpleName}.fna"), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: callcons
    tuple val(refName), val(sampleId), path("varcall/ref/${ref.simpleName}.fna"), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: contaminated
    tuple val(refName), path("varcall/vcf/${ref.simpleName}/${sampleId}.vcf"), emit: contaminant
    
    script:
    sampleId = (bam =~ /(.+)\.(.+)/)[0][1]
    refName = (ref =~ /(.+)\.(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}

    mkdir -p work varcall/ref/ varcall/vcf/${ref.simpleName}/ varcall/raw/${ref.simpleName}/ work/${ref.simpleName}/
    perl -i.bkp -pe 's/\\r\$//' ${ref}

    if [[ ! -s ${bam} ]]; then
        touch "varcall/vcf/${ref.simpleName}/${sampleId}.vcf" "varcall/raw/${ref.simpleName}/${sampleId}.vcf"

        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    #Variant calling of the major variants, then vcf filtering and formatting.
    freebayes --theta ${params.callvar_freebayes["theta"]} --ploidy ${params.callvar_freebayes["ploidy"]} --report-all-haplotype-alleles --pooled-continuous --use-best-n-alleles ${params.callvar_freebayes["use-best-n-alleles"]} --allele-balance-priors-off --haplotype-length ${params.callvar_freebayes["haplotype_length"]} --use-duplicate-reads --genotyping-max-iterations ${params.callvar_freebayes["genotyping-max-iterations"]} --genotyping-max-banddepth ${params.callvar_freebayes["genotyping-max-banddepth"]} --min-mapping-quality ${params.callvar_freebayes["min_mapping_quality"]} --min-base-quality ${params.callvar_freebayes["min_base_quality"]} -F 0.5 -C ${params.callvar_freebayes["min_var_depth"]} --min-coverage ${params.callvar_freebayes["min_depth"]} --f ${ref} -b ${bam} > "work/${ref.simpleName}/${sampleId}_cons.vcf"

    vt decompose -s "work/${ref.simpleName}/${sampleId}_cons.vcf" -o "work/${ref.simpleName}/${sampleId}_decomposed_snps.vcf"
    python2.7 ${params.nfpath}/script/recalc_fbvcf.py -i "work/${ref.simpleName}/${sampleId}_decomposed_snps.vcf" -o "work/${ref.simpleName}/${sampleId}_recalc_snps.vcf" -n 0.0 -x 1.0 -t snp,del,ins,mnp,complex
    bcftools filter -i 'STRLEN(REF)==STRLEN(ALT) & AF >= 0.5 & ${params.callvar_freebayes["bcftools_filter_snps"]}' "work/${ref.simpleName}/${sampleId}_recalc_snps.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_snps.vcf"
    #bcftools filter -i 'STRLEN(REF)==STRLEN(ALT) & AF >= 0.5 & ((SAF>=${params.callvar_freebayes["saf"]} & SAR>=${params.callvar_freebayes["sar"]}) | AF >= 0.9) & INFO/QA>=${params.callvar_freebayes["qa"]}' "work/${ref.simpleName}/${sampleId}_recalc_snps.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_snps.vcf"
    vt normalize "work/${ref.simpleName}/${sampleId}_filtered_snps.vcf" -r "${ref}" -o "work/${ref.simpleName}/${sampleId}_normalized_snps.vcf"

    bgzip < "work/${ref.simpleName}/${sampleId}_normalized_snps.vcf" > "work/${ref.simpleName}/${sampleId}_snps.vcf.gz"
    bcftools index "work/${ref.simpleName}/${sampleId}_snps.vcf.gz"
    cat ${ref} | bcftools consensus "work/${ref.simpleName}/${sampleId}_snps.vcf.gz" > "work/${ref.simpleName}/${sampleId}_snps.fna"

    vt decompose -s "work/${ref.simpleName}/${sampleId}_cons.vcf" -o "work/${ref.simpleName}/${sampleId}_decomposed_complex.vcf"
    python2.7 ${params.nfpath}/script/recalc_fbvcf.py -i "work/${ref.simpleName}/${sampleId}_decomposed_complex.vcf" -o "work/${ref.simpleName}/${sampleId}_recalc_complex.vcf" -n 0.0 -x 1.0 -t snp,del,ins,mnp,complex
    bcftools filter -i 'STRLEN(REF)!=STRLEN(ALT) & AF >= 0.5 & ${params.callvar_freebayes["bcftools_filter_complex"]}' "work/${ref.simpleName}/${sampleId}_recalc_complex.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_complex.vcf"
    vt normalize "work/${ref.simpleName}/${sampleId}_filtered_complex.vcf" -r "${ref}" -o "work/${ref.simpleName}/${sampleId}_normalized_complex.vcf"

    #Variant calling of the minor variants, then vcf filtering and formatting. Permits to verify SNPs reversion as minor variants.
    freebayes --theta ${params.callvar_freebayes["theta"]} --ploidy ${params.callvar_freebayes["ploidy"]} --report-all-haplotype-alleles --pooled-continuous --use-best-n-alleles ${params.callvar_freebayes["use-best-n-alleles"]} --allele-balance-priors-off --haplotype-length ${params.callvar_freebayes["haplotype_length"]} --use-duplicate-reads --genotyping-max-iterations ${params.callvar_freebayes["genotyping-max-iterations"]} --genotyping-max-banddepth ${params.callvar_freebayes["genotyping-max-banddepth"]} --min-mapping-quality ${params.callvar_freebayes["min_mapping_quality"]} --min-base-quality ${params.callvar_freebayes["min_base_quality"]} -F 0.01 -C ${params.callvar_freebayes["min_var_depth"]} --min-coverage ${params.callvar_freebayes["min_depth"]} --f "work/${ref.simpleName}/${sampleId}_snps.fna" -b ${bam} > "work/${ref.simpleName}/${sampleId}_var.vcf"

    vt decompose -s "work/${ref.simpleName}/${sampleId}_var.vcf" -o "work/${ref.simpleName}/${sampleId}_decomposed_var.vcf"
    python2.7 ${params.nfpath}/script/recalc_fbvcf.py -i "work/${ref.simpleName}/${sampleId}_decomposed_var.vcf" -o "work/${ref.simpleName}/${sampleId}_recalc_var.vcf" -n 0.0 -x 0.499 -t snp,del,ins,mnp,complex
    bcftools filter -i 'AF < 0.5 & ${params.callvar_freebayes["bcftools_filter_minor"]}' "work/${ref.simpleName}/${sampleId}_recalc_var.vcf" > "work/${ref.simpleName}/${sampleId}_filtered_var.vcf"
    vt normalize "work/${ref.simpleName}/${sampleId}_filtered_var.vcf" -r "work/${ref.simpleName}/${sampleId}_snps.fna" -o "work/${ref.simpleName}/${sampleId}_normalized_var.vcf"

    #Prepare final filtered vcf file
    cp "work/${ref.simpleName}/${sampleId}_normalized_snps.vcf" "work/${ref.simpleName}/${sampleId}_unsorted_cons.vcf"
    grep -v '#' "work/${ref.simpleName}/${sampleId}_normalized_complex.vcf" >> "work/${ref.simpleName}/${sampleId}_unsorted_cons.vcf"
    bcftools sort -Ov -o "work/${ref.simpleName}/${sampleId}_sorted_cons.vcf" "work/${ref.simpleName}/${sampleId}_unsorted_cons.vcf"

    cp "work/${ref.simpleName}/${sampleId}_sorted_cons.vcf" "work/${ref.simpleName}/${sampleId}_sorted_concat.vcf"
    grep -v '#' "work/${ref.simpleName}/${sampleId}_normalized_var.vcf" >> "work/${ref.simpleName}/${sampleId}_sorted_concat.vcf"
    bcftools sort -Ov -o "work/${ref.simpleName}/${sampleId}_sorted.vcf" "work/${ref.simpleName}/${sampleId}_sorted_concat.vcf"
    vcfuniq "work/${ref.simpleName}/${sampleId}_sorted.vcf" > "varcall/vcf/${ref.simpleName}/${sampleId}.vcf"

    #Prepare vcf file without any filter post-calling
    vt normalize "work/${ref.simpleName}/${sampleId}_recalc_snps.vcf" -r "${ref}" -o "work/${ref.simpleName}/${sampleId}_normalized_raw_snps.vcf"
    vt normalize "work/${ref.simpleName}/${sampleId}_recalc_complex.vcf" -r "${ref}" -o "work/${ref.simpleName}/${sampleId}_normalized_raw_complex.vcf"
    vt normalize "work/${ref.simpleName}/${sampleId}_recalc_var.vcf" -r "work/${ref.simpleName}/${sampleId}_snps.fna" -o "work/${ref.simpleName}/${sampleId}_normalized_raw_var.vcf"
    cp "work/${ref.simpleName}/${sampleId}_normalized_raw_snps.vcf" "work/${ref.simpleName}/${sampleId}_unsorted_raw_cons.vcf"
    grep -v '#' "work/${ref.simpleName}/${sampleId}_normalized_raw_complex.vcf" >> "work/${ref.simpleName}/${sampleId}_unsorted_raw_cons.vcf"
    bcftools sort -Ov -o "work/${ref.simpleName}/${sampleId}_sorted_raw_cons.vcf" "work/${ref.simpleName}/${sampleId}_unsorted_raw_cons.vcf"

    cp "work/${ref.simpleName}/${sampleId}_sorted_raw_cons.vcf" "work/${ref.simpleName}/${sampleId}_sorted_raw_concat.vcf"
    grep -v '#' "work/${ref.simpleName}/${sampleId}_normalized_raw_var.vcf" >> "work/${ref.simpleName}/${sampleId}_sorted_raw_concat.vcf"
    bcftools sort -Ov -o "work/${ref.simpleName}/${sampleId}_sorted_raw.vcf" "work/${ref.simpleName}/${sampleId}_sorted_raw_concat.vcf"
    vcfuniq "work/${ref.simpleName}/${sampleId}_sorted_raw.vcf" > "varcall/raw/${ref.simpleName}/${sampleId}.vcf"

    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    #rm -rf work
    """
}

process genbed_bedtools {

    // Generate coverage files in two different format

    label 'freebayes'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.genbed_bedtools.todo == 1

    input:
    tuple(val(sampleId), path(bam), path(bai), path(ref))
    
    output:
    path "varcall/ref/${ref.simpleName}.fna"
    path "varcall/bga/${ref.simpleName}/${sampleId}.bed"
    tuple val(refName), val(sampleId), path("varcall/bga/${ref.simpleName}/${sampleId}.bed"), emit: callcons
    tuple val(refName), val(sampleId), path("varcall/bga/${ref.simpleName}/${sampleId}.bed"), emit: contaminated
    tuple val(refName), path("varcall/bga/${ref.simpleName}/${sampleId}.cov"), emit: covsummary
    
    script:
    sampleId = (bam =~ /(.+)\.(.+)/)[0][1]
    refName = (ref =~ /(.+)\.(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}

    mkdir -p work varcall/ref/ varcall/cons/${ref.simpleName}/ varcall/vcf/${ref.simpleName}/ work/${ref.simpleName}/ varcall/raw/${ref.simpleName}/ varcall/bga/${ref.simpleName}/
    perl -i.bkp -pe 's/\\r\$//' ${ref}

    if [[ ! -s ${bam} ]]; then
        touch "varcall/bga/${ref.simpleName}/${sampleId}.bed"

        length=\$(grep -v ">" ${ref} | tr -d '\n' | wc -m)
        for line in \$(grep ">" ${ref} | tr -d '>'); do echo -e "${ref.simpleName}\\tmean_depth_\${line}\\t0.00\\t${sampleId}" >> "varcall/bga/${ref.simpleName}/${sampleId}.cov"; done;

        cp "${ref}" "varcall/ref/${ref.simpleName}.fna"

        exit 0
    fi

    bedtools genomecov -ibam ${bam} -bga > "varcall/bga/${ref.simpleName}/${sampleId}.bed"
    for line in \$(samtools coverage ${bam} | tail -n +2 | tr '\\t' ','); do SEG="\${line%%,*}"; echo -e "${ref.simpleName}\\tmean_depth_\${SEG}\\t\$(echo \$line | cut -d ',' -f 7)\\t${sampleId}" >> "varcall/bga/${ref.simpleName}/${sampleId}.cov"; done;
    cp "${ref}" "varcall/ref/${ref.simpleName}.fna"
    #rm -rf work
    """
}

process genfasta_consensus {

    // Generate a consensus sequence, then mask regions with an insufficient depth and following an additional bed file

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }

    when:
    params.genfasta_consensus.todo == 1

    input:
    tuple(val(refName), val(sampleId), path(ref), path(vcf), path(bga))
    
    output:
    tuple val(refName), path("varcall/cons/${ref.simpleName}/${sampleId}.fna"), emit: consensus

    
    script:
    sampleId = (vcf =~ /(.+)\.(.+)/)[0][1]
    refName = (ref =~ /(.+)\.(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}

    mkdir -p work/${ref.simpleName}/ varcall/cons/${ref.simpleName}/
    perl -i.bkp -pe 's/\\r\$//' ${ref}

    if [[ ! -s ${vcf} ]]; then
        echo ">${sampleId}" > "varcall/cons/${ref.simpleName}/${sampleId}.fna"
        seq \$(grep -v ">" ${ref} | tr -d '\n' | wc -m) | sed "c N" | tr -d '\n' | sed -e '\$a\\' >> "varcall/cons/${ref.simpleName}/${sampleId}.fna"

        exit 0
    fi

    awk '\$4 < ${params.genfasta_consensus["min_depth"]}' "${sampleId}.bed" > "work/${ref.simpleName}/${sampleId}.bed"


    if [[ "${params.genfasta_consensus["mask"]}" != "" ]]; then
        cat "${params.genfasta_consensus["mask"]}" >> "work/${ref.simpleName}/${sampleId}.bed"
    fi

    bcftools filter -i 'AF >= 0.5' "${sampleId}.vcf" > "work/${ref.simpleName}/${sampleId}_unsorted_cons.vcf"
    bcftools sort -Ov -o "work/${ref.simpleName}/${sampleId}_sorted_cons.vcf" "work/${ref.simpleName}/${sampleId}_unsorted_cons.vcf"
    bgzip < "work/${ref.simpleName}/${sampleId}_sorted_cons.vcf" > "work/${ref.simpleName}/${sampleId}_cons.vcf.gz"
    bcftools index "work/${ref.simpleName}/${sampleId}_cons.vcf.gz"
    bedtools maskfasta -fi ${ref} -bed "work/${ref.simpleName}/${sampleId}.bed" -fo "work/${ref.simpleName}/${sampleId}_cons_masked.fna" -soft
    cat "work/${ref.simpleName}/${sampleId}_cons_masked.fna" | bcftools consensus "work/${ref.simpleName}/${sampleId}_cons.vcf.gz" > "work/${ref.simpleName}/${sampleId}_cons.fna"
    sed 's/\\(>[^|]*\\)|\\(.*\\)/>${sampleId}|\\2/' "work/${ref.simpleName}/${sampleId}_cons.fna" > "work/${ref.simpleName}/${sampleId}_rehead.fna"
    python2.7 ${params.nfpath}/script/fasta_masklowercase.py -i "work/${ref.simpleName}/${sampleId}_rehead.fna" -o "varcall/cons/${ref.simpleName}/${sampleId}.fna"

    #rm -rf work
    """
}

process seek_coinf {

    // Compare vcf files to count minor/major variants in a sample matching the profile expected for each lineage, the best match will be considered "major lineage" of the sample.
    // Then count minor variants matching the profile expected for each lineage, excluding variants matching the "major lineage", the best match will be considered "minor lineage" of the sample and indicate a potential co-infection.

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { sampleId }
    beforeScript 'ulimit -Ss unlimited'

    when:
    params.seek_coinf.todo == 1 && target.size() > 0 && refName == (bed =~ /([^_]+)_([0-9a-zA-Z_\-]+)(.+)/)[0][1]

    input:
    tuple(val(refName), val(sampleId), path(ref), path(target), path(bga), val(vcf), path(bed))
    
    output:
    tuple val(refName), val(region), path("varcall/coinf/${refName}/${sampleId}_${region}_min.tsv"), emit: min
    tuple val(refName), val(region), path("varcall/coinf/${refName}/${sampleId}_${region}_maj.tsv"), emit: maj
    tuple val(refName), val(region), path("varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv"), emit: raw
    
    script:
    region = (bed =~ /([0-9a-zA-Z_\-]+)(.+)/)[0][1]
    """
    #!/bin/bash
    export TMPDIR=/srv/scratch/iai/seqmet/tmp/
    mkdir -p work varcall/coinf/${refName} varcall/contadb/${refName} varcall/raw/${refName}
    tbn=\$(basename "${target}")

    for file in ${vcf}*.vcf; do
    bn=\$(basename "\${file}")
    if [[ \${bn} != \${tbn} && -s \${file} ]]; then
        python3 ${params.nfpath}/script/compare_vcf.py --reference \${file} --target ${target} --depth ${bga} --bed ${bed} --min_depth ${params.seek_coinf["min_depth"]} --min_freq ${params.seek_coinf["min_freq"]} --mode "maj" --output "varcall/coinf/${refName}/\${tbn%%.*}_${region}_maj.tsv"
    fi
    done

    echo -e "${refName}\\tcoinf_maj_match\\t\$(awk -F "\\t" 'BEGIN {OFS="\\t"}; \$4>=${params.seek_coinf["min_var"]} {print \$0} \$4<${params.seek_coinf["min_var"]} {\$6=0.0; print \$0}' "varcall/coinf/${refName}/\${tbn%%.*}_${region}_maj.tsv" | sort -t\$'\\t' -k6,6 -k4,4 -nr | awk -F "\\t" '{print \$1}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv"
    echo -e "${refName}\\tcoinf_maj_common\\t\$(awk -F "\\t" 'BEGIN {OFS="\\t"}; \$4>=${params.seek_coinf["min_var"]} {print \$0} \$4<${params.seek_coinf["min_var"]} {\$6=0.0; print \$0}' "varcall/coinf/${refName}/\${tbn%%.*}_${region}_maj.tsv" | sort -t\$'\\t' -k6,6 -k4,4 -nr | awk -F "\\t" '{print \$4}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv"
    echo -e "${refName}\\tcoinf_maj_ratio\\t\$(awk -F "\\t" 'BEGIN {OFS="\\t"}; \$4>=${params.seek_coinf["min_var"]} {print \$0} \$4<${params.seek_coinf["min_var"]} {\$6=0.0; print \$0}' "varcall/coinf/${refName}/\${tbn%%.*}_${region}_maj.tsv" | sort -t\$'\\t' -k6,6 -k4,4 -nr | awk -F "\\t" '{print \$6}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv"
    echo -e "${refName}\\tcoinf_maj_median\\t\$(awk -F "\\t" 'BEGIN {OFS="\\t"}; \$4>=${params.seek_coinf["min_var"]} {print \$0} \$4<${params.seek_coinf["min_var"]} {\$6=0.0; print \$0}' "varcall/coinf/${refName}/\${tbn%%.*}_${region}_maj.tsv" | sort -t\$'\\t' -k6,6 -k4,4 -nr | awk -F "\\t" '{print \$7}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv"
    echo -e "${refName}\\tcoinf_maj_iqr\\t\$(awk -F "\\t" 'BEGIN {OFS="\\t"}; \$4>=${params.seek_coinf["min_var"]} {print \$0} \$4<${params.seek_coinf["min_var"]} {\$6=0.0; print \$0}' "varcall/coinf/${refName}/\${tbn%%.*}_${region}_maj.tsv" | sort -t\$'\\t' -k6,6 -k4,4 -nr | awk -F "\\t" '{print \$8}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv"

    MAJ_MATCH=\$(cat "varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv" | grep 'coinf_maj_match' | awk -F "\\t" '{print \$3}')

    for file in ${vcf}*.vcf; do
    bn=\$(basename "\${file}")
    if [[ \${bn} != \${tbn} && -s \${file} ]]; then
        python3 ${params.nfpath}/script/compare_vcf.py --reference \${file} --target ${target} --exclusion "${vcf}\${MAJ_MATCH}.vcf" --depth ${bga} --bed ${bed} --min_depth ${params.seek_coinf["min_depth"]} --min_freq ${params.seek_coinf["min_freq"]} --mode "min" --output "varcall/coinf/${refName}/\${tbn%%.*}_${region}_min.tsv"
    elif [[ \${MAJ_MATCH} == "NA" && -s \${file} ]]; then
        python3 ${params.nfpath}/script/compare_vcf.py --reference \${file} --target ${target} --depth ${bga} --bed ${bed} --min_depth ${params.seek_coinf["min_depth"]} --min_freq ${params.seek_coinf["min_freq"]} --mode "min" --output "varcall/coinf/${refName}/\${tbn%%.*}_${region}_min.tsv"
    fi
    done

    echo -e "${refName}\\tcoinf_min_match\\t\$(awk -F "\\t" 'BEGIN {OFS="\\t"}; \$4>=${params.seek_coinf["min_var"]} {print \$0} \$4<${params.seek_coinf["min_var"]} {\$1="NA"; \$6=0.0; print \$0}' "varcall/coinf/${refName}/\${tbn%%.*}_${region}_min.tsv" | sort -t\$'\\t' -k6,6 -k4,4 -nr | awk -F "\\t" '{print \$1}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv"
    echo -e "${refName}\\tcoinf_min_common\\t\$(awk -F "\\t" 'BEGIN {OFS="\\t"}; \$4>=${params.seek_coinf["min_var"]} {print \$0} \$4<${params.seek_coinf["min_var"]} {\$6=0.0; print \$0}' "varcall/coinf/${refName}/\${tbn%%.*}_${region}_min.tsv" | sort -t\$'\\t' -k6,6 -k4,4 -nr | awk -F "\\t" '{print \$4}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv"
    echo -e "${refName}\\tcoinf_min_ratio\\t\$(awk -F "\\t" 'BEGIN {OFS="\\t"}; \$4>=${params.seek_coinf["min_var"]} {print \$0} \$4<${params.seek_coinf["min_var"]} {\$6=0.0; print \$0}' "varcall/coinf/${refName}/\${tbn%%.*}_${region}_min.tsv" | sort -t\$'\\t' -k6,6 -k4,4 -nr | awk -F "\\t" '{print \$6}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv"
    echo -e "${refName}\\tcoinf_min_median\\t\$(awk -F "\\t" 'BEGIN {OFS="\\t"}; \$4>=${params.seek_coinf["min_var"]} {print \$0} \$4<${params.seek_coinf["min_var"]} {\$7=0.0; print \$0}' "varcall/coinf/${refName}/\${tbn%%.*}_${region}_min.tsv" | sort -t\$'\\t' -k6,6 -k4,4 -nr | awk -F "\\t" '{print \$7}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv"
    echo -e "${refName}\\tcoinf_min_iqr\\t\$(awk -F "\\t" 'BEGIN {OFS="\\t"}; \$4>=${params.seek_coinf["min_var"]} {print \$0} \$4<${params.seek_coinf["min_var"]} {\$8=0.0; print \$0}' "varcall/coinf/${refName}/\${tbn%%.*}_${region}_min.tsv" | sort -t\$'\\t' -k6,6 -k4,4 -nr | awk -F "\\t" '{print \$8}' | head -n 1 | tr -d '\\n')\\t${sampleId}" >> "varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv"

    touch  varcall/coinf/${refName}/${sampleId}_${region}_maj.tsv varcall/coinf/${refName}/${sampleId}_${region}_min.tsv varcall/raw/${refName}/${sampleId}_${region}_coinf.tsv
    """
}

process merge_depth {

    // Merge depth files generated for each sample

    label 'varcall'
    storeDir params.result
    echo false

    when:
    params.merge_depth.todo == 1

    input:
    file "varcall/premap/depth/*"
    
    output:
    file "varcall/premap_depth.tsv"
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    python2.7 ${params.nfpath}/script/merge_tables.py varcall/premap/depth/*.tsv > varcall/premap_depth.tsv
    """
}

process merge_coverage {

    // Merge premapping coverage files generated for each sample

    label 'varcall'
    storeDir params.result
    echo false

    when:
    params.merge_coverage.todo == 1

    input:
    file "varcall/premap/coverage/*"
    
    output:
    file "premap_coverage.tsv"
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    python2.7 ${params.nfpath}/script/merge_tables.py varcall/premap/coverage/*.tsv > premap_coverage.tsv
    """
}

process merge_poscsummary {

    // Merge premapping positive controls coverage files generated for each sample

    label 'varcall'
    storeDir (params.result)
    echo false

    when:
    params.merge_poscsummary.todo == 1

    input:
    file 'posc/*'
    
    output:
    file "premap_posc.tsv"
    
    script:
    """
    python2.7 ${params.nfpath}/script/merge_tables.py posc/*.tsv > premap_posc.tsv
    """
}

process merge_consensus {

    // Merge consensus files generated for the same reference

    label 'varcall'
    storeDir params.result
    echo false
    tag { refName }

    when:
    params.merge_consensus.todo == 1

    input:
    tuple(val(refName), path("varcall/cons/${refName}/*"))
    
    output:
    tuple val(refName), path("${refName}_consensus.fna")
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    cat varcall/cons/${refName}/*.fna > ${refName}_consensus.fna
    """
}

process merge_covsummary {

    // Merge coverage files generated for the same reference

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { refName }

    when:
    params.merge_covsummary.todo == 1

    input:
    tuple(val(refName), path("*"))

    output:
    tuple val(refName), path("varcall/${refName}_cov.cov")

    script:
    """
    mkdir -p varcall/
    cat *.cov > varcall/${refName}_cov.cov
    """
}

process merge_coinfsummary {

    // Merge coinf files generated for the same reference

    label 'varcall'
    storeDir (params.result)
    echo false
    tag { refName }

    when:
    params.merge_coinfsummary.todo == 1

    input:
    tuple(val(refName), val(region), path("*"))

    output:
    tuple val(refName), path("varcall/${refName}_coinf.tsv")

    script:
    """
    mkdir -p varcall/
    cat *.tsv > varcall/${refName}_coinf.tsv
    """
}

process gentsv_summary {

    // Generate a summary of the run for each reference used

    label 'varcall'
    storeDir (params.result)
    echo false

    when:
    params.gentsv_summary.todo == 1

    input:
    tuple val(refName), path("*")

    output:
    tuple val(refName), path("${refName}_summary.tsv")

    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8
    mkdir -p work


    CONSENSUS=\$(find ./ -name *_consensus.fna)
    COV=\$(find ./ -name *_cov.cov)
    COINF=\$(find ./ -name ${refName}_coinf.tsv)
    ARG="--output ${refName}_summary.tsv --runname=${params.prefix} --cons=\${CONSENSUS} --coverage=\${COV}"
    if [[ -s premap_posc.tsv ]]; then ARG+=" --posc premap_posc.tsv"; fi
    if [[ -s \${COINF} && ${refName} == "MN908947" ]]; then ARG+=" --coinf-table \${COINF}"; fi
    echo \$ARG
    /usr/bin/Rscript ${params.nfpath}/script/make_summary.R \$ARG
    """
}

process gentsv_pangolin {

    // Generate a pangolin report for a given fasta file

    label 'pangolin'
    storeDir params.result
    echo false

    when:
    params.gentsv_pangolin.todo == 1

    input:
    tuple(val(refName), path(consensus))
    
    output:
    tuple val(refName), path("${refName}_pangolin.csv")
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8

    if [[ ${refName} == "CY115151" || ${refName} == "CY115183" || ${refName} == "EPIISL710475" || ${refName} == "EPIISL517733" ]]; then
        cat ${consensus} | perl -i.bkp -pe 's/\\r\$//' | awk '/^>/ {printf("%s%s\\t",(N>0?"\\n":""),\$0);N++;next;} {printf("%s",\$0);} END {printf("\\n");}' | tr "\\t" "\\n" | grep '|S4' -A 1 | sed 's/|S4//' | grep -v "^--\$" > consensus.fasta
    else
        cat ${consensus} | perl -i.bkp -pe 's/\\r\$//' | awk '/^>/ {printf("%s%s\\t",(N>0?"\\n":""),\$0);N++;next;} {printf("%s",\$0);} END {printf("\\n");}' | tr "\\t" "\\n" | grep '|WG' -A 1 | sed 's/|WG//' | grep -v "^--\$" > consensus.fasta
    fi
    
    pangolin consensus.fasta --outfile "${refName}_pangolin.csv"
    """
}

process gentsv_nextclade {

    // Generate a nextclade report for a given fasta file

    label 'nextclade'
    storeDir params.result
    echo false

    when:
    params.gentsv_nextclade.todo == 1

    input:
    tuple(val(refName), path(consensus), val(dataset))
    
    output:
    tuple val(refName), path("${refName}_nextclade.tsv")
    
    script:
    """
    #!/bin/bash
    export TMPDIR=${params.tmp}
    export LANG=en_US.utf8
    export LC_ALL=en_US.utf8
    
    if [[ ${refName} == "CY115151" || ${refName} == "CY115183" || ${refName} == "EPIISL710475" || ${refName} == "EPIISL517733" ]]; then
        cat ${consensus} | perl -i.bkp -pe 's/\\r\$//' | awk '/^>/ {printf("%s%s\\t",(N>0?"\\n":""),\$0);N++;next;} {printf("%s",\$0);} END {printf("\\n");}' | tr "\\t" "\\n" | grep '|S4' -A 1 | sed 's/|S4//' | grep -v "^--\$" > consensus.fasta
    else
        cat ${consensus} | perl -i.bkp -pe 's/\\r\$//' | awk '/^>/ {printf("%s%s\\t",(N>0?"\\n":""),\$0);N++;next;} {printf("%s",\$0);} END {printf("\\n");}' | tr "\\t" "\\n" | grep '|WG' -A 1 | sed 's/|WG//' | grep -v "^--\$" > consensus.fasta
    fi
    
    nextclade run --input-dataset ${dataset} --include-reference --in-order --input-fasta consensus.fasta --output-tsv "${refName}_nextclade.tsv"
    """
}
