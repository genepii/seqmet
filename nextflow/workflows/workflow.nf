#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// import modules

include {
    filter_sra;
    downsample_seqtk;
    trim_cutadapt;
    premap_minimap2;
    map_minimap2;
    postmap_picard;
    postmap_abra2;
    postmap_ampliconclip;
    callvar_freebayes;
    genbed_bedtools;
    genfasta_consensus;
    seek_coinf;
    merge_depth;
    merge_coverage;
    merge_consensus;
    merge_covsummary;
    merge_coinfsummary;
    merge_poscsummary;
    gentsv_summary;
    gentsv_pangolin;
    gentsv_nextclade
} from "${params.nfpath}/nextflow/modules/module.nf"

workflow varcall {
    take:
        ch_filePairs

    main:

        if ( params.filter ) {
            filter_sra(ch_filePairs)
            filter_sra.out.set{ ch_filePairs }
        }

        if ( params.downsample ) {
            downsample_seqtk(ch_filePairs)
            downsample_seqtk.out.set{ ch_filePairs }
        }

        if ( params.trim ) {
            trim_cutadapt(ch_filePairs)
            trim_cutadapt.out.set{ ch_filePairs }
        }

        if ( params.premap == 'minimap2' ) {
            premap_minimap2(ch_filePairs, Channel.fromPath(params.premap_minimap2.ref).toSortedList().flatten().collect(), Channel.fromPath(params.premap_minimap2.posc).toSortedList().flatten().collect())
            premap_minimap2.out.depthtable.set{ ch_premap_depthtable }
            premap_minimap2.out.coveragetable.set{ ch_premap_coveragetable }
            premap_minimap2.out.poscsummary.set{ ch_premap_poscsummary }
        }

        if ( params.map == 'minimap2' ) {
            map_minimap2(ch_filePairs.join(ch_premap_coveragetable, by: 0).combine(Channel.fromPath(params.map_minimap2.ref)))
            map_minimap2.out.set{ ch_bamTuple }
        }

        if ( params.rmdup ) {
            postmap_picard(ch_bamTuple)
            postmap_picard.out.set{ ch_bamTuple }
        }

        if ( params.realign ) {
            postmap_abra2(ch_bamTuple)
            postmap_abra2.out.set{ ch_bamTuple }
        }

        if ( params.clipbam == 'ampliconclip' ) {
            postmap_ampliconclip(ch_bamTuple, Channel.fromPath(params.postmap_ampliconclip.ampliconclip_bed).toSortedList().flatten().collect())
            postmap_ampliconclip.out.set{ ch_bamTuple }
        }

        callvar_freebayes(ch_bamTuple)

        genbed_bedtools(ch_bamTuple)

        genfasta_consensus(callvar_freebayes.out.callcons.join(genbed_bedtools.out.callcons, by:[0,1]))

        seek_coinf((callvar_freebayes.out.contaminated.join(genbed_bedtools.out.contaminated, by:[0,1])).combine(Channel.from(params.seek_coinf.vcf), by:[0]).combine(Channel.fromPath(params.seek_coinf.bed)))

        merge_depth(ch_premap_depthtable.map{ it[1] }.collect())

        merge_coverage(ch_premap_coveragetable.map{ it[1] }.collect())

        merge_poscsummary(ch_premap_poscsummary.collect())

        merge_consensus(genfasta_consensus.out.consensus.groupTuple())

        merge_covsummary(genbed_bedtools.out.covsummary.groupTuple())

        merge_coinfsummary(seek_coinf.out.raw.groupTuple())

        merge_consensus.out.set{ ch_summary }
        ch_summary.join(merge_covsummary.out).set{ ch_summary }
        ch_summary.combine(merge_poscsummary.out).set{ ch_summary }
        ch_summary.join(merge_coinfsummary.out).set{ ch_summary }

        gentsv_summary(ch_summary.map{ [it[0], it[1..-1]] })

        gentsv_nextclade(merge_consensus.out.groupTuple().join(Channel.from(params.gentsv_nextclade.dataset)))

        gentsv_pangolin(merge_consensus.out.groupTuple())

}