//
// MITO branch: bulk mgatk mitochondrial variant calling on the shared BAMs.
// Reproduces `mgatk call --g hg38 --keep-duplicates` from the notebook.
//
include { MGATK } from '../../../modules/local/mgatk/main'
include { MGATK_POSTPROCESS } from '../../../modules/local/mgatk_postprocess/main'

workflow BAM_MITO_MGATK {

    take:
    ch_bam          // channel: [ val(meta), path(bam), path(bai) ]
    mito_contig     // value:   'chrM'
    keep_duplicates // value:   true
    run_postprocess // value:   true
    min_af          // value:   0.05
    min_vmr         // value:   0.01
    min_sd          // value:   0.05

    main:
    ch_versions = Channel.empty()
    ch_post_summary = Channel.empty()
    ch_post_high_confidence = Channel.empty()
    ch_post_matrix = Channel.empty()
    ch_post_mqc = Channel.empty()
    ch_post_coverage_table = Channel.empty()
    ch_post_coverage_mqc = Channel.empty()
    ch_post_coverage_circular_mqc = Channel.empty()
    ch_post_coverage_circular = Channel.empty()
    ch_post_heatmap = Channel.empty()
    ch_post_pca = Channel.empty()

    MGATK (
        ch_bam.map { meta, bam, bai -> bam }.collect(),
        ch_bam.map { meta, bam, bai -> bai }.collect(),
        mito_contig,
        keep_duplicates
    )
    ch_versions = ch_versions.mix(MGATK.out.versions)

    if ( run_postprocess ) {
        MGATK_POSTPROCESS (
            MGATK.out.rds.collect(),
            MGATK.out.coverage.collect(),
            MGATK.out.variant_stats.collect(),
            MGATK.out.vmr.collect(),
            min_af,
            min_vmr,
            min_sd
        )
        ch_post_summary = MGATK_POSTPROCESS.out.summary
        ch_post_high_confidence = MGATK_POSTPROCESS.out.high_confidence
        ch_post_matrix = MGATK_POSTPROCESS.out.matrix
        ch_post_mqc = MGATK_POSTPROCESS.out.mqc
        ch_post_coverage_table = MGATK_POSTPROCESS.out.coverage_table
        ch_post_coverage_mqc = MGATK_POSTPROCESS.out.coverage_mqc
        ch_post_coverage_circular_mqc = MGATK_POSTPROCESS.out.coverage_circular_mqc
        ch_post_coverage_circular = MGATK_POSTPROCESS.out.coverage_circular
        ch_post_heatmap = MGATK_POSTPROCESS.out.heatmap
        ch_post_pca = MGATK_POSTPROCESS.out.pca
        ch_versions = ch_versions.mix(MGATK_POSTPROCESS.out.versions)
    }

    emit:
    variant_stats   = MGATK.out.variant_stats
    vmr             = MGATK.out.vmr
    coverage        = MGATK.out.coverage
    rds             = MGATK.out.rds
    final_dir       = MGATK.out.final_dir
    outdir          = MGATK.out.outdir
    post_summary    = ch_post_summary
    post_variants   = ch_post_high_confidence
    post_matrix     = ch_post_matrix
    post_mqc        = ch_post_mqc
    post_coverage   = ch_post_coverage_table
    post_coverage_mqc = ch_post_coverage_mqc
    post_coverage_circular_mqc = ch_post_coverage_circular_mqc
    post_coverage_circular = ch_post_coverage_circular
    post_heatmap    = ch_post_heatmap
    post_pca        = ch_post_pca
    versions        = ch_versions
}
