//
// MITO branch: mgatk mitochondrial variant calling on the shared BAM.
// Reproduces `mgatk call --g hg38 --keep-duplicates` from the notebook.
//
include { MGATK } from '../../../modules/local/mgatk/main'

workflow BAM_MITO_MGATK {

    take:
    ch_bam          // channel: [ val(meta), path(bam), path(bai) ]
    mito_contig     // value:   'chrM'
    keep_duplicates // value:   true

    main:
    ch_versions = Channel.empty()

    MGATK ( ch_bam, mito_contig, keep_duplicates )
    ch_versions = ch_versions.mix(MGATK.out.versions)

    emit:
    variant_stats = MGATK.out.variant_stats
    final_dir     = MGATK.out.final_dir
    outdir        = MGATK.out.outdir
    versions      = ch_versions
}
