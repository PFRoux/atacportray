//
// TELOMERE branch: TelomereHunter telomere-content estimation on the shared BAM.
// Reproduces `telomerehunter -ibt <bam> -b <cytoBand>` from the notebook.
//
include { TELOMEREHUNTER } from '../../../modules/local/telomerehunter/main'

workflow BAM_TELOMERE_TELOMEREHUNTER {

    take:
    ch_bam       // channel: [ val(meta), path(bam), path(bai) ]
    ch_cytoband  // channel: path(cytoBand.txt)

    main:
    ch_versions = Channel.empty()

    TELOMEREHUNTER ( ch_bam, ch_cytoband )
    ch_versions = ch_versions.mix(TELOMEREHUNTER.out.versions)

    emit:
    summary    = TELOMEREHUNTER.out.summary
    outdir     = TELOMEREHUNTER.out.outdir
    plots      = TELOMEREHUNTER.out.plots
    versions   = ch_versions
}
