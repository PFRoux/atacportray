//
// CNV branch: QDNAseq copy-number profiling on the shared analysis-ready BAM.
// Reproduces the exact notebook workflow (100kb bins, sqrt segmentation,
// loss=-0.4 / gain=+0.4) inside modules/local/qdnaseq.
//
include { QDNASEQ } from '../../../modules/local/qdnaseq/main'

workflow BAM_CNV_QDNASEQ {

    take:
    ch_bam        // channel: [ val(meta), path(bam), path(bai) ]
    binsize_kb    // value:   100
    loss_threshold// value:   -0.4
    gain_threshold// value:   0.4

    main:
    ch_versions = Channel.empty()

    QDNASEQ ( ch_bam, binsize_kb, loss_threshold, gain_threshold )
    ch_versions = ch_versions.mix(QDNASEQ.out.versions)

    emit:
    calls      = QDNASEQ.out.calls
    segments   = QDNASEQ.out.segments
    plots      = QDNASEQ.out.plots
    versions   = ch_versions
}
