//
// EPIGENOME branch: bowtie2 -> filter -> blacklist -> coverage -> MACS3 peaks
//                   -> ROSE super-enhancers -> TOBIAS footprinting
// Reproduces the exact steps of the atacportray reference notebook.
//
include { BOWTIE2_ALIGN        } from '../../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_VIEW        } from '../../../modules/nf-core/samtools/view/main'
include { BEDTOOLS_INTERSECT   } from '../../../modules/nf-core/bedtools/intersect/main'
include { SAMTOOLS_SORT        } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX       } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS       } from '../../../modules/nf-core/samtools/stats/main'
include { DEEPTOOLS_BAMCOVERAGE} from '../../../modules/nf-core/deeptools/bamcoverage/main'
include { MACS3_CALLPEAK       } from '../../../modules/nf-core/macs3/callpeak/main'
include { ROSE                 } from '../../../modules/local/rose/main'
include { NUCLEOATAC           } from '../../../modules/local/nucleoatac/main'
include { CONSENSUS_PEAKS      } from '../../../modules/local/consensus_peaks/main'
include { CONSENSUS_SUPER_ENHANCERS } from '../../../modules/local/consensus_super_enhancers/main'
include { BEDTOOLS_MULTICOV_COUNTS as PEAKS_MULTICOV_COUNTS } from '../../../modules/local/bedtools_multicov_counts/main'
include { BEDTOOLS_MULTICOV_COUNTS as SUPER_ENHANCERS_MULTICOV_COUNTS } from '../../../modules/local/bedtools_multicov_counts/main'
include { ATAC_DESEQ2_QC      } from '../../../modules/local/atac_deseq2_qc/main'
include { TOBIAS_ATACORRECT    } from '../../../modules/local/tobias/atacorrect/main'
include { TOBIAS_SCOREBIGWIG   } from '../../../modules/local/tobias/scorebigwig/main'
include { TOBIAS_BINDETECT     } from '../../../modules/local/tobias/bindetect/main'
include { BED_SLOP             } from '../../../modules/local/bed_slop/main'
include { ATAC_INSERT_SIZE_MQC } from '../../../modules/local/atac_insert_size_mqc/main'
include { ATAC_FRIP_MQC        } from '../../../modules/local/atac_frip_mqc/main'
include { DEEPTOOLS_TSS_PROFILE} from '../../../modules/local/deeptools_tss_profile/main'

workflow FASTQ_EPIGENOME_BOWTIE2 {

    take:
    ch_reads          // channel: [ val(meta), [ reads ] ]  (already trimmed)
    ch_bowtie2_index  // channel: [ val(meta2), path(index) ]
    ch_fasta          // channel: [ val(meta3), path(fasta) ]
    ch_fai            // channel: [ val(meta3), path(fai) ]
    ch_blacklist      // channel: [ val(meta4), path(blacklist_bed) ]
    apply_blacklist   // value:   boolean (whether a blacklist BED was provided)
    macs3_gsize       // value:   e.g. 'hs'
    ch_motifs         // channel: path(motifs.jaspar)  (HOCOMOCO v11)
    run_rose          // value:   boolean
    rose_genome       // value:   e.g. HG38
    ch_rose_annotation // channel: path(rose_refseq.ucsc) or []
    rose_stitch       // value:   e.g. 12500
    rose_tss          // value:   e.g. 2500
    run_tobias        // value:   boolean
    run_nucleoatac    // value:   boolean
    peak_filter_slop  // value:   e.g. 500
    run_tss_profile   // value:   boolean
    ch_tss_regions    // channel: path(tss.bed)
    tss_before        // value:   e.g. 2000
    tss_after         // value:   e.g. 2000
    tss_binsize       // value:   e.g. 10
    run_deseq2_qc     // value:   boolean
    deseq2_min_count  // value:   e.g. 10
    deseq2_top_features // value: e.g. 50

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_count_tables  = Channel.empty()

    //
    // Align with bowtie2 (ext.args: -N 0 --very-sensitive-local, set in modules.config)
    //
    // NB: registry modules report versions via the `versions` topic channel
    //     (collected in the main workflow); only local modules emit versions.yml.
    BOWTIE2_ALIGN ( ch_reads, ch_bowtie2_index, ch_fasta, false, true )
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect{ it[1] }.ifEmpty([]))

    //
    // Filter unmapped reads (samtools view -F 4, ext.args in modules.config)
    //
    ch_view_in = BOWTIE2_ALIGN.out.bam.map { meta, bam -> [ meta, bam, [] ] }
    SAMTOOLS_VIEW ( ch_view_in, ch_fasta.map{ m, f -> [ m, f, [] ] }, [[:], []], [[:], []], '' )

    //
    // Remove ENCODE blacklist regions (bedtools intersect -v), if a blacklist is provided
    //
    if ( apply_blacklist ) {
        ch_bl_in = SAMTOOLS_VIEW.out.bam
            .combine( ch_blacklist.map { meta, bl -> bl } )
            .map { meta, bam, bl -> [ meta, bam, bl ] }
        BEDTOOLS_INTERSECT ( ch_bl_in, [[:], []] )
        ch_filtered_bam = BEDTOOLS_INTERSECT.out.intersect
    } else {
        ch_filtered_bam = SAMTOOLS_VIEW.out.bam
    }

    //
    // Sort + index the analysis-ready epigenome BAM
    //
    SAMTOOLS_SORT ( ch_filtered_bam, ch_fasta.map{ m,f -> [m,f,[]] }, '' )
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_bam_bai = SAMTOOLS_SORT.out.bam.join( SAMTOOLS_INDEX.out.index )
    ch_stats_ref = ch_fasta.combine(ch_fai)
    ch_stats_in = ch_bam_bai.combine(ch_stats_ref)
    SAMTOOLS_STATS (
        ch_stats_in.map { meta, bam, bai, meta_fasta, fasta, meta_fai, fai -> [ meta, bam, bai ] },
        ch_stats_in.map { meta, bam, bai, meta_fasta, fasta, meta_fai, fai -> [ meta_fasta, fasta, fai ] }
    )
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect{ it[1] }.ifEmpty([]))
    ATAC_INSERT_SIZE_MQC (
        ch_bam_bai.map { meta, bam, bai -> bam }.collect(),
        ch_bam_bai.map { meta, bam, bai -> bai }.collect()
    )
    ch_multiqc_files = ch_multiqc_files.mix(ATAC_INSERT_SIZE_MQC.out.mqc)
    ch_versions = ch_versions.mix(ATAC_INSERT_SIZE_MQC.out.versions)

    //
    // Coverage track (deepTools bamCoverage: -bs 20 --smoothLength 100 --extendReads)
    //
    DEEPTOOLS_BAMCOVERAGE (
        ch_bam_bai,
        ch_fasta.map { m, f -> f },
        ch_fai.map   { m, f -> f },
        ch_blacklist
    )

    if ( run_tss_profile ) {
        DEEPTOOLS_TSS_PROFILE (
            DEEPTOOLS_BAMCOVERAGE.out.bigwig.map { meta, bigwig -> bigwig }.collect(),
            ch_tss_regions,
            tss_before,
            tss_after,
            tss_binsize
        )
        ch_multiqc_files = ch_multiqc_files.mix(DEEPTOOLS_TSS_PROFILE.out.mqc)
        ch_versions = ch_versions.mix(DEEPTOOLS_TSS_PROFILE.out.versions)
    }

    //
    // Peak calling (MACS3 callpeak --format BAMPE --gsize hs --keep-dup all --qvalue 0.1)
    //
    ch_macs_in = SAMTOOLS_SORT.out.bam.map { meta, bam -> [ meta, bam, [] ] }
    MACS3_CALLPEAK ( ch_macs_in, macs3_gsize )
    ch_multiqc_files = ch_multiqc_files.mix(MACS3_CALLPEAK.out.xls.collect{ it[1] }.ifEmpty([]))
    ch_frip_in = ch_bam_bai.join( MACS3_CALLPEAK.out.peak )
        .map { meta, bam, bai, peak -> [ meta, bam, bai, peak ] }
    ATAC_FRIP_MQC ( ch_frip_in )
    ch_multiqc_files = ch_multiqc_files.mix(ATAC_FRIP_MQC.out.mqc.map { meta, mqc -> mqc })
    ch_versions = ch_versions.mix(ATAC_FRIP_MQC.out.versions)

    //
    // Consensus peak set used by footprinting, variant filtering and CNV masking
    //
    CONSENSUS_PEAKS ( MACS3_CALLPEAK.out.peak.collect{ it[1] } )
    ch_consensus_peaks = CONSENSUS_PEAKS.out.bed.map { bed -> [ [id:'consensus'], bed ] }
    BED_SLOP ( ch_consensus_peaks, ch_fai, peak_filter_slop )
    ch_multicov_bams = ch_bam_bai
        .collect(flat: false)
        .map { rows ->
            def sorted = rows.sort { a, b -> a[0]['id'] <=> b[0]['id'] }
            [ sorted.collect { it[0]['id'] }, sorted.collect { it[1] }, sorted.collect { it[2] } ]
        }
    ch_peak_counts_in = ch_consensus_peaks
        .combine( ch_multicov_bams )
        .map { meta, bed, sample_names, bams, bais -> [ [id:'consensus_peaks'], bed, sample_names, bams, bais ] }
    PEAKS_MULTICOV_COUNTS ( ch_peak_counts_in )
    ch_count_tables = ch_count_tables.mix(PEAKS_MULTICOV_COUNTS.out.counts)
    if ( run_deseq2_qc ) {
        ch_deseq2_meta = ch_bam_bai
            .collect(flat: false)
            .map { rows ->
                def sorted = rows.sort { a, b -> a[0]['id'] <=> b[0]['id'] }
                sorted.collect { row ->
                    [
                        sample: row[0]['id'],
                        condition: row[0].containsKey('condition') ? row[0]['condition'] : 'NA',
                        patient: row[0].containsKey('patient') ? row[0]['patient'] : 'NA'
                    ]
                }
        }
        ch_deseq2_in = PEAKS_MULTICOV_COUNTS.out.counts
            .combine( ch_deseq2_meta )
            .map { row ->
                def sample_info = row[2..-1]
                [ row[0], row[1], sample_info ]
            }
        ATAC_DESEQ2_QC ( ch_deseq2_in, deseq2_min_count, deseq2_top_features )
        ch_multiqc_files = ch_multiqc_files.mix(ATAC_DESEQ2_QC.out.mqc)
        ch_versions = ch_versions.mix(ATAC_DESEQ2_QC.out.versions)
    }
    ch_versions = ch_versions
        .mix(CONSENSUS_PEAKS.out.versions)
        .mix(BED_SLOP.out.versions)
        .mix(PEAKS_MULTICOV_COUNTS.out.versions)

    //
    // Super-enhancers with ROSE (-g/--custom annotation -s 12500 -t 2500)
    //
    ch_super_enhancers = Channel.empty()
    if ( run_rose ) {
        ch_rose_in = MACS3_CALLPEAK.out.peak.join( ch_bam_bai )
            .combine( ch_rose_annotation )
            .map { row ->
                def annotation = row.size() > 4 ? row[4] : []
                [ row[0], row[1], row[2], row[3], annotation ]
            }
        ROSE ( ch_rose_in, rose_genome, rose_stitch, rose_tss )
        ch_super_enhancers = ROSE.out.super_enhancers
        ch_super_tables = ch_super_enhancers
            .map { meta, table -> table }
            .collect()
            .filter { tables -> tables.size() > 0 }
        CONSENSUS_SUPER_ENHANCERS ( ch_super_tables )
        ch_super_counts_in = CONSENSUS_SUPER_ENHANCERS.out.bed
            .combine( ch_multicov_bams )
            .map { bed, sample_names, bams, bais -> [ [id:'consensus_super_enhancers'], bed, sample_names, bams, bais ] }
        SUPER_ENHANCERS_MULTICOV_COUNTS ( ch_super_counts_in )
        ch_count_tables = ch_count_tables.mix(SUPER_ENHANCERS_MULTICOV_COUNTS.out.counts)
        ch_versions = ch_versions.mix(ROSE.out.versions)
            .mix(CONSENSUS_SUPER_ENHANCERS.out.versions)
            .mix(SUPER_ENHANCERS_MULTICOV_COUNTS.out.versions)
    }

    //
    // Nucleosome positioning with NucleoATAC (occ -> vprocess -> nuc -> merge -> nfr)
    // over the per-sample MACS3 peak regions, using the analysis-ready BAM.
    //
    ch_nucleosomes = Channel.empty()
    if ( run_nucleoatac ) {
        ch_nuc_in = ch_bam_bai.join( MACS3_CALLPEAK.out.peak )
            .map { meta, bam, bai, peak -> [ meta, bam, bai, peak ] }
        NUCLEOATAC (
            ch_nuc_in,
            ch_fasta,
            ch_fai
        )
        ch_nucleosomes = NUCLEOATAC.out.nucpos
        ch_versions = ch_versions.mix(NUCLEOATAC.out.versions)
    }

    //
    // Footprinting with TOBIAS (ATACorrect -> ScoreBigwig -> BINDetect, HOCOMOCO v11)
    //
    ch_bindetect = Channel.empty()
    if ( run_tobias ) {
        ch_atac_in = ch_bam_bai
            .combine( ch_consensus_peaks.map { meta, bed -> bed } )
            .map { meta, bam, bai, bed -> [ meta, bam, bai, bed ] }
        TOBIAS_ATACORRECT ( ch_atac_in, ch_fasta )

        ch_score_in = TOBIAS_ATACORRECT.out.corrected
            .combine( ch_consensus_peaks.map { meta, bed -> bed } )
            .map { meta, bw, bed -> [ meta, bw, bed ] }
        TOBIAS_SCOREBIGWIG ( ch_score_in )

        ch_bindetect_in = TOBIAS_SCOREBIGWIG.out.footprints
            .collect()
            .map { rows ->
                def sorted = rows.sort { a, b -> a[0].id <=> b[0].id }
                def footprints = sorted.collect { it[1] }
                def cond_names = sorted.collect { it[0].condition ?: it[0].id }.join(' ')
                [ [ id: 'all_samples' ], footprints, cond_names ]
            }

        TOBIAS_BINDETECT (
            ch_bindetect_in,
            ch_motifs,
            ch_fasta,
            ch_consensus_peaks
        )
        ch_bindetect = TOBIAS_BINDETECT.out.results
        ch_versions = ch_versions
            .mix(TOBIAS_ATACORRECT.out.versions)
            .mix(TOBIAS_SCOREBIGWIG.out.versions)
            .mix(TOBIAS_BINDETECT.out.versions)
    }

    emit:
    bam            = ch_bam_bai                    // channel: [ meta, bam, bai ]
    bigwig         = DEEPTOOLS_BAMCOVERAGE.out.bigwig
    peaks          = MACS3_CALLPEAK.out.peak
    consensus_peaks= ch_consensus_peaks
    peak_regions   = BED_SLOP.out.bed
    super_enhancers= ch_super_enhancers
    nucleosomes    = ch_nucleosomes
    footprints     = ch_bindetect
    count_tables   = ch_count_tables
    multiqc_files  = ch_multiqc_files
    versions       = ch_versions                   // channel: [ versions.yml ]
}
