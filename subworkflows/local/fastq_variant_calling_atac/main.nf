//
// VARIANTS branch: BWA-MEM -> MarkDuplicates -> BQSR -> ApplyBQSR
//                  -> multi-caller (DeepVariant, FreeBayes, GATK HaplotypeCaller,
//                     bcftools/mpileup) -> VEP -> vcf2maf -> maftools oncoplot
// Produces a shared analysis-ready BAM consumed by CNV / telomere / mito branches.
//
include { BWA_MEM                    } from '../../../modules/nf-core/bwa/mem/main'
include { GATK4_MARKDUPLICATES       } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR     } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR            } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { GATK4_HAPLOTYPECALLER      } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { DEEPVARIANT_RUNDEEPVARIANT } from '../../../modules/nf-core/deepvariant/rundeepvariant/main'
include { FREEBAYES                  } from '../../../modules/nf-core/freebayes/main'
include { BCFTOOLS_MPILEUP           } from '../../../modules/nf-core/bcftools/mpileup/main'
include { SAMTOOLS_INDEX as INDEX_ANALYSIS } from '../../../modules/nf-core/samtools/index/main'
include { ENSEMBLVEP_VEP             } from '../../../modules/nf-core/ensemblvep/vep/main'
include { VCF2MAF                    } from '../../../modules/nf-core/vcf2maf/main'
include { MAFTOOLS_ONCOPLOT          } from '../../../modules/local/maftools_oncoplot/main'

workflow FASTQ_VARIANT_CALLING_ATAC {

    take:
    ch_reads         // channel: [ val(meta), [ reads ] ]  (trimmed)
    ch_bwa_index     // channel: [ val(meta2), path(index) ]
    ch_fasta         // channel: [ val(meta3), path(fasta) ]
    ch_fai           // channel: [ val(meta3), path(fai) ]
    ch_dict          // channel: [ val(meta4), path(dict) ]
    ch_known_sites   // channel: [ val(meta5), path(vcf) ]
    ch_known_sites_tbi // channel: [ val(meta6), path(tbi) ]
    ch_vep_cache     // channel: [ val(meta7), path(cache) ]
    vep_genome       // value:   'GRCh38'
    vep_species      // value:   'homo_sapiens'
    vep_cache_version// value:   113
    callers          // value:   list, e.g. ['deepvariant','freebayes','haplotypecaller','bcftools']
    run_oncoplot     // value:   boolean

    main:
    ch_versions   = Channel.empty()
    ch_vcfs       = Channel.empty()

    //
    // Align with BWA-MEM (sort_bam = true)
    //
    BWA_MEM ( ch_reads, ch_bwa_index, ch_fasta, true )

    //
    // GATK4 MarkDuplicates
    //
    GATK4_MARKDUPLICATES (
        BWA_MEM.out.bam,
        ch_fasta.map { m, f -> f },
        ch_fai.map   { m, f -> f }
    )

    //
    // Base Quality Score Recalibration
    //
    ch_bqsr_in = GATK4_MARKDUPLICATES.out.bam
        .join( GATK4_MARKDUPLICATES.out.bai )
        .map { meta, bam, bai -> [ meta, bam, bai, [] ] }
    GATK4_BASERECALIBRATOR (
        ch_bqsr_in, ch_fasta, ch_fai, ch_dict, ch_known_sites, ch_known_sites_tbi
    )

    ch_applybqsr_in = GATK4_MARKDUPLICATES.out.bam
        .join( GATK4_MARKDUPLICATES.out.bai )
        .join( GATK4_BASERECALIBRATOR.out.table )
        .map { meta, bam, bai, table -> [ meta, bam, bai, table, [] ] }
    GATK4_APPLYBQSR (
        ch_applybqsr_in,
        ch_fasta.map { m, f -> f },
        ch_fai.map   { m, f -> f },
        ch_dict.map  { m, f -> f }
    )

    // Analysis-ready BAM (shared with CNV / telomere / mito branches)
    INDEX_ANALYSIS ( GATK4_APPLYBQSR.out.bam )
    ch_analysis_bam = GATK4_APPLYBQSR.out.bam.join( INDEX_ANALYSIS.out.index )

    //
    // Multi-caller variant calling
    //
    ch_bam_bai_iv = ch_analysis_bam.map { meta, bam, bai -> [ meta, bam, bai, [] ] }

    if ( 'haplotypecaller' in callers ) {
        ch_hc_in = ch_analysis_bam.map { meta, bam, bai -> [ meta, bam, bai, [], [] ] }
        GATK4_HAPLOTYPECALLER (
            ch_hc_in, ch_fasta, ch_fai, ch_dict, [[:], []], [[:], []]
        )
        ch_vcfs = ch_vcfs.mix(
            GATK4_HAPLOTYPECALLER.out.vcf.map { meta, vcf -> [ meta + [caller:'haplotypecaller'], vcf ] }
        )
    }

    if ( 'deepvariant' in callers ) {
        ch_dv_in = ch_analysis_bam.map { meta, bam, bai -> [ meta, bam, bai, [] ] }
        DEEPVARIANT_RUNDEEPVARIANT (
            ch_dv_in, ch_fasta, ch_fai, [[:], []], [[:], []]
        )
        ch_vcfs = ch_vcfs.mix(
            DEEPVARIANT_RUNDEEPVARIANT.out.vcf.map { meta, vcf -> [ meta + [caller:'deepvariant'], vcf ] }
        )
    }

    if ( 'freebayes' in callers ) {
        ch_fb_in = ch_analysis_bam.map { meta, bam, bai -> [ meta, bam, bai, [], [], [] ] }
        FREEBAYES (
            ch_fb_in, ch_fasta, ch_fai, [[:], []], [[:], []], [[:], []]
        )
        ch_vcfs = ch_vcfs.mix(
            FREEBAYES.out.vcf.map { meta, vcf -> [ meta + [caller:'freebayes'], vcf ] }
        )
    }

    if ( 'bcftools' in callers ) {
        ch_mp_in = ch_analysis_bam.map { meta, bam, bai -> [ meta, bam, [], [] ] }
        BCFTOOLS_MPILEUP (
            ch_mp_in, ch_fasta.map { m, f -> [ m, f, [] ] }, false
        )
        ch_vcfs = ch_vcfs.mix(
            BCFTOOLS_MPILEUP.out.vcf.map { meta, vcf -> [ meta + [caller:'bcftools'], vcf ] }
        )
    }

    //
    // Annotate with VEP
    //
    ch_vep_in = ch_vcfs.map { meta, vcf -> [ meta, vcf, [] ] }
    ENSEMBLVEP_VEP (
        ch_vep_in, vep_genome, vep_species, vep_cache_version, ch_vep_cache, ch_fasta, []
    )

    //
    // Convert to MAF (vcf2maf)
    //
    VCF2MAF (
        ENSEMBLVEP_VEP.out.vcf,
        ch_fasta.map { m, f -> f },
        ch_vep_cache.map { m, c -> c }
    )
    ch_versions = ch_versions.mix(VCF2MAF.out.versions)

    //
    // Cohort oncoplot (maftools)
    //
    ch_oncoplot = Channel.empty()
    if ( run_oncoplot ) {
        MAFTOOLS_ONCOPLOT ( VCF2MAF.out.maf.collect { it[1] } )
        ch_oncoplot = MAFTOOLS_ONCOPLOT.out.oncoplot
        ch_versions = ch_versions.mix(MAFTOOLS_ONCOPLOT.out.versions)
    }

    emit:
    bam            = ch_analysis_bam    // channel: [ meta, bam, bai ] (shared analysis-ready BAM)
    vcf            = ch_vcfs            // channel: [ meta(+caller), vcf ]
    maf            = VCF2MAF.out.maf    // channel: [ meta, maf ]
    oncoplot       = ch_oncoplot        // channel: [ *.oncoplot.pdf ]
    versions       = ch_versions
}
