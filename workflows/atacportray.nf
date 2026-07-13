/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_atacportray_pipeline'
include { getGenomeAttribute     } from '../subworkflows/local/utils_nfcore_atacportray_pipeline'

//
// MODULE: reference preparation
//
include { FASTP                         } from '../modules/nf-core/fastp/main'
include { BOWTIE2_BUILD                 } from '../modules/nf-core/bowtie2/build/main'
include { BWA_INDEX                     } from '../modules/nf-core/bwa/index/main'
include { SAMTOOLS_FAIDX                } from '../modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY} from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { GUNZIP                        } from '../modules/nf-core/gunzip/main'
include { BEDTOOLS_INTERSECT as FILTER_BAM_OUTSIDE_PEAKS } from '../modules/nf-core/bedtools/intersect/main'
include { SAMTOOLS_INDEX as INDEX_CNV_BAM                } from '../modules/nf-core/samtools/index/main'

//
// SUBWORKFLOW: analysis branches
//
include { FASTQ_EPIGENOME_BOWTIE2       } from '../subworkflows/local/fastq_epigenome_bowtie2/main'
include { FASTQ_VARIANT_CALLING_ATAC    } from '../subworkflows/local/fastq_variant_calling_atac/main'
include { BAM_CNV_QDNASEQ               } from '../subworkflows/local/bam_cnv_qdnaseq/main'
include { BAM_TELOMERE_TELOMEREHUNTER   } from '../subworkflows/local/bam_telomere_telomerehunter/main'
include { BAM_MITO_MGATK                } from '../subworkflows/local/bam_mito_mgatk/main'
include { ANANSE_BINDING                } from '../modules/local/ananse/binding/main'
include { ANANSE_NETWORK                } from '../modules/local/ananse/network/main'
include { ANANSE_INFLUENCE              } from '../modules/local/ananse/influence/main'
include { COLTRON                       } from '../modules/local/coltron/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ATACPORTRAY {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // Reference channels
    //
    ref_fasta           = params.fasta          ?: getGenomeAttribute('fasta')
    ref_bowtie2_index   = params.bowtie2_index  ?: getGenomeAttribute('bowtie2')
    ref_bwa_index       = params.bwa_index      ?: getGenomeAttribute('bwa')
    ref_blacklist       = params.blacklist      ?: getGenomeAttribute('blacklist')
    ref_known_sites     = params.known_sites    ?: getGenomeAttribute('known_sites')
    ref_known_sites_tbi = params.known_sites_tbi ?: getGenomeAttribute('known_sites_tbi')
    ref_vep_cache       = params.vep_cache      ?: getGenomeAttribute('vep_cache')

    // Decompress FASTA if gzipped. Static reference files are value channels so
    // that every sample can reuse them in multi-input processes.
    if ( ref_fasta.endsWith('.gz') ) {
        ch_fasta_raw = channel.fromPath(ref_fasta, checkIfExists: true)
            .map { fa -> [ [id:'genome'], fa ] }
        GUNZIP ( ch_fasta_raw )
        ch_fasta = GUNZIP.out.gunzip.first()
    } else {
        ch_fasta = channel.value([ [id:'genome'], file(ref_fasta) ])
    }

    // FASTA index (.fai) + chrom sizes
    if ( params.fasta_fai ) {
        ch_fai = channel.value([ [id:'genome'], file(params.fasta_fai) ])
    } else {
        SAMTOOLS_FAIDX ( ch_fasta.map { m, f -> [ m, f, [] ] }, false )
        ch_fai = SAMTOOLS_FAIDX.out.fai
    }

    // Bowtie2 index (epigenome)
    if ( params.run_epigenome ) {
        if ( ref_bowtie2_index ) {
            ch_bowtie2_index = channel.value([ [id:'genome'], file(ref_bowtie2_index) ])
        } else {
            BOWTIE2_BUILD ( ch_fasta )
            ch_bowtie2_index = BOWTIE2_BUILD.out.index
        }
    } else {
        ch_bowtie2_index = channel.value([[:], []])
    }

    // BWA index + GATK dict (variants)
    ch_dict = channel.value([[:], []])
    if ( params.run_variants ) {
        if ( ref_bwa_index ) {
            ch_bwa_index = channel.value([ [id:'genome'], file(ref_bwa_index) ])
        } else {
            BWA_INDEX ( ch_fasta )
            ch_bwa_index = BWA_INDEX.out.index
        }
        if ( params.dict ) {
            ch_dict = channel.value([ [id:'genome'], file(params.dict) ])
        } else {
            GATK4_CREATESEQUENCEDICTIONARY ( ch_fasta )
            ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
        }
    } else {
        ch_bwa_index = channel.value([[:], []])
    }

    // Blacklist / motifs / known-sites / cytoBand / VEP cache
    ch_blacklist = ref_blacklist ?
        channel.value([ [id:'blacklist'], file(ref_blacklist) ]) :
        channel.value([[:], []])
    ch_motifs = params.tobias_motifs ?
        channel.fromPath(params.tobias_motifs, checkIfExists: true).collect() :
        channel.value([])
    ch_known_sites = ref_known_sites ?
        channel.value([ [id:'known'], file(ref_known_sites) ]) :
        channel.value([[:], []])
    ch_known_sites_tbi = ref_known_sites_tbi ?
        channel.value([ [id:'known'], file(ref_known_sites_tbi) ]) :
        channel.value([[:], []])
    ch_cytoband = params.cytoband ?
        channel.fromPath(params.cytoband, checkIfExists: true).collect() :
        channel.value([])
    ch_qdnaseq_bins_rds = params.qdnaseq_bins_rds ?
        channel.fromPath(params.qdnaseq_bins_rds, checkIfExists: true).collect() :
        channel.value([])
    ch_vep_cache = ref_vep_cache ?
        channel.value([ [id:'vep'], file(ref_vep_cache) ]) :
        channel.value([[:], []])

    //
    // Read trimming (fastp)
    //
    if ( params.skip_trimming ) {
        ch_trimmed = ch_samplesheet
    } else {
        ch_fastp_in = ch_samplesheet.map { meta, reads -> [ meta, reads, [] ] }
        FASTP ( ch_fastp_in, false, false, false )
        ch_trimmed = FASTP.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] }.ifEmpty([]))
    }

    //
    // Shared analysis-ready BAM (produced by the variants branch preprocessing,
    // consumed by CNV / telomere / mito). Populated only if run_variants = true.
    //
    ch_analysis_bam = channel.empty()
    ch_peak_regions = channel.empty()
    ch_epigenome_bam = channel.empty()
    ch_super_enhancers = channel.empty()

    //
    // EPIGENOME branch
    //
    if ( params.run_epigenome ) {
        FASTQ_EPIGENOME_BOWTIE2 (
            ch_trimmed,
            ch_bowtie2_index,
            ch_fasta,
            ch_fai,
            ch_blacklist,
            ref_blacklist as boolean,
            params.macs_gsize,
            ch_motifs,
            params.run_epigenome,        // run_rose gate (ROSE always on within epigenome)
            params.rose_stitch,
            params.rose_tss,
            params.run_footprinting,
            params.run_nucleoatac,
            params.peak_filter_slop
        )
        ch_epigenome_bam = FASTQ_EPIGENOME_BOWTIE2.out.bam
        ch_peak_regions  = FASTQ_EPIGENOME_BOWTIE2.out.peak_regions
        ch_super_enhancers = FASTQ_EPIGENOME_BOWTIE2.out.super_enhancers
        ch_versions      = ch_versions.mix(FASTQ_EPIGENOME_BOWTIE2.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_EPIGENOME_BOWTIE2.out.multiqc_files)
    }

    //
    // Regulatory network branch (ANANSE / Coltron)
    //
    if ( params.run_ananse && params.run_epigenome ) {
        ch_ananse_pfm = params.ananse_pfm ?
            channel.value(file(params.ananse_pfm)) :
            channel.value([])
        ch_ananse_reference = params.ananse_reference ?
            channel.value(file(params.ananse_reference)) :
            channel.value([])
        ch_ananse_binding_in = ch_epigenome_bam
            .combine( ch_peak_regions.map { meta, bed -> bed } )
            .map { meta, bam, bai, regions -> [ meta, bam, bai, regions ] }

        ANANSE_BINDING (
            ch_ananse_binding_in,
            ch_fasta.map { meta, fasta -> fasta },
            ch_ananse_pfm,
            ch_ananse_reference,
            params.ananse_tfs ?: '',
            params.ananse_jaccard_cutoff
        )
        ch_versions = ch_versions.mix(ANANSE_BINDING.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ANANSE_BINDING.out.binding.map { meta, binding -> binding })

        if ( params.run_ananse_network ) {
            ch_ananse_expression = channel.value(file(params.ananse_expression))
            ch_ananse_annotation = params.ananse_annotation ?
                channel.value(file(params.ananse_annotation)) :
                channel.value([])

            ANANSE_NETWORK (
                ANANSE_BINDING.out.binding,
                ch_ananse_expression,
                ch_fasta.map { meta, fasta -> fasta },
                ch_ananse_annotation,
                params.ananse_expression_column ?: '',
                params.ananse_network_full_output,
                params.ananse_tfs ?: ''
            )
            ch_versions = ch_versions.mix(ANANSE_NETWORK.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(ANANSE_NETWORK.out.network.map { meta, network -> network })
        }
    }

    if ( params.run_ananse_influence ) {
        ch_ananse_influence_in = channel.value([
            [ id: 'ananse_influence' ],
            file(params.ananse_source_network),
            file(params.ananse_target_network),
            file(params.ananse_degenes)
        ])
        ch_ananse_influence_annotation = params.ananse_influence_annotation ?
            channel.value(file(params.ananse_influence_annotation)) :
            channel.value([])

        ANANSE_INFLUENCE (
            ch_ananse_influence_in,
            ch_ananse_influence_annotation,
            params.ananse_influence_edges,
            params.ananse_influence_padj,
            params.ananse_influence_full_output
        )
        ch_versions = ch_versions.mix(ANANSE_INFLUENCE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ANANSE_INFLUENCE.out.influence.map { meta, influence -> influence },
            ANANSE_INFLUENCE.out.diffnetwork.map { meta, diffnetwork -> diffnetwork }
        )
    }

    if ( params.run_coltron && params.run_epigenome ) {
        ch_coltron_subpeaks = params.coltron_subpeaks ?
            channel.value(file(params.coltron_subpeaks)) :
            channel.value([])
        ch_coltron_activity = params.coltron_activity ?
            channel.value(file(params.coltron_activity)) :
            channel.value([])
        ch_coltron_motifs = params.coltron_motifs ?
            channel.value(file(params.coltron_motifs)) :
            channel.value([])
        ch_coltron_in = ch_super_enhancers
            .join( ch_epigenome_bam )
            .map { meta, enhancers, bam, bai -> [ meta, enhancers, bam, bai ] }

        COLTRON (
            ch_coltron_in,
            params.coltron_genome,
            ch_coltron_subpeaks,
            ch_coltron_activity,
            ch_coltron_motifs,
            params.coltron_enhancer_number ?: '',
            params.coltron_exp_cutoff ?: ''
        )
        ch_versions = ch_versions.mix(COLTRON.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(COLTRON.out.outdir.map { meta, outdir -> outdir })
    }

    //
    // VARIANTS branch (also yields the shared analysis-ready BAM)
    //
    if ( params.run_variants ) {
        def callers = params.variant_callers.tokenize(',')*.trim()
        FASTQ_VARIANT_CALLING_ATAC (
            ch_trimmed,
            ch_bwa_index,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_known_sites,
            ch_known_sites_tbi,
            ch_vep_cache,
            params.skip_bqsr,
            params.skip_annotation,
            params.vep_genome,
            params.vep_species,
            params.vep_cache_version,
            callers,
            params.run_oncoplot,
            params.variants_in_peaks_only,
            ch_peak_regions
        )
        ch_analysis_bam = FASTQ_VARIANT_CALLING_ATAC.out.bam
        ch_versions     = ch_versions.mix(FASTQ_VARIANT_CALLING_ATAC.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            FASTQ_VARIANT_CALLING_ATAC.out.vcf.map { meta, vcf -> vcf },
            FASTQ_VARIANT_CALLING_ATAC.out.vcf_stats.map { meta, stats -> stats },
            FASTQ_VARIANT_CALLING_ATAC.out.maf.map { meta, maf -> maf },
            FASTQ_VARIANT_CALLING_ATAC.out.oncoplot
        )
    }

    //
    // CNV / TELOMERE / MITO branches (consume the shared analysis-ready BAM).
    // Require run_variants = true so the analysis-ready BAM exists.
    //
    if ( params.run_cnv && params.run_variants ) {
        ch_cnv_bam = ch_analysis_bam
        if ( params.qdnaseq_exclude_peaks ) {
            ch_filter_bam_in = ch_analysis_bam
                .combine( ch_peak_regions )
                .map { meta, bam, bai, peak_meta, peaks -> [ meta, bam, peaks ] }
            FILTER_BAM_OUTSIDE_PEAKS ( ch_filter_bam_in, [[:], []] )
            INDEX_CNV_BAM ( FILTER_BAM_OUTSIDE_PEAKS.out.intersect )
            ch_cnv_bam = FILTER_BAM_OUTSIDE_PEAKS.out.intersect.join( INDEX_CNV_BAM.out.index )
            ch_multiqc_files = ch_multiqc_files.mix(ch_cnv_bam.map { meta, bam, bai -> bam })
        }
        BAM_CNV_QDNASEQ (
            ch_cnv_bam,
            params.qdnaseq_binsize,
            ch_qdnaseq_bins_rds,
            params.qdnaseq_loss_threshold,
            params.qdnaseq_gain_threshold
        )
        ch_versions = ch_versions.mix(BAM_CNV_QDNASEQ.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            BAM_CNV_QDNASEQ.out.calls.map { meta, calls -> calls },
            BAM_CNV_QDNASEQ.out.segments.map { meta, segments -> segments },
            BAM_CNV_QDNASEQ.out.plots.map { meta, plots -> plots }
        )
    }

    if ( params.run_telomere && params.run_variants ) {
        BAM_TELOMERE_TELOMEREHUNTER ( ch_analysis_bam, ch_cytoband )
        ch_versions = ch_versions.mix(BAM_TELOMERE_TELOMEREHUNTER.out.versions)
    }

    if ( params.run_mito && params.run_variants ) {
        BAM_MITO_MGATK ( ch_analysis_bam, params.mito_contig, params.mgatk_keep_duplicates )
        ch_versions = ch_versions.mix(BAM_MITO_MGATK.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            BAM_MITO_MGATK.out.variant_stats.map { meta, stats -> stats },
            BAM_MITO_MGATK.out.final_dir.map { meta, final_dir -> final_dir }
        )
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'atacportray_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
