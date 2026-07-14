process MGATK_POSTPROCESS {
    tag "mgatk_postprocess"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-summarizedexperiment:1.36.0--r44hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-summarizedexperiment:1.36.0--r44hdfd78af_0' }"

    input:
    path rds
    path coverage
    path variant_stats
    path vmr
    val min_af
    val min_vmr
    val min_sd

    output:
    path "mgatk_variant_summary.tsv"         , emit: summary
    path "mgatk_high_confidence_variants.tsv", emit: high_confidence
    path "mgatk_heteroplasmy_matrix.tsv"     , emit: matrix
    path "mgatk_postprocess_mqc.tsv"         , emit: mqc
    path "mgatk_mt_coverage.tsv"             , emit: coverage_table
    path "mgatk_mt_coverage_mqc.tsv"         , emit: coverage_mqc
    path "mgatk_mt_coverage_circular_mqc.html", emit: coverage_circular_mqc
    path "mgatk_mt_coverage_circular.pdf"    , emit: coverage_circular, optional: true
    path "mgatk_heteroplasmy_heatmap.pdf"    , emit: heatmap, optional: true
    path "mgatk_heteroplasmy_pca.pdf"        , emit: pca, optional: true
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${moduleDir}/postprocess_mgatk.R \\
        --rds "${rds}" \\
        --coverage "${coverage}" \\
        --variant-stats "${variant_stats}" \\
        --vmr "${vmr}" \\
        --min-af ${min_af} \\
        --min-vmr ${min_vmr} \\
        --min-sd ${min_sd}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript --version 2>&1 | sed 's/Rscript (R) version //')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_SUMMARY > mgatk_variant_summary.tsv
    metric	value
    total_variants	0
    high_confidence_variants	0
    variable_variants	0
    END_SUMMARY

    cat <<-END_VARIANTS > mgatk_high_confidence_variants.tsv
    variant	mean_af	vmr	detected_samples	variable
    END_VARIANTS

    cat <<-END_MATRIX > mgatk_heteroplasmy_matrix.tsv
    variant
    END_MATRIX

    cat <<-END_COVERAGE > mgatk_mt_coverage.tsv
    sample	position	coverage
    END_COVERAGE

    cat <<-END_MQC > mgatk_postprocess_mqc.tsv
    # id: atacportray_mgatk_postprocess
    # section_name: mgatk heteroplasmy summary
    # description: High-confidence mitochondrial heteroplasmy calls after allele-fraction, VMR and variability filtering.
    # plot_type: table
    metric	value
    total_variants	0
    high_confidence_variants	0
    variable_variants	0
    END_MQC

    cat <<-END_COV_MQC > mgatk_mt_coverage_mqc.tsv
    # id: atacportray_mgatk_mt_coverage
    # section_name: mgatk mitochondrial coverage
    # description: Per-position mitochondrial coverage reported by mgatk.
    # plot_type: linegraph
    # pconfig:
    #   id: atacportray_mgatk_mt_coverage
    #   title: mgatk mitochondrial coverage
    #   xlab: MT position
    #   ylab: Coverage
    END_COV_MQC

    cat <<-END_CIRCULAR_MQC > mgatk_mt_coverage_circular_mqc.html
    # id: atacportray_mgatk_mt_coverage_circular
    # section_name: mgatk circular mitochondrial coverage
    # description: Circular representation of per-position mitochondrial coverage reported by mgatk.
    # plot_type: html
    <p>No mitochondrial coverage data available.</p>
    END_CIRCULAR_MQC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.4.3
    END_VERSIONS
    """
}
