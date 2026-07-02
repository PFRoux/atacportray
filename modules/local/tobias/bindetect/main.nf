process TOBIAS_BINDETECT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tobias:0.16.1--py39h24fbfe6_1' :
        'quay.io/biocontainers/tobias:0.16.1--py39h24fbfe6_1' }"

    input:
    tuple val(meta), path(footprints)
    path motifs
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(peaks)

    output:
    tuple val(meta), path("bindetect/*bindetect_results.txt"), emit: results
    tuple val(meta), path("bindetect/*bindetect_results.xlsx"), emit: results_xlsx, optional: true
    tuple val(meta), path("bindetect/*.pdf")                  , emit: figures    , optional: true
    tuple val(meta), path("bindetect/*_overview.txt")         , emit: overview   , optional: true
    tuple val(meta), path("bindetect/")                       , emit: outdir
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cond   = meta.condition ? "--cond_names ${meta.condition}" : ''
    """
    TOBIAS BINDetect \\
        --motifs $motifs \\
        --signals $footprints \\
        --genome $fasta \\
        --peaks $peaks \\
        --outdir bindetect \\
        --cores $task.cpus \\
        $cond \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tobias: \$(TOBIAS --version 2>&1 | sed 's/TOBIAS //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p bindetect
    touch bindetect/${prefix}_bindetect_results.txt
    touch bindetect/${prefix}_bindetect_results.xlsx
    touch bindetect/${prefix}_bindetect_overview.txt
    touch bindetect/${prefix}_figures.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tobias: 0.16.1
    END_VERSIONS
    """
}
