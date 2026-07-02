process TOBIAS_ATACORRECT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tobias:0.16.1--py39h24fbfe6_1' :
        'quay.io/biocontainers/tobias:0.16.1--py39h24fbfe6_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(peaks)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*_corrected.bw")  , emit: corrected
    tuple val(meta), path("*_expected.bw")   , emit: expected
    tuple val(meta), path("*_uncorrected.bw"), emit: uncorrected, optional: true
    tuple val(meta), path("*_bias.bw")       , emit: bias       , optional: true
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    TOBIAS ATACorrect \\
        --bam $bam \\
        --genome $fasta \\
        --peaks $peaks \\
        --prefix $prefix \\
        --outdir . \\
        --cores $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tobias: \$(TOBIAS --version 2>&1 | sed 's/TOBIAS //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_corrected.bw
    touch ${prefix}_expected.bw
    touch ${prefix}_uncorrected.bw
    touch ${prefix}_bias.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tobias: 0.16.1
    END_VERSIONS
    """
}
