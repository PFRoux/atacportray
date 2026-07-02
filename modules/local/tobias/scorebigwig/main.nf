process TOBIAS_SCOREBIGWIG {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tobias:0.16.1--py39h24fbfe6_1' :
        'quay.io/biocontainers/tobias:0.16.1--py39h24fbfe6_1' }"

    input:
    tuple val(meta), path(corrected), path(peaks)

    output:
    tuple val(meta), path("*_footprints.bw"), emit: footprints
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    TOBIAS ScoreBigwig \\
        --signal $corrected \\
        --regions $peaks \\
        --output ${prefix}_footprints.bw \\
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
    touch ${prefix}_footprints.bw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tobias: 0.16.1
    END_VERSIONS
    """
}
