process ANANSE_NETWORK {
    tag "$meta.id"
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ananse:0.5.1--pyhdfd78af_0' :
        'quay.io/biocontainers/ananse:0.5.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(binding)
    path  expression
    path  genome
    path  annotation
    val   expression_column
    val   full_output
    val   tfs

    output:
    tuple val(meta), path("*.network.tsv") , emit: network
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.ananse"
    def annotation_arg = annotation ? "--annotation ${annotation}" : ''
    def column = expression_column ?: meta.id
    def full_arg = full_output ? "--full-output" : ''
    def tfs_arg = tfs ? "-t ${tfs}" : ''
    """
    ananse network \\
        ${binding} \\
        --expression ${expression} \\
        --genome ${genome} \\
        --outfile ${prefix}.network.tsv \\
        --columns ${column} \\
        --ncore ${task.cpus} \\
        ${annotation_arg} \\
        ${full_arg} \\
        ${tfs_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ananse: \$(ananse --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.ananse"
    """
    touch ${prefix}.network.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ananse: 0.5.1
    END_VERSIONS
    """
}
