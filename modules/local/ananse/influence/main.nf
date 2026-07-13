process ANANSE_INFLUENCE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ananse:0.5.1--pyhdfd78af_0' :
        'quay.io/biocontainers/ananse:0.5.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(source_network), path(target_network), path(degenes)
    path  annotation
    val   edges
    val   padj
    val   full_output

    output:
    tuple val(meta), path("*.influence.tsv")   , emit: influence
    tuple val(meta), path("*diffnetwork*.tsv") , emit: diffnetwork, optional: true
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.ananse"
    def annotation_arg = annotation ? "--annotation ${annotation}" : ''
    def full_arg = full_output ? "--full-output" : ''
    """
    ananse influence \\
        --source ${source_network} \\
        --target ${target_network} \\
        --degenes ${degenes} \\
        --outfile ${prefix}.influence.tsv \\
        --edges ${edges} \\
        --padj ${padj} \\
        --ncore ${task.cpus} \\
        ${annotation_arg} \\
        ${full_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ananse: \$(ananse --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.ananse"
    """
    touch ${prefix}.influence.tsv
    touch ${prefix}_diffnetwork.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ananse: 0.5.1
    END_VERSIONS
    """
}
