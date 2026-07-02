process TELOMEREHUNTER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/pfroux/atacportray-telomerehunter:1.1.0' :
        'ghcr.io/pfroux/atacportray-telomerehunter:1.1.0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  cytoband

    output:
    tuple val(meta), path("${meta.id}/*_summary.tsv")                       , emit: summary
    tuple val(meta), path("${meta.id}/**/*_TVR_context_summary.tsv")        , emit: tvr_context, optional: true
    tuple val(meta), path("${meta.id}/**/*_singletons.tsv")                 , emit: singletons , optional: true
    tuple val(meta), path("${meta.id}/**/*.png")                            , emit: plots      , optional: true
    tuple val(meta), path("${meta.id}/")                                     , emit: outdir
    path "versions.yml"                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    telomerehunter \\
        -ibt $bam \\
        -o . \\
        -p $prefix \\
        -b $cytoband \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: \$(telomerehunter --version 2>&1 | sed 's/.*version //; s/ .*//' || echo 1.1.0)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/plots
    touch ${prefix}/${prefix}_summary.tsv
    touch ${prefix}/plots/${prefix}_TVR_context_summary.tsv
    touch ${prefix}/plots/${prefix}_singletons.tsv
    touch ${prefix}/plots/${prefix}_telomere_content.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telomerehunter: 1.1.0
    END_VERSIONS
    """
}
