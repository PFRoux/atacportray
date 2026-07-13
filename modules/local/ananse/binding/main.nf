process ANANSE_BINDING {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ananse:0.5.1--pyhdfd78af_0' :
        'quay.io/biocontainers/ananse:0.5.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(regions)
    path  genome
    path  pfm
    path  reference
    val   tfs
    val   jaccard_cutoff

    output:
    tuple val(meta), path("*/binding.h5") , emit: binding
    tuple val(meta), path("*ananse_binding") , emit: outdir
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}.ananse_binding"
    def pfm_arg = pfm ? "-p ${pfm}" : ''
    def reference_arg = reference ? "-R ${reference}" : ''
    def tfs_arg = tfs ? "-t ${tfs}" : ''
    """
    ananse binding \\
        --atac-bams ${bam} \\
        --genome ${genome} \\
        --regions ${regions} \\
        --outdir ${prefix} \\
        --jaccard-cutoff ${jaccard_cutoff} \\
        --ncore ${task.cpus} \\
        ${pfm_arg} \\
        ${reference_arg} \\
        ${tfs_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ananse: \$(ananse --version 2>&1 | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}.ananse_binding"
    """
    mkdir -p ${prefix}
    touch ${prefix}/binding.h5

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ananse: 0.5.1
    END_VERSIONS
    """
}
