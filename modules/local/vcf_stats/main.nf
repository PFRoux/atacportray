process VCF_STATS {
    tag "${meta.id}:${meta.caller}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.bcftools_stats.txt"), emit: stats
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def caller = meta.caller ?: 'variants'
    def prefix = task.ext.prefix ?: "${meta.id}.${caller}"
    """
    bcftools stats \\
        ${args} \\
        ${vcf} \\
        > ${prefix}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | sed '1!d; s/^.*bcftools //')
    END_VERSIONS
    """

    stub:
    def caller = meta.caller ?: 'variants'
    def prefix = task.ext.prefix ?: "${meta.id}.${caller}"
    """
    touch ${prefix}.bcftools_stats.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.23.1
    END_VERSIONS
    """
}
