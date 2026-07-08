process VCF_FILTER_REGIONS {
    tag "${meta.id}:${meta.caller}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf), path(regions)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def caller = meta.caller ?: 'variants'
    def prefix = task.ext.prefix ?: "${meta.id}.${caller}.peaks"
    """
    bcftools view \\
        --targets-file ${regions} \\
        --output-type z \\
        --output-file ${prefix}.vcf.gz \\
        ${args} \\
        ${vcf}

    tabix -p vcf ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | sed '1!d; s/^.*bcftools //')
        tabix: \$(tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+')
    END_VERSIONS
    """

    stub:
    def caller = meta.caller ?: 'variants'
    def prefix = task.ext.prefix ?: "${meta.id}.${caller}.peaks"
    """
    echo "stub" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.23.1
        tabix: 1.23.1
    END_VERSIONS
    """
}
