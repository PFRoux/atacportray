process BED_SLOP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta2), path(fai)
    val slop_bp

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.slop${slop_bp}"
    """
    awk 'BEGIN{OFS="\\t"} {print \$1,\$2}' ${fai} > genome.sizes

    bedtools slop \\
        -i ${bed} \\
        -g genome.sizes \\
        -b ${slop_bp} \\
        | sort -k1,1 -k2,2n \\
        | bedtools merge -i - \\
        > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e 's/bedtools v//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.slop${slop_bp}"
    """
    touch ${prefix}.bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
    END_VERSIONS
    """
}
