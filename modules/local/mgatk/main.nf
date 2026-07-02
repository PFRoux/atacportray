process MGATK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgatk:0.7.0--pyhdfd78af_2' :
        'quay.io/biocontainers/mgatk:0.7.0--pyhdfd78af_2' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val   mito_contig
    val   keep_duplicates

    output:
    tuple val(meta), path("*.mgatk/final/*.variant_stats.tsv.gz"), emit: variant_stats, optional: true
    tuple val(meta), path("*.mgatk/final/*.vmr_strand.tsv.gz")   , emit: vmr          , optional: true
    tuple val(meta), path("*.mgatk/final")                       , emit: final_dir
    tuple val(meta), path("*.mgatk")                             , emit: outdir
    path "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def keepdup = keep_duplicates ? "--keep-duplicates" : ""
    """
    mgatk call \\
        --input $bam \\
        --output ${prefix}.mgatk \\
        --name $prefix \\
        --mito-genome $mito_contig \\
        --ncores $task.cpus \\
        $keepdup \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgatk: \$(mgatk --version 2>&1 | sed 's/.*version //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.mgatk/final
    echo "" | gzip > ${prefix}.mgatk/final/${prefix}.variant_stats.tsv.gz
    echo "" | gzip > ${prefix}.mgatk/final/${prefix}.vmr_strand.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgatk: 0.7.0
    END_VERSIONS
    """
}
