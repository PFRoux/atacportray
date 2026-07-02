process QDNASEQ {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-qdnaseq:1.46.0--r45hdfd78af_0' :
        'quay.io/biocontainers/bioconductor-qdnaseq:1.46.0--r45hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val   binsize_kb
    val   loss_threshold
    val   gain_threshold

    output:
    tuple val(meta), path("*.calls.txt")     , emit: calls
    tuple val(meta), path("*.segments.txt")  , emit: segments
    tuple val(meta), path("*.copynumbers.txt"), emit: copynumbers
    tuple val(meta), path("*.igv")           , emit: igv       , optional: true
    tuple val(meta), path("*.png")           , emit: plots     , optional: true
    tuple val(meta), path("*.rds")           , emit: rds
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    qdnaseq.R \\
        --bam $bam \\
        --sample $prefix \\
        --binsize $binsize_kb \\
        --loss $loss_threshold \\
        --gain $gain_threshold \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(strsplit(R.version.string, " ")[[1]][3])')
        qdnaseq: \$(Rscript -e 'cat(as.character(packageVersion("QDNAseq")))')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.calls.txt
    touch ${prefix}.segments.txt
    touch ${prefix}.copynumbers.txt
    touch ${prefix}.igv
    touch ${prefix}.png
    touch ${prefix}.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-qdnaseq: 1.46.0
    END_VERSIONS
    """
}
