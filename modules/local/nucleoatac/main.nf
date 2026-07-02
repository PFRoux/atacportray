process NUCLEOATAC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/pfroux/atacportray-nucleoatac:1.0.0' :
        'ghcr.io/pfroux/atacportray-nucleoatac:1.0.0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(peaks)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.nucpos.bed.gz")            , emit: nucpos       , optional: true
    tuple val(meta), path("*.nucmap_combined.bed.gz")   , emit: nucmap       , optional: true
    tuple val(meta), path("*.nfrpos.bed.gz")            , emit: nfr          , optional: true
    tuple val(meta), path("*.occpeaks.bed.gz")          , emit: occpeaks     , optional: true
    tuple val(meta), path("*.occ.bedgraph.gz")          , emit: occ          , optional: true
    tuple val(meta), path("*.nucleoatac_signal.bedgraph.gz"), emit: signal   , optional: true
    tuple val(meta), path("*.nucleoatac_signal.smooth.bedgraph.gz"), emit: signal_smooth, optional: true
    tuple val(meta), path("*")                          , emit: results
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # NucleoATAC needs an indexed FASTA (pysam FastaFile). The .fai is staged
    # next to the FASTA; refresh it only if missing so a provided index is used.
    if [ ! -s ${fasta}.fai ]; then
        samtools faidx $fasta
    fi

    # `nucleoatac run` chains occ -> vprocess -> nuc -> merge -> nfr, the full
    # nucleosome-calling pipeline over the supplied peak regions.
    nucleoatac run \\
        --bed $peaks \\
        --bam $bam \\
        --fasta $fasta \\
        --out $prefix \\
        --cores $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nucleoatac: \$(nucleoatac --version 2>&1 | sed 's/.*[ ]//' || echo 1.0.0)
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.occ.bedgraph.gz
    echo | gzip > ${prefix}.occpeaks.bed.gz
    echo | gzip > ${prefix}.nucpos.bed.gz
    echo | gzip > ${prefix}.nucmap_combined.bed.gz
    echo | gzip > ${prefix}.nfrpos.bed.gz
    echo | gzip > ${prefix}.nucleoatac_signal.bedgraph.gz
    echo | gzip > ${prefix}.nucleoatac_signal.smooth.bedgraph.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nucleoatac: 1.0.0
        samtools: 1.21
    END_VERSIONS
    """
}
