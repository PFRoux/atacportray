process CONSENSUS_PEAKS {
    tag "consensus"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--h13024bc_3' :
        'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3' }"

    input:
    path peaks, stageAs: "peaks/*"

    output:
    path "consensus_peaks.bed"    , emit: bed
    path "consensus_peaks.saf"    , emit: saf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''   // e.g. min sample overlap via awk downstream
    """
    # Merge peak intervals across all samples into a single consensus set
    cat peaks/* \\
        | cut -f1-3 \\
        | sort -k1,1 -k2,2n \\
        | bedtools merge -i - $args \\
        > consensus_peaks.bed

    # SAF format (for featureCounts-style quantification downstream)
    awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"} \\
         {print \$1":"\$2"-"\$3, \$1, \$2+1, \$3, "."}' consensus_peaks.bed \\
        > consensus_peaks.saf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e 's/bedtools v//g')
    END_VERSIONS
    """

    stub:
    """
    touch consensus_peaks.bed
    touch consensus_peaks.saf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
    END_VERSIONS
    """
}
