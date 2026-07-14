process ATAC_FRIP_MQC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(bam), path(bai), path(peaks)

    output:
    tuple val(meta), path("*_frip_mqc.tsv"), emit: mqc
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    total_mapped=\$(samtools view -@ ${task.cpus - 1} -c -F 4 ${bam})
    reads_in_peaks=\$(samtools view -@ ${task.cpus - 1} -c -F 4 -L ${peaks} ${bam})
    frip=\$(awk -v in_peaks="\${reads_in_peaks}" -v total="\${total_mapped}" 'BEGIN { if (total > 0) printf "%.6f", in_peaks / total; else printf "0.000000" }')

    cat <<-END_MQC > ${prefix}_frip_mqc.tsv
    # id: atacportray_frip
    # section_name: ATAC FRiP
    # description: Fraction of mapped reads overlapping sample MACS3 peak regions.
    # plot_type: table
    Sample	Total mapped reads	Reads in peaks	FRiP
    ${meta.id}	\${total_mapped}	\${reads_in_peaks}	\${frip}
    END_MQC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | sed '1!d;s/.* //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-END_MQC > ${prefix}_frip_mqc.tsv
    # id: atacportray_frip
    # section_name: ATAC FRiP
    # description: Fraction of mapped reads overlapping sample MACS3 peak regions.
    # plot_type: table
    Sample	Total mapped reads	Reads in peaks	FRiP
    ${meta.id}	1000	500	0.500000
    END_MQC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """
}
