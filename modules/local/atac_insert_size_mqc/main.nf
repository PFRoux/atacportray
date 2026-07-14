process ATAC_INSERT_SIZE_MQC {
    tag "atac_insert_sizes"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    path bams
    path bais

    output:
    path "atac_fragment_sizes.tsv"   , emit: sizes
    path "atac_insert_sizes_mqc.tsv" , emit: mqc
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def max_size = task.ext.max_size ?: 2000
    """
    printf "sample\\tinsert_size\\tfragments\\tfragments_per_million\\n" > atac_fragment_sizes.tsv

    cat <<-END_MQC > atac_insert_sizes_mqc.tsv
    # id: atacportray_insert_sizes
    # section_name: ATAC insert size distribution
    # description: Fragment-size distribution from properly paired aligned ATAC-seq reads, counted once per pair using positive TLEN and normalized per million fragments.
    # plot_type: linegraph
    # pconfig:
    #   id: atacportray_insert_sizes
    #   title: ATAC insert size distribution
    #   xlab: Insert size (bp)
    #   ylab: Fragments per million
    #   xmax: ${max_size}
    END_MQC

    for bam in ${bams}; do
        sample=\$(basename "\${bam}")
        sample=\${sample%.filtered.sorted.bam}
        sample=\${sample%.sorted.bam}
        sample=\${sample%.bam}

        samtools view \\
            -@ ${task.cpus - 1} \\
            -f 2 \\
            -F 3844 \\
            "\${bam}" \\
            | awk -v sample="\${sample}" -v max_size=${max_size} '
            {
                size = \$9
                if (size > 0) {
                    total++
                    if (size <= max_size) {
                        counts[size]++
                    }
                }
            }

            END {
                printf "%s", sample >> "atac_insert_sizes_mqc.tsv"
                for (size = 1; size <= max_size; size++) {
                    count = counts[size] + 0
                    fpm = (total > 0) ? (count / total) * 1000000 : 0
                    printf "%s\\t%d\\t%d\\t%.6f\\n", sample, size, count, fpm >> "atac_fragment_sizes.tsv"
                    printf "\\t%.6f", fpm >> "atac_insert_sizes_mqc.tsv"
                }
                printf "\\n" >> "atac_insert_sizes_mqc.tsv"
            }
        '
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools version | sed '1!d;s/.* //')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_SIZES > atac_fragment_sizes.tsv
    sample	insert_size	fragments	fragments_per_million
    S1	1	0	0.000000
    S1	2	0	0.000000
    S2	1	0	0.000000
    S2	2	0	0.000000
    END_SIZES

    cat <<-END_MQC > atac_insert_sizes_mqc.tsv
    # id: atacportray_insert_sizes
    # section_name: ATAC insert size distribution
    # description: Fragment-size distribution from properly paired aligned ATAC-seq reads, counted once per pair using positive TLEN and normalized per million fragments.
    # plot_type: linegraph
    # pconfig:
    #   id: atacportray_insert_sizes
    #   title: ATAC insert size distribution
    #   xlab: Insert size (bp)
    #   ylab: Fragments per million
    #   xmax: 2000
    S1	0.000000	0.000000
    S2	0.000000	0.000000
    END_MQC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """
}
