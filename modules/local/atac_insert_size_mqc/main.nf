process ATAC_INSERT_SIZE_MQC {
    tag "atac_insert_sizes"
    label 'process_single'

    input:
    path stats_files

    output:
    path "atac_insert_sizes_mqc.json", emit: mqc
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def max_size = task.ext.max_size ?: 2000
    """
    awk -v max_size=${max_size} '
        FNR == 1 {
            sample = FILENAME
            sub(/^.*\\//, "", sample)
            sub(/\\.filtered\\.sorted\\.stats\$/, "", sample)
            sub(/\\.stats\$/, "", sample)
            samples[++ns] = sample
        }
        /^IS/ {
            size = \$2
            if (size <= max_size) {
                count[sample, size] = \$3
                seen[size] = 1
            }
        }
        END {
            print "{"
            print "  \\"id\\": \\"atacportray_insert_sizes\\","
            print "  \\"section_name\\": \\"ATAC insert sizes\\","
            print "  \\"plot_type\\": \\"linegraph\\","
            print "  \\"description\\": \\"Insert-size distribution from aligned, filtered ATAC-seq BAM files.\\","
            print "  \\"pconfig\\": {"
            print "    \\"id\\": \\"atacportray_insert_sizes\\","
            print "    \\"title\\": \\"ATAC insert size distribution\\","
            print "    \\"xlab\\": \\"Insert size (bp)\\","
            print "    \\"ylab\\": \\"Read pairs\\","
            print "    \\"xmax\\": " max_size
            print "  },"
            print "  \\"data\\": {"
            for (i = 1; i <= ns; i++) {
                sample = samples[i]
                printf "    \\"%s\\": {", sample
                for (size = 1; size <= max_size; size++) {
                    value = ((sample, size) in count) ? count[sample, size] : 0
                    printf "\\"%d\\": %d", size, value
                    if (size < max_size) {
                        printf ", "
                    }
                }
                printf "}"
                if (i < ns) {
                    printf ","
                }
                printf "\\n"
            }
            print "  }"
            print "}"
        }
    ' ${stats_files} > atac_insert_sizes_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -W version 2>&1 | sed -n '1s/.*Awk //p;1s/,.*//p')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_MQC > atac_insert_sizes_mqc.json
    {
      "id": "atacportray_insert_sizes",
      "section_name": "ATAC insert sizes",
      "plot_type": "linegraph",
      "description": "Insert-size distribution from aligned, filtered ATAC-seq BAM files.",
      "pconfig": {
        "id": "atacportray_insert_sizes",
        "title": "ATAC insert size distribution",
        "xlab": "Insert size (bp)",
        "ylab": "Read pairs",
        "xmax": 2000
      },
      "data": {
        "S1": {"1": 0, "2": 0},
        "S2": {"1": 0, "2": 0}
      }
    }
    END_MQC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: unknown
    END_VERSIONS
    """
}
