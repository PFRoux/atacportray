process SAMTOOLS_MT_COVERAGE {
    tag "mt_coverage"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    path bams
    path bais
    val  mito_contig

    output:
    path "samtools_mt_coverage.tsv", emit: coverage
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    set -euo pipefail

    mito_region="${mito_contig}"
    if [[ -f "\${mito_region}" ]]; then
        mito_region=\$(grep -m 1 '^>' "\${mito_region}" | sed 's/^>//; s/[[:space:]].*//')
    fi

    printf 'sample\\tposition\\tcoverage\\n' > samtools_mt_coverage.tsv

    for bam in ${bams}; do
        sample=\$(basename "\${bam}" .bam)
        sample=\${sample%.md}
        samtools depth \\
            -a \\
            -@ ${task.cpus - 1} \\
            -r "\${mito_region}" \\
            "\${bam}" \\
            | awk -v sample="\${sample}" 'BEGIN{OFS="\\t"} {print sample,\$2,\$3}' \\
            >> samtools_mt_coverage.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1s/^samtools //p')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_COVERAGE > samtools_mt_coverage.tsv
    sample	position	coverage
    sample1	1	1
    END_COVERAGE

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: 1.23.1
    END_VERSIONS
    """
}
