process BEDTOOLS_MULTICOV_COUNTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(regions), val(sample_names), path(bams), path(bais)

    output:
    tuple val(meta), path("*.raw_counts.tsv"), emit: counts
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def bam_args = bams.collect { it.toString() }.join(' ')
    def samples = sample_names.join('\t')
    """
    printf "region_id\\tchrom\\tstart\\tend\\t${samples}\\n" > ${prefix}.raw_counts.tsv

    bedtools multicov \\
        -bams ${bam_args} \\
        -bed ${regions} \\
        | awk 'BEGIN{OFS="\\t"} {print \$1":"\$2"-"\$3,\$0}' \\
        >> ${prefix}.raw_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e 's/bedtools v//g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    def samples = sample_names.join('\t')
    """
    printf "region_id\\tchrom\\tstart\\tend\\t${samples}\\n" > ${prefix}.raw_counts.tsv
    printf "chr1:100-200\\tchr1\\t100\\t200" >> ${prefix}.raw_counts.tsv
    for sample in ${sample_names.join(' ')}; do
        printf "\\t0" >> ${prefix}.raw_counts.tsv
    done
    printf "\\n" >> ${prefix}.raw_counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: 2.31.1
    END_VERSIONS
    """
}
