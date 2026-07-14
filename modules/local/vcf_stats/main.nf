process VCF_STATS {
    tag "${meta.id}:${meta.caller}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data' :
        'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.bcftools_stats.txt"), emit: stats
    tuple val(meta), path("*_mqc.tsv")           , emit: mqc
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def caller = meta.caller ?: 'variants'
    def prefix = task.ext.prefix ?: "${meta.id}.${caller}"
    """
    bcftools stats \\
        ${args} \\
        ${vcf} \\
        > ${prefix}.bcftools_stats.txt

    awk -v sample="${meta.id}" -v caller="${caller}" '
        BEGIN {
            records = snps = mnps = indels = others = multiallelic = 0
            sample_caller = sample "." caller
        }
        /^SN/ {
            key = \$3
            for (i = 4; i < NF; i++) {
                key = key " " \$i
            }
            sub(":\$", "", key)
            value = \$NF
            if (key == "number of records") records = value
            else if (key == "number of SNPs") snps = value
            else if (key == "number of MNPs") mnps = value
            else if (key == "number of indels") indels = value
            else if (key == "number of others") others = value
            else if (key == "number of multiallelic sites") multiallelic = value
        }
        END {
            print "# id: atacportray_variant_calls"
            print "# section_name: Variant calls"
            print "# description: Number and type of variants reported by each caller after atacportray filtering."
            print "# plot_type: table"
            print "Sample_caller\\tSample\\tCaller\\tRecords\\tSNPs\\tMNPs\\tIndels\\tOther\\tMultiallelic"
            print sample_caller "\\t" sample "\\t" caller "\\t" records "\\t" snps "\\t" mnps "\\t" indels "\\t" others "\\t" multiallelic
        }
    ' ${prefix}.bcftools_stats.txt > ${prefix}.variant_counts_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | sed '1!d; s/^.*bcftools //')
    END_VERSIONS
    """

    stub:
    def caller = meta.caller ?: 'variants'
    def prefix = task.ext.prefix ?: "${meta.id}.${caller}"
    """
    touch ${prefix}.bcftools_stats.txt
    cat <<-END_MQC > ${prefix}.variant_counts_mqc.tsv
    # id: atacportray_variant_calls
    # section_name: Variant calls
    # description: Number and type of variants reported by each caller after atacportray filtering.
    # plot_type: table
    Sample_caller	Sample	Caller	Records	SNPs	MNPs	Indels	Other	Multiallelic
    ${meta.id}.${caller}	${meta.id}	${caller}	0	0	0	0	0	0
    END_MQC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: 1.23.1
    END_VERSIONS
    """
}
