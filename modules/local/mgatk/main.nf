process MGATK {
    tag "mgatk_bulk"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgatk:0.7.0--pyhdfd78af_2' :
        'quay.io/biocontainers/mgatk:0.7.0--pyhdfd78af_2' }"

    input:
    path bams
    path bais
    val   mito_contig
    val   keep_duplicates

    output:
    path "*.mgatk/final/*.variant_stats.tsv.gz", emit: variant_stats, optional: true
    path "*.mgatk/final/*.vmr_strand.tsv.gz"   , emit: vmr          , optional: true
    path "*.mgatk/final/*.coverage.txt.gz"      , emit: coverage     , optional: true
    path "*.mgatk/final/*.rds"                  , emit: rds          , optional: true
    path "*.mgatk/final/*_refAllele.txt"        , emit: ref_allele   , optional: true
    path "*.mgatk/final"                       , emit: final_dir
    path "*.mgatk"                             , emit: outdir
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "mgatk_bulk"
    def keepdup = keep_duplicates ? "--keep-duplicates" : ""
    """
    mkdir -p mgatk_input
    for bam in ${bams}; do
        sample=\$(basename "\${bam}" .bam)
        ln -s ../\${bam} mgatk_input/\${sample}.bam
    done
    for bai in ${bais}; do
        bai_base=\$(basename "\${bai}")
        sample=\${bai_base%.bam.bai}
        sample=\${sample%.bai}
        ln -s ../\${bai} mgatk_input/\${sample}.bam.bai
    done

    mgatk call \\
        --input mgatk_input \\
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
    def prefix = task.ext.prefix ?: "mgatk_bulk"
    """
    mkdir -p ${prefix}.mgatk/final
    echo "" | gzip > ${prefix}.mgatk/final/${prefix}.variant_stats.tsv.gz
    echo "" | gzip > ${prefix}.mgatk/final/${prefix}.vmr_strand.tsv.gz
    echo "" | gzip > ${prefix}.mgatk/final/${prefix}.coverage.txt.gz
    touch ${prefix}.mgatk/final/${prefix}.rds
    touch ${prefix}.mgatk/final/MT_refAllele.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mgatk: 0.7.0
    END_VERSIONS
    """
}
