process MAFTOOLS_ONCOPLOT {
    tag "oncoplot"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-maftools:2.26.0--r45h01b2380_0' :
        'quay.io/biocontainers/bioconductor-maftools:2.26.0--r45h01b2380_0' }"

    input:
    path mafs, stageAs: "mafs/*"

    output:
    path "*.oncoplot.pdf"       , emit: oncoplot
    path "*.summary.pdf"        , emit: summary   , optional: true
    path "*.maf_summary.txt"    , emit: summary_txt, optional: true
    path "merged.maf.gz"        , emit: merged_maf , optional: true
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    maftools_oncoplot.R --input_dir mafs --out_prefix atacportray $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e 'cat(strsplit(R.version.string, " ")[[1]][3])')
        maftools: \$(Rscript -e 'cat(as.character(packageVersion("maftools")))')
    END_VERSIONS
    """

    stub:
    """
    touch atacportray.oncoplot.pdf
    touch atacportray.summary.pdf
    touch atacportray.maf_summary.txt
    echo "" | gzip > merged.maf.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: 4.4.0
        maftools: 2.26.0
    END_VERSIONS
    """
}
