process ROSE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/atacportray-rose:1.0.0' :
        'community.wave.seqera.io/library/atacportray-rose:1.0.0' }"

    input:
    tuple val(meta), path(peaks), path(bam), path(bai)
    val   genome_build
    val   stitch
    val   tss_exclusion

    output:
    tuple val(meta), path("*_SuperEnhancers.table.txt")        , emit: super_enhancers
    tuple val(meta), path("*_AllEnhancers.table.txt")          , emit: all_enhancers
    tuple val(meta), path("*_Enhancers_withSuper.bed")         , emit: enhancers_bed , optional: true
    tuple val(meta), path("*_SuperStitched.table.txt")         , emit: stitched      , optional: true
    tuple val(meta), path("*.pdf")                             , emit: plot          , optional: true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Convert narrowPeak/BED to the GFF ROSE expects
    rose_peaks_to_gff.py $peaks ${prefix}_peaks.gff

    ROSE_main.py \\
        -g ${genome_build} \\
        -i ${prefix}_peaks.gff \\
        -r $bam \\
        -o . \\
        -s ${stitch} \\
        -t ${tss_exclusion} \\
        $args

    # Normalize output names to the sample prefix
    for f in *_SuperEnhancers.table.txt *_AllEnhancers.table.txt; do
        [ -e "\$f" ] && mv "\$f" "${prefix}_\${f#*_}" 2>/dev/null || true
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rose: 1.0.0
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_SuperEnhancers.table.txt
    touch ${prefix}_AllEnhancers.table.txt
    touch ${prefix}_Enhancers_withSuper.bed
    touch ${prefix}_plot.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rose: 1.0.0
        bamliquidator: 1.0
    END_VERSIONS
    """
}
