process ROSE {
    tag "$meta.id"
    label 'process_medium'

    // St. Jude python3-compatible ROSE (github.com/stjude/ROSE). The official
    // image bundles ROSE_main.py + the genome annotation tables, so no local
    // build is required. Digest of v1.3.2: sha256:b46fadd8dfeaa0fe080a31aa216d4f34de9c68644dc053ae736cc612244e969b
    conda "${moduleDir}/environment.yml"
    container "ghcr.io/stjude/abralab/rose:v1.3.2"

    input:
    tuple val(meta), path(peaks), path(bam), path(bai)
    val   genome_build
    path  annotation
    val   stitch
    val   tss_exclusion

    output:
    tuple val(meta), path("*_SuperStitched.table.txt")           , emit: super_enhancers, optional: true
    tuple val(meta), path("*_AllStitched.table.txt")             , emit: all_enhancers  , optional: true
    tuple val(meta), path("*_Stitched_withSuper.bed")            , emit: enhancers_bed , optional: true
    tuple val(meta), path("*_SuperStitched.table_withGENES.txt") , emit: super_genes   , optional: true
    tuple val(meta), path("*_Plot_points.png")                   , emit: plot          , optional: true
    path "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # ROSE consumes the constituent-enhancer peaks as a BED/GFF. MACS3 narrowPeak
    # is BED6+4 with a unique peak name in column 4, which ROSE's internal
    # bedToGFF() uses as the region ID — so staging it with a .bed suffix lets
    # ROSE convert it internally (no external helper needed). ROSE derives the
    # output prefix from this filename up to the first dot.
    cp -L $peaks ${prefix}_peaks.bed

    # Prefer the user-supplied GTF converted to ROSE's UCSC/refSeq-like format.
    # Fall back to the annotation bundled with the ROSE image only when no GTF
    # annotation was provided by the pipeline.
    ROSE_HOME=\$(dirname \$(dirname \$(readlink -f \$(command -v ROSE_main.py))))
    GENOME_LC=\$(echo "${genome_build}" | tr '[:upper:]' '[:lower:]')
    if [ -s "$annotation" ]; then
        ANNOT="$annotation"
    else
        ANNOT="\${ROSE_HOME}/annotation/\${GENOME_LC}_refseq.ucsc"
    fi
    if [ ! -s "\$ANNOT" ]; then
        echo "ERROR: ROSE annotation not found at \$ANNOT" >&2
        echo "Available: \$(ls \${ROSE_HOME}/annotation/ 2>/dev/null | tr '\\n' ' ')" >&2
        exit 1
    fi

    ROSE_main.py \\
        --custom "\$ANNOT" \\
        -i ${prefix}_peaks.bed \\
        -r $bam \\
        -o . \\
        -s ${stitch} \\
        -t ${tss_exclusion} \\
        $args

    # ROSE's super-enhancer cutoff (ROSE_callSuper.R) fits a curve across the
    # ranked stitched enhancers. When a sample yields too few stitched regions
    # (e.g. very shallow data or a tiny reference), the fit cannot converge and
    # the *_SuperStitched / *_AllStitched tables are not written — yet
    # ROSE_main.py still exits 0. Surface that clearly instead of failing the
    # whole run; the corresponding outputs are declared optional.
    if [ ! -s ${prefix}_peaks_SuperStitched.table.txt ]; then
        echo "WARNING: ROSE produced no super-enhancer table for '${prefix}'. This usually means too few stitched enhancers to fit a cutoff (low peak count). Downstream super-enhancer steps will be skipped for this sample." >&2
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rose: \$(echo "\${ROSE_VERSION:-1.3.2}")
        python: \$(python3 --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_peaks_SuperStitched.table.txt
    touch ${prefix}_peaks_AllStitched.table.txt
    touch ${prefix}_peaks_Stitched_withSuper.bed
    touch ${prefix}_peaks_SuperStitched.table_withGENES.txt
    touch ${prefix}_peaks_Plot_points.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rose: 1.3.2
        python: 3.10.0
    END_VERSIONS
    """
}
