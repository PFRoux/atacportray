process TOBIAS_BINDETECT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tobias:0.16.1--py39h24fbfe6_1' :
        'quay.io/biocontainers/tobias:0.16.1--py39h24fbfe6_1' }"

    input:
    tuple val(meta), path(footprints), val(cond_names)
    path motifs
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(peaks)

    output:
    tuple val(meta), path("bindetect/*bindetect_results.txt"), emit: results
    tuple val(meta), path("bindetect/*bindetect_results.xlsx"), emit: results_xlsx, optional: true
    tuple val(meta), path("bindetect/*.pdf")                  , emit: figures    , optional: true
    tuple val(meta), path("bindetect/*_overview.txt")         , emit: overview   , optional: true
    tuple val(meta), path("bindetect/")                       , emit: outdir
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cond   = cond_names ? "--cond_names ${cond_names}" : ''
    """
    mkdir -p .matplotlib
    export MPLCONFIGDIR="\${PWD}/.matplotlib"
    export PYTHONPATH="\${PWD}:\${PYTHONPATH:-}"

    python - <<'PY'
    import pyBigWig
    from pathlib import Path

    peaks = Path("${peaks}")
    filtered = Path("bindetect_peaks.filtered.bed")
    footprint_paths = "${footprints}".split()
    chrom_sizes_by_file = {}
    for footprint in footprint_paths:
        bw = pyBigWig.open(footprint)
        chrom_sizes_by_file[footprint] = bw.chroms()
        bw.close()

    kept = 0
    dropped = 0
    with peaks.open() as handle, filtered.open("w") as out:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 3:
                dropped += 1
                continue
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                dropped += 1
                continue
            if start < 0 or end <= start:
                dropped += 1
                continue
            valid = True
            for footprint, chrom_sizes in chrom_sizes_by_file.items():
                chrom_size = chrom_sizes.get(chrom)
                if chrom_size is None or end > chrom_size:
                    valid = False
                    break
            if not valid:
                dropped += 1
                continue
            out.write(line)
            kept += 1

    if kept == 0:
        raise SystemExit(f"No peaks from {peaks} overlap chromosomes available in all footprint bigWigs; dropped {dropped} regions.")
    if dropped:
        print(f"Filtered {dropped} peak regions absent from or outside one or more footprint bigWigs; kept {kept} regions for TOBIAS BINDetect.")
    PY

    printf '%s\\n' \\
        'import inspect' \\
        '' \\
        'try:' \\
        '    import adjustText' \\
        '' \\
        '    _adjust_text = adjustText.adjust_text' \\
        '    _params = inspect.signature(_adjust_text).parameters' \\
        '' \\
        '    def _compat_adjust_text(*args, **kwargs):' \\
        '        if "add_objects" in kwargs and "add_objects" not in _params:' \\
        '            if "objects" not in kwargs:' \\
        '                kwargs["objects"] = kwargs.pop("add_objects")' \\
        '            else:' \\
        '                kwargs.pop("add_objects")' \\
        '        if "text_from_points" in kwargs and "text_from_points" not in _params:' \\
        '            kwargs.pop("text_from_points")' \\
        '        return _adjust_text(*args, **kwargs)' \\
        '' \\
        '    adjustText.adjust_text = _compat_adjust_text' \\
        'except Exception:' \\
        '    pass' \\
        > sitecustomize.py

    TOBIAS BINDetect \\
        --motifs $motifs \\
        --signals $footprints \\
        --genome $fasta \\
        --peaks bindetect_peaks.filtered.bed \\
        --outdir bindetect \\
        --cores $task.cpus \\
        $cond \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tobias: "\$(TOBIAS --version 2>&1 | tail -n 1 | sed 's/TOBIAS //g')"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p bindetect
    touch bindetect/${prefix}_bindetect_results.txt
    touch bindetect/${prefix}_bindetect_results.xlsx
    touch bindetect/${prefix}_bindetect_overview.txt
    touch bindetect/${prefix}_figures.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tobias: 0.16.1
    END_VERSIONS
    """
}
