process DEEPTOOLS_TSS_PROFILE {
    tag "tss_profile"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0':
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    path bigwigs
    path tss_bed
    val before
    val after
    val binsize

    output:
    path "atac_tss_matrix.gz"         , emit: matrix
    path "atac_tss_matrix.tab"        , emit: matrix_tab
    path "atac_tss_profile.pdf"       , emit: profile_pdf
    path "atac_tss_profile.tab"       , emit: profile_data
    path "atac_tss_profile_mqc.json"  , emit: mqc
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    computeMatrix reference-point \\
        --referencePoint TSS \\
        --regionsFileName ${tss_bed} \\
        --scoreFileName ${bigwigs} \\
        --beforeRegionStartLength ${before} \\
        --afterRegionStartLength ${after} \\
        --binSize ${binsize} \\
        --numberOfProcessors ${task.cpus} \\
        --skipZeros \\
        --missingDataAsZero \\
        --outFileName atac_tss_matrix.gz \\
        --outFileNameMatrix atac_tss_matrix.tab \\
        ${args}

    plotProfile \\
        --matrixFile atac_tss_matrix.gz \\
        --outFileName atac_tss_profile.pdf \\
        --outFileNameData atac_tss_profile.tab \\
        --plotTitle "ATAC signal at TSS" \\
        --perGroup

    python - <<'PY'
    from pathlib import Path

    before = int("${before}")
    after = int("${after}")
    binsize = int("${binsize}")
    nbins = int((before + after) / binsize)
    positions = [(-before + (idx * binsize) + (binsize // 2)) for idx in range(nbins)]

    def clean_label(label):
        label = Path(label).name
        for suffix in (".coverage.bigWig", ".coverage.bw", ".bigWig", ".bw", ".coverage"):
            if label.endswith(suffix):
                label = label[:-len(suffix)]
                break
        return label

    data = []

    # plotProfile already writes the exact mean profile used for the PDF. Use it
    # for MultiQC too; parsing computeMatrix rows directly is fragile because the
    # matrix layout depends on grouping/sample options.
    with open("atac_tss_profile.tab") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\\n").split("\\t")
            if fields[0] in ("bin labels", "bins"):
                continue
            if len(fields) < nbins + 2:
                continue
            label = clean_label(fields[0])
            values = []
            for value in fields[2:2 + nbins]:
                try:
                    values.append(float(value))
                except ValueError:
                    values.append(0.0)
            data.append((label, values))

    import json
    with open("atac_tss_profile_mqc.json", "w") as out:
        json.dump({
            "id": "atacportray_tss_profile",
            "section_name": "ATAC signal at TSS",
            "plot_type": "linegraph",
            "description": "Average normalized bamCoverage signal around transcription start sites, computed from deepTools computeMatrix output.",
            "pconfig": {
                "id": "atacportray_tss_profile",
                "title": "ATAC signal at TSS",
                "xlab": "Distance from TSS (bp)",
                "ylab": "Mean coverage signal",
                "xDecimals": False
            },
            "data": {
                label: [[positions[idx], values[idx]] for idx in range(len(values))]
                for label, values in data
            }
        }, out)
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(computeMatrix --version | sed 's/computeMatrix //')
    END_VERSIONS
    """

    stub:
    """
    touch atac_tss_matrix.gz
    cat <<-END_MATRIX > atac_tss_matrix.tab
    chr1	1000	1001	gene1	0	+	1	2	3	4
    END_MATRIX
    touch atac_tss_profile.pdf
    cat <<-END_PROFILE > atac_tss_profile.tab
    sample	-1000	0	1000
    S1	1	2	1
    END_PROFILE
    cat <<-END_MQC > atac_tss_profile_mqc.json
    {
      "id": "atacportray_tss_profile",
      "section_name": "ATAC signal at TSS",
      "plot_type": "linegraph",
      "description": "Average normalized bamCoverage signal around transcription start sites, computed from deepTools computeMatrix output.",
      "pconfig": {
        "id": "atacportray_tss_profile",
        "title": "ATAC signal at TSS",
        "xlab": "Distance from TSS (bp)",
        "ylab": "Mean coverage signal",
        "xDecimals": false
      },
      "data": {
        "S1": [[-1000, 1], [0, 2], [1000, 1]]
      }
    }
    END_MQC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: 3.5.6
    END_VERSIONS
    """
}
