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

    for bw in ${bigwigs}; do
        sample=\$(basename "\$bw")
        sample=\${sample%.bigWig}
        sample=\${sample%.bw}
        sample=\${sample%.coverage}
        printf "%s\\n" "\$sample"
    done > sample_labels.txt

    python - <<'PY'
    import math
    from pathlib import Path

    before = int("${before}")
    after = int("${after}")
    binsize = int("${binsize}")
    nbins = int((before + after) / binsize)
    labels = [line.strip() for line in Path("sample_labels.txt").read_text().splitlines() if line.strip()]
    sums = [[0.0] * nbins for _ in labels]
    counts = [[0] * nbins for _ in labels]

    with open("atac_tss_matrix.tab") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("@"):
                continue
            fields = line.rstrip("\\n").split("\\t")
            raw_values = fields[6:]
            values = []
            for value in raw_values:
                try:
                    values.append(float(value))
                except ValueError:
                    values = []
                    break
            if len(values) < nbins * len(labels):
                continue
            for sample_index in range(len(labels)):
                offset = sample_index * nbins
                for bin_index, value in enumerate(values[offset:offset + nbins]):
                    if not math.isnan(value):
                        sums[sample_index][bin_index] += value
                        counts[sample_index][bin_index] += 1

    positions = [(-before + (idx * binsize) + (binsize // 2)) for idx in range(nbins)]
    data = {}
    for sample_index, label in enumerate(labels):
        values = {}
        for bin_index, position in enumerate(positions):
            denom = counts[sample_index][bin_index]
            values[str(position)] = sums[sample_index][bin_index] / denom if denom else 0.0
        data[label] = values

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
                "ylab": "Mean coverage signal"
            },
            "data": data
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
        "ylab": "Mean coverage signal"
      },
      "data": {
        "S1": {"-1000": 1, "0": 2, "1000": 1}
      }
    }
    END_MQC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: 3.5.6
    END_VERSIONS
    """
}
