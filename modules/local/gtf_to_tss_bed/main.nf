process GTF_TO_TSS_BED {
    tag "tss_bed"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0':
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    path gtf

    output:
    path "tss.bed", emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python - <<'PY'
    import gzip
    import re
    from pathlib import Path

    gtf = Path("${gtf}")
    opener = gzip.open if gtf.name.endswith(".gz") else open
    rows = set()
    with opener(gtf, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 9 or fields[2] not in {"transcript", "gene"}:
                continue
            chrom, start, end, strand, attrs = fields[0], int(fields[3]), int(fields[4]), fields[6], fields[8]
            if strand == "+":
                tss = start - 1
            elif strand == "-":
                tss = end - 1
            else:
                continue
            tss = max(tss, 0)
            match = re.search(r'gene_id "([^"]+)"', attrs) or re.search(r'transcript_id "([^"]+)"', attrs)
            name = match.group(1) if match else "."
            rows.add((chrom, tss, tss + 1, name, "0", strand))

    with open("tss.bed", "w") as out:
        for row in sorted(rows, key=lambda item: (item[0], item[1], item[5], item[3])):
            out.write("\\t".join(map(str, row)) + "\\n")
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_BED > tss.bed
    chr1	1000	1001	gene1	0	+
    chr1	2000	2001	gene2	0	-
    END_BED

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.0
    END_VERSIONS
    """
}
