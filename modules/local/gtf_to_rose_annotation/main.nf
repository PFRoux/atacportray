process GTF_TO_ROSE_ANNOTATION {
    tag "rose_annotation"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0':
        'quay.io/biocontainers/mulled-v2-eb9e7907c7a753917c1e4d7a64384c047429618a:28424fe3aec58d2b3e4e4390025d886207657d25-0' }"

    input:
    path gtf

    output:
    path "rose_refseq.ucsc", emit: annotation
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

    def attr_value(attrs, key):
        match = re.search(rf'{key} "([^"]+)"', attrs)
        return match.group(1) if match else None

    transcripts = {}
    genes = {}

    with opener(gtf, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 9:
                continue
            chrom, feature, start, end, strand, attrs = fields[0], fields[2], int(fields[3]), int(fields[4]), fields[6], fields[8]
            if strand not in {"+", "-"}:
                continue

            gene_id = attr_value(attrs, "gene_id")
            transcript_id = attr_value(attrs, "transcript_id")
            gene_name = attr_value(attrs, "gene_name") or attr_value(attrs, "gene") or gene_id or transcript_id
            if not gene_id and not transcript_id:
                continue

            if feature == "gene":
                key = transcript_id or gene_id
                genes[key] = (chrom, start - 1, end, strand, gene_name or key)
            elif feature in {"transcript", "exon"}:
                key = transcript_id or gene_id
                if key not in transcripts:
                    transcripts[key] = [chrom, start - 1, end, strand, gene_name or gene_id or key]
                else:
                    rec = transcripts[key]
                    rec[1] = min(rec[1], start - 1)
                    rec[2] = max(rec[2], end)

    records = transcripts if transcripts else genes
    with open("rose_refseq.ucsc", "w") as out:
        out.write("\\t".join([
            "bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
            "exonCount", "exonStarts", "exonEnds", "score", "name2",
            "cdsStartStat", "cdsEndStat", "exonFrames"
        ]) + "\\n")

        for idx, (name, rec) in enumerate(sorted(records.items(), key=lambda item: (item[1][0], item[1][1], item[0]))):
            chrom, start, end, strand, gene_name = rec
            start = max(0, int(start))
            end = max(start + 1, int(end))
            row = [
                "0", name, chrom, strand, str(start), str(end), str(start), str(end),
                "1", f"{start},", f"{end},", "0", gene_name or name, "unk", "unk", "-1,"
            ]
            out.write("\\t".join(row) + "\\n")
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_ANNOT > rose_refseq.ucsc
    bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
    0	tx1	chr1	+	1000	2000	1000	2000	1	1000,	2000,	0	gene1	unk	unk	-1,
    END_ANNOT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: 3.12.0
    END_VERSIONS
    """
}
