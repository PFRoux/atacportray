process ATAC_DESEQ2_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/bioconductor-deseq2:1.50.2--r45ha27e39d_0'

    input:
    tuple val(meta), path(counts), val(sample_info)
    val min_count
    val top_features

    output:
    tuple val(meta), path("atac_deseq2_results.tsv")           , emit: results
    tuple val(meta), path("atac_deseq2_normalized_counts.tsv") , emit: normalized_counts
    tuple val(meta), path("atac_deseq2_summary.tsv")           , emit: summary
    path "atac_deseq2_qc_mqc.yaml"                             , emit: mqc
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def metadata_tsv = "sample\tcondition\tpatient\n" + sample_info.collect { row ->
        "${row.sample}\t${row.condition ?: 'NA'}\t${row.patient ?: 'NA'}"
    }.join("\n")
    def metadata_b64 = metadata_tsv.bytes.encodeBase64().toString()
    """
    printf '%s' '${metadata_b64}' | base64 -d > sample_metadata.tsv

    Rscript ${moduleDir}/run_atac_deseq2_qc.R \\
        --counts ${counts} \\
        --metadata sample_metadata.tsv \\
        --min-count ${min_count} \\
        --top-features ${top_features}

    printf '%s\\n' \\
        '"${task.process}":' \\
        "    R: \$(Rscript --version 2>&1 | sed 's/Rscript (R) version //')" \\
        "    deseq2: \$(Rscript -e 'suppressPackageStartupMessages(library(DESeq2)); cat(as.character(packageVersion(\"DESeq2\")))' 2>/dev/null)" \\
        > versions.yml
    """

    stub:
    """
    cat <<-END_RESULTS > atac_deseq2_results.tsv
    contrast	region_id	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj
    High_vs_Low	chr1:100-200	100	1	0.2	5	0.001	0.01
    END_RESULTS

    cat <<-END_COUNTS > atac_deseq2_normalized_counts.tsv
    region_id	S1	S2
    chr1:100-200	10	20
    END_COUNTS

    cat <<-END_SUMMARY > atac_deseq2_summary.tsv
    status	message	samples	features_after_filter	min_count	deseq2_run
    ok	Stub	2	1	${min_count}	true
    END_SUMMARY

    cat <<-END_MQC > atac_deseq2_qc_mqc.yaml
    id: "atacportray_deseq2_qc"
    section_name: "ATAC differential accessibility QC"
    description: "DESeq2 normalization and exploratory analysis of raw counts over consensus peaks."
    plot_type: "html"
    data: |
      <div class="mqc-custom-content"><p>Stub DESeq2 QC.</p></div>
    END_MQC

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.4.0
        deseq2: 1.46.0
    END_VERSIONS
    """
}
