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
    """
    cat > sample_metadata.tsv <<'END_METADATA'
    ${metadata_tsv}
    END_METADATA

    Rscript - <<'RSCRIPT'
    suppressPackageStartupMessages({
        library(DESeq2)
    })

    counts_file <- "${counts}"
    min_count <- as.numeric("${min_count}")
    top_features <- as.integer("${top_features}")

    counts <- read.delim(counts_file, check.names = FALSE)
    meta <- read.delim("sample_metadata.tsv", check.names = FALSE)
    meta <- meta[!duplicated(meta\$sample), , drop = FALSE]

    sample_cols <- intersect(meta\$sample, colnames(counts))
    meta <- meta[match(sample_cols, meta\$sample), , drop = FALSE]
    count_mat <- as.matrix(counts[, sample_cols, drop = FALSE])
    rownames(count_mat) <- counts\$region_id
    storage.mode(count_mat) <- "integer"

    keep <- rowSums(count_mat) >= min_count
    count_mat <- count_mat[keep, , drop = FALSE]

    status <- "ok"
    message <- "DESeq2 differential analysis completed."
    did_deseq <- FALSE
    n_features <- nrow(count_mat)
    n_samples <- ncol(count_mat)

    if (n_features == 0 || n_samples < 2) {
        status <- "skipped"
        message <- "Not enough count data to run normalization, PCA or heatmap."
    }

    norm_counts <- matrix(nrow = 0, ncol = 0)
    transform_mat <- matrix(nrow = 0, ncol = 0)
    res_all <- data.frame()

    if (status == "ok") {
        meta\$condition <- factor(meta\$condition)
        dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = meta, design = ~ condition)
        dds <- estimateSizeFactors(dds, type = "poscounts")
        norm_counts <- counts(dds, normalized = TRUE)
        transform_mat <- log2(norm_counts + 1)

        condition_counts <- table(meta\$condition)
        valid_deseq <- length(condition_counts) >= 2 && min(condition_counts) >= 2

        if (valid_deseq) {
            dds <- DESeq(dds, quiet = TRUE)
            conds <- levels(meta\$condition)
            for (i in seq_len(length(conds) - 1)) {
                for (j in seq((i + 1), length(conds))) {
                    contrast <- c("condition", conds[j], conds[i])
                    res <- as.data.frame(results(dds, contrast = contrast))
                    res\$region_id <- rownames(res)
                    res\$contrast <- paste0(conds[j], "_vs_", conds[i])
                    res_all <- rbind(res_all, res[, c("contrast", "region_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")])
                }
            }
            did_deseq <- TRUE
        } else {
            status <- "partial"
            message <- "PCA and heatmap were generated, but DESeq2 differential testing was skipped because at least two replicates per condition are required."
        }
    }

    if (nrow(res_all) == 0) {
        res_all <- data.frame(
            contrast = character(),
            region_id = character(),
            baseMean = numeric(),
            log2FoldChange = numeric(),
            lfcSE = numeric(),
            stat = numeric(),
            pvalue = numeric(),
            padj = numeric()
        )
    }
    write.table(res_all, "atac_deseq2_results.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)

    if (length(norm_counts) > 0) {
        norm_df <- data.frame(region_id = rownames(norm_counts), norm_counts, check.names = FALSE)
    } else {
        norm_df <- data.frame(region_id = character())
    }
    write.table(norm_df, "atac_deseq2_normalized_counts.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)

    summary_df <- data.frame(
        status = status,
        message = message,
        samples = n_samples,
        features_after_filter = n_features,
        min_count = min_count,
        deseq2_run = did_deseq
    )
    write.table(summary_df, "atac_deseq2_summary.tsv", sep = "\\t", quote = FALSE, row.names = FALSE)

    plot_files <- character()
    if (nrow(transform_mat) >= 2 && ncol(transform_mat) >= 2) {
        vars <- apply(transform_mat, 1, var)
        selected <- names(sort(vars, decreasing = TRUE))[seq_len(min(top_features, length(vars)))]
        plot_mat <- transform_mat[selected, , drop = FALSE]

        pca <- prcomp(t(plot_mat), center = TRUE, scale. = FALSE)
        variance <- round((pca\$sdev^2 / sum(pca\$sdev^2))[seq_len(min(2, length(pca\$sdev)))] * 100, 1)
        pca_scores <- pca\$x[, seq_len(min(2, ncol(pca\$x))), drop = FALSE]
        if (ncol(pca_scores) == 1) {
            pca_scores <- cbind(pca_scores, PC2 = 0)
            variance <- c(variance, 0)
        }
        pca_conditions <- meta\$condition[match(rownames(pca_scores), meta\$sample)]
        condition_levels <- unique(as.character(pca_conditions))
        palette <- grDevices::hcl.colors(max(3, length(condition_levels)), "Dark 3")
        point_cols <- palette[match(as.character(pca_conditions), condition_levels)]
        grDevices::svg("atac_deseq2_pca.svg", width = 7, height = 5)
        plot(
            pca_scores[, 1],
            pca_scores[, 2],
            pch = 19,
            col = point_cols,
            xlab = paste0("PC1 (", variance[1], "%)"),
            ylab = paste0("PC2 (", variance[2], "%)"),
            main = "Consensus peak PCA"
        )
        text(pca_scores[, 1], pca_scores[, 2], labels = rownames(pca_scores), pos = 3, cex = 0.7)
        legend("topright", legend = condition_levels, col = palette[seq_along(condition_levels)], pch = 19, bty = "n", cex = 0.8)
        grDevices::dev.off()
        plot_files <- c(plot_files, "atac_deseq2_pca.svg")

        scaled <- t(scale(t(plot_mat)))
        scaled[is.na(scaled)] <- 0
        grDevices::svg("atac_deseq2_heatmap.svg", width = max(6, ncol(plot_mat) * 0.35), height = 6)
        op <- par(mar = c(8, 4, 3, 6))
        image(
            x = seq_len(ncol(scaled)),
            y = seq_len(nrow(scaled)),
            z = t(scaled[nrow(scaled):1, , drop = FALSE]),
            col = grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(101),
            axes = FALSE,
            xlab = "",
            ylab = "",
            main = "Top variable consensus peaks"
        )
        axis(1, at = seq_len(ncol(scaled)), labels = colnames(scaled), las = 2, cex.axis = 0.7)
        axis(2, at = c(1, nrow(scaled)), labels = c(nrow(scaled), 1), las = 1, cex.axis = 0.7)
        mtext("Consensus peaks ranked by variance", side = 2, line = 2.5)
        par(op)
        grDevices::dev.off()
        plot_files <- c(plot_files, "atac_deseq2_heatmap.svg")
    }

    html_parts <- c(
        "<div class=\\"mqc-custom-content\\">",
        sprintf("<p><strong>Status:</strong> %s. %s</p>", status, message),
        sprintf("<p><strong>Samples:</strong> %s &nbsp; <strong>Features after filtering:</strong> %s</p>", n_samples, n_features)
    )
    for (plot_file in plot_files) {
        encoded <- paste(system2("base64", plot_file, stdout = TRUE), collapse = "")
        alt <- gsub("_", " ", tools::file_path_sans_ext(plot_file), fixed = TRUE)
        html_parts <- c(html_parts, sprintf("<img alt=\\"%s\\" style=\\"width:100%%;max-width:900px;margin:0.5rem 0;\\" src=\\"data:image/svg+xml;base64,%s\\">", alt, encoded))
    }
    if (length(plot_files) == 0) {
        html_parts <- c(html_parts, "<p>No PCA/heatmap could be generated from the filtered count matrix.</p>")
    }
    html_parts <- c(html_parts, "</div>")
    html <- paste(html_parts, collapse = "\\n")

    cat('id: "atacportray_deseq2_qc"\\n', file = "atac_deseq2_qc_mqc.yaml")
    cat('section_name: "ATAC differential accessibility QC"\\n', file = "atac_deseq2_qc_mqc.yaml", append = TRUE)
    cat('description: "DESeq2 normalization and exploratory analysis of raw counts over consensus peaks."\\n', file = "atac_deseq2_qc_mqc.yaml", append = TRUE)
    cat('plot_type: "html"\\n', file = "atac_deseq2_qc_mqc.yaml", append = TRUE)
    cat('data: |\\n', file = "atac_deseq2_qc_mqc.yaml", append = TRUE)
    cat(paste0("  ", gsub("\\n", "\\n  ", html), "\\n"), file = "atac_deseq2_qc_mqc.yaml", append = TRUE)
    RSCRIPT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(Rscript --version 2>&1 | sed 's/Rscript (R) version //')
        deseq2: \$(Rscript -e 'suppressPackageStartupMessages(library(DESeq2)); cat(as.character(packageVersion("DESeq2")))' 2>/dev/null)
    END_VERSIONS
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
