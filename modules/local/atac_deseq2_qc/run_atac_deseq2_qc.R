suppressPackageStartupMessages({
    library(DESeq2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
    idx <- match(flag, args)
    if (is.na(idx) || idx == length(args)) {
        return(default)
    }
    args[[idx + 1]]
}

counts_file <- get_arg("--counts")
metadata_file <- get_arg("--metadata")
min_count <- as.numeric(get_arg("--min-count", "10"))
top_features <- as.integer(get_arg("--top-features", "50"))

if (is.null(counts_file) || is.null(metadata_file)) {
    stop("--counts and --metadata are required", call. = FALSE)
}

empty_results <- function() {
    data.frame(
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

write_mqc <- function(status, message, n_samples, n_features, plot_files) {
    html_parts <- c(
        "<div class=\"mqc-custom-content\">",
        sprintf("<p><strong>Status:</strong> %s. %s</p>", status, message),
        sprintf(
            "<p><strong>Samples:</strong> %s &nbsp; <strong>Features after filtering:</strong> %s</p>",
            n_samples,
            n_features
        )
    )

    for (plot_file in plot_files) {
        encoded <- paste(system2("base64", plot_file, stdout = TRUE), collapse = "")
        alt <- gsub("_", " ", tools::file_path_sans_ext(plot_file), fixed = TRUE)
        html_parts <- c(
            html_parts,
            sprintf(
                "<img alt=\"%s\" style=\"width:100%%;max-width:900px;margin:0.5rem 0;\" src=\"data:image/svg+xml;base64,%s\">",
                alt,
                encoded
            )
        )
    }

    if (length(plot_files) == 0) {
        html_parts <- c(html_parts, "<p>No PCA/heatmap could be generated from the filtered count matrix.</p>")
    }
    html_parts <- c(html_parts, "</div>")
    html <- paste(html_parts, collapse = "\n")

    cat('id: "atacportray_deseq2_qc"\n', file = "atac_deseq2_qc_mqc.yaml")
    cat('section_name: "ATAC differential accessibility QC"\n', file = "atac_deseq2_qc_mqc.yaml", append = TRUE)
    cat('description: "DESeq2 normalization and exploratory analysis of raw counts over consensus peaks."\n', file = "atac_deseq2_qc_mqc.yaml", append = TRUE)
    cat('plot_type: "html"\n', file = "atac_deseq2_qc_mqc.yaml", append = TRUE)
    cat('data: |\n', file = "atac_deseq2_qc_mqc.yaml", append = TRUE)
    cat(paste0("  ", gsub("\n", "\n  ", html), "\n"), file = "atac_deseq2_qc_mqc.yaml", append = TRUE)
}

counts <- read.delim(counts_file, check.names = FALSE)
meta <- read.delim(metadata_file, check.names = FALSE)
meta <- meta[!duplicated(meta$sample), , drop = FALSE]

sample_cols <- intersect(meta$sample, colnames(counts))
meta <- meta[match(sample_cols, meta$sample), , drop = FALSE]

if (length(sample_cols) > 0) {
    count_mat <- as.matrix(counts[, sample_cols, drop = FALSE])
    rownames(count_mat) <- counts$region_id
    storage.mode(count_mat) <- "integer"
    keep <- rowSums(count_mat) >= min_count
    count_mat <- count_mat[keep, , drop = FALSE]
} else {
    count_mat <- matrix(nrow = 0, ncol = 0)
}

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
res_all <- empty_results()

if (status == "ok") {
    meta$condition <- factor(meta$condition)
    dds <- DESeqDataSetFromMatrix(countData = count_mat, colData = meta, design = ~ condition)
    dds <- estimateSizeFactors(dds, type = "poscounts")
    norm_counts <- counts(dds, normalized = TRUE)
    transform_mat <- log2(norm_counts + 1)

    condition_counts <- table(meta$condition)
    valid_deseq <- length(condition_counts) >= 2 && min(condition_counts) >= 2

    if (valid_deseq) {
        dds <- DESeq(dds, quiet = TRUE)
        conds <- levels(meta$condition)
        for (i in seq_len(length(conds) - 1)) {
            for (j in seq((i + 1), length(conds))) {
                contrast <- c("condition", conds[j], conds[i])
                res <- as.data.frame(results(dds, contrast = contrast))
                res$region_id <- rownames(res)
                res$contrast <- paste0(conds[j], "_vs_", conds[i])
                res_all <- rbind(
                    res_all,
                    res[, c("contrast", "region_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
                )
            }
        }
        did_deseq <- TRUE
    } else {
        status <- "partial"
        message <- paste(
            "PCA and heatmap were generated, but DESeq2 differential testing was skipped because",
            "at least two replicates per condition are required."
        )
    }
}

write.table(res_all, "atac_deseq2_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

if (length(norm_counts) > 0) {
    norm_df <- data.frame(region_id = rownames(norm_counts), norm_counts, check.names = FALSE)
} else {
    norm_df <- data.frame(region_id = character())
}
write.table(norm_df, "atac_deseq2_normalized_counts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

summary_df <- data.frame(
    status = status,
    message = message,
    samples = n_samples,
    features_after_filter = n_features,
    min_count = min_count,
    deseq2_run = did_deseq
)
write.table(summary_df, "atac_deseq2_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

plot_files <- character()
if (nrow(transform_mat) >= 2 && ncol(transform_mat) >= 2) {
    vars <- apply(transform_mat, 1, var)
    selected <- names(sort(vars, decreasing = TRUE))[seq_len(min(top_features, length(vars)))]
    plot_mat <- transform_mat[selected, , drop = FALSE]

    pca <- prcomp(t(plot_mat), center = TRUE, scale. = FALSE)
    variance <- round((pca$sdev^2 / sum(pca$sdev^2))[seq_len(min(2, length(pca$sdev)))] * 100, 1)
    pca_scores <- pca$x[, seq_len(min(2, ncol(pca$x))), drop = FALSE]
    if (ncol(pca_scores) == 1) {
        pca_scores <- cbind(pca_scores, PC2 = 0)
        variance <- c(variance, 0)
    }
    pca_conditions <- meta$condition[match(rownames(pca_scores), meta$sample)]
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

write_mqc(status, message, n_samples, n_features, plot_files)
