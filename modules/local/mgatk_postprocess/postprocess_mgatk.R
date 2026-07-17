args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = "") {
    idx <- match(flag, args)
    if (is.na(idx) || idx == length(args)) return(default)
    args[[idx + 1]]
}

split_paths <- function(x) {
    x <- trimws(x)
    if (!nzchar(x)) return(character())
    unlist(strsplit(x, "\\s+"))
}

open_text <- function(path) {
    if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
}

read_table <- function(path, header = TRUE) {
    if (!file.exists(path) || file.info(path)$size == 0) return(NULL)
    con <- open_text(path)
    on.exit(close(con), add = TRUE)
    first_line <- tryCatch(readLines(con, n = 1, warn = FALSE), error = function(e) character())
    if (!length(first_line)) return(NULL)
    sep <- if (grepl(",", first_line) && !grepl("\t", first_line)) "," else "\t"
    tryCatch(
        read.table(
            path,
            sep = sep,
            header = header,
            check.names = FALSE,
            comment.char = "",
            quote = "",
            stringsAsFactors = FALSE
        ),
        error = function(e) NULL
    )
}

clean_sample <- function(path, suffix) {
    x <- basename(path)
    sub(suffix, "", x)
}

clean_sample_value <- function(x) {
    x <- basename(as.character(x))
    x <- sub("\\.coverage\\.txt\\.gz$", "", x)
    x <- sub("\\.coverage\\.txt$", "", x)
    x <- sub("\\.bam$", "", x)
    x <- sub("\\.md$", "", x)
    x
}

rds_files <- split_paths(get_arg("--rds"))
coverage_files <- split_paths(get_arg("--coverage"))
min_af <- as.numeric(get_arg("--min-af", "0.05"))
min_vmr <- as.numeric(get_arg("--min-vmr", "0.01"))
min_sd <- as.numeric(get_arg("--min-sd", "0.05"))

write_empty <- function() {
    write.table(data.frame(metric = c("total_variants", "high_confidence_variants", "variable_variants"), value = 0), "mgatk_variant_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(data.frame(variant = character(), mean_af = numeric(), max_af = numeric(), sd_af = numeric(), vmr = numeric(), detected_samples = integer(), variable = logical()), "mgatk_high_confidence_variants.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(data.frame(variant = character()), "mgatk_heteroplasmy_matrix.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
}

extract_af_from_se <- function(obj, sample_name) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) return(NULL)
    assays <- tryCatch(SummarizedExperiment::assays(obj), error = function(e) NULL)
    if (is.null(assays)) return(NULL)
    assay_names <- names(assays)
    cov_name <- assay_names[grepl("^coverage$|coverage", assay_names, ignore.case = TRUE)][1]
    if (is.na(cov_name)) return(NULL)
    cov <- as.matrix(assays[[cov_name]])
    if (ncol(cov) == 0) return(NULL)

    rr <- tryCatch(SummarizedExperiment::rowRanges(obj), error = function(e) NULL)
    pos <- tryCatch(as.integer(GenomicRanges::start(rr)), error = function(e) seq_len(nrow(cov)))
    ref <- tryCatch(as.character(S4Vectors::mcols(rr)$refAllele), error = function(e) rep("N", nrow(cov)))
    if (length(ref) != nrow(cov)) ref <- rep("N", nrow(cov))

    allele_mats <- list()
    for (base in c("A", "C", "G", "T")) {
        fw <- paste0(base, "_counts_fw")
        rv <- paste0(base, "_counts_rev")
        plain <- paste0(base, "_counts")
        if (fw %in% assay_names && rv %in% assay_names) {
            allele_mats[[base]] <- as.matrix(assays[[fw]]) + as.matrix(assays[[rv]])
        } else if (plain %in% assay_names) {
            allele_mats[[base]] <- as.matrix(assays[[plain]])
        }
    }
    if (!length(allele_mats)) return(NULL)

    out <- list()
    for (base in names(allele_mats)) {
        keep <- toupper(ref) != base & toupper(ref) != "N"
        if (!any(keep)) next
        af <- allele_mats[[base]][keep, , drop = FALSE] / pmax(cov[keep, , drop = FALSE], 1)
        rownames(af) <- paste0(pos[keep], toupper(ref[keep]), ">", base)
        out[[base]] <- af
    }
    if (!length(out)) return(NULL)
    mat <- do.call(rbind, out)
    if (is.null(colnames(mat))) colnames(mat) <- sample_name
    mat
}

af_mats <- list()
for (path in rds_files) {
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(obj)) next
    sample_name <- clean_sample(path, "\\.rds$")
    mat <- extract_af_from_se(obj, sample_name)
    if (!is.null(mat)) af_mats[[sample_name]] <- mat
}

merge_af <- function(mats) {
    if (!length(mats)) return(NULL)
    all_vars <- unique(unlist(lapply(mats, rownames)))
    all_samples <- unique(unlist(lapply(mats, colnames)))
    out <- matrix(0, nrow = length(all_vars), ncol = length(all_samples), dimnames = list(all_vars, all_samples))
    for (mat in mats) {
        out[rownames(mat), colnames(mat)] <- pmax(out[rownames(mat), colnames(mat)], mat, na.rm = TRUE)
    }
    out
}

af_matrix <- merge_af(af_mats)

if (is.null(af_matrix)) {
    write_empty()
} else {
    mean_af <- rowMeans(af_matrix, na.rm = TRUE)
    max_af <- apply(af_matrix, 1, max, na.rm = TRUE)
    sd_af <- apply(af_matrix, 1, stats::sd, na.rm = TRUE)
    vmr <- apply(af_matrix, 1, stats::var, na.rm = TRUE) / pmax(mean_af, 1e-8)
    detected_samples <- rowSums(af_matrix >= min_af, na.rm = TRUE)
    high_conf <- detected_samples >= 1 & vmr >= min_vmr
    variable <- sd_af >= min_sd

    matrix_df <- data.frame(variant = rownames(af_matrix), af_matrix, check.names = FALSE)
    write.table(matrix_df, "mgatk_heteroplasmy_matrix.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

    variants_df <- data.frame(
        variant = rownames(af_matrix),
        mean_af = mean_af,
        max_af = max_af,
        sd_af = sd_af,
        vmr = vmr,
        detected_samples = detected_samples,
        variable = variable,
        high_confidence = high_conf,
        check.names = FALSE
    )
    variants_df <- variants_df[order(-variants_df$high_confidence, -variants_df$variable, -variants_df$max_af), ]
    write.table(subset(variants_df, high_confidence, select = -high_confidence), "mgatk_high_confidence_variants.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

    summary_df <- data.frame(
        metric = c("total_variants", "high_confidence_variants", "variable_variants", "min_af", "min_vmr", "min_sd"),
        value = c(nrow(variants_df), sum(high_conf), sum(high_conf & variable), min_af, min_vmr, min_sd)
    )
    write.table(summary_df, "mgatk_variant_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

    plot_mat <- af_matrix[high_conf & variable, , drop = FALSE]
    if (nrow(plot_mat) > 1 && ncol(plot_mat) > 1) {
        pdf("mgatk_heteroplasmy_heatmap.pdf", width = 8, height = 8)
        heatmap(plot_mat, scale = "none", Colv = NA, margins = c(8, 8), main = "Mitochondrial heteroplasmy")
        dev.off()

        pca <- prcomp(t(plot_mat), center = TRUE, scale. = TRUE)
        pdf("mgatk_heteroplasmy_pca.pdf", width = 7, height = 5)
        plot(pca$x[, 1], pca$x[, 2], pch = 19, xlab = "PC1", ylab = "PC2", main = "PCA on mtDNA heteroplasmy")
        text(pca$x[, 1], pca$x[, 2], labels = rownames(pca$x), pos = 3, cex = 0.7)
        dev.off()
    }
}

summary_df <- read.delim("mgatk_variant_summary.tsv")
cat("# id: atacportray_mgatk_postprocess\n# section_name: mgatk heteroplasmy summary\n# description: High-confidence mitochondrial heteroplasmy calls after allele-fraction, VMR and variability filtering.\n# plot_type: table\n", file = "mgatk_postprocess_mqc.tsv")
cat(paste(names(summary_df), collapse = "\t"), "\n", file = "mgatk_postprocess_mqc.tsv", append = TRUE, sep = "")
write.table(summary_df, "mgatk_postprocess_mqc.tsv", sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE, col.names = FALSE)

coverage_list <- list()
for (path in coverage_files) {
    dat <- read_table(path, header = FALSE)
    if (is.null(dat) || !nrow(dat)) next

    first_row <- tolower(as.character(unlist(dat[1, , drop = TRUE])))
    has_header <- any(grepl("pos|position|site|sample|cell|barcode|coverage|depth|cov", first_row))
    if (has_header) {
        names(dat) <- make.names(as.character(unlist(dat[1, , drop = TRUE])), unique = TRUE)
        dat <- dat[-1, , drop = FALSE]
    }

    sample <- clean_sample(path, "\\.coverage\\.txt\\.gz$")
    sample_vec <- rep(clean_sample_value(sample), nrow(dat))

    if (ncol(dat) >= 3 && !has_header) {
        pos <- dat[[1]]
        sample_vec <- clean_sample_value(dat[[2]])
        cov <- dat[[3]]
    } else if (ncol(dat) == 1) {
        cov <- dat[[1]]
        pos <- seq_along(cov)
    } else {
        pos_col <- grep("pos|position|site", names(dat), ignore.case = TRUE, value = TRUE)[1]
        sample_col <- grep("sample|cell|barcode", names(dat), ignore.case = TRUE, value = TRUE)[1]
        cov_col <- grep("coverage|depth|cov", names(dat), ignore.case = TRUE, value = TRUE)[1]
        pos <- if (!is.na(pos_col)) dat[[pos_col]] else seq_len(nrow(dat))
        if (!is.na(sample_col)) sample_vec <- clean_sample_value(dat[[sample_col]])
        cov <- if (!is.na(cov_col)) dat[[cov_col]] else dat[[ncol(dat)]]
    }

    parsed <- data.frame(
        sample = sample_vec,
        position = suppressWarnings(as.integer(pos)),
        coverage = suppressWarnings(as.numeric(cov))
    )
    parsed <- parsed[!is.na(parsed$position) & !is.na(parsed$coverage), , drop = FALSE]
    if (nrow(parsed)) coverage_list[[path]] <- parsed
}

coverage_df <- if (length(coverage_list)) do.call(rbind, coverage_list) else data.frame(sample = character(), position = integer(), coverage = numeric())
write.table(coverage_df, "mgatk_mt_coverage.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

cat("# id: atacportray_mgatk_mt_coverage\n# section_name: samtools mitochondrial coverage\n# description: Per-position mitochondrial coverage computed from BAM files with samtools depth, normalized to each sample mean coverage.\n# plot_type: linegraph\n# pconfig:\n#   id: atacportray_mgatk_mt_coverage\n#   title: samtools mitochondrial coverage\n#   xlab: MT position\n#   ylab: Coverage / sample mean coverage\n", file = "mgatk_mt_coverage_mqc.tsv")
cat('id: "atacportray_mgatk_mt_coverage_circular"\nsection_name: "Circular mitochondrial coverage"\ndescription: "Circular representation of per-position mitochondrial coverage computed from BAM files with samtools depth, normalized to each sample mean coverage."\nplot_type: "html"\ndata: |\n', file = "mgatk_mt_coverage_circular_mqc.yaml")
if (nrow(coverage_df)) {
    samples <- unique(coverage_df$sample)
    positions <- sort(unique(coverage_df$position))
    wide <- matrix(NA_real_, nrow = length(samples), ncol = length(positions), dimnames = list(samples, positions))
    for (i in seq_len(nrow(coverage_df))) {
        wide[coverage_df$sample[i], as.character(coverage_df$position[i])] <- coverage_df$coverage[i]
    }
    norm_wide <- wide
    for (sample in rownames(norm_wide)) {
        sample_mean <- mean(norm_wide[sample, ], na.rm = TRUE)
        if (is.finite(sample_mean) && sample_mean > 0) {
            norm_wide[sample, ] <- norm_wide[sample, ] / sample_mean
        }
    }
    write.table(norm_wide, "mgatk_mt_coverage_mqc.tsv", sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE, na = "")

    pdf("mgatk_mt_coverage_circular.pdf", width = 7, height = 7)
    par(mar = c(1, 1, 3, 1))
    plot.new()
    plot.window(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), asp = 1)
    title("Circular mitochondrial coverage (sample-normalized)")
    theta <- 2 * pi * (positions - min(positions)) / max(1, diff(range(positions)))
    cols <- grDevices::rainbow(length(samples))
    max_cov <- stats::quantile(as.numeric(norm_wide), probs = 0.99, na.rm = TRUE)
    if (!is.finite(max_cov) || max_cov <= 0) max_cov <- max(norm_wide, na.rm = TRUE)
    for (i in seq_along(samples)) {
        radius <- 0.35 + 0.6 * (pmin(norm_wide[i, ], max_cov) / max(max_cov, 1))
        lines(radius * cos(theta), radius * sin(theta), col = cols[i], lwd = 1)
    }
    legend("bottom", legend = samples, col = cols, lwd = 2, bty = "n", cex = 0.7)
    dev.off()

    svg_width <- 680
    svg_height <- 680
    center <- svg_width / 2
    theta_svg <- 2 * pi * (positions - min(positions)) / max(1, diff(range(positions))) - pi / 2
    line_tags <- character()
    label_tags <- character()
    for (i in seq_along(samples)) {
        radius <- 120 + 190 * (pmin(norm_wide[i, ], max_cov) / max(max_cov, 1))
        x <- center + radius * cos(theta_svg)
        y <- center + radius * sin(theta_svg)
        valid <- is.finite(x) & is.finite(y)
        runs <- split(seq_along(valid), cumsum(c(TRUE, diff(valid) != 0)))
        sample_lines <- character()
        for (run in runs) {
            if (!valid[run[1]] || length(run) < 2) next
            points <- paste(sprintf("%.1f,%.1f", x[run], y[run]), collapse = " ")
            sample_lines <- c(
                sample_lines,
                sprintf('<polyline points="%s" fill="none" stroke="%s" stroke-width="2" stroke-opacity="0.85"/>', points, cols[i])
            )
        }
        line_tags <- c(
            line_tags,
            sample_lines
        )
        label_tags <- c(
            label_tags,
            sprintf('<span style="display:inline-flex;align-items:center;margin-right:1rem;"><span style="width:1rem;height:0.25rem;background:%s;display:inline-block;margin-right:0.35rem;"></span>%s</span>', cols[i], samples[i])
        )
    }
    guide_tags <- paste(sprintf('<circle cx="%d" cy="%d" r="%d" fill="none" stroke="#ddd" stroke-width="1"/>', center, center, c(120, 210, 300)), collapse = "\n")
    html <- paste0(
        '<div class="mqc-custom-content" style="max-width:760px;">',
        '<svg viewBox="0 0 680 680" role="img" aria-label="Circular mitochondrial coverage" style="width:100%;height:auto;">',
        '<rect width="680" height="680" fill="white"/>',
        guide_tags,
        '<circle cx="340" cy="340" r="4" fill="#555"/>',
        paste(line_tags, collapse = "\n"),
        '<text x="340" y="34" text-anchor="middle" font-size="20" font-family="sans-serif">Circular mitochondrial coverage (sample-normalized)</text>',
        '</svg>',
        '<div style="font-size:0.9rem;margin-top:0.5rem;">',
        paste(label_tags, collapse = ""),
        '</div>',
        '</div>'
    )
    cat(paste0("  ", gsub("\n", "\n  ", html), "\n"), file = "mgatk_mt_coverage_circular_mqc.yaml", append = TRUE)
} else {
    cat("  <p>No mitochondrial coverage data available.</p>\n", file = "mgatk_mt_coverage_circular_mqc.yaml", append = TRUE)
}
