#!/usr/bin/env Rscript
# QDNAseq copy-number workflow for atacportray
# Reproduces the exact steps and thresholds from the source notebook:
#   binReadCounts -> applyFilters(residual, blacklist) -> estimateCorrection ->
#   correctBins -> normalizeBins -> smoothOutlierBins -> segmentBins(sqrt/Anscombe) ->
#   normalizeSegmentedBins -> callBins(thresholds log2 loss=-0.4 / gain=+0.4)

suppressMessages({
  library(QDNAseq)
})

parse_args <- function(args) {
  opts <- list(
    bam = NULL,
    sample = "sample",
    binsize = 100,
    bins_rds = NULL,
    loss = -0.4,
    gain = 0.4,
    genome = "hg19"
  )

  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (!key %in% names(opts)) {
      stop("Unknown option: ", args[[i]], call. = FALSE)
    }
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      stop("Missing value for option: ", args[[i]], call. = FALSE)
    }
    opts[[key]] <- args[[i + 1]]
    i <- i + 2
  }

  if (is.null(opts$bam)) {
    stop("Missing required option: --bam", call. = FALSE)
  }

  opts$binsize <- as.numeric(opts$binsize)
  opts$loss <- as.numeric(opts$loss)
  opts$gain <- as.numeric(opts$gain)
  opts
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))

# Bin annotations for the requested bin size (pre-packaged hg19 bins;
# QDNAseq calls relative to bins, liftOver of segments is handled downstream if needed)
bins <- if (!is.null(opt$bins_rds) && nzchar(opt$bins_rds)) {
  readRDS(opt$bins_rds)
} else {
  getBinAnnotations(binSize = opt$binsize, genome = opt$genome)
}

# 1. Count reads in bins
readCounts <- binReadCounts(bins, bamfiles = opt$bam)

# 2. Filter (residuals + blacklisted regions), as in the notebook
readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)

# 3. Estimate GC/mappability correction and correct
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)

# 4. Normalize + smooth outliers
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

# 5. Segment (sqrt / Anscombe transform) and normalize segments
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun = "sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

# 6. Call gains/losses with explicit log2 thresholds from the notebook
copyNumbersCalled <- callBins(
  copyNumbersSegmented,
  method = "cutoff",
  cutoffs = c(loss = opt$loss, gain = opt$gain)
)

# --- Exports ---
prefix <- opt$sample
exportBins(copyNumbersCalled, file = paste0(prefix, ".calls.txt"),   type = "calls")
exportBins(copyNumbersCalled, file = paste0(prefix, ".segments.txt"), type = "segments")
exportBins(copyNumbersCalled, file = paste0(prefix, ".copynumbers.txt"), type = "copynumber")
exportBins(copyNumbersCalled, file = paste0(prefix, ".igv"), format = "igv", type = "copynumber")
saveRDS(copyNumbersCalled, file = paste0(prefix, ".rds"))

png(paste0(prefix, ".profile.png"), width = 1600, height = 600)
plot(copyNumbersCalled)
dev.off()

cat("QDNAseq complete for", prefix, "\n")
