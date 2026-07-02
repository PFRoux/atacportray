#!/usr/bin/env Rscript
# QDNAseq copy-number workflow for atacportray.
suppressMessages({
  library(QDNAseq)
  library(optparse)
})

opt_list <- list(
  make_option("--bam", type="character"),
  make_option("--sample", type="character", default="sample"),
  make_option("--binsize", type="double", default=100),
  make_option("--bins_rds", type="character", default=NULL),
  make_option("--loss", type="double", default=-0.4),
  make_option("--gain", type="double", default=0.4),
  make_option("--genome", type="character", default="hg19")
)
opt <- parse_args(OptionParser(option_list=opt_list))

bins <- if (!is.null(opt$bins_rds) && nzchar(opt$bins_rds)) {
  readRDS(opt$bins_rds)
} else {
  getBinAnnotations(binSize = opt$binsize, genome = opt$genome)
}

readCounts <- binReadCounts(bins, bamfiles = opt$bam)
readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun = "sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

copyNumbersCalled <- callBins(
  copyNumbersSegmented,
  method = "cutoff",
  cutoffs = c(loss = opt$loss, gain = opt$gain)
)

prefix <- opt$sample
exportBins(copyNumbersCalled, file = paste0(prefix, ".calls.txt"), type = "calls")
exportBins(copyNumbersCalled, file = paste0(prefix, ".segments.txt"), type = "segments")
exportBins(copyNumbersCalled, file = paste0(prefix, ".copynumbers.txt"), type = "copynumber")
exportBins(copyNumbersCalled, file = paste0(prefix, ".igv"), format = "igv", type = "copynumber")
saveRDS(copyNumbersCalled, file = paste0(prefix, ".rds"))

png(paste0(prefix, ".profile.png"), width = 1600, height = 600)
plot(copyNumbersCalled)
dev.off()

cat("QDNAseq complete for", prefix, "\n")
