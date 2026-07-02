#!/usr/bin/env Rscript
# Merge per-sample MAFs and draw an oncoplot (maftools) for atacportray.
suppressMessages({ library(maftools); library(optparse) })

opt <- parse_args(OptionParser(option_list=list(
  make_option("--input_dir", type="character"),
  make_option("--out_prefix", type="character", default="atacportray"),
  make_option("--top", type="integer", default=20)
)))

maf_files <- list.files(opt$input_dir, pattern="\\.maf(\\.gz)?$", full.names=TRUE)
stopifnot(length(maf_files) > 0)

mafs <- lapply(maf_files, function(f) tryCatch(read.maf(maf=f, verbose=FALSE),
                                               error=function(e) NULL))
mafs <- Filter(Negate(is.null), mafs)
merged <- if (length(mafs) > 1) merge_mafs(mafs) else mafs[[1]]

# Persist merged MAF
write.table(merged@data, file=gzfile("merged.maf.gz"), sep="\t",
            quote=FALSE, row.names=FALSE)

# Oncoplot
pdf(paste0(opt$out_prefix, ".oncoplot.pdf"), width=10, height=8)
oncoplot(maf=merged, top=opt$top)
dev.off()

# Summary plot + table
pdf(paste0(opt$out_prefix, ".summary.pdf"), width=10, height=8)
plotmafSummary(maf=merged, addStat="median", dashboard=TRUE)
dev.off()

sink(paste0(opt$out_prefix, ".maf_summary.txt"))
print(getSampleSummary(merged)); cat("\n"); print(getGeneSummary(merged))
sink()

cat("maftools oncoplot complete\n")
