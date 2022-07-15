# file: git_03-plot_GC_correction.R
# author: Derek Wong, Ph.D
# date: August 4th, 2021
# edited by: Ming Han
# date: July 15th, 2022

tib.list <- tib.list %>% dplyr::select(-matches("X"))

## Plot GC Correction metrics
pdf(file = file.path(outdir, paste0(id, "_GC_metrics.pdf")))
par(mfrow=c(2,2))
## short
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$short,
              main = "short",
              xlab = "frag_GC",
              ylab = "short")
## short corrected
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$short.corrected,
              main = "short corrected",
              xlab = "frag_GC",
              ylab = "short_corrected")
## short vs short predicted
smoothScatter(x = tib.list$short.predicted,
              y = tib.list$short,
              main = "short predicted vs actual",
              xlab = "short_predicted",
              ylab = "short")
## short corrected vs short predicted
smoothScatter(x = tib.list$short.predicted,
              y = tib.list$short.corrected,
              main = "short predicted vs corrected",
              xlab = "short_predicted",
              ylab = "short_corrected")
## long
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$long,
              main = "long",
              xlab = "frag_GC",
              ylab = "long")
## long corrected
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$long.corrected,
              main = "corrected long",
              xlab = "frag_GC",
              ylab = "long_corrected")
## long vs long predicted
smoothScatter(x = tib.list$long.predicted,
              y = tib.list$long,
              main = "long predicted vs actual",
              xlab = "long_predicted",
              ylab = "long")
## long corrected vs long predicted
smoothScatter(x = tib.list$long.predicted,
              y = tib.list$long.corrected,
              main = "long predicted vs corrected",
              xlab = "long_predicted",
              ylab = "long_corrected")
## total fragments
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$nfrags,
              main = "nfrags",
              xlab = "frag_GC",
              ylab = "nfrags")
## corrected total fragments
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$nfrags.corrected,
              main = "corrected nfrags",
              xlab = "frag_GC",
              ylab = "nfrags_corrected")
## fragments vs predicted fragments
smoothScatter(x = tib.list$nfrags.predicted,
              y = tib.list$nfrags,
              main = "nfrags predicted vs actual",
              xlab = "nfrags_predicted",
              ylab = "nfrags")
## corrected fragments vs predicted fragments
smoothScatter(x = tib.list$nfrags.predicted,
              y = tib.list$nfrags.corrected,
              main = "nfrags predicted vs corrected",
              xlab = "nfrags_predicted",
              ylab = "nfrags_corrected")
## ratios
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$ratio,
              main = "ratios",
              xlab = "frag_gc",
              ylab = "ratio")
## corrected ratios
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$ratio.corrected,
              main = "corrected ratios",
              xlab = "frag_gc",
              ylab = "ratio_corrected")
## predicted ratios
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$ratio.predicted,
              main = "predicted ratios",
              xlab = "frag_gc",
              ylab = "ratio_predicted")
## bin GC vs frag GC
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$C.G,
              main = "GC content",
              xlab = "frag_GC",
              ylab = "bin_GC")
## coverage
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$coverage,
              main = "coverage",
              xlab = "frag_GC",
              ylab = "coverage")
## corrected coverage
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$coverage.corrected,
              main = "corrected coverage",
              xlab = "frag_GC",
              ylab = "coverage_corrected")
## predicted coverage
smoothScatter(x = tib.list$frag.gc,
              y = tib.list$coverage.predicted,
              main = "predicted coverage",
              xlab = "frag_GC",
              ylab = "coverage_predicted")
dev.off()

## EOF