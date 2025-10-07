#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript download_gisaid.R <lineage> <output_file>")
}

lineage <- args[1]
outfile <- args[2]

# install library outbreakinfo https://rdrr.io/github/outbreak-info/R-outbreak-info/
# modified from https://github.com/AlexiaNomena/SC2_VASIL/blob/main/Hands_On.R

library(outbreakinfo)

cat("Downloading mutation data for lineage:", lineage, "\n")

# Authenticate once per Snakemake job
outbreakinfo::authenticateUser()

mutations <- getMutationsByLineage(
  pangolin_lineage = lineage,
  frequency = 0.75,
  logInfo = FALSE
)

if (is.null(mutations) || nrow(mutations) == 0) {
  warning("No mutations found for ", lineage)
  quit(status = 0)
}

lin_mut <- paste0(mutations$ref_aa, mutations$codon_num, mutations$alt_aa)

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
writeLines(lin_mut, con = outfile)
cat("Saved mutations for", lineage, "to", outfile, "\n")

cat("Download completed.\n")
