#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript download_gisaid.R <lineage> <output_file>")
}

lineage <- args[1]
outfile <- args[2]

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



# # install library outbreakinfo https://rdrr.io/github/outbreak-info/R-outbreak-info/
# # modified from https://github.com/AlexiaNomena/SC2_VASIL/blob/main/Hands_On.R

# library(outbreakinfo)

# # Define the lineage of interest
# lineages_of_interest <- c("JN.1", "JN.2", "JN.3", "KP.3", "XBB.1.5")
# lineages_file_name <- gsub("\\.", "_", lineages_of_interest)

# # Loop through each lineage and download mutation data
# # Note: You need to be logged in to GISAID for this to work

# for (i in seq_along(lineages_of_interest)) {
#     lineage <- lineages_of_interest[i]
#     file_name <- lineages_file_name[i]

#     print(paste("\nDownloading mutation data for lineage:", lineage, "\n"))

#     print("Please login to gisaid.org in browser to initialize the session in this script.")
#     outbreakinfo::authenticateUser()
#     Sys.sleep(5)  # Pause for 5 seconds to allow time for login

#     # Get mutations for the specified lineage with at least 75% frequency
#     mutations <- getMutationsByLineage(pangolin_lineage=lineage, frequency=0.75, logInfo = FALSE)
#     print(head(mutations))

#     # Save mutations as txt files with each mutations in a new line
#     lin_mut <- paste0(mutations$ref_aa,as.character(mutations$codon_num),mutations$alt_aa)
#     print(lin_mut)
#     outdir <- file.path("data", "input", "mut_profiles")
#     dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
#     outfile <- file.path(outdir, paste0(file_name, ".txt"))
#     writeLines(lin_mut, con = outfile)  # one mutation per line
#     print(paste("Saved mutations for lineage", lineage, "to", outfile))
#     Sys.sleep(5)  # Pause for 5 seconds between requests to avoid overwhelming the server

# }
# print("All downloads completed.")
