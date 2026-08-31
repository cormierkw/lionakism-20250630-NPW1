#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "Usage: calculate_unnormalized_msstats_protein_abundance.R <msstats_objects.rds> <output.rds>",
    call. = FALSE
  )
}

input_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_file <- args[[2]]
objects <- readRDS(input_file)
if (!"msstats_input" %in% names(objects)) {
  stop("The input RDS does not contain msstats_input.", call. = FALSE)
}

processed_unnormalized <- MSstats::dataProcess(
  raw = objects$msstats_input,
  logTrans = 2,
  normalization = FALSE,
  summaryMethod = "TMP",
  MBimpute = FALSE,
  use_log_file = FALSE,
  verbose = FALSE,
  numberOfCores = 1
)

protein_level <- processed_unnormalized[["ProteinLevelData"]]
if (is.null(protein_level) || nrow(protein_level) == 0L) {
  stop("MSstats did not return unnormalized ProteinLevelData.", call. = FALSE)
}

saveRDS(protein_level, output_file)
cat(
  "Saved", nrow(protein_level),
  "unnormalized, non-imputed protein/run rows to", output_file, "\n"
)
