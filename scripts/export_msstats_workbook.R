#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(
    "Usage: export_msstats_workbook.R <input.xlsx> <msstats_objects.rds> <results_dir> [output.xlsx]",
    call. = FALSE
  )
}

input_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
objects_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
results_dir <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
output_file <- if (length(args) >= 4L) args[[4]] else input_file
unnormalized_file <- file.path(results_dir, "msstats_unnormalized_protein_level.rds")

objects <- readRDS(objects_file)
required_objects <- c("msstats_input", "protein_level_data", "annotation", "comparison_results")
missing_objects <- setdiff(required_objects, names(objects))
if (length(missing_objects)) {
  stop("Missing exported object(s): ", paste(missing_objects, collapse = ", "), call. = FALSE)
}
if (!file.exists(unnormalized_file)) {
  stop(
    "Missing unnormalized protein results: ", unnormalized_file,
    ". Run calculate_unnormalized_msstats_protein_abundance.R first.",
    call. = FALSE
  )
}

clean_peptide <- function(x) {
  x <- sub("NA$", "", as.character(x))
  flanked <- grepl("^\\[[^]]+\\]\\..*\\.\\[[^]]+\\]$", x)
  x[flanked] <- sub("^\\[[^]]+\\]\\.(.*)\\.\\[[^]]+\\]$", "\\1", x[flanked])
  toupper(gsub("[^A-Za-z]", "", x))
}

annotation <- as.data.table(objects$annotation)
annotation[, Run := as.character(Run)]
annotation[, SampleName := as.character(SampleName)]
run_order <- annotation$Run
sample_names <- make.unique(annotation$SampleName)
names(sample_names) <- run_order
source_analysis_description <- if ("F22" %in% run_order) {
  "Full 30-run MSstats analysis in lionakism-20250630-NPW1_MSstats_PSM.qmd."
} else {
  paste0(
    "29-run MSstats analysis excluding F22 in ",
    "lionakism-20250630-NPW1_MSstats_PSM_exclude_F22.qmd."
  )
}

# PDtoMSstatsFormat has already resolved duplicate feature/run rows according to
# the report settings. Equalize-medians is an additive run shift on log2 values.
peptides <- as.data.table(objects$msstats_input)
peptides <- peptides[is.finite(Intensity) & Intensity > 0]
peptides[, `:=`(
  Protein = as.character(ProteinName),
  Peptide = clean_peptide(PeptideModifiedSequence),
  Run = as.character(Run),
  Log2Intensity = log2(Intensity)
)]
run_medians <- peptides[, .(RunMedian = median(Log2Intensity, na.rm = TRUE)), by = Run]
global_median <- median(run_medians$RunMedian, na.rm = TRUE)
peptides[run_medians, on = "Run", RunMedian := i.RunMedian]
peptides[, NormalizedIntensity := Intensity * 2^(global_median - RunMedian)]

peptide_summary <- peptides[, .(
  CalculatedAbundance = sum(Intensity, na.rm = TRUE),
  NormalizedAbundance = sum(NormalizedIntensity, na.rm = TRUE)
), by = .(Protein, Peptide, Run)]
peptide_summary[, Key := paste(Protein, Peptide, sep = "\r")]

make_wide_matrix <- function(dt, value_col, keys, runs) {
  wide <- dcast(dt, Key ~ Run, value.var = value_col)
  missing_runs <- setdiff(runs, names(wide))
  for (run in missing_runs) wide[, (run) := NA_real_]
  wide <- wide[, c("Key", runs), with = FALSE]
  idx <- match(keys, wide$Key)
  as.matrix(wide[idx, ..runs])
}

wb <- loadWorkbook(input_file)
template_sheet_order <- names(wb)
main_sheet <- if ("All_conditions Filtered" %in% names(wb)) {
  "All_conditions Filtered"
} else {
  names(wb)[[1]]
}
meta <- readWorkbook(
  wb, main_sheet, cols = 1:4, colNames = FALSE,
  skipEmptyRows = FALSE, skipEmptyCols = FALSE
)
n_rows <- nrow(meta)

protein_rows <- which(meta[[1]] == "Master Protein")
peptide_rows <- which(
  is.na(meta[[1]]) & !is.na(meta[[2]]) & meta[[2]] != "Confidence" &
    !is.na(meta[[4]]) & nzchar(as.character(meta[[4]]))
)
parent_index <- findInterval(peptide_rows, protein_rows)
parent_accession <- as.character(meta[[3]][protein_rows[parent_index]])
peptide_keys <- paste(parent_accession, clean_peptide(meta[[4]][peptide_rows]), sep = "\r")
protein_accessions <- as.character(meta[[3]][protein_rows])

calculated_matrix <- make_wide_matrix(
  peptide_summary, "CalculatedAbundance", peptide_keys, run_order
)
normalized_peptide_matrix <- make_wide_matrix(
  peptide_summary, "NormalizedAbundance", peptide_keys, run_order
)

protein_level <- as.data.table(objects$protein_level_data)
protein_level[, `:=`(Protein = as.character(Protein), originalRUN = as.character(originalRUN))]
protein_wide <- dcast(
  protein_level, Protein ~ originalRUN, value.var = "LogIntensities",
  fun.aggregate = function(x) if (length(x)) mean(x, na.rm = TRUE) else NA_real_
)
for (run in setdiff(run_order, names(protein_wide))) protein_wide[, (run) := NA_real_]
protein_idx <- match(protein_accessions, protein_wide$Protein)
protein_sample_matrix <- as.matrix(protein_wide[protein_idx, ..run_order])

unnormalized_protein_level <- as.data.table(readRDS(unnormalized_file))
unnormalized_protein_level[, `:=`(
  Protein = as.character(Protein),
  originalRUN = as.character(originalRUN)
)]
unnormalized_wide <- dcast(
  unnormalized_protein_level, Protein ~ originalRUN,
  value.var = "LogIntensities",
  fun.aggregate = function(x) if (length(x)) mean(x, na.rm = TRUE) else NA_real_
)
for (run in setdiff(run_order, names(unnormalized_wide))) {
  unnormalized_wide[, (run) := NA_real_]
}
unnormalized_idx <- match(protein_accessions, unnormalized_wide$Protein)
unnormalized_protein_matrix <- as.matrix(
  unnormalized_wide[unnormalized_idx, ..run_order]
)

group_order <- c("Healthy", "FC", "LC")
group_summary <- protein_level[, .(
  GroupAbundanceLog2 = mean(LogIntensities, na.rm = TRUE)
), by = .(Protein, GROUP)]
group_wide <- dcast(group_summary, Protein ~ GROUP, value.var = "GroupAbundanceLog2")
for (group in setdiff(group_order, names(group_wide))) group_wide[, (group) := NA_real_]
group_idx <- match(protein_accessions, group_wide$Protein)
group_matrix <- as.matrix(group_wide[group_idx, ..group_order])

comparison_order <- c(
  "FC_vs_Healthy", "LC_vs_Healthy", "LC_vs_FC", "Infected_vs_Healthy"
)
stat_fields <- c(
  "Ratio", "log2FC", "SE", "Tvalue", "DF", "pvalue", "adj.pvalue",
  "MissingPercentage", "ImputationPercentage", "issue"
)

protein_headers <- c(
  paste0("MSstats Protein Unnormalized TMP Abundance (log2; no imputation): ", sample_names),
  paste0("MSstats Protein Normalized Abundance (log2): ", sample_names),
  paste0("MSstats Group Abundance (mean log2): ", group_order),
  unlist(lapply(comparison_order, function(label) {
    paste0("MSstats ", label, ": ", stat_fields)
  }), use.names = FALSE)
)
peptide_headers <- c(
  paste0("MSstats Peptide Calculated Abundance (linear): ", sample_names),
  paste0("MSstats Peptide Normalized Abundance (linear): ", sample_names),
  rep(NA_character_, length(group_order) + length(comparison_order) * length(stat_fields))
)

n_new_cols <- length(protein_headers)
new_data <- setNames(
  lapply(seq_len(n_new_cols), function(i) rep(NA_real_, n_rows)),
  protein_headers
)
issue_columns <- grep(": issue$", protein_headers)
for (i in issue_columns) new_data[[i]] <- rep(NA_character_, n_rows)

calc_cols <- seq_along(run_order)
norm_cols <- max(calc_cols) + seq_along(run_order)
group_cols <- max(norm_cols) + seq_along(group_order)
new_data[calc_cols] <- lapply(seq_along(calc_cols), function(j) {
  x <- rep(NA_real_, n_rows)
  x[peptide_rows] <- calculated_matrix[, j]
  x[protein_rows] <- unnormalized_protein_matrix[, j]
  x
})
new_data[norm_cols] <- lapply(seq_along(norm_cols), function(j) {
  x <- rep(NA_real_, n_rows)
  x[peptide_rows] <- normalized_peptide_matrix[, j]
  x[protein_rows] <- protein_sample_matrix[, j]
  x
})
new_data[group_cols] <- lapply(seq_along(group_cols), function(j) {
  x <- rep(NA_real_, n_rows); x[protein_rows] <- group_matrix[, j]; x
})

comparisons <- as.data.table(objects$comparison_results)
comparisons[, Ratio := 2^log2FC]
stats_start <- max(group_cols) + 1L
for (comparison_i in seq_along(comparison_order)) {
  label <- comparison_order[[comparison_i]]
  result <- comparisons[Label == label]
  result_idx <- match(protein_accessions, result$Protein)
  col_start <- stats_start + (comparison_i - 1L) * length(stat_fields)
  for (field_i in seq_along(stat_fields)) {
    field <- stat_fields[[field_i]]
    out_col <- col_start + field_i - 1L
    if (field == "issue") {
      x <- rep(NA_character_, n_rows)
      x[protein_rows] <- as.character(result[[field]][result_idx])
    } else {
      x <- rep(NA_real_, n_rows)
      x[protein_rows] <- as.numeric(result[[field]][result_idx])
    }
    new_data[[out_col]] <- x
  }
}
new_data <- as.data.frame(new_data, check.names = FALSE)

original_header <- readWorkbook(
  wb, main_sheet, rows = 1, colNames = FALSE,
  skipEmptyRows = FALSE, skipEmptyCols = FALSE
)
original_header_values <- as.character(unlist(original_header[1, ], use.names = FALSE))
existing_msstats_cols <- grep("^MSstats ", original_header_values)
start_col <- if (length(existing_msstats_cols)) {
  min(existing_msstats_cols)
} else {
  ncol(original_header) + 1L
}
if (length(existing_msstats_cols)) {
  # The full analysis has one more run than the exclude-F22 analysis. Clear the
  # complete pre-existing calculated block before writing so columns belonging
  # only to the 30-run export cannot remain at the right edge of the worksheet.
  deleteData(
    wb, main_sheet,
    cols = seq.int(min(existing_msstats_cols), max(existing_msstats_cols)),
    rows = seq_len(n_rows + 1L),
    gridExpand = TRUE
  )
}
writeData(
  wb, main_sheet, new_data[-1, , drop = FALSE], startCol = start_col,
  startRow = 2, colNames = FALSE, rowNames = FALSE, keepNA = FALSE
)
writeData(wb, main_sheet, t(protein_headers), startCol = start_col, startRow = 1, colNames = FALSE)
writeData(wb, main_sheet, t(peptide_headers), startCol = start_col, startRow = 3, colNames = FALSE)

protein_header_style <- createStyle(
  fontColour = "#FFFFFF", fgFill = "#27445C", textRotation = 90,
  halign = "center", valign = "bottom", wrapText = TRUE, border = "Bottom"
)
peptide_header_style <- createStyle(
  fontColour = "#000000", fgFill = "#D9EAF7", textRotation = 90,
  halign = "center", valign = "bottom", wrapText = TRUE, border = "Bottom"
)
addStyle(wb, main_sheet, protein_header_style, rows = 1, cols = start_col:(start_col + n_new_cols - 1), gridExpand = TRUE)
addStyle(wb, main_sheet, peptide_header_style, rows = 3, cols = start_col:(start_col + n_new_cols - 1), gridExpand = TRUE)
setColWidths(wb, main_sheet, cols = start_col:(start_col + n_new_cols - 1), widths = 13)

description <- as.character(meta[[4]][protein_rows])
gene_symbol <- sub(".* GN=([^ ]+).*", "\\1", description)
gene_symbol[gene_symbol == description] <- NA_character_
protein_annotation <- data.table(
  Protein = protein_accessions,
  GeneSymbol = gene_symbol,
  Description = description
)

header_style <- createStyle(
  fontColour = "#FFFFFF", fgFill = "#27445C", textDecoration = "bold",
  halign = "center", valign = "center", wrapText = TRUE
)
significant_fill <- createStyle(fgFill = "#E2F0D9")

replace_sheet <- function(sheet_name, data) {
  if (sheet_name %in% names(wb)) removeWorksheet(wb, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data, withFilter = nrow(data) > 0, headerStyle = header_style)
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(data)), widths = "auto")
  if (nrow(data) > 0) {
    addStyle(wb, sheet_name, significant_fill, rows = 2:(nrow(data) + 1), cols = 1:ncol(data), gridExpand = TRUE, stack = TRUE)
  }
}

for (label in comparison_order) {
  de <- merge(
    comparisons[
      Label == label & is.finite(adj.pvalue) & adj.pvalue < 0.05 &
        is.finite(log2FC) & abs(log2FC) >= 1
    ],
    protein_annotation,
    by = "Protein",
    all.x = TRUE,
    sort = FALSE
  )
  setcolorder(de, c("Protein", "GeneSymbol", "Description", setdiff(names(de), c("Protein", "GeneSymbol", "Description"))))
  de[, .abs_log2FC_sort := abs(log2FC)]
  setorder(de, adj.pvalue, -.abs_log2FC_sort)
  de[, .abs_log2FC_sort := NULL]
  replace_sheet(substr(paste("DE", gsub("_", " ", label)), 1, 31), as.data.frame(de))
}

gsea_files <- c(
  FC_vs_Healthy = "gsea_go_bp_fc_vs_healthy.csv",
  LC_vs_Healthy = "gsea_go_bp_lc_vs_healthy.csv",
  LC_vs_FC = "gsea_go_bp_lc_vs_fc.csv",
  Infected_vs_Healthy = "gsea_go_bp_infected_vs_healthy.csv"
)
for (label in names(gsea_files)) {
  path <- file.path(results_dir, gsea_files[[label]])
  gsea <- if (file.exists(path)) data.table::fread(path) else data.frame(Note = "No GSEA result file was produced.")
  replace_sheet(substr(paste("GSEA", gsub("_", " ", label)), 1, 31), as.data.frame(gsea))
}

guide <- data.frame(
  Section = c(
    "Source analysis", "Peptide calculated abundance", "Protein unnormalized abundance",
    "Peptide normalized abundance",
    "Protein normalized abundance", "Group abundance", "Comparison ratio",
    "Differential-expression sheets", "GSEA sheets"
  ),
  Definition = c(
    source_analysis_description,
    "Sum of positive PDtoMSstatsFormat feature intensities for a protein-peptide-run across precursor charge states; linear scale.",
    "MSstats ProteinLevelData LogIntensities from a second dataProcess call with normalization and model-based imputation disabled, using log2 transformation and TMP summarization; log2 scale. Values occupy protein rows directly above their peptide rows.",
    "Calculated peptide abundance after the report's equalizeMedians run adjustment; linear scale.",
    "MSstats ProteinLevelData LogIntensities after normalization, model-based imputation, and TMP summarization; log2 scale.",
    "Arithmetic mean of MSstats protein LogIntensities within Healthy, FC, or LC; log2 scale.",
    "2 raised to the MSstats log2FC. Remaining fields are copied from groupComparison, including uncertainty, p-values, missingness, imputation, and issue flags.",
    "Proteins with adjusted p-value < 0.05 and absolute log2FC >= 1, one worksheet per contrast.",
    "Complete GO Biological Process GSEA result table produced by the report, one worksheet per contrast."
  ),
  stringsAsFactors = FALSE
)
replace_sheet("MSstats Export Guide", guide)
setColWidths(wb, "MSstats Export Guide", cols = 1, widths = 30)
setColWidths(wb, "MSstats Export Guide", cols = 2, widths = 100)
setRowHeights(wb, "MSstats Export Guide", rows = 2:(nrow(guide) + 1), heights = 36)
addStyle(
  wb, "MSstats Export Guide", createStyle(wrapText = TRUE, valign = "top"),
  rows = 2:(nrow(guide) + 1), cols = 1:2, gridExpand = TRUE, stack = TRUE
)

# openxlsx intentionally resets the worksheet dimension when loading. Restore
# the source report's final row and expand only the final column for the newly
# appended fields so Excel's used range remains accurate.
main_sheet_index <- match(main_sheet, names(wb))
wb$worksheets[[main_sheet_index]]$dimension <- sprintf(
  "<dimension ref=\"A1:%s%d\"/>",
  int2col(start_col + n_new_cols - 1L),
  n_rows + 1L
)
current_sheet_names <- names(wb)
desired_sheet_order <- c(
  template_sheet_order[template_sheet_order %in% current_sheet_names],
  setdiff(current_sheet_names, template_sheet_order)
)
worksheetOrder(wb) <- match(desired_sheet_order, current_sheet_names)

saveWorkbook(wb, output_file, overwrite = TRUE)
cat("Saved:", normalizePath(output_file, winslash = "/", mustWork = TRUE), "\n")
cat("Protein rows:", length(protein_rows), "\n")
cat("Peptide rows:", length(peptide_rows), "\n")
cat("Appended columns:", n_new_cols, "\n")
