as_plain_tibble <- function(x) {
  tibble::as_tibble(as.data.frame(x))
}

first_existing_col <- function(data, candidates) {
  found <- intersect(candidates, names(data))
  if (length(found) == 0) NULL else found[[1]]
}

is_truthy_flag <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  vals <- tolower(trimws(as.character(x)))
  !is.na(vals) & vals %in% c("true", "t", "yes", "y", "1", "+")
}

has_modification <- function(x) {
  vals <- tolower(trimws(as.character(x)))
  !is.na(vals) & nzchar(vals) &
    !vals %in% c("none", "unmodified", "no modifications", "na")
}

add_psm_feature_id <- function(data) {
  data |>
    mutate(
      .qc_feature_id = paste(
        .data$`Master Protein Accessions`,
        .data$`Annotated Sequence`,
        .data$Charge,
        sep = " | "
      )
    )
}

coerce_log2_intensity <- function(x, value_col) {
  vals <- suppressWarnings(as.numeric(x))
  finite_vals <- vals[is.finite(vals)]
  if (length(finite_vals) == 0) return(vals)

  looks_linear <- grepl("^Intensity$|Quan Value", value_col, ignore.case = TRUE) ||
    stats::median(finite_vals, na.rm = TRUE) > 100

  if (looks_linear) {
    log2(pmax(vals, .Machine$double.xmin))
  } else {
    vals
  }
}

impute_matrix_by_row_mean <- function(mat) {
  if (nrow(mat) == 0 || ncol(mat) == 0) return(mat)

  keep <- rowSums(is.finite(mat)) >= 2
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) == 0) return(mat)

  row_means <- rowMeans(mat, na.rm = TRUE)
  for (i in seq_len(nrow(mat))) {
    mat[i, !is.finite(mat[i, ])] <- row_means[[i]]
  }

  row_vars <- apply(mat, 1, stats::var, na.rm = TRUE)
  mat[is.finite(row_vars) & row_vars > 0, , drop = FALSE]
}

make_feature_run_matrix <- function(data, feature_col, run_col, value_col,
                                    max_features = params$qc_max_features) {
  matrix_data <- data |>
    mutate(
      .qc_feature = as.character(.data[[feature_col]]),
      .qc_run = as.character(.data[[run_col]]),
      .qc_log2 = coerce_log2_intensity(.data[[value_col]], value_col)
    ) |>
    filter(!is.na(.data$.qc_feature), nzchar(.data$.qc_feature),
           !is.na(.data$.qc_run), nzchar(.data$.qc_run),
           is.finite(.data$.qc_log2)) |>
    group_by(.data$.qc_feature, .data$.qc_run) |>
    summarise(log2Intensity = mean(.data$.qc_log2, na.rm = TRUE),
              .groups = "drop")

  if (nrow(matrix_data) == 0) return(matrix(numeric(), nrow = 0, ncol = 0))

  selected_features <- matrix_data |>
    group_by(.data$.qc_feature) |>
    summarise(
      n_runs = n_distinct(.data$.qc_run),
      variance = stats::var(.data$log2Intensity, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(variance = if_else(is.finite(.data$variance), .data$variance, 0)) |>
    arrange(desc(.data$n_runs), desc(.data$variance)) |>
    slice_head(n = max_features) |>
    pull(.data$.qc_feature)

  wide <- matrix_data |>
    filter(.data$.qc_feature %in% selected_features) |>
    tidyr::pivot_wider(
      names_from = ".qc_run",
      values_from = "log2Intensity"
    )

  mat <- wide |>
    tibble::column_to_rownames(".qc_feature") |>
    as.matrix()
  storage.mode(mat) <- "numeric"
  mat
}

plot_correlation_heatmap <- function(feature_matrix, sample_annotation, title) {
  if (ncol(feature_matrix) < 2 || nrow(feature_matrix) < 2) {
    message("Correlation heatmap skipped: not enough quantified features or runs.")
    return(invisible(NULL))
  }

  mat <- impute_matrix_by_row_mean(feature_matrix)
  if (ncol(mat) < 2 || nrow(mat) < 2) {
    message("Correlation heatmap skipped: not enough complete features after imputation.")
    return(invisible(NULL))
  }

  corr <- stats::cor(mat, use = "pairwise.complete.obs")
  corr_df <- as.data.frame(as.table(corr), stringsAsFactors = FALSE)
  names(corr_df) <- c("Run1", "Run2", "Correlation")

  correlation_annotation <- sample_annotation |>
    mutate(
      CorrelationLabel = if_else(
        !is.na(.data$Condition) & nzchar(.data$Condition) &
          !is.na(.data$BioReplicate) & nzchar(as.character(.data$BioReplicate)),
        paste0(.data$Condition, "_", .data$BioReplicate),
        coalesce(.data$FigureLabel, .data$Run)
      )
    )
  label_lookup <- stats::setNames(
    correlation_annotation$CorrelationLabel,
    correlation_annotation$Run
  )
  display_levels <- unname(label_lookup[colnames(mat)])
  display_levels[is.na(display_levels)] <- colnames(mat)[is.na(display_levels)]
  corr_df <- corr_df |>
    mutate(
      Sample1 = coalesce(unname(label_lookup[as.character(.data$Run1)]), as.character(.data$Run1)),
      Sample2 = coalesce(unname(label_lookup[as.character(.data$Run2)]), as.character(.data$Run2)),
      Sample1 = factor(.data$Sample1, levels = unique(display_levels)),
      Sample2 = factor(.data$Sample2, levels = unique(display_levels))
    )

  ggplot(corr_df, aes(x = Sample1, y = Sample2, fill = Correlation)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(
      low = "#2C7BB6", mid = "white", high = "#D7191C",
      midpoint = 0.9, limits = c(min(corr_df$Correlation, na.rm = TRUE), 1)
    ) +
    labs(x = NULL, y = NULL, title = title, fill = "Pearson r") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

plot_pca_from_matrix <- function(feature_matrix, sample_annotation, title) {
  if (ncol(feature_matrix) < 3 || nrow(feature_matrix) < 2) {
    message("PCA skipped: at least 3 runs and 2 quantified features are required.")
    return(invisible(NULL))
  }

  mat <- impute_matrix_by_row_mean(feature_matrix)
  if (ncol(mat) < 3 || nrow(mat) < 2) {
    message("PCA skipped: not enough complete, variable features after imputation.")
    return(invisible(NULL))
  }

  pca <- stats::prcomp(t(mat), center = TRUE, scale. = FALSE)
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)

  pca_df <- as.data.frame(pca$x[, seq_len(min(2, ncol(pca$x))), drop = FALSE]) |>
    tibble::rownames_to_column("Run") |>
    left_join(sample_annotation, by = "Run") |>
    mutate(
      PcaLabel = if_else(
        !is.na(.data$BioReplicate) & nzchar(as.character(.data$BioReplicate)),
        paste0("BR", .data$BioReplicate),
        coalesce(.data$FigureLabel, .data$Run)
      )
    )

  ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = PcaLabel)) +
    geom_point(size = 3, alpha = 0.9) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20, show.legend = FALSE) +
    labs(
      title = title,
      x = paste0("PC1 (", scales::percent(var_explained[[1]], accuracy = 0.1), ")"),
      y = paste0("PC2 (", scales::percent(var_explained[[2]], accuracy = 0.1), ")")
    ) +
    theme_bw(base_size = 12)
}

make_protein_level_matrix <- function(processed_data, max_features = params$qc_max_features) {
  protein_level <- as_plain_tibble(processed_data[["ProteinLevelData"]])
  protein_col <- first_existing_col(protein_level, c("Protein", "PROTEIN", "ProteinName"))
  run_col <- first_existing_col(protein_level, c("Run", "RUN", "originalRUN", "SUBJECT", "BioReplicate"))
  value_col <- first_existing_col(protein_level, c("LogIntensities", "log2Intensity", "ABUNDANCE", "Intensity"))

  if (is.null(protein_col) || is.null(run_col) || is.null(value_col)) {
    message("Protein-level matrix skipped: required MSstats columns were not found.")
    return(matrix(numeric(), nrow = 0, ncol = 0))
  }

  make_feature_run_matrix(protein_level, protein_col, run_col, value_col, max_features)
}

extract_distribution_data <- function(data, stage_label) {
  tbl <- as_plain_tibble(data)
  run_col <- first_existing_col(tbl, c("Run", "RUN", "originalRUN", "SpectrumFile"))
  condition_col <- first_existing_col(tbl, c("Condition", "GROUP", "GROUP_ORIGINAL", "Condition_ORIGINAL"))
  value_col <- first_existing_col(tbl, c("LogIntensities", "log2Intensity", "ABUNDANCE", "Intensity", "Quan Value"))

  if (is.null(run_col) || is.null(value_col)) {
    return(tibble(Stage = character(), Run = character(),
                  Condition = character(), log2Intensity = numeric()))
  }

  tbl |>
    transmute(
      Stage = stage_label,
      Run = as.character(.data[[run_col]]),
      Condition = if (!is.null(condition_col)) as.character(.data[[condition_col]]) else NA_character_,
      log2Intensity = coerce_log2_intensity(.data[[value_col]], value_col)
    ) |>
    filter(is.finite(.data$log2Intensity))
}

extract_normalized_distribution_data <- function(
  data,
  normalization_method = "equalizeMedians",
  stage_label = "After normalization"
) {
  tbl <- as_plain_tibble(data)
  run_col <- first_existing_col(tbl, c("Run", "RUN", "originalRUN", "SpectrumFile"))
  condition_col <- first_existing_col(tbl, c("Condition", "GROUP", "GROUP_ORIGINAL", "Condition_ORIGINAL"))
  value_col <- first_existing_col(tbl, c("LogIntensities", "log2Intensity", "ABUNDANCE", "Intensity", "Quan Value"))
  fraction_col <- first_existing_col(tbl, c("Fraction", "FRACTION"))
  label_col <- first_existing_col(tbl, c("IsotopeLabelType", "LABEL"))

  empty_result <- tibble(
    Stage = character(), Run = character(),
    Condition = character(), log2Intensity = numeric()
  )
  if (is.null(run_col) || is.null(value_col)) return(empty_result)

  distribution <- tbl |>
    transmute(
      Run = as.character(.data[[run_col]]),
      Condition = if (!is.null(condition_col)) as.character(.data[[condition_col]]) else NA_character_,
      Fraction = if (!is.null(fraction_col)) as.character(.data[[fraction_col]]) else "1",
      Label = if (!is.null(label_col)) as.character(.data[[label_col]]) else "L",
      log2Intensity = coerce_log2_intensity(.data[[value_col]], value_col)
    ) |>
    filter(is.finite(.data$log2Intensity))

  method <- toupper(as.character(normalization_method))
  if (method %in% c("NONE", "FALSE")) {
    return(distribution |>
      transmute(Stage = stage_label, .data$Run, .data$Condition, .data$log2Intensity))
  }
  if (method != "EQUALIZEMEDIANS") {
    stop(
      "Normalization diagnostic currently supports equalizeMedians or no normalization; requested: ",
      normalization_method,
      call. = FALSE
    )
  }

  label_values <- unique(distribution$Label[!is.na(distribution$Label)])
  reference_label <- if (length(label_values) == 1L) label_values[[1]] else "H"

  run_medians <- distribution |>
    group_by(.data$Run, .data$Fraction) |>
    summarise(
      RunMedian = median(
        .data$log2Intensity[.data$Label == reference_label],
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  distribution |>
    left_join(run_medians, by = c("Run", "Fraction")) |>
    group_by(.data$Fraction) |>
    mutate(FractionMedian = median(.data$RunMedian, na.rm = TRUE)) |>
    ungroup() |>
    transmute(
      Stage = stage_label,
      .data$Run,
      .data$Condition,
      log2Intensity = .data$log2Intensity - .data$RunMedian + .data$FractionMedian
    ) |>
    filter(is.finite(.data$log2Intensity))
}

summarize_flag_columns <- function(data, patterns) {
  tbl <- as_plain_tibble(data)
  flag_cols <- names(tbl)[stringr::str_detect(tolower(names(tbl)), patterns)]
  if (length(flag_cols) == 0) {
    return(tibble(column = character(), value = character(), n = integer()))
  }

  purrr::map_dfr(flag_cols, function(col) {
    tbl |>
      mutate(.value = as.character(.data[[col]]),
             .value = if_else(is.na(.data$.value) | !nzchar(.data$.value),
                              "Missing/blank", .data$.value)) |>
      count(value = .data$.value, name = "n") |>
      mutate(column = col, .before = 1) |>
      arrange(.data$column, desc(.data$n))
  })
}

two_sided_t_power <- function(ncp, df, alpha = 0.05) {
  ncp <- abs(suppressWarnings(as.numeric(ncp)))
  df <- suppressWarnings(as.numeric(df))
  alpha <- suppressWarnings(as.numeric(alpha))

  if (!is.finite(ncp) || !is.finite(df) || df <= 0 ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    return(NA_real_)
  }

  critical_t <- stats::qt(1 - alpha / 2, df = df)
  stats::pt(-critical_t, df = df, ncp = ncp) +
    stats::pt(critical_t, df = df, ncp = ncp, lower.tail = FALSE)
}

minimum_detectable_ncp <- function(df, alpha = 0.05, target_power = 0.8) {
  df <- suppressWarnings(as.numeric(df))
  alpha <- suppressWarnings(as.numeric(alpha))
  target_power <- suppressWarnings(as.numeric(target_power))

  if (!is.finite(df) || df <= 0 ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      !is.finite(target_power) || target_power <= 0 || target_power >= 1) {
    return(NA_real_)
  }

  base_power <- two_sided_t_power(0, df = df, alpha = alpha)
  if (is.finite(base_power) && target_power <= base_power) return(0)

  upper <- 1
  upper_power <- two_sided_t_power(upper, df = df, alpha = alpha)
  while (is.finite(upper_power) && upper_power < target_power && upper < 1e4) {
    upper <- upper * 2
    upper_power <- two_sided_t_power(upper, df = df, alpha = alpha)
  }
  if (!is.finite(upper_power) || upper_power < target_power) return(NA_real_)

  stats::uniroot(
    function(x) two_sided_t_power(x, df = df, alpha = alpha) - target_power,
    interval = c(0, upper)
  )$root
}

add_power_metrics <- function(results, alpha = 0.05, target_power = 0.8,
                              effect_log2fc = 1) {
  required_cols <- c("Label", "log2FC", "SE", "DF", "pvalue", "adj.pvalue")
  missing_cols <- setdiff(required_cols, names(results))
  if (length(missing_cols) > 0) {
    stop("Power diagnostics require column(s): ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  power_results <- as_plain_tibble(results) |>
    mutate(
      .power_log2fc = suppressWarnings(as.numeric(.data$log2FC)),
      .power_se = suppressWarnings(as.numeric(.data$SE)),
      .power_df = suppressWarnings(as.numeric(.data$DF)),
      .power_pvalue = suppressWarnings(as.numeric(.data$pvalue)),
      .power_adj_pvalue = suppressWarnings(as.numeric(.data$adj.pvalue)),
      ncp_observed = if_else(
        is.finite(.data$.power_log2fc) &
          is.finite(.data$.power_se) & .data$.power_se > 0,
        abs(.data$.power_log2fc) / .data$.power_se,
        NA_real_
      ),
      ncp_target_effect = if_else(
        is.finite(.data$.power_se) & .data$.power_se > 0,
        abs(effect_log2fc) / .data$.power_se,
        NA_real_
      ),
      observed_effect_power = purrr::map2_dbl(
        .data$ncp_observed,
        .data$.power_df,
        two_sided_t_power,
        alpha = alpha
      ),
      target_log2fc_power = purrr::map2_dbl(
        .data$ncp_target_effect,
        .data$.power_df,
        two_sided_t_power,
        alpha = alpha
      )
    )

  df_power_lookup <- power_results |>
    distinct(.data$.power_df) |>
    filter(is.finite(.data$.power_df), .data$.power_df > 0) |>
    mutate(
      .min_ncp_for_target_power = purrr::map_dbl(
        .data$.power_df,
        minimum_detectable_ncp,
        alpha = alpha,
        target_power = target_power
      )
    )

  power_results |>
    left_join(df_power_lookup, by = ".power_df") |>
    mutate(
      detectable_log2fc_at_target_power = .data$.min_ncp_for_target_power * .data$.power_se,
      nominal_detected = .data$.power_pvalue < alpha,
      fdr_detected = .data$.power_adj_pvalue < alpha,
      target_effect_observed = abs(.data$.power_log2fc) >= abs(effect_log2fc),
      adequately_powered_for_target = .data$target_log2fc_power >= target_power,
      observed_effect_adequately_powered = .data$observed_effect_power >= target_power,
      power_detection_status = case_when(
        .data$fdr_detected & .data$target_effect_observed ~
          "FDR detected and at least target effect",
        .data$fdr_detected ~
          "FDR detected below target effect",
        !.data$fdr_detected & .data$adequately_powered_for_target ~
          "Not FDR detected, powered for target effect",
        !.data$fdr_detected & !.data$adequately_powered_for_target ~
          "Not FDR detected, underpowered for target effect",
        TRUE ~ "Insufficient information"
      )
    ) |>
    select(
      -any_of(c(
        ".power_log2fc",
        ".power_se",
        ".power_df",
        ".power_pvalue",
        ".power_adj_pvalue",
        ".min_ncp_for_target_power"
      ))
    )
}

percent_true <- function(x) {
  if (!any(!is.na(x))) return(NA_real_)
  100 * mean(x, na.rm = TRUE)
}

summarize_power_metrics <- function(power_data) {
  as_plain_tibble(power_data) |>
    group_by(.data$Label) |>
    summarise(
      proteins = n(),
      evaluable_power_rows = sum(is.finite(.data$target_log2fc_power)),
      median_power_for_target_log2fc = median(.data$target_log2fc_power, na.rm = TRUE),
      adequately_powered_for_target_percent = percent_true(.data$adequately_powered_for_target),
      median_detectable_log2fc_at_target_power = median(
        .data$detectable_log2fc_at_target_power,
        na.rm = TRUE
      ),
      fdr_detected = sum(.data$fdr_detected, na.rm = TRUE),
      fdr_detected_and_target_effect = sum(
        .data$fdr_detected & .data$target_effect_observed,
        na.rm = TRUE
      ),
      median_observed_effect_power = median(.data$observed_effect_power, na.rm = TRUE),
      observed_effect_adequately_powered_percent =
        percent_true(.data$observed_effect_adequately_powered),
      .groups = "drop"
    )
}
