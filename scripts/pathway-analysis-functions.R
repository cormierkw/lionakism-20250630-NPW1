empty_gsea_table <- function() {
  tibble::tibble(
    ID = character(),
    Description = character(),
    setSize = integer(),
    enrichmentScore = numeric(),
    NES = numeric(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    rank = integer(),
    leading_edge = character(),
    core_enrichment = character()
  )
}

format_comparison_label <- function(label) {
  stringr::str_replace_all(label, "_", " ")
}

normalize_uniprot_accession <- function(accession) {
  accession <- stringr::str_trim(as.character(accession))
  pipe_fields <- stringr::str_split_fixed(accession, "\\|", 3)
  has_pipe_accession <- pipe_fields[, 2] != ""
  accession[has_pipe_accession] <- pipe_fields[has_pipe_accession, 2]
  stringr::str_remove(accession, "-[0-9]+$")
}

map_uniprot_to_gene_symbol <- function(accession) {
  accession <- as.character(accession)
  long_accessions <- tibble::tibble(
    .row_id = seq_along(accession),
    Protein = accession
  ) |>
    tidyr::separate_rows(Protein, sep = ";\\s*") |>
    dplyr::mutate(UNIPROT = normalize_uniprot_accession(.data$Protein)) |>
    dplyr::filter(!is.na(.data$UNIPROT), nzchar(.data$UNIPROT))

  symbols <- rep(NA_character_, length(accession))
  if (nrow(long_accessions) == 0L) return(symbols)

  keys <- unique(long_accessions$UNIPROT)
  symbol_map <- suppressMessages(
    AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = keys,
      column = "SYMBOL",
      keytype = "UNIPROT",
      multiVals = "first"
    )
  )

  mapped_symbols <- long_accessions |>
    dplyr::mutate(GeneSymbol = unname(symbol_map[.data$UNIPROT])) |>
    dplyr::filter(!is.na(.data$GeneSymbol), nzchar(.data$GeneSymbol)) |>
    dplyr::group_by(.data$.row_id) |>
    dplyr::summarise(
      GeneSymbol = paste(unique(.data$GeneSymbol), collapse = "/"),
      .groups = "drop"
    )

  symbols[mapped_symbols$.row_id] <- mapped_symbols$GeneSymbol
  symbols
}

prepare_gsea_gene_list <- function(results, comparison_label) {
  comparison_results <- results |>
    dplyr::filter(.data$Label == comparison_label) |>
    dplyr::transmute(
      Protein = as.character(.data$Protein),
      log2FC = suppressWarnings(as.numeric(.data$log2FC)),
      pvalue = suppressWarnings(as.numeric(.data$pvalue))
    ) |>
    dplyr::filter(
      !is.na(.data$Protein),
      nzchar(.data$Protein),
      is.finite(.data$log2FC),
      is.finite(.data$pvalue),
      .data$pvalue >= 0
    ) |>
    tidyr::separate_rows(Protein, sep = ";\\s*") |>
    dplyr::mutate(
      UNIPROT = normalize_uniprot_accession(.data$Protein),
      rank_score = sign(.data$log2FC) * -log10(pmax(.data$pvalue, .Machine$double.xmin))
    ) |>
    dplyr::filter(nzchar(.data$UNIPROT), is.finite(.data$rank_score))

  input_accessions <- unique(comparison_results$UNIPROT)
  if (length(input_accessions) == 0) {
    return(list(
      gene_list = numeric(),
      input_proteins = 0L,
      input_accessions = 0L,
      mapped_accessions = 0L,
      mapped_entrez = 0L
    ))
  }

  mapping <- tryCatch(
    suppressMessages(
      clusterProfiler::bitr(
        input_accessions,
        fromType = "UNIPROT",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db::org.Hs.eg.db
      )
    ),
    error = function(e) data.frame(UNIPROT = character(), ENTREZID = character())
  ) |>
    tibble::as_tibble() |>
    dplyr::distinct(.data$UNIPROT, .data$ENTREZID)

  ranked_mapping <- comparison_results |>
    dplyr::inner_join(mapping, by = "UNIPROT", relationship = "many-to-many") |>
    dplyr::group_by(.data$ENTREZID) |>
    dplyr::arrange(dplyr::desc(abs(.data$rank_score)), .by_group = TRUE) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(.data$rank_score))

  gene_list <- ranked_mapping$rank_score
  names(gene_list) <- ranked_mapping$ENTREZID

  list(
    gene_list = gene_list,
    input_proteins = dplyr::n_distinct(comparison_results$Protein),
    input_accessions = length(input_accessions),
    mapped_accessions = dplyr::n_distinct(ranked_mapping$UNIPROT),
    mapped_entrez = length(gene_list)
  )
}

run_go_bp_gsea <- function(results, comparison_label, min_gene_set_size = 100,
                           max_gene_set_size = 500, pvalue_cutoff = 0.05,
                           dot_show_categories = 10,
                           tree_show_categories = 20,
                           tree_label_wrap = 70) {
  prepared <- prepare_gsea_gene_list(results, comparison_label)
  display_label <- format_comparison_label(comparison_label)

  empty_result <- list(
    comparison = comparison_label,
    input_proteins = prepared$input_proteins,
    input_accessions = prepared$input_accessions,
    mapped_accessions = prepared$mapped_accessions,
    mapped_entrez = prepared$mapped_entrez,
    tested_terms = 0L,
    significant_terms = 0L,
    gsea = NULL,
    readable = NULL,
    dotplot = NULL,
    treeplot = NULL,
    dotplot_status = "Not generated",
    treeplot_status = "Not generated",
    table = empty_gsea_table(),
    status = "Not run"
  )

  if (length(prepared$gene_list) < min_gene_set_size) {
    empty_result$status <- paste0(
      "Skipped: only ", length(prepared$gene_list),
      " unique human Entrez IDs were mapped; at least ",
      min_gene_set_size, " are required."
    )
    return(empty_result)
  }

  gsea_fit <- tryCatch(
    suppressMessages(
      clusterProfiler::gseGO(
        geneList = prepared$gene_list,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        minGSSize = min_gene_set_size,
        maxGSSize = max_gene_set_size,
        pvalueCutoff = pvalue_cutoff,
        pAdjustMethod = "BH",
        verbose = FALSE
      )
    ),
    error = function(e) e
  )

  if (inherits(gsea_fit, "error")) {
    empty_result$status <- paste0("GSEA failed: ", conditionMessage(gsea_fit))
    return(empty_result)
  }

  gsea_table <- tibble::as_tibble(as.data.frame(gsea_fit))
  empty_result$gsea <- gsea_fit
  empty_result$tested_terms <- if ("geneSets" %in% methods::slotNames(gsea_fit)) {
    length(methods::slot(gsea_fit, "geneSets"))
  } else {
    nrow(gsea_table)
  }

  if (nrow(gsea_table) == 0) {
    empty_result$status <- paste0(
      "No GO Biological Process terms met p-value cutoff ", pvalue_cutoff, "."
    )
    return(empty_result)
  }

  readable_fit <- tryCatch(
    suppressMessages(
      clusterProfiler::setReadable(
        gsea_fit,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "ENTREZID"
      )
    ),
    error = function(e) gsea_fit
  )
  readable_table <- tibble::as_tibble(as.data.frame(readable_fit))

  dotplot_error <- NULL
  dot_plot <- tryCatch(
    enrichplot::dotplot(
      readable_fit,
      showCategory = dot_show_categories,
      split = ".sign",
      font.size = 9
    ) +
      ggplot2::facet_grid(. ~ .sign) +
      ggplot2::ggtitle(paste0(display_label, ": GSEA GO BP top pathways")) +
      ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)),
    error = function(e) {
      dotplot_error <<- conditionMessage(e)
      NULL
    }
  )

  treeplot_error <- NULL
  tree_plot <- tryCatch({
    term_similarity <- enrichplot::pairwise_termsim(readable_fit)
    enrichplot::treeplot(
      term_similarity,
      showCategory = tree_show_categories,
      label_format = tree_label_wrap,
      fontsize_tiplab = 3,
      fontsize_cladelab = 3,
      extend = 0.5,
      hilight = FALSE,
      hexpand = 0.2
    ) +
      ggplot2::ggtitle(paste0(display_label, ": GSEA GO BP tree plot"))
  }, error = function(e) {
    treeplot_error <<- conditionMessage(e)
    NULL
  })

  empty_result$readable <- readable_fit
  empty_result$dotplot <- dot_plot
  empty_result$treeplot <- tree_plot
  empty_result$dotplot_status <- if (is.null(dotplot_error)) {
    "Generated"
  } else {
    paste0("Unavailable: ", dotplot_error)
  }
  empty_result$treeplot_status <- if (is.null(treeplot_error)) {
    "Generated"
  } else {
    paste0("Unavailable: ", treeplot_error)
  }
  empty_result$table <- readable_table
  empty_result$significant_terms <- nrow(readable_table)
  empty_result$status <- paste0(nrow(readable_table), " enriched GO BP terms returned.")
  empty_result
}

summarize_gsea_runs <- function(gsea_results) {
  purrr::imap_dfr(gsea_results, function(result, comparison_label) {
    tibble::tibble(
      Comparison = format_comparison_label(comparison_label),
      input_proteins = result$input_proteins,
      unique_accessions = result$input_accessions,
      mapped_accessions = result$mapped_accessions,
      ranked_entrez_ids = result$mapped_entrez,
      tested_terms = result$tested_terms,
      returned_terms = result$significant_terms,
      dotplot_status = result$dotplot_status,
      treeplot_status = result$treeplot_status,
      status = result$status
    )
  })
}

render_analysis_conclusion <- function(comparison_results, gsea_table,
                                       fdr_cutoff = 0.05,
                                       log2fc_cutoff = 1,
                                       top_pathways = 3,
                                       power_summary = NULL,
                                       target_power = 0.8,
                                       power_effect_log2fc = 1) {
  comparison_figure_anchors <- c(
    FC_vs_Healthy = "#fig-volcano-fc-vs-healthy",
    LC_vs_Healthy = "#fig-volcano-lc-vs-healthy",
    LC_vs_FC = "#fig-volcano-lc-vs-fc",
    Infected_vs_Healthy = "#fig-infected-volcano"
  )

  linked_comparison_label <- function(label, anchor = NULL) {
    display_label <- format_comparison_label(label)
    if (is.null(anchor)) {
      anchor <- unname(comparison_figure_anchors[as.character(label)])
    }
    if (length(anchor) == 0 || is.na(anchor) || !nzchar(anchor)) {
      return(paste0("**", display_label, "**"))
    }
    paste0("**[", display_label, "](", anchor, ")**")
  }

  differential_summary <- comparison_results |>
    dplyr::mutate(
      .conclusion_log2fc = suppressWarnings(as.numeric(.data$log2FC)),
      .conclusion_fdr = suppressWarnings(as.numeric(.data$adj.pvalue)),
      .conclusion_valid = !is.na(.data$.conclusion_log2fc) &
        is.finite(.data$.conclusion_fdr),
      .conclusion_fdr_hit = .data$.conclusion_valid &
        .data$.conclusion_fdr < fdr_cutoff,
      .conclusion_large_hit = .data$.conclusion_fdr_hit &
        abs(.data$.conclusion_log2fc) >= abs(log2fc_cutoff)
    ) |>
    dplyr::group_by(.data$Label) |>
    dplyr::summarise(
      tested_proteins = sum(.data$.conclusion_valid),
      fdr_significant = sum(.data$.conclusion_fdr_hit),
      fdr_large_effect = sum(.data$.conclusion_large_hit),
      higher_abundance = sum(.data$.conclusion_large_hit & .data$.conclusion_log2fc > 0),
      lower_abundance = sum(.data$.conclusion_large_hit & .data$.conclusion_log2fc < 0),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$fdr_large_effect))

  cat("### Differential protein analysis\n\n")
  cat(
    "Differential proteins were summarized at FDR < ",
    format(fdr_cutoff, trim = TRUE),
    " and a large-effect threshold of |log2FC| >= ",
    format(abs(log2fc_cutoff), trim = TRUE),
    ".\n\n",
    sep = ""
  )

  if (nrow(differential_summary) == 0) {
    cat("No valid differential-comparison results were available.\n\n")
  } else {
    purrr::pwalk(differential_summary, function(Label, tested_proteins,
                                                fdr_significant,
                                                fdr_large_effect,
                                                higher_abundance,
                                                lower_abundance) {
      cat(
        "- ", linked_comparison_label(Label), ": ",
        format(tested_proteins, big.mark = ","), " proteins tested; ",
        format(fdr_significant, big.mark = ","), " met the FDR threshold, and ",
        format(fdr_large_effect, big.mark = ","),
        " also met the effect-size threshold (",
        format(higher_abundance, big.mark = ","), " higher and ",
        format(lower_abundance, big.mark = ","), " lower in the numerator group).\n",
        sep = ""
      )
    })
    cat("\n")

    if (nrow(differential_summary) > 1) {
      strongest <- differential_summary |>
        dplyr::slice_max(.data$fdr_large_effect, n = 1, with_ties = FALSE)
      weakest <- differential_summary |>
        dplyr::slice_min(.data$fdr_large_effect, n = 1, with_ties = FALSE)
      cat(
        "The largest large-effect differential signal occurred for **",
        "[", format_comparison_label(strongest$Label[[1]]), "](",
        unname(comparison_figure_anchors[strongest$Label[[1]]]), ")** (",
        format(strongest$fdr_large_effect[[1]], big.mark = ","),
        " proteins), whereas ", linked_comparison_label(weakest$Label[[1]]),
        " had the smallest (", format(weakest$fdr_large_effect[[1]], big.mark = ","),
        ").\n\n",
        sep = ""
      )
    }
  }

  if (!is.null(power_summary) && nrow(power_summary) > 0) {
    power_conclusion <- power_summary |>
      dplyr::mutate(
        median_power_for_target_log2fc = suppressWarnings(
          as.numeric(.data$median_power_for_target_log2fc)
        ),
        adequately_powered_for_target_percent = suppressWarnings(
          as.numeric(.data$adequately_powered_for_target_percent)
        ),
        median_detectable_log2fc_at_target_power = suppressWarnings(
          as.numeric(.data$median_detectable_log2fc_at_target_power)
        )
      ) |>
      dplyr::filter(
        !is.na(.data$Label),
        is.finite(.data$median_power_for_target_log2fc),
        is.finite(.data$adequately_powered_for_target_percent),
        is.finite(.data$median_detectable_log2fc_at_target_power)
      )

    cat("### Power analysis\n\n")
    cat(
      "[Power was evaluated](#fig-power-target-effect) for a target |log2FC| of ",
      format(abs(power_effect_log2fc), trim = TRUE),
      " (", sprintf("%.2f", 2^abs(power_effect_log2fc)),
      "-fold) at a target power of ",
      scales::percent(target_power, accuracy = 1), ".\n\n",
      sep = ""
    )

    if (nrow(power_conclusion) == 0) {
      cat("No finite comparison-level power summaries were available.\n\n")
    } else {
      purrr::pwalk(
        power_conclusion,
        function(Label, median_power_for_target_log2fc,
                 adequately_powered_for_target_percent,
                 median_detectable_log2fc_at_target_power, ...) {
          coverage_text <- dplyr::case_when(
            adequately_powered_for_target_percent >= 80 ~ "broad power coverage",
            adequately_powered_for_target_percent >= 50 ~ "mixed power coverage",
            TRUE ~ "limited power coverage"
          )
          cat(
            "- ", linked_comparison_label(Label, "#fig-power-target-effect"),
            ": median power was ",
            scales::percent(median_power_for_target_log2fc, accuracy = 1),
            "; ", sprintf("%.1f", adequately_powered_for_target_percent),
            "% of proteins met the power target, indicating ", coverage_text,
            ". The [median minimum detectable |log2FC|](#fig-minimum-detectable-effect) was ",
            sprintf("%.2f", median_detectable_log2fc_at_target_power), ".\n",
            sep = ""
          )
        }
      )
      cat("\n")

      limited_comparisons <- power_conclusion |>
        dplyr::filter(.data$adequately_powered_for_target_percent < 50) |>
        dplyr::pull(.data$Label) |>
        purrr::map_chr(format_comparison_label)
      if (length(limited_comparisons) > 0) {
        cat(
          "Non-significant results in ",
          paste0("**", limited_comparisons, "**", collapse = ", "),
          " should be interpreted cautiously because fewer than half of proteins were adequately powered for the target effect.\n\n",
          sep = ""
        )
      } else {
        cat(
          "At least half of evaluated proteins were adequately powered for the target effect in every comparison, although protein-specific uncertainty still varies.\n\n"
        )
      }
    }

    cat(
      "These calculations describe sensitivity under the fitted model and observed standard errors. [Post-hoc power](#fig-post-hoc-power) is effect-dependent and does not turn a non-significant result into evidence of no biological difference.\n\n"
    )
  }

  cat("### Pathway analysis\n\n")
  usable_gsea <- gsea_table |>
    dplyr::mutate(
      p.adjust = suppressWarnings(as.numeric(.data$p.adjust)),
      NES = suppressWarnings(as.numeric(.data$NES))
    ) |>
    dplyr::filter(
      !is.na(.data$Comparison),
      !is.na(.data$Description),
      nzchar(.data$Description),
      is.finite(.data$p.adjust),
      is.finite(.data$NES)
    )

  if (nrow(usable_gsea) == 0) {
    cat("No interpretable GO Biological Process GSEA results were available.\n\n")
  } else {
    for (comparison_label in unique(usable_gsea$Comparison)) {
      comparison_gsea <- usable_gsea |>
        dplyr::filter(.data$Comparison == comparison_label)
      adjusted_hits <- sum(comparison_gsea$p.adjust < fdr_cutoff)
      top_terms <- comparison_gsea |>
        dplyr::arrange(.data$p.adjust, dplyr::desc(abs(.data$NES))) |>
        dplyr::distinct(.data$Description, .keep_all = TRUE) |>
        dplyr::slice_head(n = top_pathways)
      term_text <- paste0(
        top_terms$Description,
        " (NES ", sprintf("%.2f", top_terms$NES), ")",
        collapse = "; "
      )
      cat(
        "- ", linked_comparison_label(comparison_label, "#fig-gsea-dotplots"), ": ",
        format(adjusted_hits, big.mark = ","),
        " GO terms met the adjusted-p-value threshold. Leading terms were ",
        term_text, ".\n",
        sep = ""
      )
    }
    cat("\n")
  }

  cat(
    "GSEA uses the full ranked protein list, so pathway-level enrichment may be ",
    "present even when individual proteins do not pass both differential-analysis ",
    "thresholds. Closely related GO terms are biologically redundant and should be ",
    "interpreted together with the [similarity trees](#fig-gsea-treeplots) and ",
    "[protein-level volcano plots](#fig-volcano-fc-vs-healthy).\n"
  )

  invisible(list(
    differential_summary = differential_summary,
    power_summary = power_summary,
    pathway_results = usable_gsea
  ))
}
