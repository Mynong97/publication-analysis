#!/usr/bin/env Rscript

## ============================================================
## Downstream microbiome and virome analysis script
## ============================================================
##
## Generalized public version for manuscript reproducibility.
##
## Analyses included:
##   1. Alpha diversity
##   2. Beta diversity, PCoA, PERMANOVA, ANOSIM
##   3. Taxon filtering
##   4. MaAsLin2 differential abundance analysis
##   5. Selected-taxon abundance plots
##   6. S-plot using log2 fold change and CLR-Spearman correlation
##   7. CLR-based Spearman taxon-taxon network
##   8. CLR-based Spearman taxon-clinical heatmap
##
## Input tables should be prepared as:
##   - Feature table: taxa/features as rows and samples as columns,
##                    or samples as rows and taxa/features as columns.
##   - Metadata: one row per sample, containing SampleID and group.
##
## No raw sequencing data or local server paths are included.
##
## ============================================================


## ============================================================
## 0. User configuration
## ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

CONFIG <- list(
  ## Project directory.
  ## Use "." when this script is run from the project folder.
  project_dir = ".",

  ## Output folder
  output_dir = "results_downstream",

  ## Metadata and feature table names.
  ## Change these names according to your input files.
  feature_table_all          = "otu_table.csv",
  feature_table_longitudinal = "otu_table_longitudinal.csv",
  feature_table_endpoint     = "otu_table_endpoint.csv",

  metadata_all          = "metadata.csv",
  metadata_longitudinal = "metadata_longitudinal.csv",
  metadata_endpoint     = "metadata_endpoint.csv",

  ## Column names in metadata
  sample_col  = "SampleID",
  subject_col = "SubjectID",
  group_col   = "group",

  ## Candidate column names for the feature/taxon ID column
  taxon_col_candidates = c(
    "Taxon", "taxon", "FeatureID", "feature", "Species",
    "species", "Genus", "genus", "X.OTU.ID", "OTU", "OTU_ID"
  ),

  ## Analysis settings
  pseudocount = 1e-6,
  pseudo_count_scale = 10000,
  prevalence_threshold = 0.50,

  ## Filtering thresholds
  mean_abundance_threshold_bacteriome = 0.01,   # 1%
  mean_abundance_threshold_virome     = 0.001,  # 0.1%

  ## S-plot settings
  splot_log2fc_threshold = 1.0,
  splot_pcorr_threshold = 0.4,
  splot_point_size = 5,
  splot_label_size = 4,
  splot_use_log2fc = TRUE,

  ## Spearman network settings
  network_r_threshold = 0.40,
  network_p_threshold = 0.01,

  ## Heatmap settings
  use_q_for_heatmap_stars = FALSE,

  ## Plot/export settings
  save_csv  = TRUE,
  save_png  = TRUE,
  save_pptx = TRUE,
  save_emf  = FALSE,

  ## Optional selected taxa for abundance plots.
  ## Replace with taxa used in the manuscript, or set character(0) to skip.
  selected_taxa = c(
    "Akkermansia_muciniphila",
    "Xylanibacter_rodentium",
    "Prevotella_sp_MGM1",
    "Duncaniella_muris",
    "Duncaniella_dubosii",
    "Muribaculum_gordoncarteri"
  )
)


## ============================================================
## 1. Package loading
## ============================================================

load_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf(
        "Package '%s' is required but not installed. Please install it before running this script.",
        pkg
      )
    )
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

required_packages <- c(
  "dplyr", "tidyr", "readr", "tibble", "ggplot2",
  "vegan", "ape", "Hmisc", "igraph", "pheatmap",
  "grid", "gridExtra"
)

invisible(lapply(required_packages, load_or_stop))

optional_package_available <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

if (CONFIG$save_pptx) {
  if (!optional_package_available("officer") ||
      !optional_package_available("rvg")) {
    warning("Packages 'officer' and/or 'rvg' are not installed. PPTX export will be skipped.")
    CONFIG$save_pptx <- FALSE
  } else {
    suppressPackageStartupMessages({
      library(officer)
      library(rvg)
    })
  }
}

if (CONFIG$save_emf) {
  if (!optional_package_available("devEMF")) {
    warning("Package 'devEMF' is not installed. EMF export will be skipped.")
    CONFIG$save_emf <- FALSE
  }
}


## ============================================================
## 2. Helper functions
## ============================================================

setwd(CONFIG$project_dir)
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

read_table_auto <- function(fp, check.names = FALSE) {
  if (!file.exists(fp)) {
    stop("Input file not found: ", fp)
  }

  ext <- tolower(tools::file_ext(fp))

  if (ext %in% c("csv")) {
    out <- read.csv(fp, check.names = check.names, stringsAsFactors = FALSE)
  } else {
    out <- read.table(
      fp,
      header = TRUE,
      sep = "\t",
      check.names = check.names,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )
  }

  out
}

write_csv_safe <- function(x, fp, row.names = FALSE) {
  dir.create(dirname(fp), showWarnings = FALSE, recursive = TRUE)
  write.csv(x, fp, row.names = row.names)
}

p_to_label <- function(p) {
  if (is.na(p)) {
    "ns"
  } else if (p < 0.0001) {
    "****"
  } else if (p < 0.001) {
    "***"
  } else if (p < 0.01) {
    "**"
  } else if (p < 0.05) {
    "*"
  } else {
    "ns"
  }
}

clr_transform <- function(X, pseudocount = 1e-6) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"

  X2 <- X + pseudocount
  X2 <- sweep(X2, 2, colSums(X2, na.rm = TRUE), "/")
  L <- log(X2)
  L - matrix(
    colMeans(L, na.rm = TRUE),
    nrow = nrow(L),
    ncol = ncol(L),
    byrow = TRUE
  )
}

prepare_metadata <- function(meta, sample_col, group_col) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  names(meta) <- trimws(names(meta))

  if (!sample_col %in% colnames(meta)) {
    colnames(meta)[1] <- sample_col
  }

  if (!group_col %in% colnames(meta)) {
    group_candidates <- c(
      "Group", "group", "Condition", "condition",
      "Status", "status", "Class", "class"
    )

    detected <- intersect(group_candidates, colnames(meta))

    if (length(detected) == 0) {
      stop("Group column was not found in metadata.")
    }

    group_col <- detected[1]
    message("Detected group column: ", group_col)
  }

  rownames(meta) <- meta[[sample_col]]
  meta[[group_col]] <- as.factor(meta[[group_col]])

  meta
}

prepare_feature_table <- function(feature_df, meta, sample_col, taxon_col_candidates) {
  df <- as.data.frame(feature_df, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- trimws(names(df))

  taxon_col <- intersect(taxon_col_candidates, colnames(df))

  if (length(taxon_col) > 0) {
    taxon_col <- taxon_col[1]
    rownames(df) <- make.unique(as.character(df[[taxon_col]]))
    df[[taxon_col]] <- NULL
  } else {
    rownames(df) <- make.unique(as.character(df[[1]]))
    df[[1]] <- NULL
  }

  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  mat <- as.matrix(df)

  sample_ids <- rownames(meta)

  n_col_match <- length(intersect(colnames(mat), sample_ids))
  n_row_match <- length(intersect(rownames(mat), sample_ids))

  if (n_row_match > n_col_match) {
    mat <- t(mat)
  }

  common_samples <- intersect(colnames(mat), sample_ids)

  if (length(common_samples) < 2) {
    stop("Too few matched samples between feature table and metadata.")
  }

  mat <- mat[, common_samples, drop = FALSE]
  meta2 <- meta[common_samples, , drop = FALSE]

  list(
    taxa_by_samples = mat,
    samples_by_taxa = t(mat),
    metadata = meta2
  )
}

load_feature_and_metadata <- function(feature_file, metadata_file) {
  meta <- read_table_auto(metadata_file)
  meta <- prepare_metadata(
    meta = meta,
    sample_col = CONFIG$sample_col,
    group_col = CONFIG$group_col
  )

  feature_df <- read_table_auto(feature_file)

  prepare_feature_table(
    feature_df = feature_df,
    meta = meta,
    sample_col = CONFIG$sample_col,
    taxon_col_candidates = CONFIG$taxon_col_candidates
  )
}

save_plot_png <- function(plot_obj, fp, width = 7, height = 5, dpi = 300) {
  if (!CONFIG$save_png) return(invisible(NULL))

  dir.create(dirname(fp), showWarnings = FALSE, recursive = TRUE)

  if (optional_package_available("ragg")) {
    ragg::agg_png(
      filename = fp,
      width = width,
      height = height,
      units = "in",
      res = dpi
    )
    print(plot_obj)
    dev.off()
  } else {
    ggplot2::ggsave(fp, plot_obj, width = width, height = height, dpi = dpi)
  }
}

save_plot_pptx <- function(plot_obj, fp, width = 9, height = 6) {
  if (!CONFIG$save_pptx) return(invisible(NULL))

  doc <- officer::read_pptx()
  doc <- officer::add_slide(doc, layout = "Blank", master = "Office Theme")
  doc <- officer::ph_with(
    doc,
    rvg::dml(ggobj = plot_obj),
    location = officer::ph_location(
      left = 0.5,
      top = 0.5,
      width = width,
      height = height
    )
  )
  print(doc, target = fp)
}

save_plot_emf <- function(plot_obj, fp, width = 8, height = 6) {
  if (!CONFIG$save_emf) return(invisible(NULL))

  devEMF::emf(file = fp, width = width, height = height)
  print(plot_obj)
  dev.off()
}


## ============================================================
## 3. Alpha diversity
## ============================================================

run_alpha_diversity <- function(feature_file, metadata_file, label = "all") {
  if (!file.exists(feature_file) || !file.exists(metadata_file)) {
    message("Alpha diversity skipped: missing input files for ", label)
    return(NULL)
  }

  obj <- load_feature_and_metadata(feature_file, metadata_file)
  X <- obj$samples_by_taxa
  meta <- obj$metadata

  X_count <- round(X * CONFIG$pseudo_count_scale)

  alpha_df <- data.frame(
    SampleID = rownames(X_count),
    Shannon = vegan::diversity(X_count, index = "shannon"),
    Simpson = vegan::diversity(X_count, index = "simpson"),
    stringsAsFactors = FALSE
  )

  alpha_df <- alpha_df %>%
    left_join(
      meta %>%
        rownames_to_column("SampleID_meta"),
      by = c("SampleID" = "SampleID_meta")
    )

  out_dir <- file.path(CONFIG$output_dir, "alpha_diversity")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (CONFIG$save_csv) {
    write_csv_safe(
      alpha_df,
      file.path(out_dir, paste0("alpha_diversity_", label, ".csv")),
      row.names = FALSE
    )
  }

  for (metric in c("Shannon", "Simpson")) {
    p <- ggplot(alpha_df, aes(x = .data[[CONFIG$group_col]], y = .data[[metric]])) +
      geom_boxplot(outlier.shape = NA, color = "black", fill = "white") +
      geom_jitter(width = 0.12, size = 2, alpha = 0.9) +
      theme_classic(base_size = 12) +
      labs(
        title = paste0(metric, " diversity"),
        x = CONFIG$group_col,
        y = metric
      )

    save_plot_png(
      p,
      file.path(out_dir, paste0("alpha_diversity_", metric, "_", label, ".png")),
      width = 5,
      height = 4
    )

    save_plot_pptx(
      p,
      file.path(out_dir, paste0("alpha_diversity_", metric, "_", label, ".pptx")),
      width = 5.5,
      height = 4.5
    )
  }

  stat_list <- list()

  for (metric in c("Shannon", "Simpson")) {
    df <- alpha_df %>%
      filter(!is.na(.data[[metric]]), !is.na(.data[[CONFIG$group_col]]))

    if (nlevels(as.factor(df[[CONFIG$group_col]])) >= 2) {
      kw <- kruskal.test(df[[metric]] ~ df[[CONFIG$group_col]])

      stat_list[[paste0(metric, "_global")]] <- data.frame(
        metric = metric,
        comparison = "Global",
        test = "Kruskal-Wallis",
        p_value = kw$p.value,
        stringsAsFactors = FALSE
      )

      pw <- pairwise.wilcox.test(
        x = df[[metric]],
        g = df[[CONFIG$group_col]],
        p.adjust.method = "BH",
        exact = FALSE
      )

      pw_df <- as.data.frame(as.table(pw$p.value)) %>%
        filter(!is.na(Freq)) %>%
        transmute(
          metric = metric,
          comparison = paste(Var1, "vs", Var2),
          test = "Pairwise Wilcoxon with BH correction",
          p_value = Freq
        )

      stat_list[[paste0(metric, "_pairwise")]] <- pw_df
    }
  }

  stat_df <- bind_rows(stat_list)

  if (nrow(stat_df) > 0 && CONFIG$save_csv) {
    write_csv_safe(
      stat_df,
      file.path(out_dir, paste0("alpha_diversity_statistics_", label, ".csv")),
      row.names = FALSE
    )
  }

  alpha_df
}


## ============================================================
## 4. Beta diversity, PCoA, PERMANOVA, ANOSIM
## ============================================================

run_beta_diversity <- function(feature_file, metadata_file, label = "all") {
  if (!file.exists(feature_file) || !file.exists(metadata_file)) {
    message("Beta diversity skipped: missing input files for ", label)
    return(NULL)
  }

  obj <- load_feature_and_metadata(feature_file, metadata_file)
  X <- obj$samples_by_taxa
  meta <- obj$metadata

  out_dir <- file.path(CONFIG$output_dir, "beta_diversity")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  dist_bc <- vegan::vegdist(X, method = "bray")
  dist_mat <- as.matrix(dist_bc)

  pcoa_res <- ape::pcoa(dist_bc)
  explained_var <- pcoa_res$values$Relative_eig * 100

  pcoa_df <- data.frame(
    SampleID = rownames(pcoa_res$vectors),
    PC1 = pcoa_res$vectors[, 1],
    PC2 = pcoa_res$vectors[, 2],
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      meta %>% rownames_to_column("SampleID_meta"),
      by = c("SampleID" = "SampleID_meta")
    )

  if (CONFIG$save_csv) {
    write_csv_safe(
      dist_mat,
      file.path(out_dir, paste0("bray_curtis_matrix_", label, ".csv")),
      row.names = TRUE
    )

    write_csv_safe(
      pcoa_df,
      file.path(out_dir, paste0("pcoa_coordinates_", label, ".csv")),
      row.names = FALSE
    )
  }

  pc1_label <- sprintf("PC1: %.2f%%", explained_var[1])
  pc2_label <- sprintf("PC2: %.2f%%", explained_var[2])

  p <- ggplot(
    pcoa_df,
    aes(x = PC1, y = PC2, color = .data[[CONFIG$group_col]])
  ) +
    geom_point(size = 2.5) +
    stat_ellipse(type = "norm", level = 0.65, show.legend = FALSE) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("Bray-Curtis PCoA — ", label),
      x = pc1_label,
      y = pc2_label,
      color = CONFIG$group_col
    )

  save_plot_png(
    p,
    file.path(out_dir, paste0("pcoa_", label, ".png")),
    width = 5.5,
    height = 4.5
  )

  save_plot_pptx(
    p,
    file.path(out_dir, paste0("pcoa_", label, ".pptx")),
    width = 6.5,
    height = 5
  )

  set.seed(123)
  permanova_res <- vegan::adonis2(
    dist_bc ~ meta[[CONFIG$group_col]],
    permutations = 999
  )

  set.seed(123)
  anosim_res <- vegan::anosim(
    x = dist_bc,
    grouping = meta[[CONFIG$group_col]],
    permutations = 999
  )

  permanova_df <- as.data.frame(permanova_res)
  anosim_df <- data.frame(
    statistic_R = anosim_res$statistic,
    significance_p = anosim_res$signif,
    permutations = 999
  )

  if (CONFIG$save_csv) {
    write_csv_safe(
      permanova_df,
      file.path(out_dir, paste0("PERMANOVA_", label, ".csv")),
      row.names = TRUE
    )

    write_csv_safe(
      anosim_df,
      file.path(out_dir, paste0("ANOSIM_", label, ".csv")),
      row.names = FALSE
    )
  }

  list(
    distance = dist_bc,
    pcoa = pcoa_df,
    permanova = permanova_res,
    anosim = anosim_res
  )
}


## ============================================================
## 5. Taxon filtering
## ============================================================

filter_taxa_by_prevalence_and_abundance <- function(
    feature_file,
    metadata_file,
    label = "all",
    mean_abundance_threshold = CONFIG$mean_abundance_threshold_bacteriome
) {
  if (!file.exists(feature_file) || !file.exists(metadata_file)) {
    message("Taxon filtering skipped: missing input files for ", label)
    return(NULL)
  }

  obj <- load_feature_and_metadata(feature_file, metadata_file)
  taxa_by_samples <- obj$taxa_by_samples
  meta <- obj$metadata

  groups <- unique(meta[[CONFIG$group_col]])

  taxa_keep_by_prevalence <- lapply(groups, function(g) {
    samples_g <- rownames(meta)[meta[[CONFIG$group_col]] == g]
    Xg <- taxa_by_samples[, samples_g, drop = FALSE]

    prevalence <- rowMeans(Xg > 0, na.rm = TRUE)
    names(prevalence)[prevalence > CONFIG$prevalence_threshold]
  })

  taxa_keep_by_prevalence <- unique(unlist(taxa_keep_by_prevalence))

  X_prev <- taxa_by_samples[taxa_keep_by_prevalence, , drop = FALSE]

  rel <- sweep(X_prev, 2, colSums(X_prev, na.rm = TRUE), "/")
  mean_abund <- rowMeans(rel, na.rm = TRUE)

  taxa_keep_final <- names(mean_abund)[mean_abund >= mean_abundance_threshold]
  X_filtered <- X_prev[taxa_keep_final, , drop = FALSE]

  out_samples_by_taxa <- as.data.frame(t(X_filtered))
  out_samples_by_taxa <- rownames_to_column(out_samples_by_taxa, "SampleID")

  out_dir <- file.path(CONFIG$output_dir, "filtered_feature_tables")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  out_file <- file.path(
    out_dir,
    paste0("filtered_feature_table_", label, ".csv")
  )

  if (CONFIG$save_csv) {
    write_csv_safe(out_samples_by_taxa, out_file, row.names = FALSE)
  }

  message(
    "Taxon filtering completed for ",
    label,
    ": retained ",
    length(taxa_keep_final),
    " taxa."
  )

  list(
    taxa_by_samples = X_filtered,
    samples_by_taxa = t(X_filtered),
    metadata = meta,
    output_file = out_file
  )
}


## ============================================================
## 6. MaAsLin2 differential abundance analysis
## ============================================================

run_maaslin2_analysis <- function(
    filtered_table_file,
    metadata_file,
    label = "all",
    random_effect = NULL
) {
  if (!optional_package_available("Maaslin2")) {
    warning("Package 'Maaslin2' is not installed. MaAsLin2 analysis skipped.")
    return(NULL)
  }

  if (!file.exists(filtered_table_file) || !file.exists(metadata_file)) {
    message("MaAsLin2 skipped: missing input files for ", label)
    return(NULL)
  }

  out_dir <- file.path(CONFIG$output_dir, "maaslin2", label)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  metadata <- read_table_auto(metadata_file)

  fixed_effects <- CONFIG$group_col

  if (!is.null(random_effect) && random_effect %in% colnames(metadata)) {
    random_effects <- random_effect
  } else {
    random_effects <- NULL
  }

  fit <- Maaslin2::Maaslin2(
    input_data = filtered_table_file,
    input_metadata = metadata_file,
    output = out_dir,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    normalization = "CLR",
    transform = "NONE",
    analysis_method = "LM",
    max_significance = 1.0
  )

  all_results_file <- file.path(out_dir, "all_results.tsv")

  if (file.exists(all_results_file)) {
    res <- readr::read_tsv(all_results_file, show_col_types = FALSE) %>%
      filter(metadata == CONFIG$group_col) %>%
      select(feature, coef, stderr, pval, qval) %>%
      rename(
        Taxon = feature,
        Beta = coef,
        SE = stderr,
        p_value = pval,
        q_value = qval
      )

    write_csv_safe(
      res,
      file.path(out_dir, paste0("MaAsLin2_group_results_", label, ".csv")),
      row.names = FALSE
    )
  }

  fit
}


## ============================================================
## 7. Selected-taxon abundance plots
## ============================================================

run_selected_taxa_plots <- function(feature_file, metadata_file, label = "all") {
  if (length(CONFIG$selected_taxa) == 0) {
    message("Selected-taxon plots skipped: no selected taxa specified.")
    return(NULL)
  }

  if (!file.exists(feature_file) || !file.exists(metadata_file)) {
    message("Selected-taxon plots skipped: missing input files for ", label)
    return(NULL)
  }

  obj <- load_feature_and_metadata(feature_file, metadata_file)
  taxa_by_samples <- obj$taxa_by_samples
  meta <- obj$metadata

  present_taxa <- CONFIG$selected_taxa[CONFIG$selected_taxa %in% rownames(taxa_by_samples)]
  missing_taxa <- setdiff(CONFIG$selected_taxa, rownames(taxa_by_samples))

  if (length(missing_taxa) > 0) {
    warning("Missing selected taxa: ", paste(missing_taxa, collapse = ", "))
  }

  if (length(present_taxa) == 0) {
    message("No selected taxa were present in the feature table.")
    return(NULL)
  }

  otu_long <- as.data.frame(taxa_by_samples[present_taxa, , drop = FALSE]) %>%
    rownames_to_column("Taxon") %>%
    pivot_longer(
      cols = -Taxon,
      names_to = "SampleID",
      values_to = "RelAbund"
    ) %>%
    left_join(
      meta %>% rownames_to_column("SampleID"),
      by = "SampleID"
    ) %>%
    mutate(
      RelAbund = as.numeric(RelAbund),
      RelAbund_percent = RelAbund * 100,
      Taxon_label = gsub("_", " ", Taxon)
    )

  out_dir <- file.path(CONFIG$output_dir, "selected_taxa_plots")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  p_all <- ggplot(
    otu_long,
    aes(x = .data[[CONFIG$group_col]], y = RelAbund_percent)
  ) +
    geom_boxplot(outlier.shape = NA, color = "black", fill = "white") +
    geom_jitter(width = 0.12, size = 1.8, alpha = 0.9) +
    facet_wrap(~ Taxon_label, scales = "free_y") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = paste0("Selected taxa — ", label),
      x = CONFIG$group_col,
      y = "Relative abundance (%)"
    )

  save_plot_png(
    p_all,
    file.path(out_dir, paste0("selected_taxa_boxplots_", label, ".png")),
    width = 9,
    height = 6
  )

  save_plot_pptx(
    p_all,
    file.path(out_dir, paste0("selected_taxa_boxplots_", label, ".pptx")),
    width = 9,
    height = 6
  )

  if (CONFIG$save_csv) {
    summary_table <- otu_long %>%
      group_by(Taxon, .data[[CONFIG$group_col]]) %>%
      summarise(
        mean = mean(RelAbund, na.rm = TRUE),
        sd = sd(RelAbund, na.rm = TRUE),
        median = median(RelAbund, na.rm = TRUE),
        .groups = "drop"
      )

    write_csv_safe(
      summary_table,
      file.path(out_dir, paste0("selected_taxa_group_summary_", label, ".csv")),
      row.names = FALSE
    )
  }

  invisible(p_all)
}


## ============================================================
## 8. S-plot: log2 fold change + CLR-Spearman correlation
## ============================================================

run_splot_analysis <- function(feature_file, metadata_file, label = "endpoint") {
  if (!file.exists(feature_file) || !file.exists(metadata_file)) {
    message("S-plot skipped: missing input files for ", label)
    return(NULL)
  }

  obj <- load_feature_and_metadata(feature_file, metadata_file)
  taxa_by_samples <- obj$taxa_by_samples
  meta <- obj$metadata

  group_values <- as.character(meta[[CONFIG$group_col]])
  group_levels <- unique(group_values)

  if (length(group_levels) != 2) {
    warning("S-plot requires exactly two groups. Skipping: ", label)
    return(NULL)
  }

  group_A <- group_levels[1]
  group_B <- group_levels[2]

  message(
    sprintf(
      "S-plot group A = '%s' reference; group B = '%s' comparison",
      group_A,
      group_B
    )
  )

  grp_num <- ifelse(group_values == group_B, 1, 0)

  mean_A <- rowMeans(taxa_by_samples[, group_values == group_A, drop = FALSE], na.rm = TRUE)
  mean_B <- rowMeans(taxa_by_samples[, group_values == group_B, drop = FALSE], na.rm = TRUE)

  if (CONFIG$splot_use_log2fc) {
    log2FC <- log2((mean_B + CONFIG$pseudocount) / (mean_A + CONFIG$pseudocount))
    fc_label <- "log2 Fold Change (B vs A)"
  } else {
    log2FC <- (mean_B + CONFIG$pseudocount) / (mean_A + CONFIG$pseudocount)
    fc_label <- "Fold Change (B vs A)"
  }

  if (optional_package_available("compositions")) {
    otu_clr <- apply(
      taxa_by_samples + CONFIG$pseudocount,
      2,
      compositions::clr
    )
  } else {
    otu_clr <- clr_transform(taxa_by_samples, CONFIG$pseudocount)
  }

  p_corr <- apply(otu_clr, 1, function(x) {
    suppressWarnings(
      cor(as.numeric(x), grp_num, method = "spearman")
    )
  })

  res <- tibble(
    Taxon = rownames(taxa_by_samples),
    mean_A = mean_A,
    mean_B = mean_B,
    log2FC = log2FC,
    p_corr = p_corr
  ) %>%
    mutate(
      color_group = case_when(
        log2FC > CONFIG$splot_log2fc_threshold  ~ "B_enriched",
        log2FC < -CONFIG$splot_log2fc_threshold ~ "A_enriched",
        TRUE ~ "No_diff"
      ),
      label = ifelse(
        abs(p_corr) >= CONFIG$splot_pcorr_threshold &
          abs(log2FC) >= CONFIG$splot_log2fc_threshold,
        Taxon,
        ""
      )
    )

  color_map <- c(
    "B_enriched" = "green",
    "A_enriched" = "red",
    "No_diff" = "black"
  )

  out_dir <- file.path(CONFIG$output_dir, "splot")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (CONFIG$save_csv) {
    write_csv_safe(
      res,
      file.path(out_dir, paste0("splot_CLR_spearman_results_", label, ".csv")),
      row.names = FALSE
    )
  }

  p_splot_colored <- ggplot(
    res,
    aes(x = log2FC, y = p_corr, color = color_group)
  ) +
    geom_point(
      size = CONFIG$splot_point_size,
      alpha = 1,
      shape = 16,
      stroke = 0
    ) +
    scale_color_manual(values = color_map) +
    geom_vline(xintercept = 0, color = "gray", linewidth = 0.4) +
    geom_hline(yintercept = 0, color = "gray", linewidth = 0.4) +
    geom_vline(
      xintercept = c(CONFIG$splot_log2fc_threshold, -CONFIG$splot_log2fc_threshold),
      linetype = "dotted"
    ) +
    geom_hline(
      yintercept = c(CONFIG$splot_pcorr_threshold, -CONFIG$splot_pcorr_threshold),
      linetype = "dotted"
    ) +
    labs(
      title = sprintf(
        "S-plot (CLR + Spearman) — A = '%s', B = '%s'",
        group_A,
        group_B
      ),
      x = fc_label,
      y = "Spearman r (CLR abundance vs group)",
      color = "Enrichment"
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme_minimal(base_size = 13)

  if (optional_package_available("ggrepel")) {
    p_splot_labeled <- p_splot_colored +
      ggrepel::geom_text_repel(
        aes(label = label),
        color = "black",
        max.overlaps = Inf,
        box.padding = 0.4,
        point.padding = 0.2,
        seed = 7,
        min.segment.length = 0,
        segment.alpha = 0.6,
        size = CONFIG$splot_label_size
      ) +
      labs(title = "S-plot (CLR + Spearman) — highlighted taxa")
  } else {
    p_splot_labeled <- p_splot_colored +
      geom_text(
        aes(label = label),
        color = "black",
        size = CONFIG$splot_label_size,
        check_overlap = TRUE
      ) +
      labs(title = "S-plot (CLR + Spearman) — highlighted taxa")
  }

  print(p_splot_colored)
  print(p_splot_labeled)

  save_plot_png(
    p_splot_colored,
    file.path(out_dir, paste0("splot_CLR_spearman_colored_", label, ".png")),
    width = 7,
    height = 5
  )

  save_plot_png(
    p_splot_labeled,
    file.path(out_dir, paste0("splot_CLR_spearman_labeled_", label, ".png")),
    width = 7,
    height = 5
  )

  save_plot_emf(
    p_splot_colored,
    file.path(out_dir, paste0("splot_CLR_spearman_colored_", label, ".emf")),
    width = 7,
    height = 5
  )

  save_plot_emf(
    p_splot_labeled,
    file.path(out_dir, paste0("splot_CLR_spearman_labeled_", label, ".emf")),
    width = 7,
    height = 5
  )

  if (CONFIG$save_pptx) {
    pptx_file <- file.path(out_dir, paste0("Splot_CLR_spearman_", label, ".pptx"))

    doc <- officer::read_pptx()

    doc <- officer::add_slide(
      doc,
      layout = "Title and Content",
      master = "Office Theme"
    )

    doc <- officer::ph_with(
      doc,
      value = "S-plot: CLR + Spearman",
      location = officer::ph_location_type(type = "title")
    )

    doc <- officer::ph_with(
      doc,
      rvg::dml(ggobj = p_splot_colored),
      location = officer::ph_location_type(type = "body")
    )

    doc <- officer::add_slide(
      doc,
      layout = "Title and Content",
      master = "Office Theme"
    )

    doc <- officer::ph_with(
      doc,
      value = "S-plot: highlighted taxa",
      location = officer::ph_location_type(type = "title")
    )

    doc <- officer::ph_with(
      doc,
      rvg::dml(ggobj = p_splot_labeled),
      location = officer::ph_location_type(type = "body")
    )

    print(doc, target = pptx_file)
  }

  list(
    results = res,
    plot_colored = p_splot_colored,
    plot_labeled = p_splot_labeled
  )
}


## ============================================================
## 9. CLR-Spearman taxon-taxon network
## ============================================================

run_spearman_network <- function(feature_file, metadata_file, label = "endpoint") {
  if (!file.exists(feature_file) || !file.exists(metadata_file)) {
    message("Spearman network skipped: missing input files for ", label)
    return(NULL)
  }

  obj <- load_feature_and_metadata(feature_file, metadata_file)
  taxa_by_samples <- obj$taxa_by_samples
  meta <- obj$metadata

  if (ncol(taxa_by_samples) < 4) {
    warning("Too few samples for network analysis.")
    return(NULL)
  }

  group_levels <- levels(as.factor(meta[[CONFIG$group_col]]))

  if (length(group_levels) < 2) {
    warning("At least two groups are required for group-specific network analysis.")
    return(NULL)
  }

  out_dir <- file.path(CONFIG$output_dir, "spearman_network")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  taxa_by_samples_clr <- clr_transform(taxa_by_samples, CONFIG$pseudocount)

  all_cor <- suppressWarnings(
    Hmisc::rcorr(t(taxa_by_samples_clr), type = "spearman")$r
  )

  all_cor[is.na(all_cor)] <- 0

  dist_mat <- matrix(
    0,
    nrow = nrow(all_cor),
    ncol = ncol(all_cor),
    dimnames = dimnames(all_cor)
  )

  dist_mat[all_cor >= 0] <- 1 - abs(all_cor[all_cor >= 0])
  dist_mat[all_cor < 0]  <- 1 + abs(all_cor[all_cor < 0])

  g_all <- igraph::graph_from_adjacency_matrix(
    1 / (dist_mat + 1e-6),
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  set.seed(123)
  base_layout <- igraph::layout_with_fr(g_all, weights = igraph::E(g_all)$weight)

  build_edges_for_group <- function(group_label) {
    samp_idx <- which(meta[[CONFIG$group_col]] == group_label)

    if (length(samp_idx) < 3) {
      warning("Group ", group_label, " has fewer than 3 samples. Skipping.")
      return(data.frame())
    }

    Xg <- taxa_by_samples[, samp_idx, drop = FALSE]
    Xg_clr <- clr_transform(Xg, CONFIG$pseudocount)

    cor_out <- suppressWarnings(
      Hmisc::rcorr(t(Xg_clr), type = "spearman")
    )

    cor_mat <- cor_out$r
    p_mat <- cor_out$P

    upper_idx <- which(upper.tri(cor_mat), arr.ind = TRUE)
    rho_vals <- cor_mat[upper_idx]
    p_vals <- p_mat[upper_idx]

    sel <- which(
      abs(rho_vals) >= CONFIG$network_r_threshold &
        p_vals < CONFIG$network_p_threshold
    )

    if (length(sel) == 0) {
      return(data.frame())
    }

    data.frame(
      from = rownames(cor_mat)[upper_idx[sel, 1]],
      to = colnames(cor_mat)[upper_idx[sel, 2]],
      rho = rho_vals[sel],
      p_value = p_vals[sel],
      sign = ifelse(rho_vals[sel] > 0, "positive", "negative"),
      group = group_label,
      stringsAsFactors = FALSE
    )
  }

  edge_list <- lapply(group_levels, build_edges_for_group)
  names(edge_list) <- group_levels

  all_edges <- bind_rows(edge_list)

  if (CONFIG$save_csv) {
    write_csv_safe(
      all_edges,
      file.path(out_dir, paste0("spearman_network_edges_", label, ".csv")),
      row.names = FALSE
    )
  }

  if (nrow(all_edges) == 0) {
    message("No network edges passed the threshold.")
    return(NULL)
  }

  node_df <- data.frame(
    name = rownames(taxa_by_samples),
    stringsAsFactors = FALSE
  )

  graph_list <- list()

  for (g_label in group_levels) {
    edges_g <- edge_list[[g_label]]

    graph_g <- igraph::graph_from_data_frame(
      edges_g,
      vertices = node_df,
      directed = FALSE
    )

    deg <- igraph::degree(graph_g)
    igraph::V(graph_g)$size <- 10 + deg * 2
    igraph::V(graph_g)$color <- "grey80"

    if (igraph::ecount(graph_g) > 0) {
      igraph::E(graph_g)$color <- ifelse(
        igraph::E(graph_g)$sign == "positive",
        "#B22222",
        "steelblue"
      )
      igraph::E(graph_g)$width <- pmax(0.3, abs(igraph::E(graph_g)$rho) * 1.2)
    }

    graph_list[[g_label]] <- graph_g
  }

  if (CONFIG$save_pptx) {
    pptx_file <- file.path(out_dir, paste0("spearman_network_", label, ".pptx"))
    doc <- officer::read_pptx()

    for (g_label in names(graph_list)) {
      graph_g <- graph_list[[g_label]]

      doc <- officer::add_slide(
        doc,
        layout = "Title and Content",
        master = "Office Theme"
      )

      doc <- officer::ph_with(
        doc,
        value = sprintf(
          "Spearman network — %s (|rho| >= %.2f, p < %.2f)",
          g_label,
          CONFIG$network_r_threshold,
          CONFIG$network_p_threshold
        ),
        location = officer::ph_location_type(type = "title")
      )

      doc <- officer::ph_with(
        doc,
        rvg::dml(code = {
          plot(
            graph_g,
            layout = base_layout,
            vertex.label = igraph::V(graph_g)$name,
            vertex.label.cex = 0.6,
            vertex.label.color = "black",
            vertex.label.dist = 0.5,
            main = NULL
          )
        }),
        location = officer::ph_location_type(type = "body")
      )
    }

    print(doc, target = pptx_file)
  }

  graph_list
}


## ============================================================
## 10. CLR-Spearman taxon-clinical heatmap
## ============================================================

run_spearman_heatmap <- function(feature_file, metadata_file, label = "endpoint") {
  if (!file.exists(feature_file) || !file.exists(metadata_file)) {
    message("Spearman heatmap skipped: missing input files for ", label)
    return(NULL)
  }

  obj <- load_feature_and_metadata(feature_file, metadata_file)
  taxa_by_samples <- obj$taxa_by_samples
  meta <- obj$metadata

  id_cols <- c(
    CONFIG$sample_col,
    CONFIG$subject_col,
    CONFIG$group_col
  )

  clinical_candidates <- setdiff(colnames(meta), id_cols)
  clinical_keep <- clinical_candidates[
    sapply(meta[clinical_candidates], is.numeric)
  ]

  if (length(clinical_keep) < 1) {
    message("No numeric clinical variables found. Heatmap skipped.")
    return(NULL)
  }

  group_levels <- levels(as.factor(meta[[CONFIG$group_col]]))

  if (length(group_levels) < 1) {
    message("No group level found. Heatmap skipped.")
    return(NULL)
  }

  out_dir <- file.path(CONFIG$output_dir, "spearman_heatmap")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  taxa_order <- rownames(taxa_by_samples)
  otu_clr <- clr_transform(taxa_by_samples, CONFIG$pseudocount)

  make_q_matrix_rect <- function(p_mat, method = "BH") {
    q_mat <- matrix(
      NA_real_,
      nrow = nrow(p_mat),
      ncol = ncol(p_mat),
      dimnames = dimnames(p_mat)
    )

    p_vec <- as.vector(p_mat)
    valid <- !is.na(p_vec) & is.finite(p_vec)

    q_vec <- rep(NA_real_, length(p_vec))
    q_vec[valid] <- p.adjust(p_vec[valid], method = method)

    q_mat[] <- q_vec
    q_mat
  }

  make_long_table <- function(rho_mat, p_mat, q_mat, group_label) {
    out <- expand.grid(
      Taxon = rownames(rho_mat),
      Clinical = colnames(rho_mat),
      stringsAsFactors = FALSE
    )

    out$group <- group_label
    out$rho <- as.vector(rho_mat)
    out$p_value <- as.vector(p_mat)
    out$q_value <- as.vector(q_mat)

    out %>%
      select(group, Taxon, Clinical, rho, p_value, q_value) %>%
      arrange(q_value, p_value)
  }

  make_group_heatmap <- function(gval) {
    meta_g <- meta %>% filter(.data[[CONFIG$group_col]] == gval)
    common_samples <- intersect(rownames(meta_g), colnames(otu_clr))

    if (length(common_samples) < 3) {
      warning("Group ", gval, " has fewer than 3 samples. Skipping.")
      return(NULL)
    }

    meta_g <- meta_g[common_samples, , drop = FALSE]
    otu_g <- otu_clr[taxa_order, common_samples, drop = FALSE]
    clin_df <- meta_g[, clinical_keep, drop = FALSE]

    res_rho <- matrix(
      NA_real_,
      nrow = length(taxa_order),
      ncol = length(clinical_keep),
      dimnames = list(taxa_order, clinical_keep)
    )

    res_p <- matrix(
      NA_real_,
      nrow = length(taxa_order),
      ncol = length(clinical_keep),
      dimnames = list(taxa_order, clinical_keep)
    )

    for (i in seq_along(taxa_order)) {
      x <- as.numeric(otu_g[taxa_order[i], ])

      for (j in seq_along(clinical_keep)) {
        y <- as.numeric(clin_df[[clinical_keep[j]]])
        ok <- is.finite(x) & is.finite(y)

        if (sum(ok) >= 3 && sd(x[ok]) > 0 && sd(y[ok]) > 0) {
          tmp <- Hmisc::rcorr(x[ok], y[ok], type = "spearman")
          res_rho[i, j] <- tmp$r[1, 2]
          res_p[i, j] <- tmp$P[1, 2]
        }
      }
    }

    res_q <- make_q_matrix_rect(res_p, method = "BH")
    star_mat <- if (CONFIG$use_q_for_heatmap_stars) res_q else res_p

    stars <- ifelse(
      star_mat < 0.01,
      "**",
      ifelse(star_mat < 0.05, "*", "")
    )

    stars[is.na(stars)] <- ""

    bk <- seq(-1, 1, length.out = 101)
    cols <- colorRampPalette(c("blue", "white", "red"))(length(bk) - 1)

    ph <- pheatmap::pheatmap(
      res_rho,
      color = cols,
      breaks = bk,
      display_numbers = stars,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      fontsize_number = 10,
      angle_col = 45,
      na_col = "grey90",
      main = paste0("CLR + Spearman — group = ", gval),
      silent = TRUE
    )

    list(
      gtable = ph$gtable,
      rho = res_rho,
      p = res_p,
      q = res_q,
      long = make_long_table(res_rho, res_p, res_q, gval)
    )
  }

  heatmap_list <- lapply(group_levels, make_group_heatmap)
  names(heatmap_list) <- group_levels
  heatmap_list <- Filter(Negate(is.null), heatmap_list)

  if (length(heatmap_list) == 0) {
    message("No heatmap could be generated.")
    return(NULL)
  }

  if (CONFIG$save_csv) {
    all_long <- list()

    for (gval in names(heatmap_list)) {
      gsafe <- safe_filename(gval)

      write_csv_safe(
        heatmap_list[[gval]]$rho,
        file.path(out_dir, paste0("rho_matrix_group_", gsafe, "_", label, ".csv")),
        row.names = TRUE
      )

      write_csv_safe(
        heatmap_list[[gval]]$p,
        file.path(out_dir, paste0("p_matrix_group_", gsafe, "_", label, ".csv")),
        row.names = TRUE
      )

      write_csv_safe(
        heatmap_list[[gval]]$q,
        file.path(out_dir, paste0("q_matrix_group_", gsafe, "_", label, ".csv")),
        row.names = TRUE
      )

      write_csv_safe(
        heatmap_list[[gval]]$long,
        file.path(out_dir, paste0("rho_p_q_long_group_", gsafe, "_", label, ".csv")),
        row.names = FALSE
      )

      all_long[[gval]] <- heatmap_list[[gval]]$long
    }

    write_csv_safe(
      bind_rows(all_long),
      file.path(out_dir, paste0("rho_p_q_long_all_groups_", label, ".csv")),
      row.names = FALSE
    )
  }

  if (CONFIG$save_pptx) {
    pptx_file <- file.path(out_dir, paste0("CLR_Spearman_heatmaps_", label, ".pptx"))
    doc <- officer::read_pptx()

    for (gval in names(heatmap_list)) {
      gt <- heatmap_list[[gval]]$gtable

      doc <- officer::add_slide(
        doc,
        layout = "Title and Content",
        master = "Office Theme"
      )

      doc <- officer::ph_with(
        doc,
        value = sprintf("CLR + Spearman heatmap — group = %s", gval),
        location = officer::ph_location_type(type = "title")
      )

      doc <- officer::ph_with(
        doc,
        rvg::dml(code = {
          grid::grid.newpage()
          grid::grid.draw(gt)
        }),
        location = officer::ph_location_type(type = "body")
      )
    }

    if (length(heatmap_list) >= 2) {
      doc <- officer::add_slide(
        doc,
        layout = "Blank",
        master = "Office Theme"
      )

      left_gt <- heatmap_list[[1]]$gtable
      right_gt <- heatmap_list[[2]]$gtable

      doc <- officer::ph_with(
        doc,
        value = "CLR + Spearman heatmaps — side-by-side",
        location = officer::ph_location(
          left = 0.5,
          top = 0.3,
          width = 9,
          height = 0.6
        )
      )

      doc <- officer::ph_with(
        doc,
        rvg::dml(code = {
          grid::grid.newpage()
          grid::grid.draw(left_gt)
        }),
        location = officer::ph_location(
          left = 0.5,
          top = 1.1,
          width = 4.5,
          height = 6.0
        )
      )

      doc <- officer::ph_with(
        doc,
        rvg::dml(code = {
          grid::grid.newpage()
          grid::grid.draw(right_gt)
        }),
        location = officer::ph_location(
          left = 5.3,
          top = 1.1,
          width = 4.5,
          height = 6.0
        )
      )
    }

    print(doc, target = pptx_file)
  }

  heatmap_list
}


## ============================================================
## 11. Run all analyses
## ============================================================

run_all <- function() {
  message("============================================================")
  message("Starting downstream microbiome/virome analysis")
  message("============================================================")

  ## Longitudinal data
  if (
    file.exists(CONFIG$feature_table_longitudinal) &&
    file.exists(CONFIG$metadata_longitudinal)
  ) {
    message("\n[1] Longitudinal alpha diversity")
    run_alpha_diversity(
      CONFIG$feature_table_longitudinal,
      CONFIG$metadata_longitudinal,
      label = "longitudinal"
    )

    message("\n[2] Longitudinal beta diversity")
    run_beta_diversity(
      CONFIG$feature_table_longitudinal,
      CONFIG$metadata_longitudinal,
      label = "longitudinal"
    )

    message("\n[3] Longitudinal taxon filtering")
    filtered_longitudinal <- filter_taxa_by_prevalence_and_abundance(
      CONFIG$feature_table_longitudinal,
      CONFIG$metadata_longitudinal,
      label = "longitudinal",
      mean_abundance_threshold = CONFIG$mean_abundance_threshold_bacteriome
    )

    if (!is.null(filtered_longitudinal)) {
      message("\n[4] Longitudinal MaAsLin2")
      run_maaslin2_analysis(
        filtered_longitudinal$output_file,
        CONFIG$metadata_longitudinal,
        label = "longitudinal",
        random_effect = CONFIG$subject_col
      )
    }

    message("\n[5] Longitudinal selected-taxon plots")
    run_selected_taxa_plots(
      CONFIG$feature_table_longitudinal,
      CONFIG$metadata_longitudinal,
      label = "longitudinal"
    )
  } else {
    message("Longitudinal input files not found. Longitudinal analyses skipped.")
  }

  ## Endpoint data
  if (
    file.exists(CONFIG$feature_table_endpoint) &&
    file.exists(CONFIG$metadata_endpoint)
  ) {
    message("\n[6] Endpoint alpha diversity")
    run_alpha_diversity(
      CONFIG$feature_table_endpoint,
      CONFIG$metadata_endpoint,
      label = "endpoint"
    )

    message("\n[7] Endpoint beta diversity")
    run_beta_diversity(
      CONFIG$feature_table_endpoint,
      CONFIG$metadata_endpoint,
      label = "endpoint"
    )

    message("\n[8] Endpoint taxon filtering")
    filtered_endpoint <- filter_taxa_by_prevalence_and_abundance(
      CONFIG$feature_table_endpoint,
      CONFIG$metadata_endpoint,
      label = "endpoint",
      mean_abundance_threshold = CONFIG$mean_abundance_threshold_bacteriome
    )

    if (!is.null(filtered_endpoint)) {
      message("\n[9] Endpoint MaAsLin2")
      run_maaslin2_analysis(
        filtered_endpoint$output_file,
        CONFIG$metadata_endpoint,
        label = "endpoint",
        random_effect = NULL
      )
    }

    message("\n[10] Endpoint selected-taxon plots")
    run_selected_taxa_plots(
      CONFIG$feature_table_endpoint,
      CONFIG$metadata_endpoint,
      label = "endpoint"
    )

    message("\n[11] Endpoint S-plot")
    run_splot_analysis(
      CONFIG$feature_table_endpoint,
      CONFIG$metadata_endpoint,
      label = "endpoint"
    )

    message("\n[12] Endpoint CLR-Spearman network")
    run_spearman_network(
      CONFIG$feature_table_endpoint,
      CONFIG$metadata_endpoint,
      label = "endpoint"
    )

    message("\n[13] Endpoint CLR-Spearman heatmap")
    run_spearman_heatmap(
      CONFIG$feature_table_endpoint,
      CONFIG$metadata_endpoint,
      label = "endpoint"
    )
  } else {
    message("Endpoint input files not found. Endpoint analyses skipped.")
  }

  ## All data, optional
  if (
    file.exists(CONFIG$feature_table_all) &&
    file.exists(CONFIG$metadata_all)
  ) {
    message("\n[14] All-sample selected-taxon plots")
    run_selected_taxa_plots(
      CONFIG$feature_table_all,
      CONFIG$metadata_all,
      label = "all"
    )
  } else {
    message("All-sample input files not found. Optional all-sample analyses skipped.")
  }

  message("\n============================================================")
  message("Analysis finished.")
  message("Results saved in: ", normalizePath(CONFIG$output_dir))
  message("============================================================")
}

run_all()
