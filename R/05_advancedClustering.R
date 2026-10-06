# ==============================================================================
# 05_ADVANCED_CLUSTERING_COMPARISONS.R - Klastry K-means Profili Ekspresji
# ==============================================================================
cat("\n=================================================================\n")
cat(" >>> KROK 05: KLASTROWANIE PROFILI EKSPRESJI (K-MEANS) <<<\n")
cat("=================================================================\n")

if (!exists("PATHS") || !exists("BIO") || !exists("STATS")) stop("Brak załadowanego 01_config.R!")
utils_path <- file.path(PATHS$scripts_dir, PATHS$utils_file)
if (file.exists(utils_path)) source(utils_path)

rds_path <- file.path(PATHS$out_base_dir, "txi_salmon_data.rds")
data_bundle  <- readRDS(rds_path)
txi_full     <- data_bundle$txi
colData_Full <- data_bundle$colData
if (!exists("anno_map") || !is.data.frame(anno_map) || nrow(anno_map) == 0) {
  cat("[ADNOTACJE] Budowanie słownika genów z GTF...\n")
  anno_map <- build_custom_annotation_map(PATHS$merged_gtf)
} else {
  cat(sprintf("[ADNOTACJE] Wykorzystano istniejący słownik anno_map z pamięci (%d genów).\n", nrow(anno_map)))
}

subset_txi <- function(txi_obj, sample_names) {
  list(
    abundance           = txi_obj$abundance[, sample_names, drop = FALSE],
    counts              = txi_obj$counts[, sample_names, drop = FALSE],
    length              = txi_obj$length[, sample_names, drop = FALSE],
    countsFromAbundance = txi_obj$countsFromAbundance
  )
}

# --- FUNKCJA GŁÓWNA KLASTROWANIA ---
generate_diet_expression_clusters <- function(vsd_clean, all_sig_genes, group_dir, group_name, condition_type, k = STATS$cluster_k) {
  all_sig_genes <- unique(all_sig_genes)
  cat(sprintf("  -> Analiza klastrów dla %d unikalnych DEG (Wariant: %s, k = %d)...\n", length(all_sig_genes), condition_type, k))
  
  if (length(all_sig_genes) <= k * 2) {
    cat("     [!] Zbyt mało unikalnych genów do klastrowania K-means.\n")
    return(NULL)
  }
  
  # Wyciągnięcie macierzy VST dla DEG
  mat <- assay(vsd_clean)[rownames(vsd_clean) %in% all_sig_genes, , drop = FALSE]
  syms <- anno_map$symbol[match(rownames(mat), anno_map$gene_id)]
  valid_idx <- !grepl(BIO$noise_filter_regex, syms) & !is.na(syms)
  mat <- mat[valid_idx, , drop = FALSE]
  rownames(mat) <- syms[valid_idx]
  
  if (nrow(mat) <= k * 2) return(NULL)
  
  # Standaryzacja wierszowa (Z-score dla każdego genu)
  mat_scaled <- t(scale(t(mat)))
  
  # Algorytm K-means
  set.seed(42)
  km <- kmeans(mat_scaled, centers = k, nstart = 25)
  
  col_df <- as.data.frame(colData(vsd_clean))
  df_long <- as.data.frame(mat_scaled)
  df_long$Gene <- rownames(df_long)
  df_long$Cluster <- paste("Klaster", km$cluster)
  
  df_long <- tidyr::pivot_longer(df_long, cols = c(-Gene, -Cluster), names_to = "Sample", values_to = "Z_Score")
  df_long$Condition <- col_df$Condition[match(df_long$Sample, rownames(col_df))]
  df_long$Condition <- factor(df_long$Condition, levels = names(PALETTES$condition_colors))
  
  # Rysowanie wykresu spaghetti ze średnią trajektorią
  p <- ggplot(df_long, aes(x = Condition, y = Z_Score, group = Gene)) +
    geom_line(alpha = 0.15, color = "#2B83BA") +
    stat_summary(aes(group = 1), fun = mean, geom = "line", color = "#D7191C", linewidth = 1.8) +
    stat_summary(aes(group = 1), fun = mean, geom = "point", color = "#D7191C", size = 3) +
    facet_wrap(~ Cluster, scales = "free_y", ncol = 2) +
    labs(
      title = paste("Profile Ekspresji (K-means Z-score) -", toupper(group_name)),
      subtitle = paste("Wariant:", condition_type, "| Liczba genów:", nrow(mat)),
      x = BIO$condition_label,
      y = "Z-score Ekspresji"
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
      strip.background = element_rect(fill = "#E0F3F8")
    )
  
  out_dir_clust <- file.path(group_dir, "Visualizations", "Clusters", condition_type)
  if (!dir.exists(out_dir_clust)) dir.create(out_dir_clust, recursive = TRUE)
  
  out_png <- file.path(out_dir_clust, paste0("Clusters_K", k, "_", group_name, ".png"))
  ggsave(out_png, plot = p, width = 12, height = 8, dpi = 300)
  cat(sprintf("     [+] Zapisano wykres klastrów w: %s\n", out_png))
  
  # Zapis przypisania genów do klastrów do CSV
  cluster_assignments <- data.frame(
    Symbol  = names(km$cluster),
    Cluster = paste("Klaster", km$cluster),
    stringsAsFactors = FALSE
  )
  write.csv(cluster_assignments, file.path(out_dir_clust, paste0("Gene_Cluster_Assignments_K", k, ".csv")), row.names = FALSE)
  
  return(p)
}

# --- PĘTLA PO TKANKACH / GRUPACH ---
for (group in BIO$target_tissues) {
  group_dir <- file.path(PATHS$out_base_dir, group)
  samples_group <- if (group == "Whole_Dataset") rownames(colData_Full) else rownames(colData_Full)[colData_Full$Tissue == group]
  colData_group <- droplevels(colData_Full[samples_group, ])
  txi_group     <- subset_txi(txi_full, samples_group)
  design_formula <- ~ Condition
  
  # Dynamiczny dobór wariantu zależnie od flagi enable_filtering
  analysis_variants <- if (isTRUE(STATS$enable_filtering)) c("Unfiltered", "Filtered") else c("Standard")
  
  for (cond_type in analysis_variants) {
    if (cond_type == "Filtered") {
      best_samples  <- select_best_n_samples(txi_group, colData_group, BIO$samples_per_grp, design_formula)
      colData_run   <- colData_group[best_samples, ]
      txi_run       <- subset_txi(txi_group, best_samples)
    } else {
      colData_run   <- colData_group
      txi_run       <- txi_group
    }
    
    dds <- DESeqDataSetFromTximport(txi_run, colData = colData_run, design = design_formula)
    dds <- estimateSizeFactors(dds)
    vsd <- vst(dds, blind = FALSE)
    
    sig_dir <- file.path(group_dir, "significant", cond_type)
    sig_files <- list.files(sig_dir, pattern = "_shrink_sig\\.csv$", full.names = TRUE)
    
    if (length(sig_files) == 0) {
      cat(sprintf("[!] Brak plików _shrink_sig.csv w folderze: %s\n", sig_dir))
      next
    }
    
    # Zbieramy wszystkie unikalne ID genów ze wszystkich kontrastów
    all_sig_ids <- character()
    for (f in sig_files) {
      df_tmp <- read.csv(f)
      if (nrow(df_tmp) > 0 && "gene_id" %in% colnames(df_tmp)) {
        all_sig_ids <- c(all_sig_ids, df_tmp$gene_id)
      }
    }
    
    generate_diet_expression_clusters(vsd, unique(all_sig_ids), group_dir, group, cond_type, k = STATS$cluster_k)
  }
}
cat("\n[SUKCES] Krok 05 zakończony pomyślnie.\n")
