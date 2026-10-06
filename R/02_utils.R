# ==============================================================================
# 02_UTILS.R - Silnik Funkcji QC, Statystyki, Wizualizacji, Adnotacji i Multi-Omics
# PROJEKT: Danio rerio (SSOT Architecture)
# ==============================================================================
cat("\n[UTILS] Ładowanie silnika analitycznego...\n")

# ==============================================================================
# BLOK 1: ZAAWANSOWANA BUDOWA SŁOWNIKA ADNOTACJI I BIOTYPÓW (GTF + DB)
# ==============================================================================
build_custom_annotation_map <- function(gtf_path, ref_gtf_path = NULL) {
  cat(sprintf("[ADNOTACJE] Budowanie słownika genów z GTF: %s...\n", basename(gtf_path)))
  
  if (is.null(ref_gtf_path) && exists("PATHS") && !is.null(PATHS$ref_gtf)) {
    ref_gtf_path <- PATHS$ref_gtf
  }
  
  gtf_data <- rtracklayer::import(gtf_path)
  gtf_tx <- gtf_data[gtf_data$type %in% c("transcript", "mRNA"), ]
  mcols_tx <- mcols(gtf_tx)
  gene_ids <- as.character(gtf_tx$gene_id)
  
  if ("ref_gene_id" %in% colnames(mcols_tx)) {
    has_ref <- !is.na(mcols_tx$ref_gene_id) & mcols_tx$ref_gene_id != ""
    gene_ids[has_ref] <- as.character(mcols_tx$ref_gene_id[has_ref])
  }
  
  symbols <- if ("gene_name" %in% colnames(mcols_tx)) as.character(mcols_tx$gene_name) else rep(NA, length(gene_ids))
  symbols <- ifelse(is.na(symbols) | symbols == "", gene_ids, symbols)
  
  biotypes <- rep("unknown", length(gene_ids))
  for (attr_name in c("gene_biotype", "gene_type", "transcript_biotype", "biotype")) {
    if (attr_name %in% colnames(mcols_tx)) {
      found_bt <- as.character(mcols_tx[[attr_name]])
      valid_bt <- !is.na(found_bt) & found_bt != ""
      biotypes[valid_bt] <- found_bt[valid_bt]
      break
    }
  }
  
  df_map <- data.frame(
    gene_id     = gene_ids,
    symbol      = symbols,
    biotype     = biotypes,
    description = "No description available",
    stringsAsFactors = FALSE
  )
  
  df_map <- df_map[order(df_map$symbol == df_map$gene_id), ]
  df_map <- df_map[!duplicated(df_map$gene_id) & !is.na(df_map$gene_id), ]
  rownames(df_map) <- df_map$gene_id
  
  if (!is.null(ref_gtf_path) && file.exists(ref_gtf_path)) {
    cat(sprintf("  -> Dociąganie biotypów z referencji Ensembl: %s...\n", basename(ref_gtf_path)))
    tryCatch({
      ref_data <- rtracklayer::import(ref_gtf_path)
      ref_genes <- ref_data[ref_data$type == "gene", ]
      ref_mcols <- mcols(ref_genes)
      
      bt_col <- intersect(c("gene_biotype", "gene_type", "biotype"), colnames(ref_mcols))[1]
      if (!is.na(bt_col)) {
        ref_df <- data.frame(
          gene_id = as.character(ref_genes$gene_id),
          biotype = as.character(ref_mcols[[bt_col]]),
          stringsAsFactors = FALSE
        )
        ref_df <- ref_df[!duplicated(ref_df$gene_id) & !is.na(ref_df$biotype), ]
        
        m_ref <- match(df_map$gene_id, ref_df$gene_id)
        has_ref_bt <- !is.na(m_ref)
        df_map$biotype[has_ref_bt] <- ref_df$biotype[m_ref[has_ref_bt]]
      }
    }, error = function(e) cat("     [!] Błąd parsowania referencyjnego GTF:", conditionMessage(e), "\n"))
  }
  
  df_map$biotype[df_map$biotype == "unknown" & grepl("^MSTRG", df_map$gene_id)] <- "novel_isoform"
  
  if (requireNamespace("org.Dr.eg.db", quietly = TRUE)) {
    cat("  -> Dociąganie pełnych opisów genów z bazy org.Dr.eg.db...\n")
    tryCatch({
      all_db <- suppressMessages(
        AnnotationDbi::select(
          org.Dr.eg.db::org.Dr.eg.db,
          keys = AnnotationDbi::keys(org.Dr.eg.db::org.Dr.eg.db, keytype = "ENTREZID"),
          columns = c("SYMBOL", "GENENAME", "ENSEMBL"),
          keytype = "ENTREZID"
        )
      )
      all_db <- all_db[!is.na(all_db$GENENAME) & all_db$GENENAME != "", ]
      
      db_ens  <- all_db[!is.na(all_db$ENSEMBL), ]
      m_ens   <- match(df_map$gene_id, db_ens$ENSEMBL)
      has_ens <- !is.na(m_ens)
      df_map$description[has_ens] <- db_ens$GENENAME[m_ens[has_ens]]
      
      missing_desc <- which(df_map$description == "No description available" & !grepl(BIO$noise_filter_regex, df_map$symbol))
      db_sym_clean <- all_db[!duplicated(toupper(all_db$SYMBOL)), ]
      m_sym <- match(toupper(df_map$symbol[missing_desc]), toupper(db_sym_clean$SYMBOL))
      has_sym <- !is.na(m_sym)
      df_map$description[missing_desc[has_sym]] <- db_sym_clean$GENENAME[m_sym[has_sym]]
      
      cat(sprintf("     [+] Przypisano opisy biologiczne dla %d genów.\n", sum(df_map$description != "No description available")))
    }, error = function(e) cat("     [!] Ostrzeżenie przy adnotacji:", conditionMessage(e), "\n"))
  }
  
  cat(sprintf("[ADNOTACJE] Zbudowano słownik dla %d genów.\n", nrow(df_map)))
  return(df_map)
}

# ==============================================================================
# BLOK 2: ZARZĄDZANIE STRUKTURĄ KATALOGÓW I FILTRACJA PRÓBEK
# ==============================================================================
setup_group_folders <- function(base_dir, group_name) {
  g_dir <- file.path(base_dir, group_name)
  active_variants <- if (isTRUE(STATS$enable_filtering)) c("Filtered", "Unfiltered") else c("Standard")
  
  for (cond in active_variants) {
    dir.create(file.path(g_dir, "significant", cond), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(g_dir, "Visualizations", "Volcano", cond), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(g_dir, "Visualizations", "Heatmaps", cond), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(g_dir, "Visualizations", "QC", cond, "MA_Plots"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(g_dir, "Visualizations", "UpSet", cond), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(g_dir, "Visualizations", "Trajectories", cond), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(g_dir, "Visualizations", "Clusters", cond), showWarnings = FALSE, recursive = TRUE)
  }
  dir.create(file.path(g_dir, "Visualizations", "gProfiler"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(g_dir, "Visualizations", "Sankey", "Standard"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(g_dir, "Excel_Reports"), showWarnings = FALSE, recursive = TRUE)
  
  return(g_dir)
}

# Selekcja najbardziej spójnych próbek (gdy włączony tryb Filtered)
select_best_n_samples <- function(txi_sub, coldata_sub, n_target, design_formula) {
  selected <- character()
  for (grp in unique(coldata_sub$Condition)) {
    s_grp <- rownames(coldata_sub)[coldata_sub$Condition == grp]
    if (length(s_grp) <= n_target) {
      selected <- c(selected, s_grp)
    } else {
      counts_sub <- txi_sub$counts[, s_grp, drop = FALSE]
      cor_mat <- cor(counts_sub, method = "spearman")
      mean_cors <- rowMeans(cor_mat)
      best_s <- names(sort(mean_cors, decreasing = TRUE))[1:n_target]
      selected <- c(selected, best_s)
    }
  }
  return(selected)
}

# ==============================================================================
# BLOK 3: RAPORTY QC (SEQUENCING, SALMON, MACIERZ ODLEGŁOŚCI, PCA)
# ==============================================================================
export_comprehensive_qc_report <- function(txi, colData, salmon_dir, out_dir, qc_tr_dir = "02_qc/tr") {
  cat("\n[QC REPORT] Generowanie rozbudowanego raportu kontroli jakości sekwencjonowania i mapowania...\n")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  samples <- rownames(colData)
  salmon_stats <- list()
  
  for (s in samples) {
    json_file <- file.path(salmon_dir, s, "aux_info", "meta_info.json")
    if (file.exists(json_file)) {
      j_data <- tryCatch(jsonlite::fromJSON(json_file), error = function(e) NULL)
      if (!is.null(j_data)) {
        salmon_stats[[s]] <- data.frame(
          Processed_Reads  = if (!is.null(j_data$num_processed)) j_data$num_processed else NA,
          Mapped_Reads     = if (!is.null(j_data$num_mapped)) j_data$num_mapped else NA,
          Mapping_Rate_Pct = if (!is.null(j_data$percent_mapped)) round(j_data$percent_mapped, 2) else NA,
          Library_Type     = if (!is.null(j_data$library_types)) paste(j_data$library_types, collapse = ";") else "Unknown",
          stringsAsFactors = FALSE
        )
      }
    }
    if (!s %in% names(salmon_stats)) {
      salmon_stats[[s]] <- data.frame(Processed_Reads = NA, Mapped_Reads = NA, Mapping_Rate_Pct = NA, Library_Type = "Unknown", stringsAsFactors = FALSE)
    }
  }
  df_salmon <- do.call(rbind, salmon_stats)
  
  df_multiqc <- data.frame(Sample = samples, Raw_Reads = NA, Q30_Rate_Pct = NA, GC_Content_Pct = NA, stringsAsFactors = FALSE)
  mqc_file <- file.path(qc_tr_dir, "multiqc_data", "multiqc_general_stats.txt")
  if (!file.exists(mqc_file)) mqc_file <- list.files(qc_tr_dir, pattern = "multiqc_general_stats\\.txt$", full.names = TRUE, recursive = TRUE)[1]
  
  if (!is.na(mqc_file) && file.exists(mqc_file)) {
    mqc_raw <- tryCatch(read.delim(mqc_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(mqc_raw) && "Sample" %in% colnames(mqc_raw)) {
      for (i in seq_along(samples)) {
        s <- samples[i]
        match_idx <- grep(s, mqc_raw$Sample)
        if (length(match_idx) > 0) {
          row_m <- mqc_raw[match_idx[1], ]
          q30_col <- grep("percent_q30|q30_rate", colnames(row_m), ignore.case = TRUE, value = TRUE)[1]
          gc_col  <- grep("percent_gc|gc_content", colnames(row_m), ignore.case = TRUE, value = TRUE)[1]
          raw_col <- grep("total_reads|raw_reads", colnames(row_m), ignore.case = TRUE, value = TRUE)[1]
          
          if (!is.na(q30_col)) df_multiqc$Q30_Rate_Pct[i] <- round(as.numeric(row_m[[q30_col]]) * (if (max(as.numeric(row_m[[q30_col]]), na.rm=T) <= 1) 100 else 1), 2)
          if (!is.na(gc_col))  df_multiqc$GC_Content_Pct[i] <- round(as.numeric(row_m[[gc_col]]) * (if (max(as.numeric(row_m[[gc_col]]), na.rm=T) <= 1) 100 else 1), 2)
          if (!is.na(raw_col)) df_multiqc$Raw_Reads[i] <- as.numeric(row_m[[raw_col]])
        }
      }
    }
  }
  
  counts_mat <- txi$counts[, samples, drop = FALSE]
  tpm_mat    <- txi$abundance[, samples, drop = FALSE]
  
  qc_summary <- data.frame(
    Sample_ID                 = samples,
    Condition                 = colData$Condition,
    Raw_Sequencing_Reads      = df_multiqc$Raw_Reads,
    Trimmed_Input_Reads       = df_salmon$Processed_Reads,
    Salmon_Mapped_Reads       = df_salmon$Mapped_Reads,
    Mapping_Rate_Percent      = df_salmon$Mapping_Rate_Pct,
    Q30_Bases_Percent         = df_multiqc$Q30_Rate_Pct,
    GC_Content_Percent        = df_multiqc$GC_Content_Pct,
    Library_Strandness        = df_salmon$Library_Type,
    Total_Estimated_Counts    = round(colSums(counts_mat)),
    Detected_Genes_Count      = colSums(counts_mat > 0),
    Active_Genes_TPM_ge_1     = colSums(tpm_mat >= 1.0),
    High_Expressing_TPM_ge_10 = colSums(tpm_mat >= 10.0),
    stringsAsFactors          = FALSE
  )
  
  write.csv(qc_summary, file.path(out_dir, "Comprehensive_Mapping_and_QC_Summary.csv"), row.names = FALSE)
  
  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet("Sequencing_QC_Summary")$add_data("Sequencing_QC_Summary", x = qc_summary)
  h_dims <- openxlsx2::wb_dims(rows = 1, cols = 1:ncol(qc_summary))
  wb$add_fill("Sequencing_QC_Summary", dims = h_dims, color = openxlsx2::wb_color(hex = "1F497D"))
  wb$add_font("Sequencing_QC_Summary", dims = h_dims, color = openxlsx2::wb_color(hex = "FFFFFF"), bold = TRUE)
  wb$add_cell_style("Sequencing_QC_Summary", dims = h_dims, horizontal = "center")
  wb$set_col_widths("Sequencing_QC_Summary", cols = 1:ncol(qc_summary), widths = "auto")
  wb$save(file.path(out_dir, "Comprehensive_Mapping_and_QC_Summary.xlsx"), overwrite = TRUE)
  
  cat("  [+] Zapisano: Comprehensive_Mapping_and_QC_Summary.xlsx / .csv\n")
  return(qc_summary)
}

export_mapping_stats <- function(txi, colData, out_dir, file_name = NULL) {
  salmon_dir <- if (exists("PATHS") && !is.null(PATHS$salmon_dir)) PATHS$salmon_dir else "03_mapping_quant_GRCz11/salmon"
  export_comprehensive_qc_report(txi, colData, salmon_dir = salmon_dir, out_dir = out_dir)
}

export_coldata <- function(colData, out_dir, file_name = "colData_Summary.csv") {
  cd_df <- as.data.frame(colData)
  cd_df$Sample_ID <- rownames(cd_df)
  cd_df <- cd_df[, c("Sample_ID", setdiff(colnames(cd_df), "Sample_ID")), drop = FALSE]
  write.csv(cd_df, file.path(out_dir, file_name), row.names = FALSE)
}

generate_QC_plots <- function(dds, vsd_clean, group_dir, group_name, condition_type, do_dispersion = TRUE, do_pca = TRUE, do_dist = TRUE) {
  qc_dir <- file.path(group_dir, "Visualizations", "QC", condition_type)
  if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (do_dispersion) {
    disp_title <- if (exists("LABELS") && !is.null(LABELS$qc_dispersion_title)) sprintf(LABELS$qc_dispersion_title, group_name, condition_type) else sprintf("Dispersion Estimates - %s (%s)", group_name, condition_type)
    png(file.path(qc_dir, "Dispersion_Plot.png"), width = 8, height = 6, units = "in", res = 300)
    plotDispEsts(dds, main = disp_title)
    dev.off()
  }
  
  if (do_pca) {
    pca_data  <- plotPCA(vsd_clean, intgroup = c("Condition"), returnData = TRUE)
    pca_title <- if (exists("LABELS") && !is.null(LABELS$qc_pca_title)) sprintf(LABELS$qc_pca_title, group_name) else sprintf("PCA - %s", group_name)
    pca_sub   <- if (exists("LABELS") && !is.null(LABELS$qc_pca_sub)) sprintf(LABELS$qc_pca_sub, condition_type) else condition_type
    cond_lbl  <- if (exists("BIO") && !is.null(BIO$condition_label)) BIO$condition_label else "Condition"
    
    p_pca <- ggplot(pca_data, aes(PC1, PC2, fill = Condition)) + 
      geom_point(size = 4, alpha = 0.9, shape = 21, color = "black") + 
      ggrepel::geom_text_repel(aes(label = name), size = 3, show.legend = FALSE, color = "black") +
      scale_fill_manual(values = PALETTES$condition_colors) + 
      labs(title = pca_title, subtitle = pca_sub, fill = cond_lbl) +
      theme_bw(base_size = 11)
    
    ggsave(file.path(qc_dir, "PCA_Plot.png"), plot = p_pca, width = 9, height = 7, dpi = 300)
  }
  
  if (do_dist) {
    sampleDists <- dist(t(assay(vsd_clean)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- colnames(vsd_clean)
    colnames(sampleDistMatrix) <- colnames(vsd_clean)
    
    ann_df <- data.frame(
      Condition = colData(vsd_clean)$Condition,
      row.names = colnames(vsd_clean)
    )
    
    present_conds <- levels(factor(ann_df$Condition))
    ann_colors <- list(
      Condition = PALETTES$condition_colors[names(PALETTES$condition_colors) %in% present_conds]
    )
    
    dist_title <- if (exists("LABELS") && !is.null(LABELS$qc_dist_title)) sprintf(LABELS$qc_dist_title, group_name, condition_type) else sprintf("Sample Distance Matrix - %s (%s)", group_name, condition_type)
    
    png(file.path(qc_dir, "Sample_Distance_Heatmap.png"), width = 10, height = 8, units = "in", res = 300)
    pheatmap::pheatmap(
      sampleDistMatrix,
      clustering_distance_rows = sampleDists,
      clustering_distance_cols = sampleDists,
      col                      = PALETTES$dist_colors,
      annotation_col           = ann_df,
      annotation_colors        = ann_colors,
      main                     = dist_title,
      fontsize                 = 9,
      fontsize_row             = 8,
      fontsize_col             = 8
    )
    dev.off()
    cat(sprintf("     [+] Zapisano: %s/Sample_Distance_Heatmap.png\n", qc_dir))
  }
}

# ==============================================================================
# BLOK 4: KONTRASTY I PODSUMOWANIE ANALIZY RÓŻNICOWEJ (DESEQ2)
# ==============================================================================
run_contrast_analysis <- function(dds, contrast_vec, group_name, group_dir, condition_type, anno_map, alpha_val = STATS$alpha_thr, lfc_val = STATS$lfc_thr) {
  num <- contrast_vec[2]
  denom <- contrast_vec[3]
  comp_name <- paste0(num, "_vs_", denom)
  cat(sprintf("  -> Kontrast: %s (Alpha=%.2f, LFC=%.2f)...\n", comp_name, alpha_val, lfc_val))
  
  # 1. DESeq2 Results & Shrinkage
  res_raw <- if (isTRUE(STATS$use_ihw)) {
    results(dds, contrast = contrast_vec, alpha = alpha_val, filterFun = ihw)
  } else {
    results(dds, contrast = contrast_vec, alpha = alpha_val)
  }
  
  res_sh <- lfcShrink(dds, contrast = contrast_vec, type = STATS$lfc_shrink, res = res_raw, quiet = TRUE)
  
  df_sh <- as.data.frame(res_sh)
  df_sh$gene_id     <- rownames(df_sh)
  df_sh$symbol      <- anno_map$symbol[match(df_sh$gene_id, anno_map$gene_id)]
  df_sh$biotype     <- anno_map$biotype[match(df_sh$gene_id, anno_map$gene_id)]
  df_sh$description <- anno_map$description[match(df_sh$gene_id, anno_map$gene_id)]
  df_sh$symbol      <- ifelse(is.na(df_sh$symbol) | df_sh$symbol == "", df_sh$gene_id, df_sh$symbol)
  
  # 2. Selekcja istotnych genów (DEG)
  out_dir_sig <- file.path(group_dir, "significant", condition_type)
  df_sig <- df_sh[!is.na(df_sh$padj) & df_sh$padj <= alpha_val & abs(df_sh$log2FoldChange) >= lfc_val, ]
  
  write.csv(df_sh[order(df_sh$padj, na.last = TRUE), ], file.path(out_dir_sig, paste0(group_name, "_", comp_name, "_shrink_full.csv")), row.names = FALSE)
  write.csv(df_sig[order(df_sig$padj, na.last = TRUE), ], file.path(out_dir_sig, paste0(group_name, "_", comp_name, "_shrink_sig.csv")), row.names = FALSE)
  
  # 3. MA Plot (zabezpieczony graphics device)
  ma_dir <- file.path(group_dir, "Visualizations", "QC", condition_type, "MA_Plots")
  if (!dir.exists(ma_dir)) dir.create(ma_dir, recursive = TRUE)
  ma_title <- sprintf("MA Plot - %s (%s)", comp_name, condition_type)
  
  png(file.path(ma_dir, paste0("MA_Plot_", comp_name, ".png")), width = 7, height = 6, units = "in", res = 300)
  DESeq2::plotMA(res_sh, main = ma_title, ylim = c(-5, 5))
  dev.off()
  
  # 4. NATYWNY, KULOODPORNY VOLCANO PLOT (Czysty ggplot2)
  vol_dir <- file.path(group_dir, "Visualizations", "Volcano", condition_type)
  if (!dir.exists(vol_dir)) dir.create(vol_dir, recursive = TRUE)
  
  df_plot <- df_sh[!is.na(df_sh$padj) & !is.na(df_sh$log2FoldChange) & is.finite(df_sh$log2FoldChange), ]
  
  # Zabezpieczenie przed padj == 0 (-log10(0) = Inf), co ubijało silnik C++ ggrepel
  df_plot$neg_log10_padj <- -log10(pmax(df_plot$padj, 1e-300))
  
  df_plot$Regulation <- "NS"
  df_plot$Regulation[df_plot$padj <= alpha_val & df_plot$log2FoldChange >= lfc_val]  <- "Up"
  df_plot$Regulation[df_plot$padj <= alpha_val & df_plot$log2FoldChange <= -lfc_val] <- "Down"
  df_plot$Regulation <- factor(df_plot$Regulation, levels = c("Up", "Down", "NS"))
  
  # Etykietowanie wyłącznie Top 20 najważniejszych genów
  top_labels_df <- df_plot[df_plot$Regulation %in% c("Up", "Down") & !grepl(BIO$noise_filter_regex, df_plot$symbol), ]
  if (nrow(top_labels_df) > 0) {
    top_labels_df <- head(top_labels_df[order(top_labels_df$padj, -abs(top_labels_df$log2FoldChange)), ], 20)
  }
  
  p_vol <- ggplot(df_plot, aes(x = log2FoldChange, y = neg_log10_padj, color = Regulation)) +
    geom_point(alpha = 0.55, size = 1.6) +
    scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey70")) +
    geom_vline(xintercept = c(-lfc_val, lfc_val), linetype = "dashed", color = "grey30", linewidth = 0.6) +
    geom_hline(yintercept = -log10(alpha_val), linetype = "dashed", color = "grey30", linewidth = 0.6) +
    labs(
      title    = sprintf("Volcano Plot: %s", comp_name),
      subtitle = sprintf("Cutoffs: FDR <= %.2f, |Log2FC| >= %.2f | Annotated Top 20 DEGs", alpha_val, lfc_val),
      x        = "log2(Fold Change)",
      y        = "-log10(Adjusted P-Value)",
      color    = "Expression"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "right",
      plot.title      = element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle   = element_text(size = 9.5, hjust = 0.5, color = "grey30")
    )
  
  if (nrow(top_labels_df) > 0) {
    p_vol <- p_vol + ggrepel::geom_text_repel(
      data          = top_labels_df,
      aes(label     = symbol),
      size          = 3.5,
      color         = "black",
      fontface      = "bold.italic",
      max.overlaps  = 15,
      box.padding   = 0.35,
      point.padding = 0.3
    )
  }
  
  ggsave(file.path(vol_dir, paste0("Volcano_", comp_name, ".png")), plot = p_vol, width = 8.5, height = 6.5, dpi = 300)
  rm(p_vol, df_plot)
  
  # 5. Podsumowanie do tabeli
  sig_up <- df_sig[df_sig$log2FoldChange > 0, ]
  sig_dn <- df_sig[df_sig$log2FoldChange < 0, ]
  
  top_up_str <- if (nrow(sig_up) > 0) {
    paste(head(sprintf("%s (+%.2f)", sig_up[order(sig_up$log2FoldChange, decreasing = TRUE), "symbol"], 
                       sig_up[order(sig_up$log2FoldChange, decreasing = TRUE), "log2FoldChange"]), 3), collapse = "; ")
  } else "None"
  
  top_dn_str <- if (nrow(sig_dn) > 0) {
    paste(head(sprintf("%s (%.2f)", sig_dn[order(sig_dn$log2FoldChange, decreasing = FALSE), "symbol"], 
                       sig_dn[order(sig_dn$log2FoldChange, decreasing = FALSE), "log2FoldChange"]), 3), collapse = "; ")
  } else "None"
  
  summary_row <- data.frame(
    Group                 = group_name,
    Contrast              = comp_name,
    Baseline_Condition    = denom,
    Target_Condition      = num,
    FDR_Cutoff            = alpha_val,
    LFC_Cutoff            = lfc_val,
    Total_Genes_Tested    = nrow(df_sh),
    Total_DEGs            = nrow(df_sig),
    DEGs_Percent_Tested   = round(100 * (nrow(df_sig) / max(1, nrow(df_sh))), 2),
    Up_Regulated_Count    = nrow(sig_up),
    Down_Regulated_Count  = nrow(sig_dn),
    Top3_Up_Regulated     = top_up_str,
    Top3_Down_Regulated   = top_dn_str,
    stringsAsFactors      = FALSE
  )
  
  df_sig$contrast <- comp_name
  df_sig$status   <- ifelse(df_sig$log2FoldChange > 0, "Up", "Down")
  
  return(list(
    summary   = summary_row,
    sig_data  = df_sig[, c("gene_id", "symbol", "status", "contrast")],
    sig_df    = df_sig,
    full_df   = df_sh,
    comp_name = comp_name
  ))
}

# ==============================================================================
# BLOK 5: EKSPORT DO MASTER EXCELA (openxlsx2)
# ==============================================================================
export_to_master_excel <- function(res_list, group_dir, group_name, condition_type) {
  cat(sprintf("[EXCEL] Generowanie raportu .xlsx dla: %s (%s)...\n", group_name, condition_type))
  wb <- openxlsx2::wb_workbook()
  
  for (comp_name in names(res_list)) {
    df <- res_list[[comp_name]]$full_df
    df <- df[order(df$padj, na.last = TRUE), ]
    sheet_name <- substr(comp_name, 1, 31)
    
    wb$add_worksheet(sheet_name)$add_data(sheet = sheet_name, x = df)
    
    header_dims <- openxlsx2::wb_dims(rows = 1, cols = 1:ncol(df))
    wb$add_fill(sheet = sheet_name, dims = header_dims, color = openxlsx2::wb_color(hex = "4F81BD"))
    wb$add_font(sheet = sheet_name, dims = header_dims, color = openxlsx2::wb_color(hex = "FFFFFF"), bold = TRUE)
    wb$add_cell_style(sheet = sheet_name, dims = header_dims, horizontal = "center")
    
    padj_col_idx <- which(colnames(df) == "padj")
    if (length(padj_col_idx) > 0) {
      sig_rows <- which(!is.na(df$padj) & df$padj <= STATS$alpha_thr) + 1
      if (length(sig_rows) > 0) {
        sig_dims <- openxlsx2::wb_dims(rows = sig_rows, cols = padj_col_idx)
        wb$add_fill(sheet = sheet_name, dims = sig_dims, color = openxlsx2::wb_color(hex = "C6EFCE"))
        wb$add_font(sheet = sheet_name, dims = sig_dims, color = openxlsx2::wb_color(hex = "006100"))
      }
    }
    wb$set_col_widths(sheet = sheet_name, cols = 1:ncol(df), widths = "auto")
  }
  
  out_path <- file.path(group_dir, "Excel_Reports", paste0("Master_DE_Report_", group_name, "_", condition_type, ".xlsx"))
  wb$save(out_path, overwrite = TRUE)
}

# ==============================================================================
# BLOK 6: WIZUALIZACJE ZBIORCZE (HEATMAPY, UPSET, TRAJEKTORIE)
# ==============================================================================

# --- 6.1. Wielowariantowe Heatmapy DEGs (Top 50, 75, 100, 150, All) ---
generate_master_heatmaps <- function(vsd_clean, all_sig_genes, group_name, group_dir, condition_type, anno_map, top_n_vec = c(50, 75, 100, 150, "all")) {
  all_sig_genes <- unique(all_sig_genes)
  if (length(all_sig_genes) < 3) {
    cat("     [!] Zbyt mało istotnych DEG do wygenerowania heatmap.\n")
    return(NULL)
  }
  
  mat_sig <- assay(vsd_clean)[rownames(vsd_clean) %in% all_sig_genes, , drop = FALSE]
  syms <- anno_map$symbol[match(rownames(mat_sig), anno_map$gene_id)]
  clean_syms <- ifelse(is.na(syms) | syms == "", rownames(mat_sig), syms)
  
  valid_rows <- !grepl(BIO$noise_filter_regex, clean_syms)
  if (sum(valid_rows) >= 3) {
    mat_sig    <- mat_sig[valid_rows, , drop = FALSE]
    clean_syms <- clean_syms[valid_rows]
  }
  
  rownames(mat_sig) <- make.unique(clean_syms)
  gene_vars <- matrixStats::rowVars(mat_sig)
  ordered_idx <- order(gene_vars, decreasing = TRUE)
  
  out_hm_dir <- file.path(group_dir, "Visualizations", "Heatmaps", condition_type)
  if (!dir.exists(out_hm_dir)) dir.create(out_hm_dir, recursive = TRUE, showWarnings = FALSE)
  
  ann_df <- data.frame(
    Condition = colData(vsd_clean)$Condition,
    row.names = colnames(vsd_clean)
  )
  present_conds <- levels(factor(ann_df$Condition))
  ann_colors <- list(Condition = PALETTES$condition_colors[names(PALETTES$condition_colors) %in% present_conds])
  
  for (n_val in top_n_vec) {
    target_n <- if (n_val == "all") nrow(mat_sig) else min(as.numeric(n_val), nrow(mat_sig))
    sub_mat  <- mat_sig[ordered_idx[1:target_n], , drop = FALSE]
    
    file_label <- if (n_val == "all") "Heatmap_All_SigDEGs" else paste0("Heatmap_Top", target_n)
    png_path   <- file.path(out_hm_dir, paste0(file_label, "_", condition_type, ".png"))
    plot_height <- max(8, min(24, 0.18 * target_n + 3))
    
    hm_title <- if (exists("LABELS") && !is.null(LABELS$heatmap_top_title)) {
      if (n_val == "all") {
        sprintf("All Significant DEGs (%d) - %s", nrow(mat_sig), condition_type)
      } else {
        sprintf(gsub("75", as.character(target_n), LABELS$heatmap_top_title), condition_type)
      }
    } else {
      sprintf("Top %d DEGs - %s", target_n, condition_type)
    }
    
    png(png_path, width = 11, height = plot_height, units = "in", res = 300)
    pheatmap::pheatmap(
      sub_mat,
      scale             = "row",
      annotation_col    = ann_df,
      annotation_colors = ann_colors,
      col               = PALETTES$heatmap_colors,
      show_rownames     = (target_n <= 100),
      show_colnames     = TRUE,
      cluster_cols      = TRUE,
      cluster_rows      = TRUE,
      fontsize_row      = if (target_n <= 50) 8 else 6.5,
      fontsize_col      = 8,
      angle_col         = 45,
      main              = hm_title
    )
    dev.off()
    cat(sprintf("     [+] Zapisano: %s\n", basename(png_path)))
  }
}

# --- 6.2. Zaawansowany UpSet (ComplexUpset: Omnibus + Per Group) ---
generate_upset_and_contrast_breakdowns <- function(sig_data_list, res_list, group_dir, group_name, condition_type) {
  cat("\n  -> [UpSet] Generowanie wykresów UpSet (Omnibus + Per Condition Group)...\n")
  out_upset_dir <- file.path(group_dir, "Visualizations", "UpSet", condition_type)
  dir.create(out_upset_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
    install.packages("ComplexUpset", repos = "https://cloud.r-project.org")
  }
  
  all_genes <- unique(unlist(lapply(res_list, function(x) x$sig_df$gene_id)))
  all_genes <- all_genes[!grepl(BIO$noise_filter_regex, all_genes)]
  
  if (length(all_genes) < 2) {
    cat("     [!] Zbyt mało DEG do wygenerowania wykresów UpSet.\n")
    return(NULL)
  }
  
  df_binary_all  <- data.frame(gene_id = all_genes, stringsAsFactors = FALSE)
  df_binary_up   <- data.frame(gene_id = all_genes, stringsAsFactors = FALSE)
  df_binary_down <- data.frame(gene_id = all_genes, stringsAsFactors = FALSE)
  
  contrast_names <- names(res_list)
  
  for (comp_name in contrast_names) {
    df_sig   <- res_list[[comp_name]]$sig_df
    sig_ids  <- df_sig$gene_id[df_sig$gene_id %in% all_genes]
    up_ids   <- df_sig$gene_id[df_sig$log2FoldChange > 0 & df_sig$gene_id %in% all_genes]
    down_ids <- df_sig$gene_id[df_sig$log2FoldChange < 0 & df_sig$gene_id %in% all_genes]
    
    df_binary_all[[comp_name]]  <- df_binary_all$gene_id %in% sig_ids
    df_binary_up[[comp_name]]   <- df_binary_up$gene_id %in% up_ids
    df_binary_down[[comp_name]] <- df_binary_down$gene_id %in% down_ids
  }
  
  upset_ylab <- if (exists("LABELS") && !is.null(LABELS$upset_ylabel)) LABELS$upset_ylabel else "Shared DEGs"
  
  build_custom_upset <- function(data_bin, cols_to_plot, title_text, ylab_text, bar_color = "#333333") {
    ComplexUpset::upset(
      data_bin,
      cols_to_plot,
      name = "Contrasts",
      width_ratio = 0.28,
      min_size = 1,
      set_sizes = (
        ComplexUpset::upset_set_size(
          geom = geom_bar(fill = bar_color, width = 0.65)
        ) +
          geom_text(
            aes(label = after_stat(count)),
            stat = "count",
            position = position_stack(vjust = 0.5),
            color = "white",
            size = 3.2,
            fontface = "bold"
          ) +
          ylab("Total DEGs") +
          theme(
            axis.text.x  = element_text(angle = 45, hjust = 1, size = 8.5),
            axis.title.x = element_text(face = "bold", size = 9.5)
          )
      ),
      base_annotations = list(
        'Intersection' = (
          ComplexUpset::intersection_size(
            text = list(size = 3.2, vjust = -0.5, fontface = "bold"),
            bar_number_threshold = 1,
            fill = bar_color
          ) +
            ylab(ylab_text) +
            theme(axis.title.y = element_text(face = "bold", size = 10))
        )
      )
    )
  }
  
  # Omnibus: All / Up / Down
  p_all <- build_custom_upset(df_binary_all, contrast_names, "All DEGs", upset_ylab, bar_color = "#333333")
  ggsave(file.path(out_upset_dir, paste0("UpSet_Omnibus_All_DEGs_", condition_type, ".png")),
         plot = p_all, width = 14, height = 8.5, dpi = 300)
  
  sub_up <- df_binary_up[rowSums(df_binary_up[, contrast_names]) > 0, ]
  if (nrow(sub_up) > 0) {
    p_up <- build_custom_upset(sub_up, contrast_names, "Up-regulated", paste(upset_ylab, "(Up-regulated)"), bar_color = "#B2182B")
    ggsave(file.path(out_upset_dir, paste0("UpSet_Omnibus_UpRegulated_", condition_type, ".png")),
           plot = p_up, width = 14, height = 8.5, dpi = 300)
  }
  
  sub_down <- df_binary_down[rowSums(df_binary_down[, contrast_names]) > 0, ]
  if (nrow(sub_down) > 0) {
    p_down <- build_custom_upset(sub_down, contrast_names, "Down-regulated", paste(upset_ylab, "(Down-regulated)"), bar_color = "#2166AC")
    ggsave(file.path(out_upset_dir, paste0("UpSet_Omnibus_DownRegulated_", condition_type, ".png")),
           plot = p_down, width = 14, height = 8.5, dpi = 300)
  }
  
  # Wykresy per warunek eksperymentalny
  group_subset_dir <- file.path(out_upset_dir, "Grouped_By_Condition")
  dir.create(group_subset_dir, showWarnings = FALSE)
  
  unique_conditions <- if (exists("PALETTES") && !is.null(names(PALETTES$condition_colors))) {
    names(PALETTES$condition_colors)
  } else {
    unique(unlist(strsplit(contrast_names, "_vs_")))
  }
  
  for (cond in unique_conditions) {
    matching_contrasts <- grep(paste0("(^|_)", cond, "($|_)"), contrast_names, value = TRUE)
    
    if (length(matching_contrasts) >= 2) {
      sub_cond_df <- df_binary_all[rowSums(df_binary_all[, matching_contrasts, drop = FALSE]) > 0, c("gene_id", matching_contrasts)]
      
      if (nrow(sub_cond_df) >= 2) {
        bar_col <- if (!is.null(PALETTES$condition_colors[cond])) PALETTES$condition_colors[cond] else "#404040"
        p_cond <- build_custom_upset(
          data_bin     = sub_cond_df,
          cols_to_plot = matching_contrasts,
          title_text   = sprintf("Contrasts involving %s", cond),
          ylab_text    = sprintf("Shared DEGs (%s)", cond),
          bar_color    = bar_col
        )
        
        out_file <- file.path(group_subset_dir, paste0("UpSet_Grouped_", cond, "_", condition_type, ".png"))
        ggsave(out_file, plot = p_cond, width = 10, height = 6.2, dpi = 300)
        cat(sprintf("     [+] Zapisano UpSet dla diety [%s]: %s\n", cond, basename(out_file)))
      }
    }
  }
}

generate_upset_plot <- function(...) invisible(NULL)

# --- 6.3. Trajektorie Ekspresji Top DEG (Ścisłe SSOT z LABELS) ---
generate_top_genes_trajectories <- function(vsd_clean, res_list, group_dir, group_name, condition_type, anno_map, top_n_genes = 16) {
  cat("  -> Generowanie profili ekspresji dla czołowych DEG...\n")
  out_traj_dir <- file.path(group_dir, "Visualizations", "Trajectories", condition_type)
  if (!dir.exists(out_traj_dir)) dir.create(out_traj_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Wybór najbardziej istotnych genów
  top_gene_ids <- character()
  for (comp_name in names(res_list)) {
    df_sig <- res_list[[comp_name]]$sig_df
    if (!is.null(df_sig) && nrow(df_sig) > 0) {
      df_clean <- df_sig[!grepl(BIO$noise_filter_regex, df_sig$symbol) & !is.na(df_sig$padj), ]
      if (nrow(df_clean) > 0) {
        df_sorted <- df_clean[order(df_clean$padj, -abs(df_clean$log2FoldChange)), ]
        top_gene_ids <- c(top_gene_ids, head(df_sorted$gene_id, 4))
      }
    }
  }
  
  top_gene_ids <- unique(top_gene_ids)
  top_gene_ids <- intersect(top_gene_ids, rownames(vsd_clean))
  
  if (length(top_gene_ids) == 0) {
    rv <- matrixStats::rowVars(assay(vsd_clean))
    top_gene_ids <- head(rownames(vsd_clean)[order(rv, decreasing = TRUE)], top_n_genes)
  } else {
    top_gene_ids <- head(top_gene_ids, top_n_genes)
  }
  
  mat <- assay(vsd_clean)[top_gene_ids, , drop = FALSE]
  syms <- anno_map$symbol[match(rownames(mat), anno_map$gene_id)]
  clean_syms <- ifelse(is.na(syms) | syms == "", rownames(mat), syms)
  rownames(mat) <- make.unique(clean_syms)
  
  # check.names = FALSE zapobiega psuciu nazw próbek (brak NA na osi X)
  df_mat <- data.frame(mat, Gene = rownames(mat), check.names = FALSE, stringsAsFactors = FALSE)
  df_long <- tidyr::pivot_longer(df_mat, cols = -Gene, names_to = "Sample", values_to = "Expression")
  
  col_df <- as.data.frame(colData(vsd_clean))
  df_long$Condition <- col_df$Condition[match(df_long$Sample, rownames(col_df))]
  
  # Usunięcie ewentualnych wierszy bez przypisania
  df_long <- df_long[!is.na(df_long$Condition), ]
  
  # Zachowanie kolejności poziomów z palety SSOT
  if (exists("PALETTES") && !is.null(names(PALETTES$condition_colors))) {
    present_conds <- unique(as.character(df_long$Condition))
    ordered_levels <- intersect(names(PALETTES$condition_colors), present_conds)
    df_long$Condition <- factor(df_long$Condition, levels = ordered_levels)
  }
  
  traj_title <- sprintf("Top DEGs Normalized Expression - %s", condition_type)
  traj_ylab  <- "Normalized VST Expression"
  cond_lbl   <- if (exists("BIO") && !is.null(BIO$condition_label)) BIO$condition_label else "Condition"
  
  p <- ggplot(df_long, aes(x = Condition, y = Expression, fill = Condition)) + 
    geom_boxplot(alpha = 0.65, outlier.shape = NA, width = 0.55, color = "black") + 
    geom_jitter(width = 0.2, size = 2.0, alpha = 0.9, shape = 21, color = "black") + 
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3.0, fill = "red", color = "black") +
    facet_wrap(~ Gene, scales = "free_y", ncol = 4) + 
    scale_fill_manual(values = PALETTES$condition_colors) + 
    labs(
      title    = traj_title,
      subtitle = "Diamonds indicate group mean | Individual points represent biological replicates",
      x        = cond_lbl,
      y        = traj_ylab
    ) + 
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "none",
      axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text       = element_text(face = "bold.italic")
    )
  
  out_png <- file.path(out_traj_dir, paste0("Top_Trajectories_", condition_type, ".png"))
  ggsave(out_png, plot = p, width = 12, height = 9, dpi = 300)
  cat(sprintf("     [+] Zapisano poprawione profile ekspresji: %s\n", out_png))
  return(p)
}

# ==============================================================================
# BLOK 7: ANALIZA FUNKCJONALNA (g:Profiler), SANKEY ORAZ ZBIORCZY DOTPLOT
# ==============================================================================
run_gprofiler_analysis <- function(df_sig, group_name, comp_name, gp_dir, alpha_thr = STATS$alpha_thr, lfc_thr = STATS$lfc_thr) {
  cat(sprintf("  -> [g:Profiler] Analiza szlaków dla: %s...\n", comp_name))
  exclude_iea_val <- if (!is.null(BIO$gprofiler_exclude_iea)) BIO$gprofiler_exclude_iea else FALSE
  
  clean_ids <- df_sig$gene_id[!is.na(df_sig$gene_id) & !grepl(BIO$noise_filter_regex, df_sig$gene_id)]
  clean_ids <- unique(clean_ids)
  
  clean_symbols <- df_sig$symbol[
    !is.na(df_sig$symbol) & 
      df_sig$symbol != "" & 
      !grepl(BIO$noise_filter_regex, df_sig$symbol) & 
      !grepl("^ENSDARG", df_sig$symbol)
  ]
  clean_symbols <- unique(clean_symbols)
  
  if (length(clean_ids) < 5 && length(clean_symbols) < 5) {
    cat(sprintf("     [!] Zbyt mała liczba genów (< 5) do analizy dla: %s\n", comp_name))
    return(NULL)
  }
  
  call_gost <- function(query_vec) {
    tryCatch(
      suppressMessages(
        gprofiler2::gost(
          query             = query_vec,
          organism          = BIO$gprofiler_org,
          ordered_query     = FALSE,
          multi_query       = FALSE,
          exclude_iea       = exclude_iea_val,
          sources           = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP"),
          evcodes           = TRUE,
          correction_method = "fdr",
          user_threshold    = alpha_thr
        )
      ),
      error = function(e) NULL
    )
  }
  
  gp <- NULL
  query_mode <- "ENSEMBL_ID"
  
  if (length(clean_ids) >= 5) {
    gp_try_id <- call_gost(clean_ids)
    if (!is.null(gp_try_id) && !is.null(gp_try_id$result) && nrow(gp_try_id$result) > 0) {
      gp <- gp_try_id
      query_mode <- "ENSEMBL_ID"
      cat(sprintf("     [+] Zidentyfikowano szlaki po Ensembl ID (%d terminów).\n", nrow(gp$result)))
    }
  }
  
  if (is.null(gp) && length(clean_symbols) >= 5) {
    gp_try_sym <- call_gost(clean_symbols)
    if (!is.null(gp_try_sym) && !is.null(gp_try_sym$result) && nrow(gp_try_sym$result) > 0) {
      gp <- gp_try_sym
      query_mode <- "GENE_SYMBOL"
      cat(sprintf("     [+] Zidentyfikowano szlaki po Symbolach (%d terminów).\n", nrow(gp$result)))
    }
  }
  
  if (is.null(gp) || is.null(gp$result) || nrow(gp$result) == 0) {
    cat(sprintf("     [!] Brak istotnych szlaków (FDR < %.2f) dla %s.\n", alpha_thr, comp_name))
    return(NULL)
  }
  
  res_tbl <- gp$result[order(gp$result$p_value), ]
  write.csv(apply(res_tbl, 2, as.character), file.path(gp_dir, paste0("gProfiler_", comp_name, ".csv")), row.names = FALSE)
  
  res_tbl$Recall <- res_tbl$intersection_size / res_tbl$term_size
  top_terms <- res_tbl %>% dplyr::group_by(source) %>% dplyr::top_n(10, wt = -p_value) %>% dplyr::ungroup()
  top_terms$term_name <- stringr::str_wrap(top_terms$term_name, width = 45)
  
  gp_title <- if (exists("LABELS") && !is.null(LABELS$gprofiler_title)) sprintf(LABELS$gprofiler_title, comp_name) else sprintf("Top Pathways - %s", comp_name)
  gp_sub   <- if (exists("LABELS") && !is.null(LABELS$gprofiler_sub)) sprintf(LABELS$gprofiler_sub, BIO$gprofiler_org, query_mode, exclude_iea_val) else sprintf("Organism: %s", BIO$gprofiler_org)
  gp_xlab  <- if (exists("LABELS") && !is.null(LABELS$gprofiler_xlab)) LABELS$gprofiler_xlab else "Gene Ratio (Recall)"
  
  p_dot <- ggplot(top_terms, aes(x = Recall, y = reorder(term_name, -p_value), size = intersection_size, color = p_value)) + 
    geom_point(alpha = 0.8) + 
    scale_color_viridis_c(direction = -1, guide = guide_colorbar(reverse = TRUE)) + 
    facet_grid(source ~ ., scales = "free_y", space = "free_y") + 
    labs(title = gp_title, subtitle = gp_sub, x = gp_xlab, y = "", size = "Gene Count", color = "FDR") + 
    theme_bw(base_size = 11) + 
    theme(strip.text.y = element_text(angle = 0, face = "bold"), axis.text.y = element_text(size = 9))
  
  ggsave(file.path(gp_dir, paste0("DotPlot_", comp_name, ".png")), plot = p_dot, width = 11, height = 9, dpi = 300)
  return(res_tbl)
}

generate_interactive_sankey <- function(res_tbl, group_name, comp_name, out_dir, top_n_terms = 30) {
  if (is.null(res_tbl) || nrow(res_tbl) == 0 || !"intersection" %in% names(res_tbl)) return(NULL)
  top_terms <- head(res_tbl[order(res_tbl$p_value), ], top_n_terms)
  
  links_list <- list()
  for (i in seq_len(nrow(top_terms))) {
    term <- top_terms$term_name[i]
    genes <- unlist(strsplit(as.character(top_terms$intersection[i]), "[,;]+"))
    genes <- trimws(genes[genes != ""])
    if (length(genes) > 0) links_list[[i]] <- data.frame(source = genes, target = term, stringsAsFactors = FALSE)
  }
  link <- do.call(rbind, links_list)
  if (is.null(link) || nrow(link) == 0) return(NULL)
  
  link$value <- 1
  node_names <- unique(c(link$source, link$target))
  node <- data.frame(name = node_names, stringsAsFactors = FALSE)
  link$IDsource <- match(link$source, node$name) - 1
  link$IDtarget <- match(link$target, node$name) - 1
  node$color <- rainbow(nrow(node), s = 0.6, v = 0.9)
  
  src_colors <- node$color[link$IDsource + 1]
  link$color <- sapply(src_colors, function(hex) {
    rgb_val <- col2rgb(hex)
    sprintf("rgba(%d, %d, %d, 0.4)", rgb_val[1], rgb_val[2], rgb_val[3])
  })
  
  sankey_title <- if (exists("LABELS") && !is.null(LABELS$sankey_title)) sprintf(LABELS$sankey_title, comp_name) else sprintf("Gene-Pathway Flow (%s)", comp_name)
  plot_height  <- max(800, length(unique(link$source)) * 18)
  
  fig <- plotly::plot_ly(
    type = "sankey", orientation = "h", height = plot_height,
    node = list(label = node$name, color = node$color, pad = 10, thickness = 20, line = list(color = "black", width = 0.5)),
    link = list(source = link$IDsource, target = link$IDtarget, value = link$value, color = link$color)
  )
  fig <- plotly::layout(fig, title = sankey_title, font = list(size = 11))
  htmlwidgets::saveWidget(plotly::as_widget(fig), file.path(out_dir, paste0("Sankey_", comp_name, ".html")), selfcontained = TRUE)
}

generate_combined_gprofiler_dotplot <- function(all_gp_results_list, gp_dir, top_n_per_contrast = 5) {
  cat("\n  -> [g:Profiler] Generowanie zbiorczego wykresu szlaków dla wszystkich kontrastów...\n")
  if (!dir.exists(gp_dir)) dir.create(gp_dir, recursive = TRUE, showWarnings = FALSE)
  
  combined_list <- list()
  for (comp_name in names(all_gp_results_list)) {
    res_tbl <- all_gp_results_list[[comp_name]]
    if (!is.null(res_tbl) && is.data.frame(res_tbl) && nrow(res_tbl) > 0) {
      res_tbl$Contrast <- comp_name
      top_sub <- res_tbl %>% 
        dplyr::group_by(source) %>% 
        dplyr::top_n(top_n_per_contrast, wt = -p_value) %>% 
        dplyr::ungroup()
      combined_list[[comp_name]] <- top_sub
    }
  }
  
  if (length(combined_list) == 0) return(NULL)
  
  df_comb <- do.call(rbind, combined_list)
  df_comb$Recall <- df_comb$intersection_size / df_comb$term_size
  df_comb$term_name_wrapped <- stringr::str_wrap(df_comb$term_name, width = 42)
  
  p_comb <- ggplot(df_comb, aes(x = Contrast, y = reorder(term_name_wrapped, -p_value))) +
    geom_point(aes(size = intersection_size, color = p_value), alpha = 0.85) +
    scale_color_viridis_c(direction = -1, guide = guide_colorbar(reverse = TRUE)) +
    facet_grid(source ~ ., scales = "free_y", space = "free_y") +
    labs(
      title    = "Comparative Biological Pathways Across All Contrasts",
      subtitle = sprintf("Top %d Enriched Terms per Source | g:Profiler (FDR < %.2f)", top_n_per_contrast, STATS$alpha_thr),
      x        = "Experimental Contrast",
      y        = "",
      size     = "Gene Count",
      color    = "FDR (p-adj)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
      axis.text.y      = element_text(size = 8),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text.y     = element_text(angle = 0, face = "bold")
    )
  
  out_png <- file.path(gp_dir, "Combined_Pathways_All_Contrasts_DotPlot.png")
  ggsave(out_png, plot = p_comb, width = 13, height = max(9, 0.28 * length(unique(df_comb$term_name))), dpi = 300)
  cat(sprintf("     [+] Zapisano zbiorczy DotPlot szlaków: %s\n", out_png))
  return(p_comb)
}

# ==============================================================================
# BLOK 8: ANALIZA MIKROBIOMU (BRACKEN/KRAKEN) I INTEGRACJA MULTI-OMICS
# ==============================================================================

# --- 8.1. Kompozycja Taksonomiczna (Stacked Barplot) ---
plot_microbiome_composition <- function(microbiome_mat, colData, out_dir, top_n = 15) {
  cat("  -> [Microbiome] Generating Taxonomic Composition Plot (TSS %)...\n")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  sample_totals <- rowSums(microbiome_mat)
  valid_samples <- sample_totals > 0
  if (sum(valid_samples) < 2) {
    cat("     [!] Insufficient non-zero samples for composition plot.\n")
    return(NULL)
  }
  
  rel_mat <- microbiome_mat[valid_samples, , drop = FALSE] / sample_totals[valid_samples]
  mean_abund <- colMeans(rel_mat)
  top_n_taxa <- min(top_n, ncol(rel_mat))
  top_taxa <- names(sort(mean_abund, decreasing = TRUE))[1:top_n_taxa]
  
  other_taxa_cols <- setdiff(colnames(rel_mat), top_taxa)
  other_abund <- if (length(other_taxa_cols) > 0) rowSums(rel_mat[, other_taxa_cols, drop = FALSE]) else rep(0, nrow(rel_mat))
  
  other_name <- if (exists("LABELS") && !is.null(LABELS$micro_comp_other)) LABELS$micro_comp_other else "Other Taxa"
  plot_mat <- cbind(rel_mat[, top_taxa, drop = FALSE], setNames(data.frame(other_abund), other_name))
  
  df_long <- tidyr::pivot_longer(
    data.frame(plot_mat, Sample = rownames(plot_mat), check.names = FALSE),
    cols = -Sample, names_to = "Species", values_to = "Relative_Abundance"
  )
  
  df_long$Condition <- colData$Condition[match(df_long$Sample, rownames(colData))]
  if (exists("PALETTES") && !is.null(names(PALETTES$condition_colors))) {
    present_conds <- levels(factor(df_long$Condition))
    ordered_levels <- intersect(names(PALETTES$condition_colors), present_conds)
    df_long$Condition <- factor(df_long$Condition, levels = c(ordered_levels, setdiff(present_conds, ordered_levels)))
  }
  df_long$Species <- factor(df_long$Species, levels = c(top_taxa, other_name))
  
  n_tax <- length(top_taxa)
  base_cols <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", 
                 "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#AEC7E8", "#FFBB78")
  tax_colors <- if (n_tax <= length(base_cols)) {
    base_cols[1:n_tax]
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_tax)
  }
  palette_colors <- c(tax_colors, "#D3D3D3")
  names(palette_colors) <- c(top_taxa, other_name)
  
  title_str <- if (exists("LABELS") && !is.null(LABELS$micro_comp_title)) LABELS$micro_comp_title else "Microbiome Taxonomic Composition"
  sub_str   <- if (exists("LABELS") && !is.null(LABELS$micro_comp_sub)) sprintf(LABELS$micro_comp_sub, top_n_taxa) else sprintf("Top %d Taxa", top_n_taxa)
  xlab_str  <- if (exists("LABELS") && !is.null(LABELS$micro_comp_xlab)) LABELS$micro_comp_xlab else "Samples"
  ylab_str  <- if (exists("LABELS") && !is.null(LABELS$micro_comp_ylab)) LABELS$micro_comp_ylab else "Relative Abundance (%)"
  fill_str  <- if (exists("LABELS") && !is.null(LABELS$micro_comp_fill)) LABELS$micro_comp_fill else "Species"
  
  p <- ggplot(df_long, aes(x = Sample, y = Relative_Abundance * 100, fill = Species)) +
    geom_bar(stat = "identity", width = 0.85, color = "black", linewidth = 0.2) +
    facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = palette_colors) +
    labs(title = title_str, subtitle = sub_str, x = xlab_str, y = ylab_str, fill = fill_str) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text       = element_text(face = "bold"),
      legend.text      = element_text(size = 8, face = "italic"),
      legend.position  = "right"
    )
  
  out_png <- file.path(out_dir, "Microbiome_Composition_Barplot.png")
  ggsave(out_png, plot = p, width = 12, height = 7, dpi = 300)
  cat(sprintf("     [+] Saved composition plot: %s\n", out_png))
  return(p)
}

# --- 8.2. Alfa-Różnorodność ze statystyką (Kruskal-Wallis) ---
calculate_and_plot_alpha_diversity <- function(microbiome_mat, colData, out_dir) {
  cat("  -> [Microbiome] Calculating Alpha Diversity Metrics & Statistical Testing...\n")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  keep <- rowSums(microbiome_mat) > 0
  mat_k <- microbiome_mat[keep, , drop = FALSE]
  cd_k  <- colData[keep, , drop = FALSE]
  if (nrow(mat_k) < 2) return(NULL)
  
  richness <- rowSums(mat_k > 0)
  rel_mat  <- mat_k / rowSums(mat_k)
  
  shannon <- apply(rel_mat, 1, function(p) {
    p <- p[p > 0]
    -sum(p * log(p))
  })
  simpson  <- apply(rel_mat, 1, function(p) 1 - sum(p^2))
  evenness <- ifelse(richness > 1, shannon / log(richness), 0)
  
  alpha_df <- data.frame(
    Sample            = rownames(cd_k),
    Condition         = cd_k$Condition,
    Observed_Richness = richness,
    Shannon_Index     = shannon,
    Simpson_Index     = simpson,
    Pielous_Evenness  = evenness,
    stringsAsFactors  = FALSE
  )
  
  metric_cols <- c("Observed_Richness", "Shannon_Index", "Simpson_Index", "Pielous_Evenness")
  kw_results <- list()
  
  for (m in metric_cols) {
    kw_test <- tryCatch(kruskal.test(alpha_df[[m]] ~ alpha_df$Condition), error = function(e) NULL)
    kw_results[[m]] <- if (!is.null(kw_test)) kw_test$p.value else NA
  }
  
  stats_summary <- data.frame(
    Metric        = metric_cols,
    KW_Chi2       = sapply(metric_cols, function(m) {
      kw <- tryCatch(kruskal.test(alpha_df[[m]] ~ alpha_df$Condition), error = function(e) NULL)
      if (!is.null(kw)) unname(kw$statistic) else NA
    }),
    KW_P_Value    = sapply(metric_cols, function(m) kw_results[[m]]),
    KW_P_Adj_BH   = p.adjust(sapply(metric_cols, function(m) kw_results[[m]]), method = "BH"),
    stringsAsFactors = FALSE
  )
  
  write.csv(alpha_df, file.path(out_dir, "Microbiome_Alpha_Diversity_Metrics.csv"), row.names = FALSE)
  write.csv(stats_summary, file.path(out_dir, "Microbiome_Alpha_Diversity_Kruskal_Stats.csv"), row.names = FALSE)
  
  alpha_long <- tidyr::pivot_longer(alpha_df, cols = dplyr::all_of(metric_cols), names_to = "Metric", values_to = "Value")
  
  metric_dict <- if (exists("LABELS") && !is.null(LABELS$alpha_metrics)) {
    LABELS$alpha_metrics
  } else {
    c("Observed_Richness" = "Observed Richness (S)", "Shannon_Index" = "Shannon Index (H')",
      "Simpson_Index" = "Simpson Index (1 - D)", "Pielous_Evenness" = "Pielou's Evenness (J')")
  }
  
  facet_labels <- sapply(metric_cols, function(m) {
    base_lbl <- metric_dict[m]
    pval <- kw_results[[m]]
    if (!is.na(pval)) sprintf("%s\n(KW p = %.3f)", base_lbl, pval) else base_lbl
  })
  
  alpha_long$Metric <- factor(alpha_long$Metric, levels = metric_cols, labels = facet_labels)
  cond_palette <- if (exists("PALETTES") && !is.null(PALETTES$condition_colors)) PALETTES$condition_colors else NULL
  
  p <- ggplot(alpha_long, aes(x = Condition, y = Value, fill = Condition)) +
    geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.2, size = 2.5, shape = 21, color = "black", alpha = 0.9) +
    facet_wrap(~ Metric, scales = "free_y", ncol = 4) +
    (if (!is.null(cond_palette)) scale_fill_manual(values = cond_palette) else list()) +
    labs(
      title    = if (exists("LABELS") && !is.null(LABELS$alpha_div_title)) LABELS$alpha_div_title else "Alpha Diversity",
      subtitle = if (exists("LABELS") && !is.null(LABELS$alpha_div_sub)) LABELS$alpha_div_sub else "Across Groups",
      x        = BIO$condition_label,
      y        = if (exists("LABELS") && !is.null(LABELS$alpha_div_ylab)) LABELS$alpha_div_ylab else "Metric Value"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position  = "none",
      axis.text.x      = element_text(angle = 45, hjust = 1, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text       = element_text(face = "bold")
    )
  
  out_png <- file.path(out_dir, "Microbiome_Alpha_Diversity_Boxplots.png")
  ggsave(out_png, plot = p, width = 12, height = 5, dpi = 300)
  cat(sprintf("     [+] Saved Alpha Diversity plot: %s\n", out_png))
  return(list(metrics = alpha_df, stats = stats_summary))
}

# --- 8.3. Beta-Różnorodność (PCoA Bray-Curtis + PERMANOVA) ---
run_beta_diversity_pcoa <- function(microbiome_mat, colData, out_dir) {
  cat("  -> [Microbiome] Analyzing Beta Diversity (Bray-Curtis PCoA + PERMANOVA)...\n")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  keep <- rowSums(microbiome_mat) > 0
  mat_k <- microbiome_mat[keep, , drop = FALSE]
  cd_k  <- colData[keep, , drop = FALSE]
  if (nrow(mat_k) < 3) return(NULL)
  
  rel_mat <- mat_k / rowSums(mat_k)
  permanova_res <- "Bray-Curtis Distances"
  
  if (requireNamespace("vegan", quietly = TRUE)) {
    bc_dist <- vegan::vegdist(rel_mat, method = "bray")
    tryCatch({
      ad <- vegan::adonis2(bc_dist ~ Condition, data = as.data.frame(cd_k), permutations = 999)
      permanova_res <- sprintf("PERMANOVA: R² = %.3f, p = %.4f", ad$R2[1], ad$`Pr(>F)`[1])
    }, error = function(e) { permanova_res <- "PERMANOVA: N/A" })
  } else {
    n_samp <- nrow(rel_mat)
    d_mat <- matrix(0, n_samp, n_samp, dimnames = list(rownames(rel_mat), rownames(rel_mat)))
    for (i in 1:n_samp) {
      for (j in 1:n_samp) {
        sum_abs <- sum(abs(rel_mat[i, ] - rel_mat[j, ]))
        sum_tot <- sum(rel_mat[i, ] + rel_mat[j, ])
        d_mat[i, j] <- if (sum_tot > 0) sum_abs / sum_tot else 0
      }
    }
    bc_dist <- as.dist(d_mat)
  }
  
  pcoa <- cmdscale(bc_dist, k = 2, eig = TRUE)
  eig_vals <- pcoa$eig[pcoa$eig > 0]
  var_explained <- if (length(eig_vals) >= 2) round(100 * eig_vals[1:2] / sum(eig_vals), 1) else c(0, 0)
  
  present_conds <- as.character(cd_k$Condition)
  ordered_levels <- if (exists("PALETTES") && !is.null(names(PALETTES$condition_colors))) {
    intersect(names(PALETTES$condition_colors), unique(present_conds))
  } else {
    unique(present_conds)
  }
  
  pcoa_df <- data.frame(
    Sample    = rownames(cd_k),
    Condition = factor(present_conds, levels = ordered_levels),
    PCoA1     = pcoa$points[, 1],
    PCoA2     = pcoa$points[, 2],
    stringsAsFactors = FALSE
  )
  
  cond_palette <- if (exists("PALETTES") && !is.null(names(PALETTES$condition_colors))) {
    PALETTES$condition_colors[names(PALETTES$condition_colors) %in% levels(pcoa_df$Condition)]
  } else {
    NULL
  }
  
  p <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, fill = Condition)) +
    stat_ellipse(geom = "polygon", alpha = 0.15, level = 0.8, color = NA) +
    geom_point(size = 4, alpha = 0.9, shape = 21, color = "black") +
    ggrepel::geom_text_repel(aes(label = Sample), size = 3.5, show.legend = FALSE, color = "black") +
    (if (!is.null(cond_palette)) scale_fill_manual(values = cond_palette) else list()) +
    labs(
      title    = if (exists("LABELS") && !is.null(LABELS$beta_pcoa_title)) LABELS$beta_pcoa_title else "Beta Diversity (PCoA)",
      subtitle = permanova_res,
      x        = if (exists("LABELS") && !is.null(LABELS$beta_pcoa_xlab)) sprintf(LABELS$beta_pcoa_xlab, var_explained[1]) else sprintf("PCoA 1 [%.1f%%]", var_explained[1]),
      y        = if (exists("LABELS") && !is.null(LABELS$beta_pcoa_ylab)) sprintf(LABELS$beta_pcoa_ylab, var_explained[2]) else sprintf("PCoA 2 [%.1f%%]", var_explained[2]),
      fill     = BIO$condition_label
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right", panel.grid.minor = element_blank())
  
  out_png <- file.path(out_dir, "Microbiome_Beta_Diversity_PCoA.png")
  ggsave(out_png, plot = p, width = 8.5, height = 6.5, dpi = 300)
  cat(sprintf("     [+] Saved Beta Diversity PCoA: %s\n", out_png))
  return(list(pcoa_df = pcoa_df, permanova = permanova_res))
}

# --- 8.4. Korelacje Krzyżowe Multi-Omics z FDR (BH) ---
run_host_microbiome_cross_correlation <- function(vsd_clean, microbiome_mat, top_deg_ids, anno_map, colData, out_dir, top_n_taxa = 20) {
  cat("  -> [Multi-Omics] Analyzing Cross-Correlations (Host DEGs vs Microbiome with FDR BH)...\n")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  common_samples <- intersect(colnames(vsd_clean), rownames(microbiome_mat))
  if (length(common_samples) < 3) {
    cat("     [!] Insufficient common samples (< 3) for correlation.\n")
    return(NULL)
  }
  
  valid_gene_ids <- intersect(top_deg_ids, rownames(vsd_clean))
  if (length(valid_gene_ids) < 3) {
    rv <- matrixStats::rowVars(assay(vsd_clean)[, common_samples, drop = FALSE])
    valid_gene_ids <- head(rownames(vsd_clean)[order(rv, decreasing = TRUE)], 30)
  }
  
  host_mat <- t(assay(vsd_clean)[valid_gene_ids, common_samples, drop = FALSE])
  gene_symbols <- anno_map$symbol[match(colnames(host_mat), anno_map$gene_id)]
  clean_syms <- ifelse(is.na(gene_symbols) | gene_symbols == "", colnames(host_mat), gene_symbols)
  colnames(host_mat) <- make.unique(clean_syms)
  
  gene_vars <- matrixStats::colVars(host_mat)
  host_mat  <- host_mat[, gene_vars > 1e-6, drop = FALSE]
  
  micro_clean <- microbiome_mat[common_samples, , drop = FALSE]
  sample_totals <- rowSums(micro_clean)
  non_zero_rows <- sample_totals > 0
  
  if (sum(non_zero_rows) < 3 || ncol(host_mat) < 2) return(NULL)
  
  micro_clean <- micro_clean[non_zero_rows, , drop = FALSE]
  host_mat    <- host_mat[non_zero_rows, , drop = FALSE]
  
  micro_rel <- micro_clean / rowSums(micro_clean)
  n_tax <- min(top_n_taxa, ncol(micro_rel))
  top_taxa <- names(sort(colMeans(micro_rel), decreasing = TRUE))[1:n_tax]
  micro_mat <- micro_rel[, top_taxa, drop = FALSE]
  
  tax_vars  <- matrixStats::colVars(micro_mat)
  micro_mat <- micro_mat[, tax_vars > 1e-8, drop = FALSE]
  
  n_genes <- ncol(host_mat)
  n_micro <- ncol(micro_mat)
  
  cor_matrix  <- matrix(NA, nrow = n_genes, ncol = n_micro, dimnames = list(colnames(host_mat), colnames(micro_mat)))
  pval_matrix <- matrix(NA, nrow = n_genes, ncol = n_micro, dimnames = list(colnames(host_mat), colnames(micro_mat)))
  
  for (i in 1:n_genes) {
    for (j in 1:n_micro) {
      test <- suppressWarnings(cor.test(host_mat[, i], micro_mat[, j], method = "spearman", exact = FALSE))
      cor_matrix[i, j]  <- test$estimate
      pval_matrix[i, j] <- test$p.value
    }
  }
  
  qval_matrix <- matrix(
    p.adjust(as.vector(pval_matrix), method = "BH"),
    nrow = n_genes, ncol = n_micro,
    dimnames = list(colnames(host_mat), colnames(micro_mat))
  )
  
  sig_stars <- matrix("", nrow = n_genes, ncol = n_micro)
  sig_stars[qval_matrix < 0.05]  <- "*"
  sig_stars[qval_matrix < 0.01]  <- "**"
  sig_stars[qval_matrix < 0.001] <- "***"
  
  df_cor  <- data.frame(Gene = rownames(cor_matrix), cor_matrix, check.names = FALSE)
  df_pval <- data.frame(Gene = rownames(pval_matrix), pval_matrix, check.names = FALSE)
  df_qval <- data.frame(Gene = rownames(qval_matrix), qval_matrix, check.names = FALSE)
  
  wb_cor <- openxlsx2::wb_workbook()
  
  wb_cor$add_worksheet("Spearman_r")$add_data("Spearman_r", x = df_cor)
  h1 <- openxlsx2::wb_dims(rows = 1, cols = 1:ncol(df_cor))
  wb_cor$add_fill("Spearman_r", dims = h1, color = openxlsx2::wb_color(hex = "4F81BD"))
  wb_cor$add_font("Spearman_r", dims = h1, color = openxlsx2::wb_color(hex = "FFFFFF"), bold = TRUE)
  wb_cor$set_col_widths("Spearman_r", cols = 1:ncol(df_cor), widths = "auto")
  
  wb_cor$add_worksheet("Raw_P_Values")$add_data("Raw_P_Values", x = df_pval)
  h2 <- openxlsx2::wb_dims(rows = 1, cols = 1:ncol(df_pval))
  wb_cor$add_fill("Raw_P_Values", dims = h2, color = openxlsx2::wb_color(hex = "4F81BD"))
  wb_cor$add_font("Raw_P_Values", dims = h2, color = openxlsx2::wb_color(hex = "FFFFFF"), bold = TRUE)
  wb_cor$set_col_widths("Raw_P_Values", cols = 1:ncol(df_pval), widths = "auto")
  
  wb_cor$add_worksheet("FDR_Q_Values")$add_data("FDR_Q_Values", x = df_qval)
  h3 <- openxlsx2::wb_dims(rows = 1, cols = 1:ncol(df_qval))
  wb_cor$add_fill("FDR_Q_Values", dims = h3, color = openxlsx2::wb_color(hex = "4F81BD"))
  wb_cor$add_font("FDR_Q_Values", dims = h3, color = openxlsx2::wb_color(hex = "FFFFFF"), bold = TRUE)
  wb_cor$set_col_widths("FDR_Q_Values", cols = 1:ncol(df_qval), widths = "auto")
  
  xlsx_out <- file.path(out_dir, "Host_Microbiome_Spearman_Correlations.xlsx")
  wb_cor$save(xlsx_out, overwrite = TRUE)
  
  main_title <- if (exists("LABELS") && !is.null(LABELS$cross_cor_title)) LABELS$cross_cor_title else "Host Expression vs Microbiome Abundance (FDR adjusted)"
  png_out <- file.path(out_dir, "Host_Microbiome_Correlation_Heatmap.png")
  
  png(png_out, width = 11, height = max(8, 0.22 * n_genes), units = "in", res = 300)
  pheatmap::pheatmap(
    cor_matrix,
    display_numbers = sig_stars,
    fontsize_number = 9,
    number_color    = "black",
    color           = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
    breaks          = seq(-1, 1, length.out = 101),
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    main            = main_title,
    fontsize_row    = 8,
    fontsize_col    = 8,
    angle_col       = 45
  )
  dev.off()
  
  cat(sprintf("     [+] Saved correlation matrices: %s and heatmap: %s\n", xlsx_out, png_out))
  return(list(cor = cor_matrix, pval = pval_matrix, qval = qval_matrix))
}


###################################

##############################
#' Dynamiczny Generator Głównego Manifestu Multiomicznego Próbek
#' Kompatybilny z modułami 03-08 potoku RNA-seq / Multi-Omics
generate_master_coldata <- function(dds           = NULL, 
                                    df_qc         = NULL, 
                                    df_microbiome = NULL, 
                                    df_dtu        = NULL,
                                    df_variants   = NULL, 
                                    df_adar       = NULL,
                                    out_dir       = PATHS$out_base_dir) {
  
  # 1. Inicjalizacja bazy próbek (Sample_ID i Condition)
  if (!is.null(dds)) {
    cd_df <- as.data.frame(SummarizedExperiment::colData(dds))
    sample_ids <- rownames(cd_df)
    cond_vec   <- as.character(cd_df$Condition)
  } else if (exists("colData_Full")) {
    sample_ids <- rownames(colData_Full)
    cond_vec   <- as.character(colData_Full$Condition)
  } else {
    stop("Brak obiektu dds lub colData_Full do zainicjalizowania listy próbek!")
  }
  
  master_df <- data.frame(
    Sample_ID = sample_ids,
    Condition = cond_vec,
    stringsAsFactors = FALSE
  )
  
  # 2. Moduł DESeq2 / Ekspresji
  if (!is.null(dds)) {
    sf <- DESeq2::sizeFactors(dds)
    if (!is.null(sf)) {
      master_df$DESeq2_Size_Factor <- round(sf[master_df$Sample_ID], 3)
    }
    raw_counts <- tryCatch(DESeq2::counts(dds, normalized = FALSE), error = function(e) NULL)
    if (!is.null(raw_counts)) {
      master_df$Active_Genes_Count <- colSums(raw_counts[, master_df$Sample_ID, drop = FALSE] >= 10)
    }
  }
  
  # 3. Moduł Mapowania & QC (Krok 03 - Comprehensive_Mapping_and_QC_Summary)
  if (!is.null(df_qc) && is.data.frame(df_qc) && nrow(df_qc) > 0) {
    qc_clean <- df_qc
    if (!"Sample_ID" %in% colnames(qc_clean) && "Sample" %in% colnames(qc_clean)) {
      qc_clean$Sample_ID <- qc_clean$Sample
    }
    
    # Obsługa przeliczeń na miliony odczytów
    if ("Salmon_Mapped_Reads" %in% colnames(qc_clean)) {
      qc_clean$Mapped_Reads_M <- round(qc_clean$Salmon_Mapped_Reads / 1e6, 2)
    } else if ("Mapped_Reads" %in% colnames(qc_clean)) {
      qc_clean$Mapped_Reads_M <- round(qc_clean$Mapped_Reads / 1e6, 2)
    }
    
    if ("Mapping_Rate_Percent" %in% colnames(qc_clean)) {
      qc_clean$Mapping_Rate_Pct <- qc_clean$Mapping_Rate_Percent
    }
    
    target_cols <- intersect(c("Mapped_Reads_M", "Mapping_Rate_Pct", "GC_Content_Percent", "Q30_Bases_Percent"), colnames(qc_clean))
    if (length(target_cols) > 0 && "Sample_ID" %in% colnames(qc_clean)) {
      master_df <- dplyr::left_join(master_df, qc_clean[, c("Sample_ID", target_cols)], by = "Sample_ID")
    }
  }
  
  # 4. Moduł Mikrobiomu (Krok 06 - Alpha Diversity)
  if (!is.null(df_microbiome) && is.data.frame(df_microbiome) && nrow(df_microbiome) > 0) {
    micro_clean <- df_microbiome
    if (!"Sample_ID" %in% colnames(micro_clean) && "Sample" %in% colnames(micro_clean)) {
      micro_clean$Sample_ID <- micro_clean$Sample
    }
    if ("Shannon_Index" %in% colnames(micro_clean)) {
      micro_clean$Shannon_Diversity <- round(micro_clean$Shannon_Index, 2)
    }
    
    target_cols <- intersect(c("Shannon_Diversity", "Observed_Richness", "Simpson_Index"), colnames(micro_clean))
    if (length(target_cols) > 0 && "Sample_ID" %in% colnames(micro_clean)) {
      master_df <- dplyr::left_join(master_df, micro_clean[, c("Sample_ID", target_cols)], by = "Sample_ID")
    }
  }
  
  # 5. Moduł Splicingu / DTU (Krok 07 - opcjonalny na przyszłość)
  if (!is.null(df_dtu) && is.data.frame(df_dtu) && nrow(df_dtu) > 0) {
    dtu_clean <- df_dtu
    if (!"Sample_ID" %in% colnames(dtu_clean) && "Sample" %in% colnames(dtu_clean)) {
      dtu_clean$Sample_ID <- dtu_clean$Sample
    }
    target_cols <- intersect(c("DTU_Genes_Count", "Isoform_Switches_Count"), colnames(dtu_clean))
    if (length(target_cols) > 0 && "Sample_ID" %in% colnames(dtu_clean)) {
      master_df <- dplyr::left_join(master_df, dtu_clean[, c("Sample_ID", target_cols)], by = "Sample_ID")
    }
  }
  
  # 6. Moduł Wariantów i ADAR (Krok 08 - opcjonalny na przyszłość)
  if (!is.null(df_variants) && is.data.frame(df_variants) && nrow(df_variants) > 0) {
    var_clean <- df_variants
    if (!"Sample_ID" %in% colnames(var_clean) && "Sample" %in% colnames(var_clean)) {
      var_clean$Sample_ID <- var_clean$Sample
    }
    target_cols <- intersect(c("Total_PASS_Vars", "Functional_Vars", "Ti_Tv_Ratio"), colnames(var_clean))
    if (length(target_cols) > 0 && "Sample_ID" %in% colnames(var_clean)) {
      master_df <- dplyr::left_join(master_df, var_clean[, c("Sample_ID", target_cols)], by = "Sample_ID")
    }
  }
  
  # 6.1. ADAR A>G % (z tabeli spektrum mutacji w kroku 08)
  if (!is.null(df_adar) && is.data.frame(df_adar) && nrow(df_adar) > 0) {
    adar_clean <- df_adar %>% 
      dplyr::filter(grepl("ADAR", Type)) %>% 
      dplyr::select(Sample_ID, ADAR_Rate_Pct = Percentage)
    if (nrow(adar_clean) > 0) {
      master_df <- dplyr::left_join(master_df, adar_clean, by = "Sample_ID")
    }
  }
  
  # 7. Wielopoziomowa Flaga Jakości (QC Status)
  flags <- rep("PASS", nrow(master_df))
  
  if ("DESeq2_Size_Factor" %in% colnames(master_df)) {
    sf_bad <- master_df$DESeq2_Size_Factor < 0.3 | master_df$DESeq2_Size_Factor > 3.0
    flags[sf_bad] <- paste0(flags[sf_bad], "; WARN(SizeFactor)")
  }
  if ("Mapping_Rate_Pct" %in% colnames(master_df)) {
    map_bad <- master_df$Mapping_Rate_Pct < 60.0
    flags[map_bad] <- paste0(flags[map_bad], "; WARN(LowMapping)")
  }
  
  master_df$QC_Status <- gsub("^PASS; ", "", flags)
  
  # 8. Eksport do CSV oraz formatowanego Excela
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_csv  <- file.path(out_dir, "Comprehensive_Master_Sample_Metadata.csv")
  out_xlsx <- file.path(out_dir, "Comprehensive_Master_Sample_Metadata.xlsx")
  
  write.csv(master_df, out_csv, row.names = FALSE)
  
  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet("Sample_Master_Manifest")$add_data(sheet = "Sample_Master_Manifest", x = master_df)
  h_dim <- openxlsx2::wb_dims(rows = 1, cols = 1:ncol(master_df))
  wb$add_fill(sheet = "Sample_Master_Manifest", dims = h_dim, color = openxlsx2::wb_color(hex = "1F497D"))
  wb$add_font(sheet = "Sample_Master_Manifest", dims = h_dim, color = openxlsx2::wb_color(hex = "FFFFFF"), bold = TRUE)
  wb$set_col_widths(sheet = "Sample_Master_Manifest", cols = 1:ncol(master_df), widths = "auto")
  wb$save(out_xlsx, overwrite = TRUE)
  
  cat(sprintf("  [+] Wygenerowano Master Manifest (%d próbek, %d kolumn): %s\n", 
              nrow(master_df), ncol(master_df), basename(out_xlsx)))
  return(master_df)
}
