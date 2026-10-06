# ==============================================================================
# 04_DESEQ2_ANALYSIS.R - Silnik różnic ekspresji, Wizualizacje, Szlaki i inne cuda
# ==============================================================================
cat("\n=================================================================\n")
cat(" >>> KROK 04: ANALIZA RÓŻNICOWA DESEQ2, WIZUALIZACJE I SZLAKI <<<\n")
cat("=================================================================\n")

# --- 1. Weryfikacja środowiska i SSOT ---
if (!exists("PATHS") || !exists("BIO") || !exists("STATS") || !exists("PALETTES") || !exists("LABELS")) {
  if (file.exists("01_config.R")) {
    source("01_config.R")
  } else if (file.exists("01_configGRCz11.R")) {
    source("01_configGRCz11.R")
  } else {
    stop("CRITICAL ERROR: Brak załadowanego pliku 01_config.R!")
  }
}

utils_path <- file.path(PATHS$scripts_dir, PATHS$utils_file)
if (file.exists(utils_path)) {
  source(utils_path)
} else if (file.exists("02_utils.R")) {
  source("02_utils.R")
} else {
  stop("CRITICAL ERROR: Brak pliku 02_utils.R!")
}

# --- 2. Wczytanie danych z kwantyfikacji Salmona i słownika adnotacji ---
rds_path <- if (!is.null(PATHS$txi_rds)) PATHS$txi_rds else file.path(PATHS$out_base_dir, "txi_salmon_data.rds")
if (!file.exists(rds_path)) {
  stop(sprintf("CRITICAL ERROR: Brak pliku %s! Uruchom najpierw krok 03_salmonLoad.R.", rds_path))
}

data_bundle  <- readRDS(rds_path)
txi_full     <- data_bundle$txi
colData_Full <- data_bundle$colData

if (!exists("anno_map") || !is.data.frame(anno_map) || nrow(anno_map) == 0) {
  cat("[ADNOTACJE] Budowanie słownika genów z GTF...\n")
  anno_map <- build_custom_annotation_map(PATHS$merged_gtf, ref_gtf_path = PATHS$ref_gtf)
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

# --- 3. Główna pętla analityczna po zdefiniowanych tkankach/zbiorach ---
for (group in BIO$target_tissues) {
  cat(sprintf("\n >>> ANALIZA DLA GRUPY: [%s] <<<\n", toupper(group)))
  group_dir <- setup_group_folders(PATHS$out_base_dir, group)
  
  samples_group <- if (group == "Whole_Dataset") rownames(colData_Full) else rownames(colData_Full)[colData_Full$Tissue == group]
  if (length(samples_group) == 0) next
  
  colData_group <- droplevels(colData_Full[samples_group, ])
  txi_group     <- subset_txi(txi_full, samples_group)
  
  # --- Dynamiczny dobór pre-filteringu ekspresji (SSOT STATS) ---
  if (exists("STATS") && identical(STATS$prefilter_mode, "strict")) {
    min_cnt <- if (!is.null(STATS$min_count_thr)) STATS$min_count_thr else 10
    min_smp <- if (!is.null(STATS$min_samples_thr) && is.numeric(STATS$min_samples_thr)) {
      STATS$min_samples_thr
    } else {
      min(table(colData_group$Condition)) # Rozmiar najmniejszej grupy badanej
    }
    
    keep <- rowSums(txi_group$counts >= min_cnt) >= min_smp
    cat(sprintf("  [FILTER] Tryb: STRICT (>= %d counts w >= %d probkach) | Pozostalo: %d genow.\n", 
                min_cnt, min_smp, sum(keep)))
  } else {
    # Tryb klasyczny (legacy)
    keep <- rowSums(txi_group$counts) > 0
    cat(sprintf("  [FILTER] Tryb: LEGACY (rowSums > 0) | Pozostalo: %d genow.\n", sum(keep)))
  }
  
  txi_group$abundance <- txi_group$abundance[keep, , drop = FALSE]
  txi_group$counts    <- txi_group$counts[keep, , drop = FALSE]
  txi_group$length    <- txi_group$length[keep, , drop = FALSE]
  
  design_formula  <- ~ Condition
  group_contrasts <- CONTRASTS_CONFIG[[group]]
  
  # Dynamiczny dobór wariantów analizy (SSOT)
  analysis_variants <- if (isTRUE(STATS$enable_filtering)) c("Unfiltered", "Filtered") else c("Standard")
  last_res_list <- list()
  
  for (cond_variant in analysis_variants) {
    cat(sprintf("\n--- Przebieg: [%s] ---\n", cond_variant))
    
    if (cond_variant == "Filtered") {
      best_samples  <- select_best_n_samples(txi_group, colData_group, BIO$samples_per_grp, design_formula)
      colData_run   <- colData_group[best_samples, ]
      txi_run       <- subset_txi(txi_group, best_samples)
    } else {
      colData_run   <- colData_group
      txi_run       <- txi_group
    }
    
    # 3.1. Uruchomienie modelu DESeq2 i transformacji VST
    dds <- DESeqDataSetFromTximport(txi_run, colData = colData_run, design = design_formula)
    dds <- DESeq(dds, quiet = TRUE)
    vsd <- vst(dds, blind = FALSE)
    
    # 3.2. Wykresy Kontroli Jakości (QC)
    generate_QC_plots(dds, vsd, group_dir, group, condition_type = cond_variant)
    
    # 3.3. Pętla analizy poszczególnych kontrastów
    summary_df <- data.frame()
    sig_data_list <- list()
    res_list <- list()
    
    for (comp_name in names(group_contrasts)) {
      cfg <- group_contrasts[[comp_name]]
      res_out <- run_contrast_analysis(
        dds            = dds, 
        contrast_vec   = cfg$vec, 
        group_name     = group, 
        group_dir      = group_dir, 
        condition_type = cond_variant, 
        anno_map       = anno_map, 
        alpha_val      = cfg$alpha, 
        lfc_val        = cfg$lfc
      )
      summary_df <- rbind(summary_df, res_out$summary)
      sig_data_list[[length(sig_data_list) + 1]] <- res_out$sig_data
      res_list[[comp_name]] <- res_out
    }
    
    all_sig_ids <- unique(do.call(rbind, sig_data_list)$gene_id)
    
    # 3.4. Kompleks wielowariantowych heatmap (Top 50, 75, 100, 150, All)
    generate_master_heatmaps(vsd, all_sig_ids, group, group_dir, condition_type = cond_variant, anno_map = anno_map)
    
    # 3.5. Kompleks wykresów UpSet (All / Up / Down) + wykres słupkowy z liczbami
    generate_upset_and_contrast_breakdowns(sig_data_list, res_list, group_dir, group, cond_variant)
    
    # 3.6. Wykresy trajektorii ekspresji Top DEG
    generate_top_genes_trajectories(vsd, res_list, group_dir, group, cond_variant, anno_map = anno_map, top_n_genes = 16)
    
    # 3.7. Zapis raportów podsumowujących (CSV + stylizowany XLSX przez openxlsx2)
    summary_out_path <- file.path(group_dir, "significant", cond_variant, paste0("Summary_DE_", group, "_", cond_variant, ".csv"))
    write.csv(summary_df, summary_out_path, row.names = FALSE)
    
    wb_sum <- openxlsx2::wb_workbook()
    wb_sum$add_worksheet("DE_Summary")
    wb_sum$add_data(sheet = "DE_Summary", x = summary_df)
    h_dims <- openxlsx2::wb_dims(rows = 1, cols = 1:ncol(summary_df))
    wb_sum$add_fill(sheet = "DE_Summary", dims = h_dims, color = openxlsx2::wb_color(hex = "4F81BD"))
    wb_sum$add_font(sheet = "DE_Summary", dims = h_dims, color = openxlsx2::wb_color(hex = "FFFFFF"), bold = TRUE)
    wb_sum$set_col_widths(sheet = "DE_Summary", cols = 1:ncol(summary_df), widths = "auto")
    wb_sum$save(file.path(group_dir, "Excel_Reports", paste0("Summary_DE_Report_", group, "_", cond_variant, ".xlsx")), overwrite = TRUE)
    
    export_to_master_excel(res_list, group_dir, group, condition_type = cond_variant)
    last_res_list <- res_list
  }
  
  # ==========================================================================
  # 3.7b. SCALENIE I ZAPIS WSZYSTKICH KONTRASTÓW DO JEDNYCH PLIKÓW (FULL & SIG)
  # ==========================================================================
  cat("  -> Generowanie zbiorczych tabel ze wszystkich kontrastów (FULL & SIG)...\n")
  
  # 1. Złożenie tabeli wszystkich genów (FULL)
  all_full_list <- lapply(names(res_list), function(comp) {
    df <- res_list[[comp]]$full_df
    parts <- strsplit(comp, "_vs_")[[1]]
    df$Contrast           <- comp
    df$Target_Condition   <- parts[1]
    df$Baseline_Condition <- parts[2]
    df$Status <- ifelse(is.na(df$padj) | df$padj > STATS$alpha_thr | abs(df$log2FoldChange) < STATS$lfc_thr, 
                        "Not_Sig", 
                        ifelse(df$log2FoldChange > 0, "Up", "Down"))
    # Uporządkowanie kolejności kolumn (metadane kontrastu z przodu)
    cols_front <- c("Contrast", "Target_Condition", "Baseline_Condition", "Status", "gene_id", "symbol")
    df[, c(cols_front, setdiff(colnames(df), cols_front))]
  })
  master_full_df <- do.call(rbind, all_full_list)
  
  # 2. Złożenie tabeli wyłącznie istotnych DEG (SIG)
  all_sig_list <- lapply(names(res_list), function(comp) {
    df <- res_list[[comp]]$sig_df
    if (nrow(df) > 0) {
      parts <- strsplit(comp, "_vs_")[[1]]
      df$Contrast           <- comp
      df$Target_Condition   <- parts[1]
      df$Baseline_Condition <- parts[2]
      df$Status             <- ifelse(df$log2FoldChange > 0, "Up", "Down")
      cols_front <- c("Contrast", "Target_Condition", "Baseline_Condition", "Status", "gene_id", "symbol")
      df[, c(cols_front, setdiff(colnames(df), cols_front))]
    } else {
      NULL
    }
  })
  master_sig_df <- do.call(rbind, all_sig_list)
  
  # 3. Zapis do plików CSV w folderze significant/[cond_variant]
  out_dir_sig <- file.path(group_dir, "significant", cond_variant)
  
  csv_master_full <- file.path(out_dir_sig, paste0("ALL_Contrasts_", group, "_", cond_variant, "_shrink_full.csv"))
  csv_master_sig  <- file.path(out_dir_sig, paste0("ALL_Contrasts_", group, "_", cond_variant, "_shrink_sig.csv"))
  
  write.csv(master_full_df, csv_master_full, row.names = FALSE)
  write.csv(master_sig_df,  csv_master_sig,  row.names = FALSE)
  
  cat(sprintf("     [+] Zapisano master FULL: %s (%d wierszy)\n", basename(csv_master_full), nrow(master_full_df)))
  cat(sprintf("     [+] Zapisano master SIG:  %s (%d wierszy)\n", basename(csv_master_sig),  nrow(master_sig_df)))
  
  # --- ETAP III: G:PROFILER, SANKEY ORAZ ZBIORCZY DOTPLOT DLA WSZYSTKICH KONTRASTÓW ---
  cat(sprintf("\n >>> [%s] ANALIZA SZLAKÓW BIOLOGICZNYCH (g:Profiler) <<<\n", toupper(group)))
  gp_dir     <- file.path(group_dir, "Visualizations", "gProfiler")
  sankey_dir <- file.path(group_dir, "Visualizations", "Sankey", "Standard")
  
  all_gp_results <- list()
  for (comp_name in names(last_res_list)) {
    cfg    <- group_contrasts[[comp_name]]
    df_sig <- last_res_list[[comp_name]]$sig_df
    
    if (!is.null(df_sig) && nrow(df_sig) > 0) {
      res_tbl <- run_gprofiler_analysis(df_sig, group, comp_name, gp_dir, alpha_thr = cfg$alpha, lfc_thr = cfg$lfc)
      if (!is.null(res_tbl)) {
        all_gp_results[[comp_name]] <- res_tbl
        generate_interactive_sankey(res_tbl, group, comp_name, sankey_dir, top_n_terms = 30)
      }
    }
  }
  
  # Generowanie zbiorczego DotPlota porównującego wzbogacenie szlaków we wszystkich kontrastach
  generate_combined_gprofiler_dotplot(all_gp_results, gp_dir, top_n_per_contrast = 5)
}

cat("\n[SUKCES] Krok 04 ukończony pomyślnie.\n")
