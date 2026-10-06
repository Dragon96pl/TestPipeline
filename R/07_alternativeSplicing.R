# ==============================================================================
# 07_ALTERNATIVE_SPLICING.R - Analiza Splicingu i Przełączania Izoform
# ==============================================================================
cat("\n=================================================================\n")
cat(" >>> KROK 07: ALTERNATYWNY SPLICING (ISOFORMSWITCHANALYZER) <<<\n")
cat("=================================================================\n")

# --- 1. Weryfikacja środowiska i SSOT ---
if (!exists("PATHS") || !exists("BIO") || !exists("STATS") || !exists("PALETTES")) {
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

# --- 2. Wczytanie metadanych i słownika adnotacji (z SSOT PATHS) ---
rds_path <- if (!is.null(PATHS$txi_rds)) PATHS$txi_rds else file.path(PATHS$out_base_dir, "txi_salmon_data.rds")
if (!file.exists(rds_path)) {
  stop(sprintf("CRITICAL ERROR: Brak pliku %s! Uruchom najpierw krok 03_salmonLoad.R.", rds_path))
}

data_bundle  <- readRDS(rds_path)
colData_Full <- data_bundle$colData

if (!exists("anno_map") || !is.data.frame(anno_map) || nrow(anno_map) == 0) {
  cat("[ADNOTACJE] Budowanie słownika genów z GTF...\n")
  anno_map <- build_custom_annotation_map(PATHS$merged_gtf)
} else {
  cat(sprintf("[ADNOTACJE] Wykorzystano istniejący słownik anno_map (%d genów).\n", nrow(anno_map)))
}

# --- 3. Wyszukanie plików kwantyfikacji Salmona ---
quant_mode <- tolower(if (!is.null(PATHS$quant_type)) PATHS$quant_type else "salmon")
if (quant_mode == "salmon") {
  quant_files <- list.files(PATHS$salmon_dir, pattern = "^quant\\.sf$", recursive = TRUE, full.names = TRUE)
  if (length(quant_files) == 0) stop(sprintf("Nie znaleziono plików quant.sf w: %s", PATHS$salmon_dir))
  names(quant_files) <- basename(dirname(quant_files))
} else {
  stop(sprintf("Nieobsługiwany typ quant_type: %s", quant_mode))
}

common_samples <- intersect(rownames(colData_Full), names(quant_files))
quant_files    <- quant_files[common_samples]

# Progi statystyczne z SSOT
dif_cutoff <- if (exists("DTU_STATS") && !is.null(DTU_STATS$dIF_thr)) DTU_STATS$dIF_thr else 0.10
fdr_cutoff <- if (exists("DTU_STATS") && !is.null(DTU_STATS$fdr_thr)) DTU_STATS$fdr_thr else STATS$alpha_thr

# --- 4. Główna pętla analityczna ---
for (group in BIO$target_tissues) {
  cat(sprintf("\n>>> ANALIZA SPLICINGU DLA: [%s] <<<\n", toupper(group)))
  group_dir <- file.path(PATHS$out_base_dir, group)
  
  samples_group <- if (group == "Whole_Dataset") {
    rownames(colData_Full)
  } else {
    rownames(colData_Full)[colData_Full$Tissue == group]
  }
  samples_group <- intersect(samples_group, names(quant_files))
  
  if (length(samples_group) < 4) {
    cat(sprintf("  [!] Zbyt mało próbek w grupie %s do analizy splicingu.\n", group))
    next
  }
  
  colData_group <- droplevels(colData_Full[samples_group, ])
  quant_group   <- quant_files[samples_group]
  
  analysis_variants <- if (isTRUE(STATS$enable_filtering)) c("Unfiltered", "Filtered") else c("Standard")
  
  for (cond_variant in analysis_variants) {
    cat(sprintf("\n--- Przebieg: [%s] ---\n", cond_variant))
    
    out_dir_as <- file.path(group_dir, "Visualizations", "Alternative_Splicing", cond_variant)
    plots_dir  <- file.path(out_dir_as, "SwitchPlots")
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
    
    colData_run <- colData_group
    quant_run   <- quant_group
    
    # Dynamiczne dopasowanie kolejności grup na podstawie konfiguracji palety
    cond_vector    <- as.character(colData_run$Condition)
    all_levels     <- if (!is.null(names(PALETTES$condition_colors))) names(PALETTES$condition_colors) else unique(cond_vector)
    present_levels <- all_levels[all_levels %in% unique(cond_vector)]
    
    design_df <- data.frame(
      sampleID  = rownames(colData_run),
      condition = factor(cond_vector, levels = present_levels),
      stringsAsFactors = FALSE
    )
    
    # 4.1. Import kwantyfikacji
    cat("  -> [1/5] Import kwantyfikacji transkryptowej...\n")
    salmon_quant <- suppressMessages(
      importIsoformExpression(sampleVector = quant_run, calculateCountsFromAbundance = TRUE, showProgress = FALSE)
    )
    
    # 4.2. Budowa switchAnalyzeRlist
    cat("  -> [2/5] Budowanie struktury danych z GTF...\n")
    switchList <- suppressWarnings(
      importRdata(
        isoformCountMatrix            = salmon_quant$counts,
        isoformRepExpression          = salmon_quant$abundance,
        designMatrix                  = design_df,
        isoformExonAnnoation          = PATHS$merged_gtf,
        estimateDifferentialGeneRange = FALSE,
        showProgress                  = TRUE
      )
    )
    
    # 4.3. Filtracja szumu
    cat("  -> [3/5] Filtracja niskiej ekspresji (isoCount = 10, IFcutoff = 0.01)...\n")
    switchList <- IsoformSwitchAnalyzeR::preFilter(
      switchAnalyzeRlist       = switchList,
      isoCount                 = 10,
      IFcutoff                 = 0.01,
      removeSingleIsoformGenes = TRUE,
      quiet                    = FALSE
    )
    
    # 4.4. Test DEXSeq (sekwencyjny, stabilny przebieg jednowątkowy)
    cat("  -> [4/5] Test statystyczny DEXSeq...\n")
    switchList <- tryCatch(
      isoformSwitchTestDEXSeq(
        switchAnalyzeRlist     = switchList,
        alpha                  = fdr_cutoff,
        dIFcutoff              = dif_cutoff,
        ncores                 = 1,
        reduceToSwitchingGenes = FALSE,
        quiet                  = FALSE
      ),
      error = function(e) {
        cat(sprintf("     [!] Błąd DEXSeq: %s\n", e$message))
        return(NULL)
      }
    )
    
    if (is.null(switchList)) next
    
    # 4.5. Klasyfikacja zdarzeń alternatywnego splicingu
    cat("  -> [5/5] Klasyfikacja zdarzeń alternatywnego splicingu...\n")
    switchList <- analyzeAlternativeSplicing(switchAnalyzeRlist = switchList, onlySwitchingGenes = FALSE, quiet = TRUE)
    
    # --- 5. Eksport wyników i raportów ---
    switch_summary <- extractSwitchSummary(switchList, dIFcutoff = dif_cutoff, alpha = fdr_cutoff)
    write.csv(switch_summary, file.path(out_dir_as, "Summary_Isoform_Switches_Counts.csv"), row.names = FALSE)
    
    top_switches <- extractTopSwitches(switchList, filterForConsequences = FALSE, n = Inf, sortByQvals = TRUE)
    
    if (nrow(top_switches) > 0) {
      top_switches$symbol <- anno_map$symbol[match(top_switches$gene_id, anno_map$gene_id)]
      top_switches$symbol <- ifelse(is.na(top_switches$symbol) | top_switches$symbol == "", top_switches$gene_id, top_switches$symbol)
      top_switches$description <- anno_map$description[match(top_switches$gene_id, anno_map$gene_id)]
      
      write.csv(top_switches, file.path(out_dir_as, "Isoform_Switch_Significant_Transcripts.csv"), row.names = FALSE)
      
      splicing_summary <- extractSplicingSummary(switchList, asFractionTotal = FALSE, returnSummary = TRUE)
      write.csv(splicing_summary, file.path(out_dir_as, "Summary_Alternative_Splicing_Events.csv"), row.names = FALSE)
      
      png(file.path(out_dir_as, "Plot_Splicing_Events_Distribution.png"), width = 2400, height = 1800, res = 300)
      print(extractSplicingSummary(switchList, plot = TRUE))
      dev.off()
      
      # Generowanie wykresów SwitchPlot dla Top genów
      unique_top_genes <- unique(top_switches$gene_id)[1:min(15, length(unique(top_switches$gene_id)))]
      cat(sprintf("  -> Generowanie %d wykresów SwitchPlots...\n", length(unique_top_genes)))
      
      for (gid in unique_top_genes) {
        gene_sym <- top_switches$symbol[top_switches$gene_id == gid][1]
        safe_gene_name <- gsub("[^A-Za-z0-9_-]", "_", gene_sym)
        pdf_out <- file.path(plots_dir, paste0("SwitchPlot_", safe_gene_name, "_", gid, ".pdf"))
        
        tryCatch({
          pdf(pdf_out, width = 10, height = 7)
          IsoformSwitchAnalyzeR::switchPlot(switchList, gene = gid)
          dev.off()
        }, error = function(e) {
          if (file.exists(pdf_out)) file.remove(pdf_out)
        })
      }
      cat(sprintf("     [+] Zapisano SwitchPlots w: %s\n", plots_dir))
    } else {
      cat("     [!] Brak istotnych przełączeń izoform przy zadanych progach.\n")
    }
  }
}

cat("\n[SUKCES] Krok 07: Analiza alternatywnego splicingu zakończona pomyślnie.\n")