# ==============================================================================
# 06_MULTIOMICS_INTEGRATION.R - Integracja Mikrobiomu (Bracken) i Transkryptomu
# ==============================================================================
cat("\n=================================================================\n")
cat(" >>> KROK 06: POTOK MIKROBIOMU ORAZ INTEGRACJA MULTI-OMICS <<<\n")
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
  stop("CRITICAL ERROR: Nie odnaleziono pliku 02_utils.R!")
}

# --- 2. Wczytanie danych z kwantyfikacji Salmona i słownika genów ---
rds_path <- file.path(PATHS$out_base_dir, "txi_salmon_data.rds")
if (!file.exists(rds_path)) {
  stop(sprintf("CRITICAL ERROR: Brak pliku %s! Uruchom najpierw krok 03_salmonLoad.R.", rds_path))
}

data_bundle  <- readRDS(rds_path)
txi_full     <- data_bundle$txi
colData_Full <- data_bundle$colData

if (!exists("anno_map") || !is.data.frame(anno_map) || nrow(anno_map) == 0) {
  cat("[ADNOTACJE] Budowanie słownika genów z GTF...\n")
  anno_map <- build_custom_annotation_map(PATHS$merged_gtf)
} else {
  cat(sprintf("[ADNOTACJE] Wykorzystano istniejący słownik anno_map (%d genów).\n", nrow(anno_map)))
}

# --- 3. Wyszukanie i parsowanie raportów Bracken (SSOT) ---
bracken_dir <- if (!is.null(PATHS$bracken_dir) && dir.exists(PATHS$bracken_dir)) {
  PATHS$bracken_dir
} else if (dir.exists("04_advanced_omics/bracken")) {
  "04_advanced_omics/bracken"
} else if (dir.exists(file.path("..", "04_advanced_omics", "bracken"))) {
  file.path("..", "04_advanced_omics", "bracken")
} else {
  stop("CRITICAL ERROR: Nie znaleziono katalogu z raportami Brackena (PATHS$bracken_dir)!")
}

cat(sprintf("[MIKROBIOM] Odnaleziono katalog Brackena: %s\n", bracken_dir))
b_files <- list.files(bracken_dir, pattern = "_bracken_species\\.txt$", full.names = TRUE)

if (length(b_files) == 0) {
  stop(sprintf("CRITICAL ERROR: Brak plików *_bracken_species.txt w katalogu %s!", bracken_dir))
}

filter_pattern <- if (!is.null(BIO$microbiome_contaminants)) {
  BIO$microbiome_contaminants
} else {
  "Danio|rerio|Homo sapiens|Synthetic|artificial|vector|unclassified|Plasmid"
}

b_list <- list()
all_species <- character()

for (bf in b_files) {
  s_name <- sub("_bracken_species\\.txt$", "", basename(bf))
  df_b <- tryCatch(read.table(bf, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = ""), error = function(e) NULL)
  
  if (!is.null(df_b) && "name" %in% colnames(df_b) && "new_est_reads" %in% colnames(df_b)) {
    # Usunięcie zanieczyszczeń gospodarza, człowieka i wektorów
    df_clean <- df_b[!grepl(filter_pattern, df_b$name, ignore.case = TRUE), ]
    df_clean <- df_clean[df_clean$new_est_reads > 0, ]
    
    if (nrow(df_clean) > 0) {
      b_list[[s_name]] <- setNames(df_clean$new_est_reads, df_clean$name)
      all_species <- unique(c(all_species, df_clean$name))
    }
  }
}

cat(sprintf("  -> Sparsowano %d próbek. Zidentyfikowano %d unikalnych gatunków bakteryjnych po filtracji.\n", 
            length(b_list), length(all_species)))

# Budowa macierzy liczebności mikrobiomu (Próbki x Gatunki)
sample_names <- rownames(colData_Full)
microbiome_matrix <- matrix(0, nrow = length(sample_names), ncol = length(all_species),
                            dimnames = list(sample_names, all_species))

for (s_name in names(b_list)) {
  matched_row <- intersect(s_name, sample_names)
  if (length(matched_row) > 0) {
    spec_data <- b_list[[s_name]]
    microbiome_matrix[matched_row, names(spec_data)] <- spec_data
  }
}

microbiome_matrix <- microbiome_matrix[, colSums(microbiome_matrix) > 0, drop = FALSE]

# --- 4. Główna pętla analityczna po podzbiorach (Tissue/Dataset) ---
for (group in BIO$target_tissues) {
  cat(sprintf("\n=================================================================\n"))
  cat(sprintf(" >>> INTEGRACJA MULTI-OMICS DLA GRUPY: [%s] <<<\n", toupper(group)))
  cat(sprintf("=================================================================\n"))
  
  group_dir <- if (group == "Whole_Dataset") {
    file.path(PATHS$out_base_dir, "Microbiome_and_Multiomics")
  } else {
    file.path(PATHS$out_base_dir, group, "Microbiome_and_Multiomics")
  }
  dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
  
  samples_group <- if (group == "Whole_Dataset") {
    rownames(colData_Full)
  } else {
    rownames(colData_Full)[colData_Full$Tissue == group]
  }
  
  if (length(samples_group) < 3) {
    cat(sprintf("  [!] Zbyt mało próbek w grupie %s do przeprowadzenia analizy.\n", group))
    next
  }
  
  colData_group <- droplevels(colData_Full[samples_group, , drop = FALSE])
  micro_group   <- microbiome_matrix[samples_group, , drop = FALSE]
  
  # --- 4.1. Eksport surowej i znormalizowanej macierzy mikrobiomu (openxlsx2) ---
  cat("  -> Eksport macierzy mikrobiomu do raportu Excel...\n")
  rel_group <- (micro_group / pmax(1, rowSums(micro_group))) * 100
  
  df_raw <- data.frame(Sample = rownames(micro_group), micro_group, check.names = FALSE)
  df_rel <- data.frame(Sample = rownames(rel_group), rel_group, check.names = FALSE)
  
  summary_xlsx <- file.path(group_dir, paste0("Microbiome_Comprehensive_Report_", group, ".xlsx"))
  wb <- openxlsx2::wb_workbook()
  
  wb$add_worksheet("Raw_Counts")$add_data("Raw_Counts", x = df_raw)
  h_raw <- openxlsx2::wb_dims(rows = 1, cols = 1:ncol(df_raw))
  wb$add_fill("Raw_Counts", dims = h_raw, color = openxlsx2::wb_color(hex = "4F81BD"))
  wb$add_font("Raw_Counts", dims = h_raw, color = openxlsx2::wb_color(hex = "FFFFFF"), bold = TRUE)
  wb$set_col_widths("Raw_Counts", cols = 1:ncol(df_raw), widths = "auto")
  
  wb$add_worksheet("Relative_Abundance_Percent")$add_data("Relative_Abundance_Percent", x = df_rel)
  h_rel <- openxlsx2::wb_dims(rows = 1, cols = 1:ncol(df_rel))
  wb$add_fill("Relative_Abundance_Percent", dims = h_rel, color = openxlsx2::wb_color(hex = "4F81BD"))
  wb$add_font("Relative_Abundance_Percent", dims = h_rel, color = openxlsx2::wb_color(hex = "FFFFFF"), bold = TRUE)
  wb$set_col_widths("Relative_Abundance_Percent", cols = 1:ncol(df_rel), widths = "auto")
  
  wb$save(summary_xlsx, overwrite = TRUE)
  
  # --- 4.2. Wykresy ekologiczne (Kompozycja, Alfa- i Beta-Różnorodność) ---
  plot_microbiome_composition(micro_group, colData_group, out_dir = group_dir, top_n = 15)
  calculate_and_plot_alpha_diversity(micro_group, colData_group, out_dir = group_dir)
  run_beta_diversity_pcoa(micro_group, colData_group, out_dir = group_dir)
  
  # --- 4.3. Normalizacja VST transkryptomu gospodarza ---
  cat("  -> Obliczanie normalizacji transkryptomu gospodarza (DESeq2 VST)...\n")
  txi_group <- list(
    abundance           = txi_full$abundance[, samples_group, drop = FALSE],
    counts              = txi_full$counts[, samples_group, drop = FALSE],
    length              = txi_full$length[, samples_group, drop = FALSE],
    countsFromAbundance = txi_full$countsFromAbundance
  )
  
  keep_genes <- rowSums(txi_group$counts) > 0
  txi_group$abundance <- txi_group$abundance[keep_genes, , drop = FALSE]
  txi_group$counts    <- txi_group$counts[keep_genes, , drop = FALSE]
  txi_group$length    <- txi_group$length[keep_genes, , drop = FALSE]
  
  dds_host <- DESeqDataSetFromTximport(txi_group, colData = colData_group, design = ~ Condition)
  dds_host <- estimateSizeFactors(dds_host)
  vsd_host <- vst(dds_host, blind = FALSE)
  
  # --- 4.4. PUBLIKACYJNY, DETERMINISTYCZNY DOBÓR GENÓW DO MULTI-OMICS ---
  sig_dir <- file.path(PATHS$out_base_dir, group, "significant", "Standard")
  if (!dir.exists(sig_dir)) sig_dir <- file.path(PATHS$out_base_dir, group, "significant", "Unfiltered")
  
  deg_ranks_list <- list()
  if (dir.exists(sig_dir)) {
    sig_files <- list.files(sig_dir, pattern = "_shrink_sig\\.csv$", full.names = TRUE)
    for (f in sig_files) {
      df_s <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(df_s) && all(c("gene_id", "padj", "symbol") %in% colnames(df_s)) && nrow(df_s) > 0) {
        deg_ranks_list[[basename(f)]] <- df_s[, c("gene_id", "symbol", "padj", "log2FoldChange")]
      }
    }
  }
  
  if (length(deg_ranks_list) > 0) {
    df_all_degs <- do.call(rbind, deg_ranks_list)
    
    # Czyszczenie szumu i odrzucenie nienazwanych loci jeśli to możliwe
    df_clean_degs <- df_all_degs[!grepl(BIO$noise_filter_regex, df_all_degs$symbol) & 
                                   !grepl("^ENSDARG", df_all_degs$symbol) & 
                                   !is.na(df_all_degs$padj), ]
    
    # Jeśli po odrzuceniu ENSDARG zostało za mało genów, bierzemy wszystkie
    if (length(unique(df_clean_degs$gene_id)) < 20) {
      df_clean_degs <- df_all_degs[!grepl(BIO$noise_filter_regex, df_all_degs$symbol) & !is.na(df_all_degs$padj), ]
    }
    
    # Wybór Top 40 unikalnych genów o najniższym globalnym padj
    top_deg_ranking <- df_clean_degs %>%
      dplyr::group_by(gene_id, symbol) %>%
      dplyr::summarise(min_padj = min(padj, na.rm = TRUE), 
                       max_lfc  = max(abs(log2FoldChange), na.rm = TRUE), 
                       .groups  = "drop") %>%
      dplyr::arrange(min_padj, desc(max_lfc))
    
    top_target_ids <- head(top_deg_ranking$gene_id, 40)
    cat(sprintf("  [MULTI-OMICS] Wyselekcjonowano %d unikalnych DEG (metoda: min padj).\n", length(top_target_ids)))
  } else {
    # Fallback dla braku DEG: Top 40 genów o najwyższej wariancji VST
    rv <- matrixStats::rowVars(assay(vsd_host))
    clean_idx <- !grepl(BIO$noise_filter_regex, rownames(vsd_host))
    top_target_ids <- head(rownames(vsd_host)[clean_idx][order(rv[clean_idx], decreasing = TRUE)], 40)
  }
  
  # --- 4.5. Wywołanie integracji Multi-Omics ---
  run_host_microbiome_cross_correlation(
    vsd_clean      = vsd_host,
    microbiome_mat = micro_group,
    top_deg_ids    = top_target_ids,
    anno_map       = anno_map,
    colData        = colData_group,
    out_dir        = group_dir,
    top_n_taxa     = 20
  )
}

cat("\n=================================================================\n")
cat(" [SUKCES] Krok 06: Integracja Multi-Omics została pomyślnie ukończona.\n")
cat("=================================================================\n")
