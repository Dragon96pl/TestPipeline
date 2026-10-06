# ==============================================================================
# 08_VARIANT_ANALYSIS.R - Zaawansowana Analiza, Wizualizacja i Eksport Danych (SSOT)
# ==============================================================================
cat("\n=================================================================\n")
cat(" >>> KROK 08: ZAAWANSOWANA ANALIZA WARIANTÓW RNA-SEQ (GATK + VEP) <<<\n")
cat("=================================================================\n")

# --- 1. Weryfikacja środowiska SSOT ---
if (!exists("PATHS") || !exists("BIO") || !exists("PALETTES") || !exists("LABELS") || !exists("STATS")) {
  if (file.exists("01_config.R")) source("01_config.R") else stop("CRITICAL ERROR: Brak pliku 01_config.R!")
}

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(openxlsx2)
  library(pheatmap)
  library(scales)
})

# --- 2. Ścieżki i metadane ---
vcf_dir <- if (!is.null(PATHS$variants_dir) && dir.exists(PATHS$variants_dir)) PATHS$variants_dir else "04_advanced_omics/variants/05_vep_annotated"
out_var_dir <- file.path(PATHS$out_base_dir, "Variant_Analysis_GATK_VEP")
dir.create(out_var_dir, recursive = TRUE, showWarnings = FALSE)

rds_path <- if (!is.null(PATHS$txi_rds)) PATHS$txi_rds else file.path(PATHS$out_base_dir, "txi_salmon_data.rds")
if (!file.exists(rds_path)) stop("Brak pliku RDS: ", rds_path)
data_bundle  <- readRDS(rds_path)
colData_Full <- data_bundle$colData

vcf_files <- list.files(vcf_dir, pattern = "\\.vep\\.vcf\\.gz$", full.names = TRUE)
if (length(vcf_files) == 0) stop("Brak plików *.vep.vcf.gz w katalogu: ", vcf_dir)
cat(sprintf("[INIT] Odnaleziono %d plików VCF z adnotacjami VEP.\n", length(vcf_files)))

# --- 3. Ekstrakcja danych ---
sample_stats_list <- list()
detailed_consequences_list <- list()
gene_variant_list <- list()
substitution_spectrum_list <- list()

for (vf in vcf_files) {
  s_raw <- sub("\\.vep\\.vcf\\.gz$", "", basename(vf))
  s_name <- gsub("(\\.|_)(aligned|sorted|sortedbycoord|out|bam|raw|filtered)+.*$", "", s_raw, ignore.case = TRUE)
  
  matched_sample <- NULL
  if (s_name %in% rownames(colData_Full)) {
    matched_sample <- s_name
  } else {
    hit <- grep(paste0("^", s_name, "$"), rownames(colData_Full), ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) matched_sample <- hit[1]
  }
  
  cond <- if (!is.null(matched_sample)) {
    as.character(colData_Full[matched_sample, "Condition"])
  } else {
    raw_grp <- toupper(gsub("[^A-Za-z].*", "", s_name))
    if (raw_grp == "DC") "DCI" else if (raw_grp == "M") "MI" else if (raw_grp == "DIA") "DIA" else if (raw_grp == "E") "E3" else raw_grp
  }
  
  clean_id <- if (!is.null(matched_sample)) matched_sample else s_name
  cat(sprintf("  -> Przetwarzanie: [%s] (Grupa: %s)...\n", clean_id, cond))
  
  vcf <- tryCatch(suppressWarnings(readVcf(vf, genome = "GRCz11")), error = function(e) NULL)
  if (is.null(vcf) || length(vcf) == 0) next
  
  pass_mask <- filt(vcf) == "PASS" | is.na(filt(vcf))
  vcf_pass  <- vcf[pass_mask]
  if (length(vcf_pass) == 0) next
  
  # 3.1. Spektrum podstawień (ADAR Editing)
  vcf_snv_bi <- vcf_pass[isSNV(vcf_pass, singleAltOnly = TRUE)]
  if (length(vcf_snv_bi) > 0) {
    ref_s <- as.character(ref(vcf_snv_bi))
    alt_s <- as.character(unlist(alt(vcf_snv_bi)))
    subst_pairs <- paste0(ref_s, ">", alt_s)
    
    c_ag_tc <- sum(subst_pairs %in% c("A>G", "T>C"))
    c_ct_ga <- sum(subst_pairs %in% c("C>T", "G>A"))
    c_ac_tg <- sum(subst_pairs %in% c("A>C", "T>G"))
    c_at_ta <- sum(subst_pairs %in% c("A>T", "T>A"))
    c_ca_gt <- sum(subst_pairs %in% c("C>A", "G>T"))
    c_cg_gc <- sum(subst_pairs %in% c("C>G", "G>C"))
    
    total_snv <- length(vcf_snv_bi)
    substitution_spectrum_list[[clean_id]] <- data.frame(
      Sample_ID  = clean_id, 
      Condition  = cond,
      Type       = c("A>G / T>C (ADAR)", "C>T / G>A (APOBEC/Deam)", "A>C / T>G", "A>T / T>A", "C>A / G>T", "C>G / G>C"),
      Count      = c(c_ag_tc, c_ct_ga, c_ac_tg, c_at_ta, c_ca_gt, c_cg_gc),
      Percentage = round(100 * c(c_ag_tc, c_ct_ga, c_ac_tg, c_at_ta, c_ca_gt, c_cg_gc) / max(1, total_snv), 2),
      stringsAsFactors = FALSE
    )
    
    transversions <- total_snv - (c_ag_tc + c_ct_ga)
    ti_tv_ratio   <- if (transversions > 0) round((c_ag_tc + c_ct_ga) / transversions, 2) else NA
  } else {
    ti_tv_ratio <- NA
  }
  
  # 3.2. Parsowanie VEP CSQ
  high_impact <- mod_impact <- low_impact <- modifier <- 0
  
  hdr_info <- info(header(vcf_pass))
  if ("CSQ" %in% rownames(hdr_info)) {
    csq_desc   <- hdr_info["CSQ", "Description"]
    csq_fields <- unlist(strsplit(sub(".*Format: ", "", csq_desc), "\\|"))
    n_fields   <- length(csq_fields)
    
    csq_raw <- info(vcf_pass)$CSQ
    has_csq <- lengths(csq_raw) > 0
    
    if (any(has_csq)) {
      first_csq  <- sapply(csq_raw[has_csq], `[`, 1)
      split_csq  <- strsplit(first_csq, "\\|")
      csq_padded <- lapply(split_csq, function(v) {
        if (length(v) < n_fields) c(v, rep("", n_fields - length(v))) else v[1:n_fields]
      })
      
      csq_mat <- do.call(rbind, csq_padded)
      colnames(csq_mat) <- csq_fields
      csq_df  <- as.data.frame(csq_mat, stringsAsFactors = FALSE)
      
      impact_counts <- table(csq_df$IMPACT)
      high_impact   <- if ("HIGH" %in% names(impact_counts)) impact_counts[["HIGH"]] else 0
      mod_impact    <- if ("MODERATE" %in% names(impact_counts)) impact_counts[["MODERATE"]] else 0
      low_impact    <- if ("LOW" %in% names(impact_counts)) impact_counts[["LOW"]] else 0
      modifier      <- if ("MODIFIER" %in% names(impact_counts)) impact_counts[["MODIFIER"]] else 0
      
      func_csq <- csq_df[csq_df$IMPACT %in% c("HIGH", "MODERATE", "LOW"), ]
      if (nrow(func_csq) > 0) {
        main_csq <- sapply(strsplit(func_csq$Consequence, "&"), `[`, 1)
        detailed_consequences_list[[clean_id]] <- data.frame(
          Sample_ID   = clean_id,
          Condition   = cond,
          Consequence = main_csq,
          stringsAsFactors = FALSE
        )
      }
      
      prot_vars <- csq_df[csq_df$IMPACT %in% c("HIGH", "MODERATE") & csq_df$SYMBOL != "", ]
      if (nrow(prot_vars) > 0) {
        main_prot_csq <- sapply(strsplit(prot_vars$Consequence, "&"), `[`, 1)
        gene_variant_list[[clean_id]] <- data.frame(
          Sample      = clean_id,
          Condition   = cond,
          Gene        = prot_vars$SYMBOL,
          Impact      = prot_vars$IMPACT,
          Consequence = main_prot_csq,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  sample_stats_list[[clean_id]] <- data.frame(
    Sample_ID       = clean_id,
    Condition       = cond,
    Total_PASS_Vars = length(vcf_pass),
    SNPs_Count      = sum(isSNV(vcf_pass)),
    InDels_Count    = sum(isIndel(vcf_pass)),
    Ti_Tv_Ratio     = ti_tv_ratio,
    Functional_Vars = (high_impact + mod_impact),
    Impact_HIGH     = high_impact,
    Impact_MODERATE = mod_impact,
    Impact_LOW      = low_impact,
    Impact_MODIFIER = modifier,
    stringsAsFactors = FALSE
  )
}

df_var_summary <- do.call(rbind, sample_stats_list)

# Porządkowanie poziomów grup z PALETTES (SSOT)
present_conds  <- as.character(df_var_summary$Condition)
ordered_levels <- if (!is.null(names(PALETTES$condition_colors))) intersect(names(PALETTES$condition_colors), unique(present_conds)) else unique(present_conds)
df_var_summary$Condition <- factor(present_conds, levels = ordered_levels)
cond_palette   <- PALETTES$condition_colors[names(PALETTES$condition_colors) %in% ordered_levels]

# --- 4. Budowa pełnych macierzy genów (ALL vs TOP N) ---
df_gvars_all <- if (length(gene_variant_list) > 0) do.call(rbind, gene_variant_list) else data.frame()

if (nrow(df_gvars_all) > 0) {
  sample_ids <- as.character(df_var_summary$Sample_ID)
  all_unique_genes <- sort(unique(df_gvars_all$Gene))
  
  # 4.1. Pełna macierz - Wszystkie Geny
  mat_all_genes <- matrix(0L, nrow = length(all_unique_genes), ncol = length(sample_ids),
                          dimnames = list(all_unique_genes, sample_ids))
  for (i in seq_len(nrow(df_gvars_all))) {
    g <- df_gvars_all$Gene[i]
    s <- df_gvars_all$Sample[i]
    if (g %in% all_unique_genes && s %in% colnames(mat_all_genes)) {
      mat_all_genes[g, s] <- mat_all_genes[g, s] + 1L
    }
  }
  
  df_mat_all <- data.frame(
    Gene = rownames(mat_all_genes),
    Total_Mutations_All_Samples = rowSums(mat_all_genes),
    mat_all_genes,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  df_mat_all <- df_mat_all[order(-df_mat_all$Total_Mutations_All_Samples), ]
  
  # 4.2. Macierz Top N Genów
  top_n_cfg   <- if (exists("STATS") && !is.null(STATS$top_mutated_genes)) STATS$top_mutated_genes else 50
  top_n_genes <- min(top_n_cfg, nrow(df_mat_all))
  df_mat_top  <- head(df_mat_all, top_n_genes)
  
  write.csv(df_mat_all, file.path(out_var_dir, "All_Genes_Functional_Variants_Matrix.csv"), row.names = FALSE)
  write.csv(df_mat_top, file.path(out_var_dir, sprintf("Top%d_Genes_Functional_Variants_Matrix.csv", top_n_genes)), row.names = FALSE)
  cat(sprintf("  [+] Zapisano pełną macierz %d genów oraz Top %d do plików CSV.\n", nrow(df_mat_all), top_n_genes))
}

# --- 5. Zestawienia podsumowujące dla wykresów słupkowych ---
# 5.1. Tabela VEP Impact
df_impact_long <- tidyr::pivot_longer(
  df_var_summary,
  cols = c("Impact_HIGH", "Impact_MODERATE", "Impact_LOW", "Impact_MODIFIER"),
  names_to = "Impact_Type",
  values_to = "Count"
)
df_impact_long$Impact_Type <- factor(gsub("Impact_", "", df_impact_long$Impact_Type), levels = names(PALETTES$vep_impact_colors))
df_impact_long <- df_impact_long %>%
  dplyr::group_by(Sample_ID) %>%
  dplyr::mutate(
    Total_Vars = sum(Count),
    Percentage = round(100 * Count / max(1, Total_Vars), 2)
  ) %>%
  dplyr::ungroup()
write.csv(df_impact_long, file.path(out_var_dir, "VEP_Impact_Distribution_Summary.csv"), row.names = FALSE)

# 5.2. Tabela Konsekwencji Funkcjonalnych
df_csq_summary <- data.frame()
if (length(detailed_consequences_list) > 0) {
  df_csq_all <- do.call(rbind, detailed_consequences_list)
  df_csq_all$Condition <- factor(df_csq_all$Condition, levels = ordered_levels)
  
  df_csq_summary <- df_csq_all %>%
    dplyr::group_by(Sample_ID, Condition, Consequence) %>%
    dplyr::summarise(Count = n(), .groups = "drop") %>%
    dplyr::group_by(Sample_ID) %>%
    dplyr::mutate(
      Total_Functional = sum(Count),
      Percentage = round(100 * Count / max(1, Total_Functional), 2)
    ) %>%
    dplyr::ungroup()
  write.csv(df_csq_summary, file.path(out_var_dir, "Functional_Consequences_Detailed_Summary.csv"), row.names = FALSE)
}

# 5.3. Tabela Podstawień ADAR
df_subst_summary <- if (length(substitution_spectrum_list) > 0) do.call(rbind, substitution_spectrum_list) else data.frame()
if (nrow(df_subst_summary) > 0) {
  df_subst_summary$Condition <- factor(df_subst_summary$Condition, levels = ordered_levels)
  write.csv(df_subst_summary, file.path(out_var_dir, "ADAR_Substitution_Spectrum_Summary.csv"), row.names = FALSE)
}

# --- 6. Eksport wieloarkuszowego Master Excela (openxlsx2) ---
wb_master <- openxlsx2::wb_workbook()
add_formatted_sheet <- function(wb, sheet_name, df_data, header_color = "366092") {
  if (is.null(df_data) || nrow(df_data) == 0) return(wb)
  wb$add_worksheet(sheet_name)$add_data(sheet = sheet_name, x = df_data)
  h_dim <- openxlsx2::wb_dims(rows = 1, cols = 1:ncol(df_data))
  wb$add_fill(sheet = sheet_name, dims = h_dim, color = openxlsx2::wb_color(hex = header_color))
  wb$add_font(sheet = sheet_name, dims = h_dim, color = openxlsx2::wb_color(hex = "FFFFFF"), bold = TRUE)
  wb$set_col_widths(sheet = sheet_name, cols = 1:ncol(df_data), widths = "auto")
  return(wb)
}

wb_master <- add_formatted_sheet(wb_master, "1. Sample_QC_Stats", df_var_summary, "1F497D")
if (exists("df_mat_all")) wb_master <- add_formatted_sheet(wb_master, "2. All_Genes_Matrix", df_mat_all, "366092")
if (exists("df_mat_top")) wb_master <- add_formatted_sheet(wb_master, "3. Top_Genes_Matrix", df_mat_top, "366092")
if (nrow(df_csq_summary) > 0) wb_master <- add_formatted_sheet(wb_master, "4. Functional_Consequences", df_csq_summary, "595959")
if (nrow(df_subst_summary) > 0) wb_master <- add_formatted_sheet(wb_master, "5. ADAR_Editing_Spectrum", df_subst_summary, "C00000")
wb_master <- add_formatted_sheet(wb_master, "6. VEP_Impact_Summary", df_impact_long, "ED7D31")

master_xlsx_path <- file.path(out_var_dir, "Comprehensive_Master_Variants_Report.xlsx")
wb_master$save(master_xlsx_path, overwrite = TRUE)
cat(sprintf("  [+] Zapisano wieloarkuszowy Master Excel: %s\n", basename(master_xlsx_path)))

# ==============================================================================
# 7. GENEROWANIE WYKRESÓW Z ETYKIETAMI LICZBOWYMI (PERCENTAGES)
# ==============================================================================

# --- 7.1. Profil Impaktu VEP z etykietami % ---
p_impact <- ggplot(df_impact_long, aes(x = Sample_ID, y = Count, fill = Impact_Type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.8, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = ifelse(Percentage >= 4.0, sprintf("%.1f%%", Percentage), "")),
    position = position_fill(vjust = 0.5),
    size = 3.2,
    fontface = "bold",
    color = "white"
  ) +
  facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = PALETTES$vep_impact_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title    = LABELS$vep_impact_title,
    subtitle = "Percentages displayed inside slices (>= 4%)",
    x        = "Samples",
    y        = LABELS$vep_impact_ylab,
    fill     = LABELS$vep_impact_fill
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text       = element_text(face = "bold")
  )
ggsave(file.path(out_var_dir, "Plot_VEP_Impact_Distribution.png"), plot = p_impact, width = 12, height = 6.5, dpi = 300)

# --- 7.2. Zbliżenie na konsekwencje funkcjonalne z etykietami % ---
if (nrow(df_csq_summary) > 0) {
  top_csq <- names(sort(table(df_csq_all$Consequence), decreasing = TRUE))[1:min(8, length(unique(df_csq_all$Consequence)))]
  df_csq_sub <- df_csq_summary[df_csq_summary$Consequence %in% top_csq, ]
  
  p_csq <- ggplot(df_csq_sub, aes(x = Sample_ID, y = Count, fill = Consequence)) +
    geom_bar(stat = "identity", position = "fill", width = 0.8, color = "black", linewidth = 0.2) +
    geom_text(
      aes(label = ifelse(Percentage >= 3.0, sprintf("%.1f%%", Percentage), "")),
      position = position_fill(vjust = 0.5),
      size = 3.0,
      fontface = "bold",
      color = "black"
    ) +
    facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title    = LABELS$consequences_title,
      subtitle = "Percentages displayed inside slices (>= 3%)",
      x        = "Samples",
      y        = "Proportion of Functional Variants",
      fill     = "Consequence"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text       = element_text(face = "bold")
    )
  ggsave(file.path(out_var_dir, "Plot_Functional_Consequences_Zoom.png"), plot = p_csq, width = 12, height = 6.5, dpi = 300)
}

# --- 7.3. Spektrum Edycji ADAR z etykietami % ---
if (nrow(df_subst_summary) > 0) {
  # Zabezpieczenie: automatyczne utworzenie kolumny Percentage, jeśli jej brakuje
  if (!"Percentage" %in% colnames(df_subst_summary)) {
    if ("Fraction" %in% colnames(df_subst_summary)) {
      df_subst_summary$Percentage <- df_subst_summary$Fraction * 100
    } else {
      df_subst_summary <- df_subst_summary %>%
        dplyr::group_by(Sample_ID) %>%
        dplyr::mutate(Percentage = round(100 * Count / max(1, sum(Count)), 2)) %>%
        dplyr::ungroup()
    }
  }
  
  subst_cols <- if (exists("PALETTES") && !is.null(PALETTES$substitution_colors)) PALETTES$substitution_colors else scales::hue_pal()(6)
  
  p_subst <- ggplot(df_subst_summary, aes(x = Sample_ID, y = Percentage, fill = Type)) +
    geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.2) +
    geom_text(
      aes(label = ifelse(Percentage >= 3.5, sprintf("%.1f%%", Percentage), "")),
      position = position_stack(vjust = 0.5),
      size = 3.0,
      fontface = "bold",
      color = "black"
    ) +
    facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = subst_cols) +
    labs(
      title    = LABELS$radar_adar_title,
      subtitle = "High A>G/T>C fraction indicates post-transcriptional ADAR enzymatic RNA-editing",
      x        = "Samples",
      y        = "Percentage of Total SNVs (%)",
      fill     = "Substitution Class"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text       = element_text(face = "bold")
    )
  
  ggsave(file.path(out_var_dir, "Plot_ADAR_RNA_Editing_Spectrum.png"), plot = p_subst, width = 12, height = 6.5, dpi = 300)
  cat("  [+] Zapisano wykres: Plot_ADAR_RNA_Editing_Spectrum.png\n")
}

# --- 7.4. Heatmapa Top N Genów (Czysty, klasyczny styl bez cyfr w kafelkach) ---
if (exists("df_mat_top") && nrow(df_mat_top) > 0) {
  mat_plot <- as.matrix(df_mat_top[, sample_ids, drop = FALSE])
  rownames(mat_plot) <- df_mat_top$Gene
  
  ann_col <- data.frame(Condition = df_var_summary$Condition, row.names = sample_ids)
  ann_colors <- list(Condition = cond_palette)
  
  hm_palette <- if (!is.null(PALETTES$variant_heatmap_colors)) {
    if (is.function(PALETTES$variant_heatmap_colors)) PALETTES$variant_heatmap_colors(50) else PALETTES$variant_heatmap_colors
  } else {
    colorRampPalette(c("#F7F7F7", "#FEE0D2", "#DE2D26", "#67000D"))(50)
  }
  
  max_v <- max(mat_plot, na.rm = TRUE)
  hm_breaks <- seq(0, max(1, max_v), length.out = length(hm_palette) + 1)
  
  plot_height <- max(7, min(25, 0.22 * nrow(mat_plot) + 2.5))
  hm_file <- file.path(out_var_dir, sprintf("Plot_Top%d_Mutated_Genes_Heatmap.png", nrow(mat_plot)))
  
  png(hm_file, width = 11, height = plot_height, units = "in", res = 300)
  pheatmap::pheatmap(
    mat_plot,
    annotation_col    = ann_col,
    annotation_colors = ann_colors,
    color             = hm_palette,
    breaks            = hm_breaks,
    display_numbers   = FALSE,
    main              = sprintf(LABELS$oncoprint_title, nrow(mat_plot)),
    fontsize          = 9,
    fontsize_row      = if (nrow(mat_plot) <= 30) 8 else 6.5,
    fontsize_col      = 8,
    cluster_cols      = (ncol(mat_plot) >= 2 && stats::sd(mat_plot) > 0),
    cluster_rows      = (nrow(mat_plot) >= 2 && stats::sd(mat_plot) > 0),
    show_colnames     = TRUE,
    show_rownames     = TRUE
  )
  dev.off()
  cat(sprintf("  [+] Zapisano czystą Heatmapę Top %d: %s\n", nrow(mat_plot), basename(hm_file)))
}

# --- 7.5. Integracja Multi-Omics: Mutacje w DEG (Czysty wykres bez nakładających się liczb) ---
sig_files <- list.files(file.path(PATHS$out_base_dir, "Whole_Dataset", "significant", "Standard"),
                        pattern = "_shrink_sig\\.csv$", full.names = TRUE)

if (length(sig_files) > 0 && nrow(df_gvars_all) > 0) {
  cat("\n  -> [Multi-Omics] Integracja wariantów z genami różnicującymi (DEGs)...\n")
  deg_list <- lapply(sig_files, function(f) read.csv(f, stringsAsFactors = FALSE)$symbol)
  all_degs <- unique(unlist(deg_list))
  all_degs <- all_degs[!is.na(all_degs) & all_degs != "" & !grepl(BIO$noise_filter_regex, all_degs)]
  
  deg_mutations <- df_gvars_all[df_gvars_all$Gene %in% all_degs, ]
  
  if (nrow(deg_mutations) > 0) {
    deg_mut_summary <- deg_mutations %>%
      dplyr::group_by(Gene, Condition, Consequence) %>%
      dplyr::summarise(Variant_Occurrences = n(), .groups = "drop")
    
    write.csv(deg_mut_summary, file.path(out_var_dir, "MultiOmics_All_DEGs_With_Functional_Variants.csv"), row.names = FALSE)
    
    wb_master <- add_formatted_sheet(wb_master, "7. DEGs_With_Variants", deg_mut_summary, "375623")
    wb_master$save(master_xlsx_path, overwrite = TRUE)
    
    top_n_cfg <- if (exists("STATS") && !is.null(STATS$top_mutated_genes)) STATS$top_mutated_genes else 50
    gene_rank <- deg_mut_summary %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(Total_Vars = sum(Variant_Occurrences), .groups = "drop") %>%
      dplyr::arrange(desc(Total_Vars))
    
    top_n_genes     <- min(top_n_cfg, nrow(gene_rank))
    top_genes_vec   <- head(gene_rank$Gene, top_n_genes)
    deg_mut_plot_df <- deg_mut_summary[deg_mut_summary$Gene %in% top_genes_vec, ]
    
    plot_height  <- max(6, 0.22 * top_n_genes + 2.0)
    out_png_file <- file.path(out_var_dir, sprintf("Plot_MultiOmics_Top%d_DEGs_With_Variants.png", top_n_genes))
    
    p_deg_mut <- ggplot(deg_mut_plot_df, aes(x = reorder(Gene, Variant_Occurrences, sum), y = Variant_Occurrences, fill = Condition)) +
      geom_bar(stat = "identity", position = "dodge", width = 0.75, color = "black", linewidth = 0.2) +
      coord_flip() +
      scale_fill_manual(values = cond_palette) +
      labs(
        title    = LABELS$deg_mut_overlap_title,
        subtitle = sprintf("Top %d Differentially Expressed Genes carrying HIGH/MODERATE functional variants", top_n_genes),
        x        = "Differentially Expressed Gene (DEG)",
        y        = "Number of Detected Variants"
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.y     = element_text(face = "bold.italic", size = if (top_n_genes <= 30) 8.5 else 7),
        legend.position = "right"
      )
    
    ggsave(out_png_file, plot = p_deg_mut, width = 10, height = plot_height, dpi = 300, limitsize = FALSE)
    cat(sprintf("     [+] Zapisano czysty wykres Top %d DEG: %s\n", top_n_genes, basename(out_png_file)))
  }
}

cat("\n=================================================================\n")
cat(" [SUKCES] Krok 08: Analiza wariantów zakończona pomyślnie.\n")
cat("=================================================================\n")