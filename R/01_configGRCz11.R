# ==============================================================================
# 01_CONFIG.R - Single Source of Truth (SSOT) dla Potoku RNA-seq
# PROJEKT: Danio rerio (Zebrafish) | Referencja: GRCz11
# ==============================================================================
cat("\n[INIT] Inicjalizacja środowiska i sprawdzanie zależności...\n")

# --- 1. PAKIETY CRAN I BIOCONDUCTOR ---
cran_packages <- c(
  "ggplot2", "dplyr", "tidyr", "patchwork", "RColorBrewer", "stringr", 
  "ggalluvial", "pheatmap", "UpSetR", "plotly", "htmlwidgets", 
  "openxlsx2", "matrixStats", "ggrepel", "gprofiler2", "scales"
)

bioc_packages <- c(
  "DESeq2", "tximport", "rtracklayer", "limma", "sva", "IHW", 
  "EnhancedVolcano", "ComplexHeatmap", "AnnotationDbi", "IsoformSwitchAnalyzeR",
  "org.Dr.eg.db", "ashr", "DEXSeq", "VariantAnnotation", "maftools"
)


for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  + Instalacja pakietu CRAN: %s...\n", pkg))
    install.packages(pkg, dependencies = TRUE)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  + Instalacja pakietu Bioconductor: %s...\n", pkg))
    BiocManager::install(pkg, update = TRUE, ask = FALSE)
  }
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(RColorBrewer)
  library(stringr)
  library(ggalluvial)
  library(pheatmap)
  library(UpSetR)
  library(plotly)
  library(htmlwidgets)
  library(openxlsx2)
  library(matrixStats)
  library(ggrepel)
  library(gprofiler2)
  library(tximport)
  library(rtracklayer)
  library(limma)
  library(sva)
  library(IHW)
  library(EnhancedVolcano)
  library(ComplexHeatmap)
  library(AnnotationDbi)
  library(IsoformSwitchAnalyzeR)
  library(org.Dr.eg.db)
  library(DEXSeq)
  library(VariantAnnotation)
  library(maftools)
})
cat("[INIT] Wszystkie biblioteki załadowane pomyślnie.\n")

# ==============================================================================
# 2. PARAMETRY ŚCIEŻEK (PATHS CONFIGURATION) - GRCz11
# ==============================================================================
# ==============================================================================
# 2. PARAMETRY ŚCIEŻEK (PATHS CONFIGURATION) - SINGLE SOURCE OF TRUTH (SSOT)
# ==============================================================================
PATHS <- list(
  quant_type   = "salmon",                                                                # Typ/silnik użyty do kwantyfikacji transkrypcyjnej obsługiwany przez tximport (np. "salmon", "kallisto")
  salmon_dir   = "03_mapping_quant_GRCz11/salmon",                                        # Ścieżka do katalogu z podfolderami próbek zawierającymi pliki kwantyfikacji 'quant.sf'
  bracken_dir  = "04_advanced_omics/bracken",                                             # Katalog z raportami liczebności mikrobiomu z algorytmu Bracken (*_bracken_species.txt)
  kraken_dir   = "04_advanced_omics/kraken",                                              # Katalog z surowymi raportami klasyfikacji taksonomicznej Kraken2 (*.kreport / *.report)
  ref_gtf      = "01_data/ref/Danio_rerio.GRCz11.116.gtf.gz",                             # Ścieżka do bazowego, referencyjnego pliku adnotacji genomowej (Ensembl GTF)
  merged_gtf   = "03_mapping_quant_GRCz11/gffcompare/SUPER_REFERENCE/super_merged.gtf",   # Scalony plik adnotacji (StringTie/Gffcompare) używany do budowy słownika transkrypt -> gen (tx2gene)
  out_base_dir = "05_DESeq2_AnalysisGRCz11",                                              # Główny katalog wyjściowy dla wszystkich raportów tabelarycznych, wykresów DE, QC i integracji omicznych
  txi_rds      = "05_DESeq2_AnalysisGRCz11/txi_salmon_data.rds",                          # Centralna ścieżka do binarnego pliku RDS ze stanem tximport, colData i tx2gene (współdzielony między krokami 04-07)
  scripts_dir  = "R/",                                                                    # Katalog zawierający moduły i skrypty wykonawcze języka R w potoku
  utils_file   = "02_utils.R",                                                            # Nazwa pliku z silnikiem funkcji analitycznych, statystycznych i graficznych (Bloki 1-8)
  metadata     = "samplesheet.tsv",                                                       # Ścieżka do pliku metadanych próbek (Sample Sheet) z przypisaniem grup doświadczalnych (Condition)
  variants_dir = "04_advanced_omics/05_vep_annotated"
)

META_MAP <- list(
  sample_id = "SAMPLE_ID",  
  condition = "Condition"
)

if (!dir.exists(PATHS$out_base_dir)) dir.create(PATHS$out_base_dir, recursive = TRUE)

# ==============================================================================
# 3. PARAMETRY BIOLOGICZNE I ORGANIZM
# ==============================================================================
BIO <- list(
  org_db                  = "org.Dr.eg.db",
  org_db_name             = "Zebrafish",
  gprofiler_org           = "drerio",
  gprofiler_exclude_iea   = FALSE,
  ensembl_prefix          = "^ENS",
  target_tissues          = c("Whole_Dataset"),
  samples_per_grp         = 6,
  noise_filter_regex      = "^MSTRG|^gene-LOC",
  microbiome_contaminants = "Danio|rerio|Homo sapiens|Synthetic|artificial|vector|unclassified|Plasmid",
  condition_label         = "Experimental Condition"
)

# ==============================================================================
# 4. PARAMETRY STATYSTYCZNE I ANALITYCZNE
# ==============================================================================
STATS <- list(
  alpha_thr        = 0.05,
  lfc_thr          = 0.58,
  use_ihw          = TRUE,
  lfc_shrink       = "ashr",
  cluster_k        = 4,
  remove_batch     = FALSE,
  enable_filtering = FALSE,
  top_mutated_genes= 50, # Dynamiczna liczba genów na OncoPrint / Heatmapie
  min_variant_depth= 10  # Minimalne pokrycie odczytami
)

DTU_STATS <- list(
  dIF_thr = 0.10,
  fdr_thr = STATS$alpha_thr
)

# --- DYNAMICZNE KONTRASTY ---
CONTRASTS_CONFIG <- list(
  "Whole_Dataset" = list(
    "E3_vs_DIA"  = list(vec = c("Condition", "E3", "DIA"),  alpha = 0.05, lfc = 0.58),
    "E3_vs_DCI"  = list(vec = c("Condition", "E3", "DCI"),  alpha = 0.05, lfc = 0.58),
    "E3_vs_SCI"  = list(vec = c("Condition", "E3", "SCI"),  alpha = 0.05, lfc = 0.58),
    "E3_vs_MI"   = list(vec = c("Condition", "E3", "MI"),   alpha = 0.05, lfc = 0.58),
    
    "DIA_vs_DCI" = list(vec = c("Condition", "DIA", "DCI"), alpha = 0.05, lfc = 0.58),
    "DIA_vs_SCI" = list(vec = c("Condition", "DIA", "SCI"), alpha = 0.05, lfc = 0.58),
    "DIA_vs_MI"  = list(vec = c("Condition", "DIA", "MI"),  alpha = 0.05, lfc = 0.58),
    
    "DCI_vs_SCI" = list(vec = c("Condition", "DCI", "SCI"), alpha = 0.05, lfc = 0.58),
    "DCI_vs_MI"  = list(vec = c("Condition", "DCI", "MI"),  alpha = 0.05, lfc = 0.58),
    
    "SCI_vs_MI"  = list(vec = c("Condition", "SCI", "MI"),  alpha = 0.05, lfc = 0.58)
  )
)

# ==============================================================================
# 5. ESTETYKA I PALETY KOLORÓW (SSOT - PALETTES)
# ==============================================================================
PALETTES <- list(
  # Kolory przypisane do grup eksperymentalnych (zachowywane we wszystkich wykresach)
  condition_colors = c(
    "E3"  = "#1B9E77",
    "DIA" = "#D95F02",
    "DCI" = "#7570B3",
    "SCI" = "#E7298A",
    "MI"  = "#66A61E"
  ),
  
  # Główny gradient ekspresji Z-Score (Niebieski -> Biały -> Czerwony)
  heatmap_colors         = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(255),
  
  # Gradient macierzy odległości próbek Euclidean/Bray-Curtis
  dist_colors            = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255),
  
  # Klasyfikacja funkcjonalna wariantów Ensembl VEP (od najcięższego do modyfikatora)
  vep_impact_colors      = c(
    "HIGH"     = "#D73027",  # Czerwony - frameshift, stop-gained, splice acceptor/donor
    "MODERATE" = "#FDAE61",  # Pomarańczowy - missense, inframe insertion/deletion
    "LOW"      = "#ABD9E9",  # Błękitny - synonymous, splice region
    "MODIFIER" = "#4575B4"   # Granatowy - intronic, UTR, intergenic
  ),
  
  # Gradient nasycenia liczby wariantów w genach (Biel -> Czerwień -> Czerń/Bordowy)
  variant_heatmap_colors = colorRampPalette(c("#F7F7F7", "#FEE0D2", "#DE2D26", "#67000D")),
  
  oncoprint_colors = c(
    "missense_variant"        = "#FDAE61", # Pomarańczowy
    "frameshift_variant"      = "#D73027", # Czerwony
    "stop_gained"             = "#7F0000", # Bordowy
    "splice_region_variant"   = "#7570B3", # Fioletowy
    "inframe_indel"           = "#66C2A5", # Morski
    "synonymous_variant"      = "#E0E0E0"  # Szary
  ),
  
  substitution_colors = c(
    "A>G / T>C (ADAR)"        = "#D7191C", # Czerwony (edycja RNA)
    "C>T / G>A (APOBEC/Deam)" = "#FDAE61",
    "A>C / T>G"               = "#ABD9E9",
    "A>T / T>A"               = "#2C7BB6",
    "C>A / G>T"               = "#A6D96A",
    "C>G / G>C"               = "#FFFFBF"
  )
)

theme_set(
  theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle    = element_text(size = 11, hjust = 0.5, color = "grey30"),
      axis.title       = element_text(face = "bold"),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "grey30"),
      strip.text       = element_text(face = "bold", size = 11)
    )
)

# ==============================================================================
# 6. CENTRALNA KONFIGURACJA PODPISÓW GRAFIK (SSOT - LABELS)
# ==============================================================================
LABELS <- list(
  # --- BLOK 3: QC ---
  qc_dispersion_title = "Dispersion Estimates - %s (%s)",
  qc_pca_title        = "Principal Component Analysis (PCA) - %s",
  qc_pca_sub          = "%s",
  qc_dist_title       = "Sample Distance Matrix - %s (%s)",
  
  # --- BLOK 4: KONTRASTY ---
  ma_title            = "MA Plot - %s (%s)",
  volcano_title       = "Volcano Plot: %s",
  volcano_sub         = "%s",
  
  # --- BLOK 6: WIZUALIZACJE ZBIORCZE ---
  heatmap_top_title   = "Top 75 DEGs - %s",
  upset_ylabel        = "Shared DEGs",
  trajectories_title  = "Top DEGs Expression Trajectories - %s",
  trajectories_ylab   = "Normalized VST Expression",
  
  # --- BLOK 7: SZLAKI BIOLOGICZNE I SANKEY ---
  gprofiler_title     = "Top Enriched Biological Pathways - %s",
  gprofiler_sub       = "Organism: %s | Query Mode: %s | IEA Excluded: %s",
  gprofiler_xlab      = "Gene Ratio (Recall)",
  sankey_title        = "Gene-Pathway Flow (%s)",
  
  # --- BLOK 5: KLASTROWANIE K-MEANS ---
  clustering_title    = "Expression Profiles (K-means Z-score) - %s",
  clustering_sub      = "Variant: %s | Total genes: %d",
  clustering_ylab     = "Expression Z-Score",
  
  # --- BLOK 8: MIKROBIOM I MULTI-OMICS ---
  micro_comp_title    = "Microbiome Taxonomic Composition (Species Level)",
  micro_comp_sub      = "Top %d Most Abundant Taxa | Total Sum Scaling (TSS)",
  micro_comp_xlab     = "Samples",
  micro_comp_ylab     = "Relative Abundance (%)",
  micro_comp_fill     = "Species",
  micro_comp_other    = "Other Taxa",
  
  alpha_div_title     = "Microbiome Alpha Diversity Metrics",
  alpha_div_sub       = "Comparison Across Experimental Groups",
  alpha_div_ylab      = "Diversity Metric Value",
  alpha_metrics       = c(
    "Observed_Richness" = "Observed Richness (S)",
    "Shannon_Index"     = "Shannon Index (H')",
    "Simpson_Index"     = "Simpson Index (1 - D)",
    "Pielous_Evenness"  = "Pielou's Evenness (J')"
  ),
  
  beta_pcoa_title     = "Microbiome Beta Diversity (Bray-Curtis PCoA)",
  beta_pcoa_xlab      = "PCoA 1 [%.1f%% variance]",
  beta_pcoa_ylab      = "PCoA 2 [%.1f%% variance]",
  
  cross_cor_title     = "Cross-Omics Correlation: Host Expression vs Microbiome Abundance (*FDR<0.05, **FDR<0.01, ***FDR<0.001)",
  
  # --- BLOK 9 (skrypt 08): WARIANTY GATK + VEP ---
  variant_burden_title = "RNA-seq Variant Burden Across Conditions",
  variant_burden_sub   = "GATK HaplotypeCaller (PASS Variants)",
  variant_burden_ylab  = "Total PASS Variants Count",
  vep_impact_title     = "VEP Functional Impact Distribution (Proportions)",
  vep_impact_sub       = "Normalized per sample across conditions",
  vep_impact_ylab      = "Percentage of Total Annotated Variants",
  vep_impact_fill      = "VEP Impact",
  vep_top_genes_title  = "Top %d Genes with Protein-Altering Variants (HIGH/MODERATE Impact)",
  oncoprint_title      ="Landscape of Functional Protein Variants (Top %d Genes)",
  consequences_title   ="Distribution of Functional Coding Consequences (Excl. Modifiers)",
  radar_adar_title     ="RNA-Seq Base Substitution Spectrum (ADAR Editing Fingerprint)",
  deg_mut_overlap_title="Multi-Omics: Functional Variants in Significant DEGs"
)