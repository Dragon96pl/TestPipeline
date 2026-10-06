# ==============================================================================
# 03_SALMONLOAD.R - Wczytywanie wyników Salmona, budowa tx2gene i colData
# ==============================================================================
cat("\n=================================================================\n")
cat(" >>> KROK 03: IMPORT DANYCH Z SALMONA I STRUKTURYZACJA METADANYCH <<<\n")
cat("=================================================================\n")

if (!exists("PATHS") || !exists("BIO") || !exists("STATS")) {
  stop("CRITICAL ERROR: Brak załadowanego configu! Odpal najpierw skrypt 01_config.R")
}

utils_path <- file.path(PATHS$scripts_dir, PATHS$utils_file)
if (file.exists(utils_path)) source(utils_path) else stop("Brak pliku 02_utils.R!")

# 1. BUDOWA TX2GENE
cat(sprintf("\n[TX2GENE] Wczytywanie pliku GTF:\n  -> %s\n", PATHS$merged_gtf))
if (!file.exists(PATHS$merged_gtf)) stop("CRITICAL ERROR: Nie znaleziono pliku GTF!")

gtf_data <- rtracklayer::import(PATHS$merged_gtf)
gtf_tx <- gtf_data[gtf_data$type %in% c("transcript", "mRNA")]

tx2gene <- data.frame(
  TXNAME = as.character(gtf_tx$transcript_id),
  GENEID = as.character(gtf_tx$gene_id),
  stringsAsFactors = FALSE
)

if ("ref_gene_id" %in% colnames(mcols(gtf_tx))) {
  valid_ref <- !is.na(gtf_tx$ref_gene_id) & gtf_tx$ref_gene_id != ""
  tx2gene$GENEID[valid_ref] <- as.character(gtf_tx$ref_gene_id[valid_ref])
}
tx2gene <- tx2gene[!is.na(tx2gene$TXNAME) & !is.na(tx2gene$GENEID) & tx2gene$TXNAME != "", ]
cat(sprintf("[TX2GENE] Zbudowano mapę dla %d transkryptów do %d loci.\n", nrow(tx2gene), length(unique(tx2gene$GENEID))))

# 2. METADANE I ŚCIEŻKI SALMONA
meta_df <- read.delim(PATHS$metadata, sep = "\t", stringsAsFactors = FALSE)
sample_ids <- as.character(meta_df[[META_MAP$sample_id]])
files <- file.path(PATHS$salmon_dir, sample_ids, "quant.sf")
names(files) <- sample_ids

missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  print(missing_files)
  stop("CRITICAL ERROR: Brak plików quant.sf dla wskazanych próbek!")
}

# 3. IMPORT PRZEZ TXIMPORT
cat("\n[TXIMPORT] Agregacja odczytów tximport...\n")
txi <- tximport(
  files, 
  type = PATHS$quant_type, 
  tx2gene = tx2gene, 
  ignoreAfterBar = TRUE, 
  ignoreTxVersion = FALSE
)

# 4. COLDATA
colData_Full <- data.frame(
  sample_id = sample_ids,
  Condition = factor(meta_df[[META_MAP$condition]]),
  row.names = sample_ids
)
if (length(BIO$target_tissues) == 1) {
  colData_Full$Tissue <- factor(BIO$target_tissues[1])
}

cat("\n[COLDATA] Zestawienie prób:\n")
print(table(colData_Full$Condition))

# 5. EKSPORT
export_mapping_stats(txi, colData_Full, PATHS$out_base_dir, file_name = "Salmon_Mapping_Statistics.csv")
export_coldata(colData_Full, PATHS$out_base_dir, file_name = "colData_Summary.csv")

rds_path <- file.path(PATHS$out_base_dir, "txi_salmon_data.rds")
saveRDS(list(txi = txi, colData = colData_Full, tx2gene = tx2gene), file = rds_path)
cat(sprintf("[EKSPORT] Zapisano stan do: %s\n", rds_path))
