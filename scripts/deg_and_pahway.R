```{r}
library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(MAST)

sce <- readH5AD("output/GSE183276_celltypist.h5ad")
seurat_obj <- as.Seurat(sce, counts = "counts", data = "logcounts")
saveRDS(seurat_obj, "output/GSE183276_with_celltypist.rds")
```


```{r}
colnames(seurat_obj@meta.data)
table(seurat_obj$Condition, seurat_obj$majority_voting)

# Remove _Diabetic, _Healthy from majority_voting
seurat_obj$majority_voting <- gsub(
  "_(Hypertensive CKD \\(H-CKD\\)|Normal Reference|Diabetic Kidney Disease \\(DKD\\)|Acute Kidney Injury \\(AKI\\))$",
  "",
  seurat_obj$majority_voting
)
table(seurat_obj$Condition, seurat_obj$majority_voting)
```


```{r small sanity check}
DefaultAssay(seurat_obj) <- "originalexp"
celltypes <- sort(unique(seurat_obj$majority_voting))
celltypes

celltypes <- gsub("-", "_", celltypes)
celltypes
```


```{r Function to compute DEGs}
run_mast_per_celltype_AKI_vs_Normal <- function(
    seu,
    celltype_col = "majority_voting",
    condition_col = "Condition",
    outdir = "output/DEG_MAST_AKI_vs_Normal",
    min_cells_per_group = 20
) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  meta <- seu@meta.data
  celltypes <- sort(unique(meta[[celltype_col]]))
  
  for (ct in celltypes) {
    
    cat("\n▶ Processing cell type:", ct, "\n")
    
    cells_ct <- rownames(meta)[meta[[celltype_col]] == ct]
    cond_tab <- table(meta[cells_ct, condition_col])
    
    if (!all(c("Acute Kidney Injury (AKI)", "Normal Reference") %in% names(cond_tab))) {
      cat("  AKI or Normal missing → skipping\n")
      next
    }
    
    if (any(cond_tab[c("Acute Kidney Injury (AKI)", "Normal Reference")] < min_cells_per_group)) {
      cat("  Not enough cells → skipping\n")
      next
    }
    
    seu_ct <- subset(seu, cells = cells_ct)
    DefaultAssay(seu_ct) <- "originalexp"
    seu_ct <- NormalizeData(seu_ct, verbose = FALSE)
    Idents(seu_ct) <- condition_col
    
    deg <- FindMarkers(
      seu_ct,
      ident.1 = "Acute Kidney Injury (AKI)",
      ident.2 = "Normal Reference",
      test.use = "MAST",
      logfc.threshold = 0,
      min.pct = 0.1,
      latent.vars = c("nCount_RNA", "percent.mt")
    )
    
    deg$gene <- rownames(deg)
    
    ct_clean <- gsub("[^a-zA-Z0-9]", "_", ct)
    out_file <- file.path(outdir, paste0("DEG_", ct_clean, "_AKI_vs_Normal.csv"))
    write.csv(deg, out_file, row.names = FALSE)
    cat("  ✔ Saved:", out_file, "\n")
  }
}
```


```{r Performing DEG analysis}
run_mast_per_celltype_AKI_vs_Normal(
  seu = seurat_obj,
  celltype_col = "majority_voting",
  condition_col = "Condition"
)
```
run_mast_per_celltype_DKD_vs_Normal <- function(
    seu,
    celltype_col = "majority_voting",
    condition_col = "Condition",
    outdir = "output/DEG_MAST_DKD_vs_Normal",
    min_cells_per_group = 20
) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  meta <- seu@meta.data
  celltypes <- sort(unique(meta[[celltype_col]]))
  
  for (ct in celltypes) {
    
    cat("\n▶ Processing cell type:", ct, "\n")
    
    cells_ct <- rownames(meta)[meta[[celltype_col]] == ct]
    cond_tab <- table(meta[cells_ct, condition_col])
    
    if (!all(c("Diabetic Kidney Disease (DKD)", "Normal Reference") %in% names(cond_tab))) {
      cat("  DKD or Normal missing → skipping\n")
      next
    }
    
    if (any(cond_tab[c("Diabetic Kidney Disease (DKD)", "Normal Reference")] < min_cells_per_group)) {
      cat("  Not enough cells → skipping\n")
      next
    }
    
    seu_ct <- subset(seu, cells = cells_ct)
    DefaultAssay(seu_ct) <- "originalexp"
    seu_ct <- NormalizeData(seu_ct, verbose = FALSE)
    Idents(seu_ct) <- condition_col
    
    deg <- FindMarkers(
      seu_ct,
      ident.1 = "Diabetic Kidney Disease (DKD)",
      ident.2 = "Normal Reference",
      test.use = "MAST",
      logfc.threshold = 0,
      min.pct = 0.1,
      latent.vars = c("nCount_RNA", "percent.mt")
    )
    
    deg$gene <- rownames(deg)
    
    ct_clean <- gsub("[^a-zA-Z0-9]", "_", ct)
    out_file <- file.path(outdir, paste0("DEG_", ct_clean, "_DKD_vs_Normal.csv"))
    write.csv(deg, out_file, row.names = FALSE)
    cat("  ✔ Saved:", out_file, "\n")
  }
}
run_mast_per_celltype_DKD_vs_Normal(
  seu = seurat_obj,
  celltype_col = "majority_voting",
  condition_col = "Condition"
)
run_mast_per_celltype_HCKD_vs_Normal <- function(
    seu,
    celltype_col = "majority_voting",
    condition_col = "Condition",
    outdir = "output/DEG_MAST_HCKD_vs_Normal",
    min_cells_per_group = 20
) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  meta <- seu@meta.data
  celltypes <- sort(unique(meta[[celltype_col]]))
  
  for (ct in celltypes) {
    
    cat("\n▶ Processing cell type:", ct, "\n")
    
    cells_ct <- rownames(meta)[meta[[celltype_col]] == ct]
    cond_tab <- table(meta[cells_ct, condition_col])
    
    if (!all(c("Hypertensive CKD (H-CKD)", "Normal Reference") %in% names(cond_tab))) {
      cat("  H-CKD or Normal missing → skipping\n")
      next
    }
    
    if (any(cond_tab[c("Hypertensive CKD (H-CKD)", "Normal Reference")] < min_cells_per_group)) {
      cat("  Not enough cells → skipping\n")
      next
    }
    
    seu_ct <- subset(seu, cells = cells_ct)
    DefaultAssay(seu_ct) <- "originalexp"
    seu_ct <- NormalizeData(seu_ct, verbose = FALSE)
    Idents(seu_ct) <- condition_col
    
    deg <- FindMarkers(
      seu_ct,
      ident.1 = "Hypertensive CKD (H-CKD)",
      ident.2 = "Normal Reference",
      test.use = "MAST",
      logfc.threshold = 0,
      min.pct = 0.1,
      latent.vars = c("nCount_RNA", "percent.mt")
    )
    
    deg$gene <- rownames(deg)
    
    ct_clean <- gsub("[^a-zA-Z0-9]", "_", ct)
    out_file <- file.path(outdir, paste0("DEG_", ct_clean, "_HCKD_vs_Normal.csv"))
    write.csv(deg, out_file, row.names = FALSE)
    cat("  ✔ Saved:", out_file, "\n")
  }
}
run_mast_per_celltype_HCKD_vs_Normal(
  seu = seurat_obj,
  celltype_col = "majority_voting",
  condition_col = "Condition"
)
```{r DEG summary}
BiocManager::install("ReactomePA")
library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(enrichplot)
deg_path <- "output/DEG_MAST_AKI_vs_Normal/"
deg_files <- list.files(deg_path, pattern = "^DEG_.*\\.csv$", full.names = TRUE)

deg_summary <- data.frame(
  CellType = character(),
  Upregulated = integer(),
  Downregulated = integer(),
  Total_DEGs = integer(),
  stringsAsFactors = FALSE
)

for (file in deg_files) {
  
  deg_data <- read.csv(file)
  
  if (!all(c("avg_log2FC", "p_val_adj") %in% colnames(deg_data))) {
    warning(paste("Skipping:", basename(file)))
    next
  }
  
  up   <- sum(deg_data$avg_log2FC > 0 & deg_data$p_val_adj < 0.05)
  down <- sum(deg_data$avg_log2FC < 0 & deg_data$p_val_adj < 0.05)
  
  celltype <- gsub("^DEG_(.*?)_AKI_vs_Normal\\.csv$", "\\1", basename(file))
  
  deg_summary <- rbind(
    deg_summary,
    data.frame(CellType = celltype,
               Upregulated = up,
               Downregulated = down,
               Total_DEGs = up + down)
  )
}

deg_summary <- deg_summary %>% arrange(desc(Total_DEGs))
print(deg_summary)

write.csv(deg_summary,
          "output/DEG_MAST_AKI_vs_Normal/DEGsummary_table.csv",
          row.names = FALSE)
deg_path <- "output/DEG_MAST_DKD_vs_Normal/"
deg_files <- list.files(deg_path, pattern = "^DEG_.*\\.csv$", full.names = TRUE)

deg_summary <- data.frame(
  CellType = character(),
  Upregulated = integer(),
  Downregulated = integer(),
  Total_DEGs = integer(),
  stringsAsFactors = FALSE
)

for (file in deg_files) {
  
  deg_data <- read.csv(file)
  
  if (!all(c("avg_log2FC", "p_val_adj") %in% colnames(deg_data))) {
    warning(paste("Skipping:", basename(file)))
    next
  }
  
  up   <- sum(deg_data$avg_log2FC > 0 & deg_data$p_val_adj < 0.05)
  down <- sum(deg_data$avg_log2FC < 0 & deg_data$p_val_adj < 0.05)
  
  celltype <- gsub("^DEG_(.*?)_DKD_vs_Normal\\.csv$", "\\1", basename(file))
  
  deg_summary <- rbind(
    deg_summary,
    data.frame(CellType = celltype,
               Upregulated = up,
               Downregulated = down,
               Total_DEGs = up + down)
  )
}

deg_summary <- deg_summary %>% arrange(desc(Total_DEGs))
print(deg_summary)

write.csv(deg_summary,
          "output/DEG_MAST_DKD_vs_Normal/DEGsummary_table.csv",
          row.names = FALSE)
deg_path <- "output/DEG_MAST_HCKD_vs_Normal/"
deg_files <- list.files(deg_path, pattern = "^DEG_.*\\.csv$", full.names = TRUE)

deg_summary <- data.frame(
  CellType = character(),
  Upregulated = integer(),
  Downregulated = integer(),
  Total_DEGs = integer(),
  stringsAsFactors = FALSE
)

for (file in deg_files) {
  
  deg_data <- read.csv(file)
  
  if (!all(c("avg_log2FC", "p_val_adj") %in% colnames(deg_data))) {
    warning(paste("Skipping:", basename(file)))
    next
  }
  
  up   <- sum(deg_data$avg_log2FC > 0 & deg_data$p_val_adj < 0.05)
  down <- sum(deg_data$avg_log2FC < 0 & deg_data$p_val_adj < 0.05)
  
  celltype <- gsub("^DEG_(.*?)_HCKD_vs_Normal\\.csv$", "\\1", basename(file))
  
  deg_summary <- rbind(
    deg_summary,
    data.frame(CellType = celltype,
               Upregulated = up,
               Downregulated = down,
               Total_DEGs = up + down)
  )
}

deg_summary <- deg_summary %>% arrange(desc(Total_DEGs))
print(deg_summary)

write.csv(deg_summary,
          "output/DEG_MAST_HCKD_vs_Normal/DEGsummary_table.csv",
          row.names = FALSE)
```

```{r Pathways}
library(dplyr)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)

# ================================
# Sort DEG files by smallest size
# ================================
file_info <- file.info(deg_files)
deg_files <- deg_files[order(file_info$size)]

file_sizes <- data.frame(
  File = basename(deg_files),
  Size_MB = round(file.info(deg_files)$size / 1024^2, 2)
)

cat("\n============================================\n")
cat(" DEG FILES SORTED BY SIZE (SMALLEST → LARGEST)\n")
cat("============================================\n")
print(file_sizes)

# OPTIONAL: test first few files only
# deg_files <- deg_files[1:3]

# ================================
# ================================
# Sort DEG files by smallest size
# ================================
file_info <- file.info(deg_files)
deg_files <- deg_files[order(file_info$size)]

file_sizes <- data.frame(
  File = basename(deg_files),
  Size_MB = round(file.info(deg_files)$size / 1024^2, 2)
)

cat("\n============================================\n")
cat(" DEG FILES SORTED BY SIZE (SMALLEST → LARGEST)\n")
cat("============================================\n")
print(file_sizes)

# OPTIONAL: test first few files only
# deg_files <- deg_files[1:3]

# ================================
# Initialize outputs
# ================================
all_pathways_df <- data.frame()
failed_files <- c()
skipped_files <- c()

total_files <- length(deg_files)
overall_start <- Sys.time()

cat("\n============================================\n")
cat(" STARTING REACTOME ORA FOR", total_files, "FILES\n")
cat(" Start time:", as.character(overall_start), "\n")
cat("============================================\n\n")

# ================================
# Loop over DEG files
# ================================
for (i in seq_along(deg_files)) {
  
  file_start <- Sys.time()
  file_path <- deg_files[i]
  file_name <- basename(file_path)
  
  cat("\n====================================================\n")
  cat(" FILE", i, "OF", total_files, ":", file_name, "\n")
  cat(" Started at:", as.character(Sys.time()), "\n")
  cat("====================================================\n")
  
  # ---- Step 1: Extract cell type ----
  cat("[1/9] Extracting cell type...\n")
  celltype <- sub("^DEG_(.*?)_Diabetic_vs_Healthy\\.csv$", "\\1", file_name)
  
  if (celltype == file_name) {
    warning(" Celltype extraction FAILED for file: ", file_name)
    failed_files <- c(failed_files, file_name)
    next
  }
  
  celltype_clean <- gsub("[^a-zA-Z0-9]", "_", celltype)
  condition <- "Diabetic_vs_Healthy"
  
  cat(" Cell type identified:", celltype, "\n")
  
  # ---- Step 2: Read DEG file ----
  cat("[2/9] Reading DEG file...\n")
  step_time <- Sys.time()
  
  res <- tryCatch(
    read.csv(file_path, stringsAsFactors = FALSE),
    error = function(e) {
      cat(" FAILED to read file:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(res)) {
    failed_files <- c(failed_files, file_name)
    next
  }
  
  cat(" DEG rows read:", nrow(res), "\n")
  cat(" Columns found:", paste(colnames(res), collapse = ", "), "\n")
  cat(" Time taken:", round(difftime(Sys.time(), step_time, units = "secs"), 2), "sec\n")
  
  # ---- Step 3: Check required columns ----
  cat("[3/9] Checking required columns...\n")
  if (!"gene" %in% colnames(res)) {
    warning(" Missing required 'gene' column in: ", file_name)
    skipped_files <- c(skipped_files, file_name)
    next
  }
  cat(" Required column found\n")
  
  cat(" Unique genes in file:", length(unique(res$gene)), "\n")
  
  # ---- Step 4: SYMBOL → ENTREZ mapping ----
  cat("[4/9] Mapping SYMBOL to ENTREZ...\n")
  step_time <- Sys.time()
  
  ncbi_map <- tryCatch(
    suppressMessages(
      clusterProfiler::bitr(
        unique(res$gene),
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      )
    ),
    error = function(e) {
      cat(" Gene mapping FAILED:", e$message, "\n")
      return(NULL)
    }
  )
  
  if (is.null(ncbi_map)) {
    failed_files <- c(failed_files, file_name)
    next
  }
  
  cat(" Genes mapped to ENTREZ:", nrow(ncbi_map), "\n")
  cat(" Time taken:", round(difftime(Sys.time(), step_time, units = "secs"), 2), "sec\n")
  
  # ---- Step 5: Prepare gene vector ----
  cat("[5/9] Preparing unique ENTREZ gene set...\n")
  step_time <- Sys.time()
  
  gene_df <- res %>%
    left_join(ncbi_map, by = c("gene" = "SYMBOL")) %>%
    filter(!is.na(ENTREZID)) %>%
    distinct(ENTREZID, .keep_all = TRUE)
  
  gene_ids <- unique(gene_df$ENTREZID)
  
  cat(" Genes after filtering:", length(gene_ids), "\n")
  cat(" Time taken:", round(difftime(Sys.time(), step_time, units = "secs"), 2), "sec\n")
  
  if (length(gene_ids) < 10) {
    cat(" SKIPPING ORA: too few mapped genes\n")
    skipped_files <- c(skipped_files, file_name)
    next
  }
  
  # ---- Step 6: Run Reactome ORA ----
  cat("[6/9] Running Reactome ORA...\n")
  step_time <- Sys.time()
  
  ora_res <- tryCatch(
    enrichPathway(
      gene          = gene_ids,
      organism      = "human",
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      minGSSize     = 10,
      maxGSSize     = 500,
      readable      = TRUE
    ),
    error = function(e) {
      cat(" ORA FAILED:", e$message, "\n")
      return(NULL)
    }
  )
  
  cat(" ORA step finished in:",
      round(difftime(Sys.time(), step_time, units = "secs"), 2), "sec\n")
  
  if (is.null(ora_res) || nrow(as.data.frame(ora_res)) == 0) {
    cat(" No significant pathways found\n")
    skipped_files <- c(skipped_files, file_name)
    next
  }
  
  # ---- Step 7: Extract results ----
  cat("[7/9] Extracting pathway results...\n")
  step_time <- Sys.time()
  
  pathways <- as.data.frame(ora_res) %>%
    filter(!is.na(Description), !is.na(p.adjust), !is.na(Count))
  
  cat(" Valid pathways:", nrow(pathways), "\n")
  cat(" Time taken:", round(difftime(Sys.time(), step_time, units = "secs"), 2), "sec\n")
  
  if (nrow(pathways) == 0) {
    skipped_files <- c(skipped_files, file_name)
    next
  }
  
  # ---- Step 8: Add metadata ----
  cat("[8/9] Adding metadata and combining...\n")
  step_time <- Sys.time()
  
  pathways$CellType  <- celltype
  pathways$Condition <- condition
  
  all_pathways_df <- bind_rows(all_pathways_df, pathways)
  
  cat(" Time taken:", round(difftime(Sys.time(), step_time, units = "secs"), 2), "sec\n")
  
  # ---- Step 9: Save per-cell CSV ----
  cat("[9/9] Saving outputs...\n")
  step_time <- Sys.time()
  
  out_csv <- file.path(
    deg_path,
    paste0(celltype_clean, "_Reactome_ORA_", condition, ".csv")
  )
  
  write.csv(pathways, out_csv, row.names = FALSE)
  
  cat(" ORA CSV written:", out_csv, "\n")
  cat(" Time taken:", round(difftime(Sys.time(), step_time, units = "secs"), 2), "sec\n")
  
  # ---- File summary ----
  file_end <- Sys.time()
  cat(" FILE COMPLETED:", file_name, "\n")
  cat(" Total file time:",
      round(difftime(file_end, file_start, units = "mins"), 2), "minutes\n")
  
  files_left <- total_files - i
  cat(" Files remaining:", files_left, "\n")
}

# ================================
# Save combined results
# ================================
cat("\n============================================\n")
cat(" SAVING FINAL OUTPUTS\n")
cat("============================================\n")

combined_out <- file.path(deg_path, "All_CellTypes_Reactome_ORA_Diabetic_vs_Healthy.csv")
write.csv(all_pathways_df, combined_out, row.names = FALSE)

cat(" Combined Reactome ORA file saved:\n", combined_out, "\n")

if (length(failed_files) > 0) {
  writeLines(failed_files, file.path(deg_path, "FAILED_ORA_FILES.txt"))
  cat(" Failed file log saved\n")
}

if (length(skipped_files) > 0) {
  writeLines(skipped_files, file.path(deg_path, "SKIPPED_ORA_FILES.txt"))
  cat(" Skipped file log saved\n")
}

overall_end <- Sys.time()

cat("\n============================================\n")
cat(" ALL FILES COMPLETED\n")
cat(" End time:", as.character(overall_end), "\n")
cat(" Total runtime:",
    round(difftime(overall_end, overall_start, units = "mins"), 2), "minutes\n")
cat(" Total successful pathway rows:", nrow(all_pathways_df), "\n")
cat("============================================\n")
```


```{r Visualizations for Pathways}
library(dplyr)
library(ggplot2)

ora_png_files <- c()

cat("\n============================================\n")
cat(" GENERATING REACTOME ORA PLOTS\n")
cat("============================================\n")

for (ct in unique(all_pathways_df$CellType)) {
  
  cat("\n--------------------------------------------\n")
  cat(" Plotting Cell Type:", ct, "\n")
  cat("--------------------------------------------\n")
  
  df_ct <- all_pathways_df %>%
    filter(CellType == ct)
  
  if (nrow(df_ct) == 0) {
    cat(" No pathways to plot for", ct, "\n")
    next
  }
  
  cat(" Total pathways available:", nrow(df_ct), "\n")
  
  ## ---- Select top pathways ----
  top_pathways <- df_ct %>%
    arrange(p.adjust, desc(Count)) %>%
    slice_head(n = 10)
  
  if (nrow(top_pathways) == 0) {
    cat(" No top pathways found for", ct, "\n")
    next
  }
  
  cat(" Top pathways selected:", nrow(top_pathways), "\n")
  
  ## ---- Convert GeneRatio to numeric ----
  top_pathways <- top_pathways %>%
    mutate(
      GeneRatio_num = sapply(GeneRatio, function(x) {
        vals <- strsplit(x, "/")[[1]]
        as.numeric(vals[1]) / as.numeric(vals[2])
      })
    )
  
  ## ---- Order pathways ----
  top_pathways$Description <- factor(
    top_pathways$Description,
    levels = rev(top_pathways$Description[order(top_pathways$GeneRatio_num)])
  )
  
  ## ---- Plot ----
  p <- ggplot(
    top_pathways,
    aes(x = GeneRatio_num, y = Description, color = p.adjust, size = Count)
  ) +
    geom_point() +
    theme_minimal(base_size = 11) +
    labs(
      title = paste("Reactome ORA -", ct),
      x = "Gene Ratio",
      y = "Pathway",
      color = "Adjusted p-value",
      size = "Gene Count"
    ) +
    theme(
      axis.text.y = element_text(size = 8, face = "bold"),
      axis.text.x = element_text(size = 8, face = "bold"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  
  ## ---- Save PNG ----
  ct_clean <- gsub("[^a-zA-Z0-9]", "_", ct)
  
  png_file <- file.path(
    deg_path,
    paste0(ct_clean, "_Reactome_ORA_dotplot.png")
  )
  
  ggsave(
    filename = png_file,
    plot = p,
    width = 9,
    height = 6,
    dpi = 300
  )
  
  cat(" Plot saved:", png_file, "\n")
  
  ora_png_files <- c(ora_png_files, png_file)
}

cat("\n============================================\n")
cat(" ALL REACTOME ORA PLOTS COMPLETED\n")
cat(" Total PNG files created:", length(ora_png_files), "\n")
cat("============================================\n")
```

