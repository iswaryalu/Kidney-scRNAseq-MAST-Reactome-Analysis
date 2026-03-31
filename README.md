**Single-Cell Transcriptomic Analysis of Kidney Disease**
# Human Kidney Single-Cell RNA-Seq Analysis (GSE183276)

**Overview**
This repository contains the code and analysis pipeline for single-cell RNA-seq (scRNA-seq) data from adult human kidney samples. The workflow includes:

**Data preprocessing & QC** – Creation of Seurat objects, filtering cells/genes, visualizations.

**Batch correction & clustering** – PCA, UMAP, Harmony integration, cluster identification.

**Cell type annotation** – Using pre-trained classifier (.pkl file) in Python.

**Differential expression & pathway analysis** – Identification of DEGs across conditions and enrichment analysis.

All analyses were performed in R (Seurat, MAST, Harmony, SingleCellExperiment,clusterProfiler,ReactomePA) and Python (CellTypist).

**Dataset**
Public dataset from GEO - (GSE183276)

**Conditions:**
Acute Kidney Injury (AKI)
Normal Reference
Diabetic Kidney Disease (DKD)
Hypertensive CKD

**Methods**
Preprocessing & QC using Seurat
Mitochondrial filtering (percent.mt)
Cell-type annotation (majority voting)
Cell-type specific DEG analysis using MAST
Reactome pathway over-representation analysis
Automated visualization of enriched pathways
**Workflow**

scRNA Data
      
      ↓

Quality Control
 
      ↓

Cell Type Annotation
    
      ↓

MAST DEG (per cell type)
       
       ↓

Gene Mapping (SYMBOL → ENTREZ)
    
      ↓

Reactome ORA
      
      ↓

Pathway Visualization

**Workflow Summary**

**1. Seurat Object Creation & QC**

Load individual Seurat objects, check for duplicate genes.

Merge objects and add metadata (Condition).

Visualize QC metrics (violin plots, density scatter plots).

Filter cells based on:

nFeature_RNA (min 200, max 7000)

percent.mt (<20%)

percent.rb (<20%)

Mahalanobis distance for outliers

Keep only protein-coding genes.

**Outputs:**

02_GSE279086_seurat_qc_filtered.rds

QC metrics CSV

Pre- and post-QC plots

**2. Clustering & Batch Correction**

Normalize and scale data (LogNormalize, variable features).

PCA and elbow plot for dimensionality selection.

UMAP visualization.

FindNeighbors and FindClusters (resolution = 0.8).

Batch correction using Harmony.

Visualize UMAP by sample and condition.

Compute batch mixing metrics before and after Harmony.

**Outputs:**

03_GSE279086_seurat_pca_umap.rds

04_GSE279086_harmony_corrected.rds

UMAP plots and batch mixing barplots

**3. Cell Type Annotation**

Use Adult_Human_Kidney.pkl in Python via CellTypist.

Add predicted cell types back to Seurat object.

**4. Differential Expression & Pathway Analysis**

Identify DEGs using MAST.

Saved results for pathway enrichment and visualization.


**Example Output**

<img width="2700" height="1800" alt="PC_Reactome_ORA_dotplot" src="https://github.com/user-attachments/assets/b915ff8b-de1c-4653-9429-1d2193ba169d" />


**Author**

Iswarya Umasankar
M.Tech Computational Biology
Pondicherry University

## Acknowledgments
I would like to thank the instructor Smriti Arora of the *SingleCellRNASeq Hands-On Workshop (Mar 2026)* for teaching the principles and practical steps of single-cell RNA-seq analysis. While the workshop used a different dataset, I applied the learned methods to the GSE183276 dataset for this project.
