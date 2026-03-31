**Single-Cell Transcriptomic Analysis of Kidney Disease**
**Overview**

This project performs cell-type specific differential gene expression and pathway enrichment analysis using single-cell RNA sequencing data from kidney samples.

**Dataset**
Public dataset from GEO
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
**Key Features**
Automated per-cell-type DEG export
Scalable multi-condition analysis
Reactome ORA integration
Publication-ready dot plots
**Tools Used**
R
Seurat
MAST
clusterProfiler
ReactomePA
org.Hs.eg.db
ggplot2
**Example Output**
<img width="2700" height="1800" alt="PC_Reactome_ORA_dotplot" src="https://github.com/user-attachments/assets/b915ff8b-de1c-4653-9429-1d2193ba169d" />


**Author**

Iswarya Umasankar
M.Tech Computational Biology
Pondicherry University
