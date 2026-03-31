# Human Kidney scRNA-Seq Analysis (GSE183276)

## Overview
This repository contains the analysis pipeline for **adult human kidney single-cell RNA-seq (scRNA-seq) data** (GSE183276), including:

- Data preprocessing & QC using Seurat  
- Batch correction & clustering (PCA, UMAP, Harmony)  
- Cell type annotation with Python (`CellTypist`)  
- Differential expression & pathway analysis (MAST, ReactomePA)  

## Dataset
- Source: GEO [GSE183276](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183276)  
- Conditions: AKI, Normal, DKD, Hypertensive CKD  

## Workflow
1. **QC & preprocessing** – Filter cells/genes, remove outliers.  
2. **Clustering & batch correction** – PCA, UMAP, Harmony integration.  
3. **Cell type annotation** – Using pre-trained classifier `.pkl` file.  
4. **DEG & pathway analysis** – Identify DEGs per cell type and perform Reactome ORA.  

## Outputs
- Filtered Seurat objects (`.rds`)  
- QC and UMAP plots  
- DEG tables and pathway visualizations  

## Author
**Iswarya Umasankar** – M.Tech Computational Biology, Pondicherry University  

## Acknowledgments
Thanks to **Smriti Arora**, instructor of the SingleCellRNASeq Workshop (Mar 2026), for teaching the analysis workflow.
