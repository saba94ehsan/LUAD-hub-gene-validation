# LUAD hub-gene study: analysis code

Code and instructions for the manuscript **"Integrative Transcriptomic and Network-Based
Identification of Cell-Cycle Hub Genes in Lung Adenocarcinoma"** (revised submission to
*Cell Biochemistry and Biophysics*).

**Version used for the revised manuscript: release `v1.0-revision`.**
Please use this release, not the `main` branch, to reproduce the numbers in the manuscript:
https://github.com/saba94ehsan/LUAD-hub-gene-validation/releases/tag/v1.0-revision

## What is in this repository

| Path | Content |
|---|---|
| `LUAD_hub_gene_analysis.R` | Complete R analysis in one script (Sections 1-6, see below) |
| `metadata/coldata.csv` | GSE81089 sample table (sample, condition, patient); inclusion and exclusion of the samples |
| `metadata/clinical_GSE81089.csv` | GSE81089 clinical data of the 127 samples (from GEO) |
| `metadata/TCGA_stage_GDC_release46.csv` | Stage of the 504 TCGA-LUAD patients from the GDC clinical table (GDC Data Release 46.0, accessed 7 November 2026); columns `patient12`, `stage_GDC` |
| `metadata/GDC_release_info.txt` | GDC release information |
| `metadata/Supplementary_Table_S6_cytoHubba.xlsx` | cytoHubba results (sheets `All_nodes`, `EPC_repeated_runs`, `Sensitivity_analysis`); input of Section 5 |
| `metadata/pathology.tsv` | Human Protein Atlas pathology table (version 23); input of Section 6 |
| `01_DESeq2/vsd_main.rds`, `01_DESeq2/res_sym_main.rds` | Intermediate objects of Section 1.11 from the verified run (used by Sections 2-4; rebuilt from Section 1 if `REBUILD_INTERMEDIATE <- TRUE`) |
| `05_survival/TCGA_survival_cohort_with_stage_source.csv` | **Sample list of the TCGA survival analysis**: one row per patient (504), with the stage from the data object and the completed stage; shows who is included in the 486-patient and 379-patient models |
| `03_PPI/` | Cytoscape input and output (see "Network analysis") |
| `00_session/` | `sessionInfo_full.txt`, `sessionInfo.txt` and the four `reference_check*.csv` files of the verified run (every row `OK`) |

Large inputs are not stored in the repository:

* the GSE81089 count file `GSE81089_readcounts_featurecounts.tsv.gz` (GEO series GSE81089);
* the TCGA-LUAD `SummarizedExperiment` `04_TCGA/luad_SE.rds`, a TCGAbiolinks download from the
  Genomic Data Commons (project TCGA-LUAD, Transcriptome Profiling, Gene Expression
  Quantification, STAR counts; assays `unstranded`, `stranded_first`, `stranded_second`,
  `tpm_unstrand`, `fpkm_unstrand`, `fpkm_uq_unstrand`).

## How to run

1. Create a project folder and put the script, the count file and the folders `metadata/` and
   `04_TCGA/` (with `luad_SE.rds`) in it. Copy the two `.rds` files into `01_DESeq2/`.
2. Open R 4.6.1 in this folder (`setwd("path/to/project")`) and run the script section by
   section in RStudio, or run `Rscript LUAD_hub_gene_analysis.R`. Each section reads what the
   previous section wrote to disk. The script creates all output folders itself.
3. Install the packages listed under "Software" first (Bioconductor 3.23).
4. Check `00_session/reference_check.csv` (Section 1), `reference_check_TCGA.csv` (Section 3),
   `reference_check_cytoHubba.csv` (Section 5) and `reference_check_HPA.csv` (Section 6): every
   row should say `OK`.

Notes:

* Section 1.8 (non-coding RNAs among the DEGs without a current gene symbol) needs the Ensembl
  archive (Ensembl 73, GRCh37). If the archive is not reachable, the script prints
  "1.8 skipped" and continues. This step does not change any other result; the PPI input
  `03_PPI/up_DEG_symbols.txt` is written by Section 1.9 in any case.
* The stage file `metadata/TCGA_stage_GDC_release46.csv` is read in Section 3.2, so no internet
  connection is needed for it. Section 3.8 (tumour purity, Aran et al.) uses the TCGAbiolinks
  package.

## Which part of the script produces which figure, table and number

| Manuscript item | Script section | Output file |
|---|---|---|
| GSE81089 DEGs, PCA | 1.1-1.4 | `01_DESeq2/DESeq2_all_genes.csv`, `01_DESeq2/DEGs_with_symbol.csv` |
| Matched-pair analysis | 1.5 | `07_tables_for_manuscript/Supplementary_Table_S2_LUAD_DESeq2_PAIRED_annotated_ALL_genes.csv` |
| GO (Table 5, Suppl. Table S3) | 1.6 | `02_GO_GSEA/topGO_all_terms.csv`, `07_tables_for_manuscript/Supplementary_Table_S3_GO_ALL_ontologies_ALL_terms.csv` |
| GSEA (Suppl. Table S4) | 1.7 | `02_GO_GSEA/fgsea_all_sets.csv`, `07_tables_for_manuscript/Supplementary_Table_S4_GSEA_ALL_databases_ALL_pathways.csv` |
| PPI input (up-regulated DEGs) | 1.9 | `03_PPI/up_DEG_symbols.txt` |
| Verification of the Section 1 numbers | 1.10 | `00_session/reference_check.csv` |
| Suppl. Table S8 (hub genes in the 9 matched pairs, tumour cell %) | 1.12 | `07_tables_for_manuscript/TableS8_paired_and_TCpct_hubgenes.csv` |
| Fig. 1 (PCA, volcano, MA), Fig. 2 (GO), Fig. 3 (GSEA), Suppl. Table S1 | 1.13 | `08_figures/Fig1_*`, `Fig2_*`, `Fig3_*`; `07_tables_for_manuscript/Supplementary_Table_S1_*` |
| TCGA-LUAD DEGs (Suppl. Table S7) | 2.3-2.4 | `04_TCGA/TCGA_DESeq2_annotated_ALL_genes.csv`, `TCGA_DEG_UP.csv`, `TCGA_DEG_DOWN.csv` |
| Table 8 (hub genes, GSE81089 vs TCGA) | 2.5 | `07_tables_for_manuscript/Table8_hubgenes_GEO_vs_TCGA.csv` |
| Suppl. Fig. S2 (GSE81089 vs TCGA log2FC) | 2.6 | `08_figures/FigS2_GEO_vs_TCGA_log2FC.png` |
| Table 3 (TCGA clinical data) | 3.2 | `07_tables_for_manuscript/Table3_TCGA_clinical.csv` |
| Suppl. Fig. S3 (TCGA Kaplan-Meier) | 3.3 | `08_figures/FigS3_TCGA_KM_hubgenes.png` / `.pdf` |
| Cox models, Table 9 (TCGA rows) | 3.4-3.5 | `07_tables_for_manuscript/Table9_TCGA_rows.csv`, `05_survival/` |
| Sensitivity analysis, stage from the data object only (379 patients) | 3.6 | `05_survival/TCGA_cox_adjusted_colData_stage_only_379.csv` |
| Suppl. Table S9 (included vs excluded, GDC stage completion) | 3.7 | `07_tables_for_manuscript/TableS9_included_vs_excluded.csv`, `TableS9_GDC_stage_completion.csv` |
| Suppl. Table S10 (proliferation score, purity, C-index) | 3.8 | `05_survival/TableS10_cox_proliferation_purity.csv`, `TableS10_incremental_value.csv` |
| Table 9 (GSE81089 rows), Fig. 6 (Kaplan-Meier, GSE81089) | 4.1-4.3 | `07_tables_for_manuscript/Table9_GSE81089_rows.csv`, `08_figures/Fig6_GSE81089_KM.*` |
| Fig. 7 (forest plot) | 4.4 | `08_figures/Fig7_forest_Cox.*` |
| Fig. 5 (hub-gene expression) | 4.5 | `08_figures/Fig5_hub_gene_expression.*` |
| Fig. 4B (UpSet plot), Fig. 4C (rank heatmap) | 5.2-5.3 | `08_figures/Fig4B_UpSet_cytoHubba.png`, `Fig4C_rank_heatmap.png` |
| Suppl. Table S11 and Suppl. Fig. S5A (HPA counts) | 6.1-6.2 | `07_tables_for_manuscript/Supplementary_Table_S11_HPA_lung_cancer_counts.csv`, `08_figures/FigS5A_HPA_counts.png` |

The figures of the manuscript were arranged and labelled (panel letters, fonts) from the same data. The script draws Fig. 6 and Suppl. Fig. S3 with panels A CDK1, B CCNA2, C CCNB1 and D BRCA1, as in the manuscript. The numbers in them are identical to the
numbers in the output files above.

### Key numbers and where they come from

| Number in the manuscript | Source |
|---|---|
| GSE81089: 108 LUAD and 19 normal samples; 4,840 DEGs (3,089 up, 1,751 down); 3,087 unique up-regulated symbols | Section 1, `reference_check.csv` |
| Nine matched pairs; 4,535 paired DEGs, 99.5% in the same direction | Section 1.5, `reference_check.csv` |
| TCGA-LUAD: 517 tumours and 59 normal samples | Section 2.2 |
| TCGA survival cohort 504 patients, 182 deaths | Section 3.2, `reference_check_TCGA.csv` |
| Main adjusted model 486 patients, 177 deaths (stage completed from GDC for 112 of 120 patients) | Sections 3.2 and 3.4, `TableS9_GDC_stage_completion.csv` |
| Sensitivity model 379 patients, 123 deaths; 125 excluded patients | Sections 3.6 and 3.7 |
| Proliferation score (323 genes), rho 0.85-0.91, C-index 0.672 to 0.685-0.691 | Section 3.8 |
| GSE81089 survival cohort 107 patients, 48 deaths | Section 4.1 |
| Four genes in the top 20 of all five cytoHubba algorithms; 14 genes in at least three | Section 5, `reference_check_cytoHubba.csv` |
| HPA counts (CDK1 5/12, CCNB1 8/12, CCNA2 5/11, BRCA1 3/12 with medium or high staining) | Section 6, `reference_check_HPA.csv` |

## Network analysis (Fig. 4A, Suppl. Fig. S1, Suppl. Tables S5-S6) and web tools

These steps were done outside R and are described here so that they can be repeated.

* **Input:** `03_PPI/up_DEG_symbols.txt` (up-regulated DEGs, written by Section 1.9).
* **Network:** Cytoscape 3.10.2 with stringApp 2.2.0, STRING database
  version 12.5, protein query, *Homo sapiens*, confidence score >= 0.4 (main analysis); 2,068
  nodes and 23,988 edges after removal of unconnected nodes.
* **Module:** MCODE 2.0.3 (degree cutoff 2, node score cutoff 0.2, K-core 2,
  maximum depth 100; haircut on, fluff off); highest-scoring cluster (score 106.0; 117 genes,
  6,148 edges). Node and edge tables of this cluster: `03_PPI/cluster1_nodes.csv`,
  `03_PPI/cluster1_edges.csv`.
* **Hub genes:** cytoHubba 0.1, algorithms Degree, MNC, EPC, Radiality and
  Closeness; genes in the top 20 of all five (tied scores included). MCC was not used for
  ranking because its scores were saturated in this network. EPC top-20 membership varied
  between repeated runs on the same network; the repeated runs are listed in
  `metadata/Supplementary_Table_S6_cytoHubba.xlsx`, sheet `EPC_repeated_runs`.
* **Sensitivity analyses:** confidence >= 0.70 and >= 0.90, experimental and database evidence
  only, and a network that includes down-regulated genes (4,837 gene symbols); Degree, MNC,
  Radiality and Closeness (top 20); results in sheet `Sensitivity_analysis` of the same file.
* **Session file:** `03_PPI/LUAD_up_STRING_MCODE.cys`.
* **Web tools (no code):** cBioPortal (Fig. 8), UALCAN (Fig. 9), TIMER 2.0 (Table 10, Fig. 10,
  Suppl. Fig. S4) and the Human Protein Atlas (Suppl. Fig. S5B-E). Fig. 10 is drawn from the
  purity-adjusted partial correlations of Table 10. The HPA counts (Suppl. Table S11 and
  Suppl. Fig. S5A) are computed by Section 6 from `metadata/pathology.tsv`.

## Software

Windows 11, R 4.6.1, Bioconductor 3.23, DESeq2 1.52.0, apeglm 1.34.0, topGO 2.64.0,
fgsea 1.38.0, org.Hs.eg.db 3.23.1, GO.db 3.23.1, msigdbr 26.1.1, survival 3.8-12,
survminer 0.5.2, ggplot2 4.0.4, patchwork 1.3.2, UpSetR 1.4.1, TCGAbiolinks 2.40.0.
Full session information: `00_session/sessionInfo_full.txt`.

The gene symbols and the GO annotation depend on the version of `org.Hs.eg.db`; the numbers
in the manuscript were reproduced with version 3.23.1.

## Data sources

GSE81089 (GEO); TCGA-LUAD (Genomic Data Commons, Data Release 46.0 for the clinical stage
table); MSigDB Hallmark, KEGG and Reactome gene sets (msigdbr); Human Protein Atlas
(version 23); TCGA tumour purity (Aran et al., Nat Commun 2015;6:8971, via TCGAbiolinks).
