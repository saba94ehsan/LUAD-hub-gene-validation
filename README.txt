README - analysis code
======================

Manuscript : Integrative Transcriptomic and Network-Based Identification of
             Cell-Cycle Hub Genes in Lung Adenocarcinoma
Discovery  : GSE81089 (108 LUAD + 19 normal lung), GEO
Validation : TCGA-LUAD (517 primary tumours + 59 normal), GDC
             (download date: [13-Oct-2026])

**Scope of reproducibility**

The R script reproduces the transcriptomic analyses, including GSE81089 differential gene expression analysis, PCA, paired-sample analysis, GO enrichment, GSEA, TCGA-LUAD validation, and survival analyses, as specified in the script. The protein–protein interaction (PPI) network construction, cytoHubba hub-gene ranking, MCODE cluster detection, and analyses using TIMER, UALCAN, cBioPortal, and the Human Protein Atlas (HPA) were performed separately and are not reproduced by this script. These analyses should be documented separately, with the corresponding methods, parameters, and source information provided where applicable:
    LUAD_hub_gene_analysis.R

Software : R 4.6.1, Bioconductor 3.23 (Ubuntu 22.04). Exact package versions
           are in expected_output/sessionInfo*.txt
Memory   : at least 8 GB RAM (TCGA part).


FOLDER CONTENT
--------------
LUAD_code/
  LUAD_hub_gene_analysis.R            the analysis script (4 sections)
  README.txt
  GSE81089_readcounts_featurecounts.tsv.gz   [ADD: download from GEO, series GSE81089]
  metadata/
    coldata.csv                       sample table (sample, condition, patient)
    clinical_GSE81089.csv             clinical data of the 127 samples (see Section 1.11)
  01_DESeq2/
    vsd_main.rds, res_sym_main.rds    intermediate objects of GSE81089 (see Section 1.11)
  04_TCGA/
    luad_SE.rds                       [ADD: TCGA-LUAD RNA-seq from GDC (STAR counts),
                                       SummarizedExperiment; see note below]
  expected_output/
    reference_check.csv               numbers the script must reproduce (22 rows, all OK)
    sessionInfo.txt, sessionInfo_5b_linux.txt


HOW TO RUN
----------
1. Put the files above in one folder (the zip already has this layout).
2. Open LUAD_hub_gene_analysis.R and change the line
        PROJ_DIR <- "~/LUAD_code"
   to the path of that folder (use "/" in the path, also on Windows).
3. Run the whole script:
        Rscript LUAD_hub_gene_analysis.R
   or open it in RStudio and run it section by section (Ctrl + Enter).
   Every section reads what the previous section wrote to disk.
4. Compare 00_session/reference_check.csv (written by the script) with
   expected_output/reference_check.csv: all rows must say OK.


WHAT THE SCRIPT DOES
--------------------
SECTION 1  GSE81089: DESeq2 (LUAD vs normal), PCA, paired analysis (9 pairs),
           topGO, fgsea, input list for STRING.
           Input : counts file, metadata/coldata.csv
           Output: 01_DESeq2/, 02_GO_GSEA/, 03_PPI/, 00_session/
           Check : 00_session/reference_check.csv (108 + 19 samples, 26,184
                   genes, 4,840 DEGs = 3,089 up + 1,751 down, 200 activated and
                   144 suppressed gene sets); last line of
                   00_session/sessionInfo.txt gives the tumour-cell % vs PC1
                   correlation (rho = -0.01, p = 0.91).
SECTION 2  TCGA-LUAD: DESeq2 and comparison with GSE81089 (Table 8, Fig S3).
           Input : 04_TCGA/luad_SE.rds, 01_DESeq2/res_sym_main.rds
           Output: 04_TCGA/, 07_tables_for_manuscript/Table8_*.csv, 08_figures/
SECTION 3  TCGA-LUAD survival: Kaplan-Meier (median split), Cox (per 1 SD,
           adjusted for age, sex and stage I / II / III-IV), Schoenfeld test.
           Output: 05_survival/, Table 3, Table 9 (TCGA rows), Fig S2
SECTION 4  GSE81089 survival, Table 9 (GSE81089 rows), Fig 6, Fig 7, Fig 8.
           Input : 01_DESeq2/*.rds, files written by Sections 2 and 3
           Output: 05_survival/, 07_tables_for_manuscript/, 08_figures/

Section 1.11 (intermediate objects). 01_DESeq2/vsd_main.rds and
res_sym_main.rds hold the VST matrix (26,184 genes x 127 samples) with the
sample data (age, sex, stage, survival) and the DESeq2 table of Section 1.
They are provided and used as they are. To rebuild them from Section 1 set
REBUILD_INTERMEDIATE <- TRUE at the top (or delete the two files); this needs
metadata/clinical_GSE81089.csv.

Note on Section 1: the counts file has 218 samples; the script keeps the 127
listed in coldata.csv (108 LUAD + 19 normal).

Note on TCGA data: [include luad_SE.rds, or give the download command here.]

NOT INCLUDED
------------
 - Protein-protein interaction network (STRING, Cytoscape, cytoHubba, MCODE):
   settings are in the Methods, cytoHubba scores in Supplementary Table S6.
 - Web tools (TIMER, UALCAN, cBioPortal, Human Protein Atlas): settings and
   access dates are in the Methods.

Contact: [saba94ehsan@gmail.com]
