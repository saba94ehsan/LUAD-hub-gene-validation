###############################################################################
# LUAD hub-gene study - complete analysis in one script
# Integrative Transcriptomic and Network-Based Identification of Cell-Cycle
# Hub Genes in Lung Adenocarcinoma
#
#   SECTION 1  GSE81089: DESeq2, PCA, paired analysis, topGO, fgsea
#   SECTION 2  TCGA-LUAD: DESeq2 validation and comparison with GSE81089
#   SECTION 3  TCGA-LUAD: survival (Kaplan-Meier, Cox, Schoenfeld), Tables 3 and 9
#   SECTION 4  GSE81089: survival; Figures 6, 7 and 8 (forest plot)
#
# NOT included: STRING / Cytoscape (cytoHubba, MCODE) and the web tools
# (TIMER, UALCAN, cBioPortal, Human Protein Atlas); see the Methods.
#
# HOW TO RUN
#   1. Set PROJ_DIR below to the folder that contains this script.
#   2. Input files (relative to PROJ_DIR):
#        GSE81089_readcounts_featurecounts.tsv.gz   (GEO, series GSE81089)
#        metadata/coldata.csv                       (sample, condition, patient)
#        04_TCGA/luad_SE.rds                        (TCGA-LUAD, GDC, SummarizedExperiment)
#        01_DESeq2/vsd_main.rds, 01_DESeq2/res_sym_main.rds
#                                                   (provided; see SECTION 1.11)
#   3. Run the whole script:  Rscript LUAD_hub_gene_analysis.R
#      or open it in RStudio and run it section by section (Ctrl + Enter);
#      each section reads what the previous section wrote to disk.
#   4. Check 00_session/reference_check.csv: every row should say OK.
#      Comments starting with CHECK say what to look at after each part.
#
# Software used: R 4.6.1, Bioconductor 3.23 (see 00_session/sessionInfo*.txt).
###############################################################################

## ---- 0: settings ----------------------------------------------------------
PROJ_DIR      <- "~/LUAD_code"                 # <-- EDIT: folder with this script
COUNTS_FILE   <- "GSE81089_readcounts_featurecounts.tsv.gz"
COLDATA_FILE  <- "metadata/coldata.csv"
CLINICAL_FILE <- "metadata/clinical_GSE81089.csv"
TCGA_SE_FILE  <- "04_TCGA/luad_SE.rds"
REBUILD_INTERMEDIATE <- FALSE   # TRUE = rebuild 01_DESeq2/*.rds from SECTION 1 (see 1.11)
PADJ_CUT <- 0.05
LFC_CUT  <- 1

setwd(PROJ_DIR); getwd()
for (d in c("00_session", "01_DESeq2", "02_GO_GSEA", "03_PPI", "04_TCGA",
            "05_survival", "07_tables_for_manuscript", "08_figures"))
  dir.create(d, showWarnings = FALSE)
stopifnot("counts file not found in PROJ_DIR" = file.exists(COUNTS_FILE),
          "metadata/coldata.csv not found"     = file.exists(COLDATA_FILE))

set.seed(42)                       # fixed seed (used again before fgsea)
options(stringsAsFactors = FALSE)

## ---- 0.1: packages and versions ------------------------------------------
need <- c("survival", "survminer", "ggplot2", "patchwork")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)
suppressPackageStartupMessages({
  library(SummarizedExperiment); library(DESeq2); library(apeglm)
  library(topGO); library(fgsea); library(org.Hs.eg.db); library(AnnotationDbi)
  library(msigdbr); library(BiocParallel)
  library(ggplot2); library(survival); library(survminer); library(patchwork)
})
REQUIRED <- c(DESeq2 = "1.52.0", topGO = "2.64.0", fgsea = "1.38.0",
              org.Hs.eg.db = "3.23.1", msigdbr = "26.1.1")
if (as.character(getRversion()) != "4.6.1")
  warning("R is ", getRversion(), "; results were produced with 4.6.1")
for (p in names(REQUIRED)) {
  v <- as.character(packageVersion(p))
  if (v != REQUIRED[[p]])
    warning(p, " is ", v, "; results were produced with ", REQUIRED[[p]])
}
register(SerialParam())            # no parallel workers -> no seed drift
print(sapply(c("survival", "survminer", "ggplot2", "patchwork", "SummarizedExperiment"),
             function(p) as.character(packageVersion(p))))


## =============================================================================
## SECTION 1  GSE81089: DESeq2, PCA, paired analysis, topGO, fgsea
## =============================================================================
## ---- 1.1: counts and sample table -----------------------------------------------
raw <- read.delim(COUNTS_FILE, row.names = 1, check.names = FALSE)
tc  <- as.numeric(raw["TC%", ])                 # tumour cell percentage
names(tc) <- colnames(raw)
raw <- raw[rownames(raw) != "TC%", , drop = FALSE]   # drop non-gene row
cd  <- read.csv(COLDATA_FILE, row.names = "sample")
raw <- raw[, colnames(raw) %in% rownames(cd), drop = FALSE]   # keep the 127 samples listed in coldata.csv
tc  <- tc[colnames(raw)]
cd  <- cd[colnames(raw), , drop = FALSE]             # same order as counts
cd$condition <- factor(cd$condition, levels = c("Normal", "LUAD"))
cd$patient   <- factor(cd$patient)
cnt <- as.matrix(raw); storage.mode(cnt) <- "integer"
rownames(cnt) <- sub("\\..*$", "", rownames(cnt))    # drop Ensembl version

# keep genes with >= 10 reads in >= 10 samples
keep <- rowSums(cnt >= 10) >= 10
cnt  <- cnt[keep, ]

## ---- 1.2: DESeq2 (design ~ condition, Wald test) --------------------------------
dds <- DESeqDataSetFromMatrix(cnt, cd, design = ~ condition)
dds <- DESeq(dds, quiet = TRUE)
# unshrunken Wald results (no lfcShrink); BH-adjusted p
res <- results(dds, contrast = c("condition", "LUAD", "Normal"),
               alpha = PADJ_CUT, pAdjustMethod = "BH")
res <- as.data.frame(res); res$ensembl <- rownames(res)
res$Direction <- ifelse(!is.na(res$padj) & res$padj < PADJ_CUT &
                          abs(res$log2FoldChange) > LFC_CUT,
                        ifelse(res$log2FoldChange > 0, "Up", "Down"), "NS")

## ---- 1.3: PCA (VST, 500 most variable genes) ------------------------------------
vsd <- vst(dds, blind = TRUE)
pca <- prcomp(t(assay(vsd)[order(-matrixStats::rowVars(assay(vsd)))[1:500], ]))
pc_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
tum <- cd$condition == "LUAD" & !is.na(tc)
sp  <- cor.test(tc[tum], pca$x[tum, 1], method = "spearman", exact = FALSE)

## ---- 1.4: gene symbols (org.Hs.eg.db) -------------------------------------------
res$symbol <- mapIds(org.Hs.eg.db, keys = res$ensembl, keytype = "ENSEMBL",
                     column = "SYMBOL", multiVals = "first")
withsym <- res[!is.na(res$symbol) & !is.na(res$padj), ]
sig     <- withsym[withsym$Direction != "NS", ]
write.csv(res[order(res$padj), ], "01_DESeq2/DESeq2_all_genes.csv", row.names = FALSE)
write.csv(sig[order(sig$padj), ], "01_DESeq2/DEGs_with_symbol.csv", row.names = FALSE)

## ---- 1.5: paired analysis (patients with a matched normal) ----------------------
pt_n  <- table(cd$patient, cd$condition)
pairs <- rownames(pt_n)[pt_n[, "LUAD"] >= 1 & pt_n[, "Normal"] >= 1]
cdp   <- droplevels(cd[cd$patient %in% pairs, ])
dds_p <- DESeqDataSetFromMatrix(cnt[, rownames(cdp)], cdp, design = ~ patient + condition)
dds_p <- DESeq(dds_p, quiet = TRUE)
rp <- as.data.frame(results(dds_p, contrast = c("condition", "LUAD", "Normal"),
                            alpha = PADJ_CUT)); rp$ensembl <- rownames(rp)
rp <- rp[!is.na(rp$padj) & rp$padj < PADJ_CUT & abs(rp$log2FoldChange) > LFC_CUT, ]
rp <- rp[rp$ensembl %in% withsym$ensembl, ]            # symbol-annotated genes
shared <- intersect(rp$ensembl, sig$ensembl)
same_dir <- mean(sign(rp$log2FoldChange[match(shared, rp$ensembl)]) ==
                 sign(sig$log2FoldChange[match(shared, sig$ensembl)]))

## ---- 1.6: topGO (BP/MF/CC; weight01 + Fisher; nodeSize = 10) --------------------
univ <- unique(withsym$symbol)
deg  <- unique(sig$symbol)
go_all <- list()
for (ont in c("BP", "MF", "CC")) {
  gl  <- factor(as.integer(univ %in% deg)); names(gl) <- univ
  tgd <- new("topGOdata", ontology = ont, allGenes = gl, nodeSize = 10,
             annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")
  rs  <- runTest(tgd, algorithm = "weight01", statistic = "fisher")
  tb  <- GenTable(tgd, raw = rs, topNodes = length(score(rs)), numChar = 120)
  tb$raw_p <- score(rs)[tb$GO.ID]
  tb$BH_p  <- p.adjust(score(rs), "BH")[tb$GO.ID]
  tb$ontology <- ont
  go_all[[ont]] <- tb[order(tb$raw_p), ]
}
go_all <- do.call(rbind, go_all)
write.csv(go_all, "02_GO_GSEA/topGO_all_terms.csv", row.names = FALSE)

## ---- 1.7: fgsea (Hallmark, KEGG_LEGACY, Reactome; size 15-500) ------------------
stat <- sign(withsym$log2FoldChange) * -log10(withsym$pvalue)
names(stat) <- withsym$symbol
stat <- stat[is.finite(stat)]
stat <- stat[order(-abs(stat))]; stat <- stat[!duplicated(names(stat))]  # one value per symbol
stat <- sort(stat, decreasing = TRUE)

msig <- function(coll, sub = NULL)
  split(msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)$gene_symbol,
        msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)$gs_name)
sets <- list(Hallmark = msig("H"),
             KEGG     = msig("C2", "CP:KEGG_LEGACY"),
             Reactome = msig("C2", "CP:REACTOME"))
gsea <- list()
for (nm in names(sets)) {
  set.seed(42)                                   # same seed for every collection
  f <- fgsea(pathways = sets[[nm]], stats = stat, minSize = 15, maxSize = 500,
             BPPARAM = SerialParam())
  f$database <- nm
  gsea[[nm]] <- f[order(f$padj), ]
}
gsea <- do.call(rbind, gsea)
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, "", collapse = ";")
write.csv(gsea, "02_GO_GSEA/fgsea_all_sets.csv", row.names = FALSE)
gs <- gsea[!is.na(gsea$padj) & gsea$padj < 0.05, ]

## ---- 1.8: DEGs without a current symbol (Ensembl 73 / GRCh37) -------------------
# Needs internet. Output: 01_DESeq2/DEGs_without_current_symbol_GRCh37_biotype.csv
nosym <- res[is.na(res$symbol) & !is.na(res$padj) & res$Direction != "NS", ]
if (requireNamespace("biomaRt", quietly = TRUE)) tryCatch({
  mart <- biomaRt::useEnsembl("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",
                              host = "https://feb2014.archive.ensembl.org")
  bm <- biomaRt::getBM(c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
                       "ensembl_gene_id", nosym$ensembl, mart)
  bm2 <- merge(nosym, bm, by.x = "ensembl", by.y = "ensembl_gene_id")
  write.csv(bm2[order(bm2$padj), ],
            "01_DESeq2/DEGs_without_current_symbol_GRCh37_biotype.csv", row.names = FALSE)
  up_nosym <- bm2[bm2$gene_biotype == "protein_coding" & bm2$Direction == "Up" &
                    grepl("^[A-Z0-9-]+$", bm2$external_gene_name) &
                    !grepl("^(AC|AL|AP|RP)[0-9]", bm2$external_gene_name), ]
  writeLines(unique(up_nosym$external_gene_name),
             "03_PPI/up_DEG_symbols_from_GRCh37_check.txt")
}, error = function(e)
  message("1.9 skipped (Ensembl archive not reachable): ", conditionMessage(e)))

## ---- 1.9: PPI input (up-regulated DEGs, unique symbols) -------------------------
writeLines(sort(unique(sig$symbol[sig$Direction == "Up"])), "03_PPI/up_DEG_symbols.txt")

## ---- 1.10: check against the verified run ----------------------------------------
chk <- function(what, got, ref) data.frame(item = what, got = got, verified = ref,
                                           status = ifelse(isTRUE(all.equal(got, ref)), "OK", "DIFFERENT"))
all_deg <- res[!is.na(res$padj) & res$Direction != "NS", ]
ck <- rbind(
  chk("LUAD samples", sum(cd$condition == "LUAD"), 108),
  chk("normal samples", sum(cd$condition == "Normal"), 19),
  chk("genes after filter", nrow(cnt), 26184),
  chk("PC1 %", round(pc_var[1]), 17), chk("PC2 %", round(pc_var[2]), 11),
  chk("tested genes with symbol", nrow(withsym), 20201),
  chk("DEGs with symbol", nrow(sig), 4840),
  chk("up", sum(sig$Direction == "Up"), 3089),
  chk("down", sum(sig$Direction == "Down"), 1751),
  chk("unique up symbols", length(unique(sig$symbol[sig$Direction == "Up"])), 3087),
  chk("DEGs, all 26,184 genes", nrow(all_deg), 6395),
  chk("  up", sum(all_deg$Direction == "Up"), 4269),
  chk("  down", sum(all_deg$Direction == "Down"), 2126),
  chk("DEGs without symbol", nrow(nosym), 1555),
  chk("matched pairs", length(pairs), 9),
  chk("paired DEGs", nrow(rp), 4535),
  chk("  up", sum(rp$log2FoldChange > 0), 2026),
  chk("  down", sum(rp$log2FoldChange < 0), 2509),
  chk("paired DEGs also in main", length(shared), 3216),
  chk("same direction %", round(100 * same_dir, 1), 99.5),
  chk("ranked genes in GSEA", length(stat), 20182),
  chk("GSEA activated", sum(gs$NES > 0), 200),
  chk("GSEA suppressed", sum(gs$NES < 0), 144)
)
write.csv(ck, "00_session/reference_check.csv", row.names = FALSE)
print(ck, row.names = FALSE)
writeLines(c(capture.output(sessionInfo()),
             "", sprintf("Spearman tumour-cell %% vs PC1: rho=%.2f p=%.2f", sp$estimate, sp$p.value)),
           "00_session/sessionInfo.txt")

## ---- 1.11: intermediate objects used by Sections 2-4 -----------------------
# 01_DESeq2/vsd_main.rds   VST matrix (26,184 genes x 127 samples) with the sample
#                          data (condition, age, sex, stage, survival)
# 01_DESeq2/res_sym_main.rds   DESeq2 result table of Section 1 (genes with symbol)
# Both files are provided. They are rebuilt here only if they are missing or
# REBUILD_INTERMEDIATE <- TRUE. The rebuild needs metadata/clinical_GSE81089.csv
# (columns: counts_column, Patient, Histology, Stage, Smoking, Age, Gender,
# WHO_PS, OS_days, Dead; clinical data of the same 127 samples from GEO).
INTERMEDIATE <- c("01_DESeq2/vsd_main.rds", "01_DESeq2/res_sym_main.rds")
if (REBUILD_INTERMEDIATE || !all(file.exists(INTERMEDIATE))) {
  clin_geo <- read.csv(CLINICAL_FILE, check.names = FALSE)
  stopifnot("counts_column" %in% names(clin_geo))
  mm1 <- match(colnames(vsd), clin_geo$counts_column)
  stopifnot(!anyNA(mm1))                                  # every sample has a clinical row
  cdm <- clin_geo[mm1, setdiff(names(clin_geo), c("counts_column", "sample_id")), drop = FALSE]
  rownames(cdm) <- colnames(vsd)
  cdm$Condition <- factor(ifelse(cd$condition == "Normal", "Normal", "Tumor"),
                          levels = c("Normal", "Tumor"))
  vsd_main <- vsd
  colData(vsd_main) <- S4Vectors::DataFrame(cdm)
  stopifnot(sum(vsd_main$Condition == "Tumor") == 108, sum(vsd_main$Condition == "Normal") == 19)
  res_sym_main <- res[!is.na(res$symbol), ]
  res_sym_main$Ensembl     <- res_sym_main$ensembl
  res_sym_main$log2FC      <- res_sym_main$log2FoldChange
  res_sym_main$SYMBOL      <- res_sym_main$symbol
  res_sym_main$Has_SYMBOL  <- TRUE
  res_sym_main$Significant <- res_sym_main$Direction != "NS"
  res_sym_main <- res_sym_main[order(res_sym_main$padj), ]
  saveRDS(vsd_main,     INTERMEDIATE[1])
  saveRDS(res_sym_main, INTERMEDIATE[2])
  cat("rebuilt:", dim(vsd_main), "|", nrow(res_sym_main), "genes\n")
}


## =============================================================================
## SECTION 2  TCGA-LUAD: DESeq2 validation and comparison with GSE81089
## =============================================================================
## ---- 2.0: check memory and input file --------------------------------------
getwd()
if (.Platform$OS.type == "unix") system("free -h")
stopifnot("TCGA file not found: 04_TCGA/luad_SE.rds" = file.exists(TCGA_SE_FILE))
list.files("04_TCGA")
# CHECK: getwd() should show .../Desktop/LUAD_project.
#        Total memory should be 8 GB or more.
#        luad_SE.rds should be listed in 04_TCGA.


# When was this file created (for the download date in the Methods):
file.info("04_TCGA/luad_SE.rds")[, c("size", "mtime")]
# NOTE: mtime is the time the file was copied to this computer, not the
#       download date. The download date is the day GDCdownload was run;
#       take it from your notes or from the GDCdata folder of the download.
#


## ---- 2.1: open and inspect the TCGA data -----------------------------------
se <- readRDS("04_TCGA/luad_SE.rds")
dim(se)
assayNames(se)
table(se$shortLetterCode)
table(se$sample_type)
# CHECK: assayNames should contain "unstranded" (these are the raw counts).
#        shortLetterCode should show TP (tumour) and NT (normal).
#        If TR (recurrent) also appears, it is removed below.


## ---- 2.2: select samples (1 tumour per patient) ----------------------------
bc   <- colnames(se)
code <- se$shortLetterCode
pat  <- substr(bc, 1, 12)

# Primary Tumor (TP) and Solid Tissue Normal (NT) only
use <- code %in% c("TP", "NT")
se2 <- se[, use]
bc2 <- colnames(se2); code2 <- se2$shortLetterCode; pat2 <- substr(bc2, 1, 12)

# if a patient has more than one tumour, keep only one
# (the first barcode in alphabetical order - a fixed rule, no choice involved)
is_tumor <- code2 == "TP"
tum_idx  <- which(is_tumor)
tum_idx  <- tum_idx[order(bc2[tum_idx])]
keep_tum <- tum_idx[!duplicated(pat2[tum_idx])]
keep_nor <- which(!is_tumor)
sel <- sort(c(keep_tum, keep_nor))

cat("Primary tumor samples downloaded :", sum(is_tumor), "\n")
cat("Unique tumor patients            :", length(unique(pat2[is_tumor])), "\n")
cat("Tumor kept (1 per patient)       :", length(keep_tum), "\n")
cat("Normal samples                   :", length(keep_nor), "\n")
cat("Patients with BOTH tumor+normal  :",
    length(intersect(pat2[keep_tum], pat2[keep_nor])), "\n")
cat("Duplicate patients among normals :", sum(duplicated(pat2[keep_nor])), "\n")

se3 <- se2[, sel]

cnt <- assay(se3, "unstranded")
cnt <- cnt[grepl("^ENSG", rownames(cnt)), ]
cnt <- cnt[!grepl("_PAR_Y", rownames(cnt)), ]
rownames(cnt) <- sub("\\..*$", "", rownames(cnt))
stopifnot(!anyDuplicated(rownames(cnt)))
storage.mode(cnt) <- "integer"

meta <- data.frame(
  barcode   = colnames(se3),
  Patient   = substr(colnames(se3), 1, 12),
  Condition = factor(ifelse(se3$shortLetterCode == "NT", "Normal", "Tumor"),
                     levels = c("Normal", "Tumor")),
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$barcode
stopifnot(all(colnames(cnt) == rownames(meta)))
meta$TCGA_paired <- meta$Patient %in% meta$Patient[meta$Condition == "Normal"] &
                    meta$Patient %in% meta$Patient[meta$Condition == "Tumor"]

table(meta$Condition)
dim(cnt)
head(rownames(cnt))
write.csv(meta, "04_TCGA/TCGA_sample_table.csv", row.names = FALSE)
# CHECK: the table shows the Normal and Tumor counts (whatever appears is the result).
#        head(rownames) should look like ENSG00000000003 (without ".15").
#        "Duplicate patients among normals" should be 0.
#        Compare these counts with Table 1 / Table 3 of the manuscript.


## ---- 2.3: filter + DESeq2 (takes a few minutes) ----------------------------
# Same filter as for GSE81089: >=10 reads in at least 10 samples
keep <- rowSums(cnt >= 10) >= 10
cat("Genes before filter:", nrow(cnt), " after:", sum(keep), "\n")

dds <- DESeqDataSetFromMatrix(countData = cnt[keep, ], colData = meta,
                              design = ~ Condition)
dds <- DESeq(dds)
resultsNames(dds)
saveRDS(dds, "04_TCGA/dds_TCGA_main.rds")
# CHECK: resultsNames should contain "Condition_Tumor_vs_Normal",
#        i.e. "Intercept" and this one only. No error.


## ---- 2.4: results (RAW and annotated, separately) --------------------------
ru <- results(dds, name = "Condition_Tumor_vs_Normal", alpha = 0.05)  # unshrunk
rs <- lfcShrink(dds, coef = "Condition_Tumor_vs_Normal", type = "apeglm")

# --- RAW sheet (DESeq2 output only) ---
raw <- as.data.frame(ru)
raw <- cbind(Ensembl = rownames(raw), raw)
write.csv(raw, "04_TCGA/TCGA_DESeq2_RAW_results_ALL_genes.csv", row.names = FALSE)

# --- annotated sheet ---
res <- data.frame(
  Ensembl      = rownames(ru),
  baseMean     = ru$baseMean,
  log2FC       = ru$log2FoldChange,          # unshrunk (used for the DEGs)
  lfcSE        = ru$lfcSE,
  stat         = ru$stat,
  pvalue       = ru$pvalue,
  padj         = ru$padj,
  log2FC_shrunk = rs$log2FoldChange,
  stringsAsFactors = FALSE
)
for (cc in c("SYMBOL", "GENENAME", "ENTREZID", "GENETYPE")) {
  res[[cc]] <- mapIds(org.Hs.eg.db, keys = res$Ensembl, column = cc,
                      keytype = "ENSEMBL", multiVals = "first")
}
nc   <- counts(dds, normalized = TRUE)
cond <- colData(dds)$Condition
res$mean_norm_Tumor  <- rowMeans(nc[res$Ensembl, cond == "Tumor"])
res$mean_norm_Normal <- rowMeans(nc[res$Ensembl, cond == "Normal"])
res$Tumor_frac  <- rowMeans(nc[res$Ensembl, cond == "Tumor"]  >= 10)
res$Normal_frac <- rowMeans(nc[res$Ensembl, cond == "Normal"] >= 10)
res$Has_SYMBOL  <- !is.na(res$SYMBOL)
res$Significant <- !is.na(res$padj) & res$padj < 0.05 & abs(res$log2FC) > 1
res$Direction   <- ifelse(res$Significant & res$log2FC > 0, "Up",
                   ifelse(res$Significant & res$log2FC < 0, "Down", "NS"))
res$Reliable <- (res$Direction == "Up"   & res$Tumor_frac  >= 0.5) |
                (res$Direction == "Down" & res$Normal_frac >= 0.5)
res <- res[order(res$padj), ]

write.csv(res, "04_TCGA/TCGA_DESeq2_annotated_ALL_genes.csv", row.names = FALSE)
res_sym_tcga <- res[res$Has_SYMBOL, ]
saveRDS(res_sym_tcga, "04_TCGA/res_sym_TCGA.rds")
write.csv(res_sym_tcga[res_sym_tcga$Direction == "Up", ],
          "04_TCGA/TCGA_DEG_UP.csv", row.names = FALSE)
write.csv(res_sym_tcga[res_sym_tcga$Direction == "Down", ],
          "04_TCGA/TCGA_DEG_DOWN.csv", row.names = FALSE)

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(RAW_DESeq2 = raw, Annotated = res,
                           DEG_UP = res[res$Direction == "Up", ],
                           DEG_DOWN = res[res$Direction == "Down", ]),
                      "04_TCGA/TCGA_DESeq2_all_sheets.xlsx")
}

cat("Genes tested        :", sum(!is.na(res$padj)), "\n")
cat("with symbol         :", sum(res$Has_SYMBOL), "\n")
print(table(res$Direction))
print(table(res$Direction[res$Has_SYMBOL]))
# CHECK: the files should exist in 04_TCGA (list.files("04_TCGA")).
#        The Up/Down/NS counts are the TCGA result.


## ---- 2.5: TCGA expression of the 4 hub genes -------------------------------
hub  <- c("CDK1", "CCNB1", "CCNA2", "BRCA1")
cmp4 <- c(hub, "GAPDH", "ALB")     # GAPDH, ALB for comparison only

# GSE81089 results
geo <- readRDS("01_DESeq2/res_sym_main.rds")
geo_ens <- if ("Ensembl" %in% names(geo)) geo$Ensembl else rownames(geo)

pick <- function(df, g) {
  d <- df[which(df$SYMBOL == g), ]
  d[order(d$padj), ][1, ]
}
tab8 <- do.call(rbind, lapply(cmp4, function(g) {
  a <- pick(geo, g); b <- pick(res_sym_tcga, g)
  data.frame(Gene = g,
             GEO_log2FC  = round(a$log2FC, 2),  GEO_padj  = signif(a$padj, 2),
             TCGA_log2FC = round(b$log2FC, 2),  TCGA_padj = signif(b$padj, 2),
             Same_direction = sign(a$log2FC) == sign(b$log2FC),
             TCGA_Tumor_frac = round(b$Tumor_frac, 2))
}))
print(tab8)
write.csv(tab8, "07_tables_for_manuscript/Table8_hubgenes_GEO_vs_TCGA.csv",
          row.names = FALSE)
# CHECK: 6 rows (4 hub genes + GAPDH + ALB). This table is the source of
#        Table 8 of the manuscript.


## ---- 2.6: GSE81089 vs TCGA - genome-wide comparison ------------------------
g1 <- geo[!duplicated(geo$SYMBOL), c("SYMBOL", "log2FC", "padj")]
g2 <- res_sym_tcga[!duplicated(res_sym_tcga$SYMBOL), c("SYMBOL", "log2FC", "padj")]
mm <- merge(g1, g2, by = "SYMBOL", suffixes = c("_GEO", "_TCGA"))
cat("Common genes with symbol:", nrow(mm), "\n")
cat("Spearman log2FC (GEO vs TCGA):",
    round(cor(mm$log2FC_GEO, mm$log2FC_TCGA, method = "spearman",
              use = "complete.obs"), 3), "\n")

geo_up   <- geo$SYMBOL[geo$Direction == "Up"]
geo_down <- geo$SYMBOL[geo$Direction == "Down"]
t_up     <- res_sym_tcga$SYMBOL[res_sym_tcga$Direction == "Up"]
t_down   <- res_sym_tcga$SYMBOL[res_sym_tcga$Direction == "Down"]
cat("GEO up in TCGA up    :", length(intersect(geo_up, t_up)),   "/", length(geo_up),   "\n")
cat("GEO down in TCGA down:", length(intersect(geo_down, t_down)), "/", length(geo_down), "\n")
cat("Opposite direction   :",
    length(intersect(geo_up, t_down)) + length(intersect(geo_down, t_up)), "\n")

p <- ggplot(mm, aes(log2FC_GEO, log2FC_TCGA)) +
  geom_point(alpha = 0.2, size = 0.5) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  coord_cartesian(xlim = c(-8, 12), ylim = c(-8, 12)) + theme_bw() +
  labs(x = "log2 fold change, GSE81089 (LUAD vs normal)",
       y = "log2 fold change, TCGA-LUAD (tumor vs normal)")
ggsave("08_figures/FigS_GEO_vs_TCGA_log2FC.png", p, width = 6, height = 6, dpi = 300)
# CHECK: the counts and the correlation are printed. This is descriptive only;
#        report what you obtain (there is no target value).
#        The plot should be written to 08_figures.


## ---- 2.7: prepare data for survival analysis -------------------------------
vsd <- vst(dds, blind = TRUE)
saveRDS(vsd, "04_TCGA/vsd_TCGA.rds")

clin <- as.data.frame(colData(se3))
clin <- clin[, !vapply(clin, is.list, logical(1))]     # drop list columns
clin$barcode <- rownames(clin)
write.csv(clin, "04_TCGA/TCGA_clinical_selected_samples.csv", row.names = FALSE)

names(clin)[grepl("days|vital|age|gender|stage|ajcc|race|smok",
                  names(clin), ignore.case = TRUE)]
# CHECK: the column list should include: vital_status, days_to_death,
#        days_to_last_follow_up, gender, age_at_index (or age_at_diagnosis),
#        ajcc_pathologic_stage. The survival script (Section 3) uses these names.
#



## =============================================================================
## SECTION 3  TCGA-LUAD: survival (Kaplan-Meier, Cox, Schoenfeld), Tables 3 and 9
## =============================================================================
## ---- 3.0: package versions -------------------------------------------------


vers <- sapply(c("survival", "survminer", "ggplot2", "SummarizedExperiment"),
               function(p) as.character(packageVersion(p)))
print(vers)


## ---- 3.1: open the data, extract hub-gene expression -----------------------
vsd  <- readRDS("04_TCGA/vsd_TCGA.rds")
rs   <- readRDS("04_TCGA/res_sym_TCGA.rds")
clin <- read.csv("04_TCGA/TCGA_clinical_selected_samples.csv",
                 stringsAsFactors = FALSE, check.names = FALSE)
dim(vsd)
dim(clin)
table(clin$shortLetterCode)
needcols <- c("barcode", "vital_status", "days_to_death",
              "days_to_last_follow_up", "age_at_index", "ajcc_pathologic_stage")
setdiff(needcols, names(clin))
intersect(c("sex_at_birth", "gender"), names(clin))

hub <- c("CDK1", "CCNB1", "CCNA2", "BRCA1")
ens <- vapply(hub, function(g) {
  d <- rs[which(rs$SYMBOL == g), ]
  d$Ensembl[order(d$padj)][1]          # same as pick() in Section 2
}, character(1))
stopifnot(all(ens %in% rownames(vsd)))
expr <- assay(vsd)[ens, , drop = FALSE]
rownames(expr) <- hub
ens
# CHECK: shortLetterCode should show TP and NT (517 and 59 in Section 2).
#        setdiff(...) should return "character(0)" (all columns present).
#        ens should show 4 ENSG IDs.


## ---- 3.2: build the survival table + Table 3 -------------------------------
parse_stage <- function(x) {
  st <- toupper(sub("^Stage\\s+", "", x))
  s4 <- ifelse(grepl("^IV", st), "IV",
        ifelse(grepl("^III", st), "III",
        ifelse(grepl("^II", st), "II",
        ifelse(grepl("^I", st), "I", NA))))
  s4 <- factor(s4, levels = c("I", "II", "III", "IV"))
  s3 <- factor(ifelse(s4 %in% c("III", "IV"), "III-IV", as.character(s4)),
               levels = c("I", "II", "III-IV"))
  data.frame(stage4 = s4, stage = s3)
}

tum <- clin[clin$shortLetterCode == "TP", ]
tum$patient12 <- substr(tum$barcode, 1, 12)
stopifnot(!anyDuplicated(tum$patient12))
stopifnot(all(tum$barcode %in% colnames(vsd)))

tum$event    <- as.integer(tum$vital_status == "Dead")
tum$time_day <- ifelse(tum$event == 1, tum$days_to_death, tum$days_to_last_follow_up)
tum$time_y   <- tum$time_day / 365.25
sexcol       <- intersect(c("sex_at_birth", "gender"), names(tum))[1]
tum$sex      <- factor(tolower(tum[[sexcol]]), levels = c("female", "male"))
tum$age      <- as.numeric(tum$age_at_index)
ps           <- parse_stage(tum$ajcc_pathologic_stage)
tum$stage4   <- ps$stage4
tum$stage    <- ps$stage

cat("Tumor patients (RNA-seq)  :", nrow(tum), "\n")
cat("Survival time missing     :", sum(is.na(tum$time_day)), "\n")
cat("Survival time <= 0        :", sum(tum$time_day <= 0, na.rm = TRUE), "\n")
sv <- tum[!is.na(tum$time_day) & tum$time_day > 0, ]
cat("Survival cohort           :", nrow(sv), " | deaths:", sum(sv$event), "\n")
cat("Complete age+sex+stage    :",
    sum(complete.cases(sv[, c("age", "sex", "stage")])), " | deaths:",
    sum(sv$event[complete.cases(sv[, c("age", "sex", "stage")])]), "\n")

for (g in hub) sv[[g]] <- as.numeric(expr[g, sv$barcode])

# median follow-up: reverse Kaplan-Meier
rk <- survfit(Surv(time_y, 1 - event) ~ 1, data = sv)
fu <- summary(rk)$table[c("median", "0.95LCL", "0.95UCL")]
print(round(fu, 2))

n_nt <- sum(clin$shortLetterCode == "NT")
tab3 <- data.frame(
  Characteristic = c("Tumor / normal samples with RNA-seq",
                     "Patients in survival analysis / deaths",
                     "Sex, female / male",
                     "Median age, years (n missing)",
                     "Stage I / II / III / IV (n missing)",
                     "Median follow-up, years (reverse KM, 95% CI)",
                     "Complete-case for adjusted model: patients / deaths"),
  Value = c(
    paste(nrow(tum), "/", n_nt),
    paste(nrow(sv), "/", sum(sv$event)),
    paste(sum(sv$sex == "female", na.rm = TRUE), "/", sum(sv$sex == "male", na.rm = TRUE)),
    sprintf("%.0f (%d)", median(sv$age, na.rm = TRUE), sum(is.na(sv$age))),
    sprintf("%d / %d / %d / %d (%d)", sum(sv$stage4 == "I", na.rm = TRUE),
            sum(sv$stage4 == "II", na.rm = TRUE), sum(sv$stage4 == "III", na.rm = TRUE),
            sum(sv$stage4 == "IV", na.rm = TRUE), sum(is.na(sv$stage4))),
    sprintf("%.1f (%.1f-%.1f)", fu[1], fu[2], fu[3]),
    paste(sum(complete.cases(sv[, c("age", "sex", "stage")])), "/",
          sum(sv$event[complete.cases(sv[, c("age", "sex", "stage")])]))
  ))
print(tab3, right = FALSE)
write.csv(tab3, "07_tables_for_manuscript/Table3_TCGA_clinical.csv", row.names = FALSE)
write.csv(sv[, c("barcode", "patient12", "time_day", "time_y", "event", "age", "sex",
                 "stage4", "stage", hub)],
          "05_survival/TCGA_survival_input.csv", row.names = FALSE)
# CHECK: compare these counts with Table 3 of the manuscript.
#        NOTE: stage may be missing for many patients
#        (ajcc_pathologic_stage is incomplete in colData). Take the stage
#        line of Table 3 from this output.


## ---- 3.3: Kaplan-Meier + log-rank (Supplementary Fig S2) -------------------
km_list <- list()
lr <- data.frame()
for (g in hub) {
  d <- sv
  d$grp <- factor(ifelse(d[[g]] > median(d[[g]]), "High", "Low"),
                  levels = c("Low", "High"))
  sdf <- survdiff(Surv(time_y, event) ~ grp, data = d)
  p   <- pchisq(sdf$chisq, df = 1, lower.tail = FALSE)
  lr  <- rbind(lr, data.frame(Gene = g, n_low = sum(d$grp == "Low"),
                              n_high = sum(d$grp == "High"),
                              events_low = sdf$obs[1], events_high = sdf$obs[2],
                              logrank_p = p))
  fit <- survfit(Surv(time_y, event) ~ grp, data = d)
  km_list[[g]] <- ggsurvplot(
    fit, data = d, risk.table = TRUE, conf.int = FALSE,
    pval = paste0("Log-rank p = ", signif(p, 2)), pval.size = 4,
    palette = c("#0072B2", "#D55E00"), legend.labs = c("Low", "High"),
    legend.title = g, xlab = "Time (years)", ylab = "Overall survival probability",
    title = paste0(g, " (n = ", nrow(d), ")"), break.time.by = 2,
    risk.table.height = 0.25, ggtheme = theme_bw())
}
print(lr)
write.csv(lr, "05_survival/TCGA_logrank.csv", row.names = FALSE)
pp <- arrange_ggsurvplots(unname(km_list), print = FALSE, ncol = 2, nrow = 2)
ggsave("08_figures/FigS2_TCGA_KM_hubgenes.png", pp, width = 12, height = 11, dpi = 300)
ggsave("08_figures/FigS2_TCGA_KM_hubgenes.pdf", pp, width = 12, height = 11)
# CHECK: 4 rows (CDK1, CCNB1, CCNA2, BRCA1); n_low + n_high = survival cohort.
#        The plot should be in 08_figures. Open the figure and check
#        that curves, numbers at risk and the p-value are shown.


## ---- 3.4: Cox models (per 1 SD) + Schoenfeld -------------------------------
# SD: computed over all patients of the survival cohort; the same z is used in the adjusted model.
for (g in hub) sv[[paste0(g, "_z")]] <- as.numeric(scale(sv[[g]]))

run_cox <- function(d, g, adjusted) {
  d$z <- d[[paste0(g, "_z")]]
  if (adjusted) {
    d <- d[complete.cases(d[, c("time_y", "event", "z", "age", "sex", "stage")]), ]
    m <- coxph(Surv(time_y, event) ~ z + age + sex + stage, data = d)
  } else {
    m <- coxph(Surv(time_y, event) ~ z, data = d)
  }
  co <- summary(m)$coefficients["z", ]
  ci <- unname(exp(confint(m)["z", ]))
  zp <- cox.zph(m)$table
  glob <- if ("GLOBAL" %in% rownames(zp)) zp["GLOBAL", "p"] else zp["z", "p"]
  data.frame(Gene = g, Model = ifelse(adjusted, "Adjusted", "Univariable"),
             N = m$n, Events = m$nevent,
             HR = unname(exp(co["coef"])), LCL = ci[1], UCL = ci[2],
             p = unname(co["Pr(>|z|)"]),
             PH_p_gene = unname(zp["z", "p"]), PH_p_global = unname(glob),
             row.names = NULL)
}

cox_all <- do.call(rbind, lapply(hub, function(g)
  rbind(run_cox(sv, g, FALSE), run_cox(sv, g, TRUE))))
print(cox_all, digits = 3)
write.csv(cox_all, "05_survival/TCGA_cox_results.csv", row.names = FALSE)
# CHECK: 8 rows (4 gene x 2 model). Univariable N = survival cohort,
#        Adjusted N = complete-case count (should match Part 2).


## ---- 3.5: TCGA rows of Table 9 ---------------------------------------------
fmt_hr <- function(h, l, u) sprintf("%.2f (%.2f–%.2f)", h, l, u)
fmt_p  <- function(p) ifelse(p < 0.001, "<0.001",
                             formatC(p, digits = 2, format = "g", flag = "#"))
uni <- cox_all[cox_all$Model == "Univariable", ]
adj <- cox_all[cox_all$Model == "Adjusted", ]
stopifnot(identical(uni$Gene, hub), identical(adj$Gene, hub), identical(lr$Gene, hub))
t9 <- data.frame(
  Gene = hub, Cohort = "TCGA-LUAD",
  `Patients / events` = paste(uni$N, "/", uni$Events),
  `Log-rank p (median split)` = formatC(lr$logrank_p, digits = 2, format = "g", flag = "#"),
  `Univariable HR (95% CI)` = fmt_hr(uni$HR, uni$LCL, uni$UCL),
  `Univariable p` = fmt_p(uni$p),
  `Adjusted HR* (95% CI)` = fmt_hr(adj$HR, adj$LCL, adj$UCL),
  `Adjusted p` = fmt_p(adj$p),
  `Adjusted N / events` = paste(adj$N, "/", adj$Events),
  `PH assumption (Schoenfeld p, gene term, adjusted model)` = formatC(adj$PH_p_gene, digits = 2, format = "f"),
  check.names = FALSE)
print(t9, right = FALSE)
write.csv(t9, "07_tables_for_manuscript/Table9_TCGA_rows.csv", row.names = FALSE)
# CHECK: 4 rows. These are the TCGA rows of Table 9; the Results paragraph
#        uses these numbers. The GSE81089 rows are produced in Section 4.
#        PH column: use the same model (univariable or adjusted) as in the
#        GSE table. Both PH p-values are saved in cox_all.


## ---- 3.6: (OPTIONAL) sensitivity analysis with stage filled from GDC -------
# Stage is missing for many patients in colData, so the complete-case adjusted
# model is smaller. This part fills stage from the GDC clinical
# table and re-runs the adjusted model.
# This is a SENSITIVITY analysis only; it does not replace the main result.
# Requires internet and TCGAbiolinks. Default FALSE.
FILL_STAGE_FROM_GDC <- FALSE
if (FILL_STAGE_FROM_GDC) {
  gc <- TCGAbiolinks::GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
  print(intersect(c("submitter_id", "ajcc_pathologic_stage"), names(gc)))
  stg <- gc$ajcc_pathologic_stage[match(sv$patient12, gc$submitter_id)]
  pf  <- parse_stage(stg)
  sv2 <- sv
  fill <- is.na(sv2$stage4)
  sv2$stage4[fill] <- pf$stage4[fill]
  sv2$stage[fill]  <- pf$stage[fill]
  cat("Stage missing before:", sum(fill), " after:", sum(is.na(sv2$stage)), "\n")
  cox_fill <- do.call(rbind, lapply(hub, function(g) run_cox(sv2, g, TRUE)))
  print(cox_fill, digits = 3)
  write.csv(cox_fill, "05_survival/TCGA_cox_adjusted_stage_filled_GDC.csv", row.names = FALSE)
}



## =============================================================================
## SECTION 4  GSE81089: survival; Figures 6, 7 and 8
## =============================================================================
## ---- 4.1: GSE81089 survival table ------------------------------------------
v   <- readRDS("01_DESeq2/vsd_main.rds")
geo <- readRDS("01_DESeq2/res_sym_main.rds")
cd  <- as.data.frame(colData(v))
table(cd$Condition)
hub <- c("CDK1", "CCNB1", "CCNA2", "BRCA1")
cmp <- hub                                         # Fig 6 shows the 4 hub genes only

geo_ens <- if ("Ensembl" %in% names(geo)) geo$Ensembl else rownames(geo)
pick_ens <- function(df, ens_ids, g) {
  i <- which(df$SYMBOL == g)
  ens_ids[i][order(df$padj[i])][1]
}
ens_geo <- vapply(cmp, function(g) pick_ens(geo, geo_ens, g), character(1))
stopifnot(all(ens_geo %in% rownames(v)))
ens_geo

parse_stage <- function(x) {
  st <- toupper(trimws(as.character(x)))
  s4 <- ifelse(grepl("^IV", st), "IV", ifelse(grepl("^III", st), "III",
        ifelse(grepl("^II", st), "II", ifelse(grepl("^I", st), "I", NA))))
  s4 <- factor(s4, levels = c("I", "II", "III", "IV"))
  s3 <- factor(ifelse(s4 %in% c("III", "IV"), "III-IV", as.character(s4)),
               levels = c("I", "II", "III-IV"))
  data.frame(stage4 = s4, stage = s3)
}

tum <- cd[cd$Condition == "Tumor", ]
tum$barcode <- rownames(tum)
tum$event   <- as.integer(tum$Dead)
tum$time_y  <- as.numeric(tum$OS_days) / 365.25
tum$age     <- as.numeric(tum$Age)
tum$sex     <- factor(tolower(tum$Gender), levels = c("female", "male"))
ps <- parse_stage(tum$Stage); tum$stage4 <- ps$stage4; tum$stage <- ps$stage

cat("Tumor patients            :", nrow(tum), " | deaths:", sum(tum$event, na.rm = TRUE), "\n")
cat("Survival time <= 0 / NA   :", sum(tum$time_y <= 0, na.rm = TRUE), "/",
    sum(is.na(tum$time_y)), "\n")
sv <- tum[!is.na(tum$time_y) & tum$time_y > 0, ]
cat("Survival cohort           :", nrow(sv), " | deaths:", sum(sv$event), "\n")
cat("Complete age+sex+stage    :",
    sum(complete.cases(sv[, c("age", "sex", "stage")])), "\n")
for (g in hub) sv[[g]] <- as.numeric(assay(v)[ens_geo[g], sv$barcode])
for (g in hub) sv[[paste0(g, "_z")]] <- as.numeric(scale(sv[[g]]))

rk <- survfit(Surv(time_y, 1 - event) ~ 1, data = sv)
fu <- summary(rk)$table[c("median", "0.95LCL", "0.95UCL")]
cat("Median follow-up (reverse KM), years:", round(fu, 2), "\n")
table(sv$stage4, useNA = "ifany")
# CHECK: the manuscript reports 108 tumours and 49 deaths; after removing 1 patient (0 days)
#        this gives 107 and 48. Report what you obtain.
#        Compare the stage counts with Table 2 (I = IA+IB, II = IIA+IIB).


## ---- 4.2: Cox + log-rank (GSE81089) ----------------------------------------
run_cox <- function(d, g, adjusted) {
  d$z <- d[[paste0(g, "_z")]]
  if (adjusted) {
    d <- d[complete.cases(d[, c("time_y", "event", "z", "age", "sex", "stage")]), ]
    m <- coxph(Surv(time_y, event) ~ z + age + sex + stage, data = d)
  } else {
    m <- coxph(Surv(time_y, event) ~ z, data = d)
  }
  co <- summary(m)$coefficients["z", ]
  ci <- unname(exp(confint(m)["z", ]))
  zp <- cox.zph(m)$table
  glob <- if ("GLOBAL" %in% rownames(zp)) zp["GLOBAL", "p"] else zp["z", "p"]
  data.frame(Gene = g, Model = ifelse(adjusted, "Adjusted", "Univariable"),
             N = m$n, Events = m$nevent,
             HR = unname(exp(co["coef"])), LCL = ci[1], UCL = ci[2],
             p = unname(co["Pr(>|z|)"]),
             PH_p_gene = unname(zp["z", "p"]), PH_p_global = unname(glob),
             row.names = NULL)
}
cox_geo <- do.call(rbind, lapply(hub, function(g)
  rbind(run_cox(sv, g, FALSE), run_cox(sv, g, TRUE))))
print(cox_geo, digits = 3)
write.csv(cox_geo, "05_survival/GSE81089_cox_results.csv", row.names = FALSE)

lr <- data.frame()
km_list <- list()
for (g in hub) {
  d <- sv
  d$grp <- factor(ifelse(d[[g]] > median(d[[g]]), "High", "Low"), levels = c("Low", "High"))
  sdf <- survdiff(Surv(time_y, event) ~ grp, data = d)
  p <- pchisq(sdf$chisq, df = 1, lower.tail = FALSE)
  lr <- rbind(lr, data.frame(Gene = g, n_low = sum(d$grp == "Low"),
                             n_high = sum(d$grp == "High"), logrank_p = p))
  fit <- survfit(Surv(time_y, event) ~ grp, data = d)
  km_list[[g]] <- ggsurvplot(
    fit, data = d, risk.table = TRUE, conf.int = FALSE,
    pval = paste0("Log-rank p = ", signif(p, 2)), pval.size = 4,
    palette = c("#0072B2", "#D55E00"), legend.labs = c("Low", "High"),
    legend.title = g, xlab = "Time (years)", ylab = "Overall survival probability",
    title = paste0(g, " (n = ", nrow(d), ")"), break.time.by = 2,
    risk.table.height = 0.25, ggtheme = theme_bw())
}
print(lr)
write.csv(lr, "05_survival/GSE81089_logrank.csv", row.names = FALSE)
# CHECK: cox_geo has 8 rows. Univariable N = survival cohort (Part 1).
#        Compare these numbers with the GSE81089 rows of Table 9.
#        Small differences (3rd decimal) can come from software versions;
#        report what you obtain.


## ---- 4.3: Table 9 (GSE81089 rows) + Fig 7 ----------------------------------
fmt_hr <- function(h, l, u) sprintf("%.2f (%.2f–%.2f)", h, l, u)
fmt_p  <- function(p) ifelse(p < 0.001, "<0.001",
                             formatC(p, digits = 2, format = "g", flag = "#"))
uni <- cox_geo[cox_geo$Model == "Univariable", ]
adj <- cox_geo[cox_geo$Model == "Adjusted", ]
stopifnot(identical(uni$Gene, hub), identical(adj$Gene, hub), identical(lr$Gene, hub))
t9g <- data.frame(
  Gene = hub, Cohort = "GSE81089",
  `Patients / events` = paste(uni$N, "/", uni$Events),
  `Log-rank p (median split)` = formatC(lr$logrank_p, digits = 2, format = "g", flag = "#"),
  `Univariable HR (95% CI)` = fmt_hr(uni$HR, uni$LCL, uni$UCL),
  `Univariable p` = fmt_p(uni$p),
  `Adjusted HR* (95% CI)` = fmt_hr(adj$HR, adj$LCL, adj$UCL),
  `Adjusted p` = fmt_p(adj$p),
  `Adjusted N / events` = paste(adj$N, "/", adj$Events),
  `PH assumption (Schoenfeld p, gene term, adjusted model)` = formatC(adj$PH_p_gene, digits = 2, format = "f"),
  check.names = FALSE)
print(t9g, right = FALSE)
write.csv(t9g, "07_tables_for_manuscript/Table9_GSE81089_rows.csv", row.names = FALSE)

pp <- arrange_ggsurvplots(unname(km_list), print = FALSE, ncol = 2, nrow = 2)
ggsave("08_figures/Fig7_GSE81089_KM.png", pp, width = 12, height = 11, dpi = 300)
ggsave("08_figures/Fig7_GSE81089_KM.pdf", pp, width = 12, height = 11)
# CHECK: the Table 9 csv and Fig 7 should be written. For the PH column use the same model
#        as in the TCGA table (Section 3 also uses the adjusted model).


## ---- 4.4: Fig 8, forest plot (GSE81089 + TCGA) -----------------------------
stopifnot(file.exists("05_survival/TCGA_cox_results.csv"))
cox_tcga <- read.csv("05_survival/TCGA_cox_results.csv")
fo <- rbind(cbind(Cohort = "GSE81089", cox_geo[, c("Gene","Model","N","Events","HR","LCL","UCL","p")]),
            cbind(Cohort = "TCGA-LUAD", cox_tcga[, c("Gene","Model","N","Events","HR","LCL","UCL","p")]))
fo$Gene   <- factor(fo$Gene, levels = rev(hub))
fo$Cohort <- factor(fo$Cohort, levels = c("GSE81089", "TCGA-LUAD"))
fo$Model  <- factor(fo$Model, levels = c("Univariable", "Adjusted"),
                    labels = c("Univariable", "Adjusted for age, sex and stage"))
pf <- ggplot(fo, aes(HR, Gene, colour = Cohort)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_pointrange(aes(xmin = LCL, xmax = UCL), position = position_dodge(width = 0.6)) +
  scale_x_log10(breaks = c(0.8, 1, 1.25, 1.5, 2)) +
  scale_colour_manual(values = c("GSE81089" = "#0072B2", "TCGA-LUAD" = "#D55E00")) +
  facet_wrap(~ Model) + theme_bw() +
  labs(x = "Hazard ratio per 1 SD increase in expression (95% CI)", y = NULL)
ggsave("08_figures/Fig8_forest_Cox.png", pf, width = 10, height = 4, dpi = 300)
ggsave("08_figures/Fig8_forest_Cox.pdf", pf, width = 10, height = 4)
# CHECK: the plot shows 4 genes, 2 panels (Univariable, Adjusted) and 2 colours (cohort).
#        HR > 1 means high expression = shorter survival.


## ---- 4.5: Fig 6, hub gene expression (GSE81089 + TCGA) ---------------------
vt   <- readRDS("04_TCGA/vsd_TCGA.rds")
rst  <- readRDS("04_TCGA/res_sym_TCGA.rds")
ens_tcga <- vapply(cmp, function(g) {
  i <- which(rst$SYMBOL == g); rst$Ensembl[i][order(rst$padj[i])][1]
}, character(1))
stopifnot(all(ens_tcga %in% rownames(vt)))

make_box <- function(vobj, ens, res, ens_res, title) {
  cdx <- as.data.frame(colData(vobj))
  stopifnot(all(cdx$Condition %in% c("Normal", "Tumor")))
  rows <- lapply(cmp, function(g) {
    padj <- res$padj[which(res$SYMBOL == g)]; padj <- min(padj, na.rm = TRUE)
    data.frame(Gene = sprintf("%s\nadj. p = %s", g, formatC(padj, digits = 2, format = "g")),
               Sample = colnames(vobj), VST = as.numeric(assay(vobj)[ens[g], ]),
               Condition = factor(cdx$Condition, levels = c("Normal", "Tumor")))
  })
  df <- do.call(rbind, rows)
  df$Gene <- factor(df$Gene, levels = unique(df$Gene))
  ggplot(df, aes(Condition, VST, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.5) +
    scale_fill_manual(values = c(Normal = "#0072B2", Tumor = "#D55E00")) +
    facet_wrap(~ Gene, nrow = 1, scales = "free_y") +
    theme_bw() + theme(legend.position = "none") +
    labs(title = title, x = NULL, y = "VST expression (log2 scale)")
}
pa <- make_box(v,  ens_geo,   geo, geo_ens, "GSE81089 (19 normal, 108 tumors)")
pb <- make_box(vt, ens_tcga,  rst, NULL,    "TCGA-LUAD (59 normal, 517 tumors)")
p6 <- (pa / pb) + plot_annotation(tag_levels = "A")
ggsave("08_figures/Fig6_hub_gene_expression.png", p6, width = 11, height = 8, dpi = 300)
ggsave("08_figures/Fig6_hub_gene_expression.pdf", p6, width = 11, height = 8)
# CHECK: 2 panels (A = GSE81089, B = TCGA), each with the 4 hub genes only (CDK1, CCNB1, CCNA2, BRCA1).
#        Each facet title shows the adjusted p-value (DESeq2). Compare with Table 8.
#        The titles state 19/108 and 59/517; confirm with table(cd$Condition)
#        in Part 1 that the counts are the same.



## ---- end: session information ------------------------------------------------
writeLines(capture.output(sessionInfo()), "00_session/sessionInfo_full.txt")
cat("Done. See 00_session/reference_check.csv, 07_tables_for_manuscript/ and 08_figures/\n")

