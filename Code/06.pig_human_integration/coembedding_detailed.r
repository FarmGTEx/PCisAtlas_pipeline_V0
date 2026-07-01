# ============================================================
# Pig-Human Cross-Species Integration — Detailed Downstream Analysis
# Complements: coembedding.r (basic co-embedding in this directory)
# ============================================================
# This script covers the detailed pig-human integration workflow:
#   1. Pig adult (180d) vs human adult co-embedding UMAP
#   2. Pig fetal (E43) vs human fetal co-embedding UMAP
#   3. Conserved CRE identification (1000-iteration liftover random sampling)
#   4. Pig 180d vs E43 differential peaks + human CPM comparison
#   5. Diff peak motif enrichment
# ============================================================

library(ArchR)
library(BSgenome.Sscrofa.UCSC.susScr11)
library(TxDb.Sscrofa.Ensembl.susScr11.V1)
library(org.Spig.eg.db)
library(dplyr)
library(Seurat)
library(Signac)

set.seed(1234)
addArchRThreads(threads = 88)
options(future.globals.maxSize = 50 * 1024^3)
source("~/pcisatlas/scripts/rFun.scATAC.r")
source("~/pcisatlas/scripts/R_function.r")

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Sscrofa.UCSC.susScr11)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Sscrofa.Ensembl.susScr11.V1, OrgDb = org.Spig.eg.db)
ATAC_project <- readRDS("~/pcisatlas/data/Save-ArchR-Project.rds")

# ================================================================
# WORKFLOW LOGIC (important cross-reference for manuscript methods)
# ================================================================
# 1. Extract pig tissues → match with human tissues in Excel
#    → extract common cell types (same tissue origin in both species)
# 2. Adult co-embedding: merged ALL human adult fragments;
#    pig: merged common-tissue fragments
# 3. Fetal co-embedding: extract 1000-cell subsets, merge on-the-fly
#    (fragment files too large to merge all)
# 4. Conserved CRE: liftover pig→human, match peak counts by random
#    sampling (1000 iterations), bar = mean of 1000 iterations
# 5. Diff peaks: call in pig (180d vs E43), liftover to human,
#    calculate human CPM from human adult/fetal fragments
# ================================================================

# ------------------------------------------------------------
# Part 1: Pig Adult (180d) vs Human Adult Co-embedding
# ------------------------------------------------------------
# Step 1: Extract pig 180d tissues
sub_ATAC_project <- ATAC_project[rownames(ATAC_project@cellColData[ATAC_project@cellColData$Date %in% c("180d"),]),]
write.table(unique(sub_ATAC_project@cellColData$Tissue), file = "01.pig_180d_tissue", sep = "\t", quote = F, col.names = F, row.names = F)

# Step 2: Select common pig-human tissues (from Excel sheet matching)
# 16 common tissues selected
ATAC_tissues <- c("Aorta", "LongissimusDorsi", "Lung", "Ovary", "FrontalLobe", "Liver",
  "Uterus", "RightAtrium", "LeftAtrium", "LeftVentricle", "RightVentricle",
  "Colon", "Duodenum", "Ileum", "Jejunum", "Esophagus")
sub_ATAC_project <- sub_ATAC_project[rownames(sub_ATAC_project@cellColData[sub_ATAC_project@cellColData$Tissue %in% ATAC_tissues,]),]

# Step 3: Extract pig cell types; do the same for human
write.table(unique(sub_ATAC_project@cellColData$CellType), file = "02.pig_180d_celltype", sep = "\t", quote = F, col.names = F, row.names = F)

# Step 4: Human metadata for common tissues (extract from full human metadata)
# Step 5: Select common cell types (MUST have same tissue origin)
#   Check human: cat 01.adult_metadata | awk '{if ($5=="CellType"){print $4}}' | sort | uniq -c
#   Check pig: table(sub_ATAC_project@cellColData$Tissue) for given cell type
# Step 6: Prepare human metadata with common cell type annotations
# Step 7: Run co-embedding (Signac/Seurat integration)
#
# 30 cell types selected for adult co-embedding; neuron subtypes merged to IN/EX

# ------------------------------------------------------------
# Part 2: Pig Fetal (E43) vs Human Fetal Co-embedding
# ------------------------------------------------------------
sub_ATAC_project <- ATAC_project[rownames(ATAC_project@cellColData[ATAC_project@cellColData$Date %in% c("E43"),]),]
write.table(unique(sub_ATAC_project@cellColData$Tissue), file = "01.pig_E43_tissue", sep = "\t", quote = F, col.names = F, row.names = F)

ATAC_tissues <- c("Liver", "LongissimusDorsi", "Heart", "Lung", "SmallIntestine",
  "Kidney", "Stomach", "Brain", "Cerebellum")
sub_ATAC_project <- sub_ATAC_project[rownames(sub_ATAC_project@cellColData[sub_ATAC_project@cellColData$Tissue %in% ATAC_tissues,]),]
write.table(unique(sub_ATAC_project@cellColData$CellType), file = "02.pig_E43_celltype", sep = "\t", quote = F, col.names = F, row.names = F)

# Fetal fragment handling: extract per-cell-type, not merge all (files too large)
# Human fetal: .../01.fetal_fragments/99.select_fragments/work.sh
# Pig fetal: .../02.pig_fragments/99.select_cells_fragments/work.sh
# Co-embedding uses ALL pig peaks + ALL human fetal peaks (via Signac integration)

# ------------------------------------------------------------
# Part 3: Conserved CRE — Pig Adult vs Human Adult
# ------------------------------------------------------------
sub_ATAC_project <- ATAC_project[rownames(ATAC_project@cellColData[ATAC_project@cellColData$Date == "180d",]),]
sub_ATAC_project <- sub_ATAC_project[rownames(sub_ATAC_project@cellColData[sub_ATAC_project@cellColData$Tissue %in% ATAC_tissues,]),]

ATAC_celltypes <- c("Hepatocyte", "Oligodendrocyte", "Oligodendrocyte.progenitor.cell",
  "Endocardial.cell", "Microglia", "Satellite.cell", "Type.1.myonuclei", "Type.2.myonuclei",
  "Alveolar.Type.1.cell", "Alveolar.Type.2.cell", "Enterocyte", "Astrocyte", "Goblet.cell",
  "Cardiomyocyte", "Sst", "Sst.Chodl", "Lamp5.1", "Lamp5.2", "Vip", "Sncg", "Pvalb",
  "L2.3.IT", "L5.ET", "L5.IT", "L5.6.NP", "L6.CT", "L6.IT", "L6.IT.Car3", "L6b",
  "Cholinergic.neuron")
sub_ATAC_project <- sub_ATAC_project[rownames(sub_ATAC_project@cellColData[sub_ATAC_project@cellColData$CellType %in% ATAC_celltypes,]),]

# Merge neuron subtypes to IN/EX for peak calling (matching human annotation)
original <- ATAC_celltypes
change <- c("Hepatocyte", "Oligodendrocyte", "Oligodendrocyte.progenitor.cell",
  "Endocardial.cell", "Microglia", "Satellite.cell", "Type.1.myonuclei", "Type.2.myonuclei",
  "Alveolar.Type.1.cell", "Alveolar.Type.2.cell", "Enterocyte", "Astrocyte", "Goblet.cell",
  "Cardiomyocyte", "IN", "IN", "IN", "IN", "IN", "IN", "IN",
  "EX", "EX", "EX", "EX", "EX", "EX", "EX", "EX", "EX")
sub_ATAC_project@cellColData$select_180d_CellType <- mapvalues(sub_ATAC_project@cellColData$CellType, original, change)

sub_ATAC_project <- addGroupCoverages(ArchRProj = sub_ATAC_project, groupBy = "select_180d_CellType",
  useLabels = TRUE, minCells = 50, maxCells = 500, maxFragments = 30 * 10^6,
  minReplicates = 2, maxReplicates = 10, sampleRatio = 0.8, kmerLength = 6, force = T)
sub_ATAC_project <- addReproduciblePeakSet(ArchRProj = sub_ATAC_project, groupBy = "select_180d_CellType",
  peakMethod = "Macs2", reproducibility = "2", peaksPerCell = 1000, maxPeaks = 200000,
  minCells = 40, excludeChr = c("chrM", "chrY"), genomeSize = 2367939388,
  shift = -75, extsize = 150, extendSummits = 250, promoterRegion = c(2500, 1000),
  plot = TRUE, force = T, cutOff = 0.05)

# Human cell-type CRE organization: .../02.human_merge_celltype_CRE/work.sh
# Liftover + random sampling (1000 iterations): .../03.pig_celltype_CRE/02.lift_2_hg38/work.sh
# Statistics of 1000 iterations: .../03.pig_celltype_CRE/02.lift_2_hg38/merge.sh

# ------------------------------------------------------------
# Part 4: Conserved CRE — Pig Fetal (E43) vs Human Fetal
# ------------------------------------------------------------
sub_ATAC_project <- ATAC_project[rownames(ATAC_project@cellColData[ATAC_project@cellColData$Date == c("E43"),]),]
ATAC_tissues <- c("Stomach", "Heart", "SmallIntestine", "Kidney", "Liver", "Lung",
  "Brain", "Cerebellum", "LongissimusDorsi")
sub_ATAC_project <- sub_ATAC_project[rownames(sub_ATAC_project@cellColData[sub_ATAC_project@cellColData$Tissue %in% ATAC_tissues,]),]

ATAC_celltypes <- c("Cardiomyocyte", "Endocardial.cell", "Goblet.cell", "Enteric.neuron",
  "Enterocyte", "Myocyte", "Erythroid.cell", "Hepatocyte", "Lung.epithelial.cell",
  "Metanephric.mesenchymal.cell", "Astrocyte", "Excitatory.neuron", "Cholangiocyte")
sub_ATAC_project <- sub_ATAC_project[rownames(sub_ATAC_project@cellColData[sub_ATAC_project@cellColData$CellType %in% ATAC_celltypes,]),]

original <- ATAC_celltypes
change <- c("Cardiomyocyte", "Endocardial.cell", "Goblet.cell", "Enteric.neuron",
  "Enterocyte", "Myocyte", "Erythroid.cell", "Hepatocyte", "Lung.epithelial.cell",
  "Metanephric.mesenchymal.cell", "Astrocyte", "EX", "Cholangiocyte")
sub_ATAC_project@cellColData$select_E43_CellType <- mapvalues(sub_ATAC_project@cellColData$CellType, original, change)

sub_ATAC_project <- addGroupCoverages(ArchRProj = sub_ATAC_project, groupBy = "select_E43_CellType",
  useLabels = TRUE, minCells = 50, maxCells = 500, maxFragments = 30 * 10^6,
  minReplicates = 2, maxReplicates = 10, sampleRatio = 0.8, kmerLength = 6, force = T)
sub_ATAC_project <- addReproduciblePeakSet(ArchRProj = sub_ATAC_project, groupBy = "select_E43_CellType",
  peakMethod = "Macs2", reproducibility = "2", peaksPerCell = 1000, maxPeaks = 200000,
  minCells = 40, excludeChr = c("chrM", "chrY"), genomeSize = 2367939388,
  shift = -75, extsize = 150, extendSummits = 250, promoterRegion = c(2500, 1000),
  plot = TRUE, force = T, cutOff = 0.05)

# ------------------------------------------------------------
# Part 5: Pig 180d vs E43 Differential Peaks + Human Comparison
# ------------------------------------------------------------
# Step 1: Call diff peaks for 7 cell types (180d vs E43) using monocle3/LR-test
# Step 2: Generate bigWig files → deeptools heatmap
# Step 3: Liftover pig diff peaks to hg38 → calculate human CPM (adult vs fetal)
#   Human fetal fragments not merged (files too large) — extract per cell type
#   Use 96.fetal_human.R and 97.adult_human.R to compute per-peak CPM

# ------------------------------------------------------------
# Part 6: Diff Peak Motif Enrichment
# ------------------------------------------------------------
# HOMER motif enrichment on 180d-UP and E43-UP peak sets
# See: .../04.compare_fetal_adult/10.diff_peak_motif/work1.sh and work2.sh
