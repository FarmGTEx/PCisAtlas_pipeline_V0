# ---- 02.GWAS_enrichment.r ----
# GWAS enrichment analysis using scPagwas (TRS score) and qgg Tsum
# Covers: scPagwas at cell type and cell lineage levels,
#         qgg Tsum permutation enrichment
# Extracted from Figure5 pipeline; standalone analysis script

library(ArchR)
addArchRThreads(threads = 32)
options(future.globals.maxSize = 40 * 1024^3)

set.seed(1234)
source("~/pcisatlas/scripts/rFun.scATAC.r")

# Color scheme for cell lineages
mycolors_lineage <- c(
  `Neural cell` = "#2E8B57", `Epithelial cell` = "#377eb8",
  `Stromal cell` = "#FFA500", `Muscle cell` = "#f781bf",
  `Immune cell` = "#984ea3", `Endothelial cell` = "#48D1CC",
  `Endocrine cell` = "#DC143C", `Germ layer cell` = "#CD853F",
  `Germline cell` = "#9C9C9C"
)

# ==============================================================================
# SECTION 1: scPagwas Integration (TRS score analysis)
# ==============================================================================

## ---- 1.1 Prepare single cell data ----

setwd("~/pcisatlas/GWAS/scPagwas/01.sc_data")
archrPig <- readRDS("~/pcisatlas/04.call_peaks/archrPig.Tissue_Date_CellType.peaks.rds")
saveRDS(archrPig@cellColData, "archrPig.metadata.rds")
retained_cells <- subset_cells_archr(archrPig, subset_by = 'Tissue_Date_CellType', n = 500, seed = 1234)

# Split into chunks to avoid memory issues
archrPig1 <- archrPig[retained_cells[1:250000], ]
peakMTX1 <- getMatrixFromProject(ArchRProj = archrPig1, useMatrix = "GeneScoreMatrix", binarize = FALSE)
saveRDS(peakMTX1, "peakMTX.1.rds")
archrPig2 <- archrPig[retained_cells[250001:530221], ]
peakMTX2 <- getMatrixFromProject(ArchRProj = archrPig2, useMatrix = "GeneScoreMatrix", binarize = FALSE)
saveRDS(peakMTX2, "peakMTX.2.rds")

# Only keep protein coding genes
peakMTX1 <- readRDS("peakMTX.1.rds")
peakMTX2 <- readRDS("peakMTX.2.rds")
rownames(peakMTX1@assays@data$GeneScoreMatrix) <- peakMTX1@elementMetadata@listData$name
rownames(peakMTX2@assays@data$GeneScoreMatrix) <- peakMTX2@elementMetadata@listData$name
genes <- read.table("~/pcisatlas/reference/genome/Sus_scrofa.Sscrofa11.1.gene.bed", stringsAsFactor = F, sep = "\t")
genes <- subset(genes, V5 == "protein_coding")
peakMTX1 <- peakMTX1@assays@data$GeneScoreMatrix[
  intersect(rownames(peakMTX1@assays@data$GeneScoreMatrix), genes$V4), ]
peakMTX2 <- peakMTX2@assays@data$GeneScoreMatrix[
  intersect(rownames(peakMTX2@assays@data$GeneScoreMatrix), genes$V4), ]
peakMTX <- cbind(peakMTX1, peakMTX2)
saveRDS(peakMTX, "peakMTX.pc.rds")
rm(peakMTX1, peakMTX2, archrPig1, archrPig2); gc()

# Create Seurat object
library(Seurat)
library(dplyr)
peakMTX <- readRDS("peakMTX.pc.rds")
metadata <- readRDS("archrPig.metadata.rds")
metadata <- metadata[colnames(peakMTX), c('Sample','Tissue','Date','CellType','CellLineage')]
pbmc <- CreateSeuratObject(counts = peakMTX, project = "pcisatlas",
                           min.cells = 3, min.features = 200,
                           meta.data = as.data.frame(metadata))
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
Idents(pbmc) <- "CellType"
saveRDS(pbmc, "seurat.rds")

## ---- 1.2 Prepare GWAS summary statistics ----
# Shell commands: Convert GWAS results to Z-score format (chrom/pos/rsid/beta/se/maf)

"

cd ~/pcisatlas/GWAS/scPagwas/02.gwasdata

# DDLLYY population (growth, reproduction, sperm)
for i in ~/pcisatlas/GWAS/gwas_sumstats/growth/*.assoc.txt;
do
    NAME=$(basename $i .assoc.txt);
    sed '1d' $i | awk -v FS='\t' -v OFS='\t' '{print $2,$3,$1,$7/$8,$8,$6}' | sed -e '1ichrom\tpos\trsid\tbeta\tse\tmaf' > $NAME &
done
mv NMP NMP_repro

for i in ~/pcisatlas/GWAS/gwas_sumstats/sperm/*.assoc.txt;
do
    NAME=$(basename $i .assoc.txt);
    sed '1d' $i | awk -v FS='\t' -v OFS='\t' '{print $2,$3,$1,$7/$8,$8,$6}' | sed -e '1ichrom\tpos\trsid\tbeta\tse\tmaf' > $NAME &
done
mv NMP NMP_sperm

# ZL population (amino acid, meat quality)
for i in ~/pcisatlas/GWAS/gwas_sumstats/AA/*.assoc.txt;
do
    NAME=$(basename $i .assoc.txt);
    sed '1d' $i | awk -v FS='\t' -v OFS='\t' '{print $2,$3,$1,$7/$8,$8,$6}' | sed -e '1ichrom\tpos\trsid\tbeta\tse\tmaf' > $NAME &
done
mv IMF IMF_AA

for i in ~/pcisatlas/GWAS/gwas_sumstats/meat_quality/*.assoc.txt;
do
    NAME=$(basename $i .assoc.txt);
    sed '1d' $i | awk -v FS='\t' -v OFS='\t' '{print $2,$3,$1,$7/$8,$8,$6}' | sed -e '1ichrom\tpos\trsid\tbeta\tse\tmaf' > $NAME &
done
mv IMF IMF_meat_quality

# CT population (CT body composition)
for i in ~/pcisatlas/GWAS/gwas_sumstats/CT/*.assoc.txt;
do
    NAME=$(basename $i .assoc.txt);
    sed '1d' $i | awk -v FS='\t' -v OFS='\t' '{print $2,$3,$1,$7/$8,$8,$6}' | sed -e '1ichrom\tpos\trsid\tbeta\tse\tmaf' > $NAME &
done
"

## ---- 1.3 Prepare KEGG pathway data ----
# Download human KEGG pathways and map to pig genes via biomaRt
setwd("~/pcisatlas/GWAS/scPagwas/03.kegg_pathway")

# Pre-computed outputs are saved as RDS files.

library(KEGGREST)
pathways.list <- keggList("pathway", "hsa")
pathway.codes <- sub("path:", "", names(pathways.list))
genes.by.pathway_kegg <- sapply(pathway.codes, function(pwid) {
  pw <- keggGet(pwid)
  if (is.null(pw[[1]]$GENE)) return(NA)
  pw2 <- pw[[1]]$GENE[c(FALSE, TRUE)]
  pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x) x[1]))
  return(pw2)
})
names(genes.by.pathway_kegg) <- names(pathways.list)
genes.by.pathway_kegg <- genes.by.pathway_kegg[!is.na(genes.by.pathway_kegg)]
saveRDS(genes.by.pathway_kegg, "genes.by.pathway_kegg.human.rds")

# Map human genes to pig orthologs using biomaRt (Ensembl 101 archive)
library(biomaRt)
host_101 <- "https://aug2020.archive.ensembl.org"
human_mart_101 <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host_101)
pig_mart_101 <- useMart("ensembl", dataset = "sscrofa_gene_ensembl", host = host_101)

all_h_genes <- unique(unlist(genes.by.pathway_kegg))
homo_map_raw <- getLDS(
  attributes = c("external_gene_name"),
  filters = "external_gene_name",
  values = all_h_genes,
  mart = human_mart_101,
  attributesL = c("ensembl_gene_id"),
  martL = pig_mart_101,
  uniqueRows = TRUE
)
colnames(homo_map_raw) <- c("Human_Symbol", "Pig_EnsemblID")

pig_genes_table <- read.table("~/pcisatlas/reference/genome/Sus_scrofa.Sscrofa11.1.gene.bed",
                              stringsAsFactor = F, sep = "\t")
pig_genes_table <- pig_genes_table[, c('V4', 'V8')]
colnames(pig_genes_table) <- c('Custom_Name', 'ensembl_gene_id')

homo_map_custom <- merge(homo_map_raw, pig_genes_table,
                         by.x = "Pig_EnsemblID", by.y = "ensembl_gene_id")
map_dict <- split(homo_map_custom$Custom_Name, homo_map_custom$Human_Symbol)

genes.by.pathway_pig <- lapply(genes.by.pathway_kegg, function(h_genes) {
  p_genes <- unname(unlist(map_dict[h_genes]))
  return(unique(p_genes[!is.na(p_genes)]))
})
genes.by.pathway_pig <- genes.by.pathway_pig[sapply(genes.by.pathway_pig, length) > 0]
saveRDS(genes.by.pathway_pig, "genes.by.pathway_kegg.pig.rds")

## ---- 1.4 Prepare block annotation ----
setwd("~/pcisatlas/GWAS/scPagwas/04.block_annotation")

pig_block_annotation <- read.table("~/pcisatlas/reference/genome/Sus_scrofa.Sscrofa11.1.gene.bed",
                                   stringsAsFactor = F, sep = "\t")
pig_block_annotation <- pig_block_annotation[, 1:4]
colnames(pig_block_annotation) <- c("chrom", "start", "end", "label")
pig_block_annotation$chrom <- gsub("^chr", "", pig_block_annotation$chrom)
pig_block_annotation$start <- as.numeric(pig_block_annotation$start)
pig_block_annotation$end <- as.numeric(pig_block_annotation$end)
pig_block_annotation <- pig_block_annotation[!duplicated(pig_block_annotation$label), ]
pig_block_annotation <- pig_block_annotation[
  pig_block_annotation$chrom %in% c(as.character(1:18), "X", "Y", "MT"), ]
saveRDS(pig_block_annotation, "pig_block_annotation_custom.RDS")

## ---- 1.5 Prepare LD data ----
# Shell commands: Calculate LD using PLINK for each population
"

cd ~/pcisatlas/GWAS/scPagwas/05.LD

plink --bfile ~/pcisatlas/GWAS/gwas_sumstats/merge_DDLLYY \
  --r2 --ld-window-kb 1000 --ld-window-r2 0.2 --autosome --threads 64 --out DDLLYY2
plink --bfile ~/pcisatlas/GWAS/gwas_sumstats/merge_ZL \
  --r2 --ld-window-kb 1000 --ld-window-r2 0.2 --autosome --threads 64 --out ZL2
plink --bfile ~/pcisatlas/GWAS/gwas_sumstats/sort_seq_795ind_imp_maf0.05 \
  --r2 --ld-window-kb 1000 --ld-window-r2 0.2 --autosome --threads 64 --out CT2
"

# Process LD into per-chromosome lists for scPagwas
library(data.table)

all_ld <- fread("DDLLYY2.ld")
all_ld[, R := sqrt(R2)]
chrom_ld <- lapply(1:18, function(i) {
  tmp <- all_ld[CHR_A == i, .(SNP_A, SNP_B, R)]
  return(tmp)
})
saveRDS(chrom_ld, file = "DDLLYY2.LD.RDS")
# Repeat for ZL2 and CT2 LD files

## ---- 1.6 Run scPagwas (Cell Type level, 500 cells/type) ----

### 1.6.1 Prepare single cell data (subset by CellType)
setwd("~/pcisatlas/GWAS/scPagwas/12.sc_data3/01.sc_data")
archrPig <- readRDS("~/pcisatlas/04.call_peaks/archrPig.Tissue_Date_CellType.peaks.rds")
retained_cells <- subset_cells_archr(archrPig, subset_by = 'CellType', n = 500, seed = 1234)
archrPig <- archrPig[retained_cells, ]
saveRDS(archrPig@cellColData, "archrPig.metadata.rds")
peakMTX <- getMatrixFromProject(ArchRProj = archrPig, useMatrix = "GeneScoreMatrix", binarize = FALSE)
saveRDS(peakMTX, "peakMTX.rds")

library(Seurat)
library(dplyr)
peakMTX <- readRDS("peakMTX.rds")
rownames(peakMTX@assays@data$GeneScoreMatrix) <- peakMTX@elementMetadata@listData$name
peakMTX <- peakMTX@assays@data$GeneScoreMatrix
metadata <- readRDS("archrPig.metadata.rds")
metadata <- metadata[colnames(peakMTX), c('Sample','Tissue','Date','CellType','CellLineage')]
pbmc <- CreateSeuratObject(counts = peakMTX, project = "pcisatlas",
                           min.cells = 3, min.features = 200,
                           meta.data = as.data.frame(metadata))
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
Idents(pbmc) <- "CellType"
saveRDS(pbmc, "seurat.rds")

### 1.6.2 Run scPagwas (Single_data_input + Pathway_pcascore_run)
setwd("~/pcisatlas/GWAS/scPagwas/12.sc_data3/02.run")
library(scPagwas)
Pagwas <- list()
Single_data <- readRDS('~/pcisatlas/GWAS/12.sc_data3/01.sc_data/seurat.rds')
Genes_by_pathway_kegg <- readRDS("~/pcisatlas/GWAS/03.kegg_pathway/genes.by.pathway_kegg.pig.rds")
Pagwas <- Single_data_input(Pagwas = Pagwas, assay = "RNA",
                            Single_data = Single_data, Pathway_list = Genes_by_pathway_kegg)
Pagwas <- Pathway_pcascore_run(Pagwas = Pagwas, Pathway_list = Genes_by_pathway_kegg,
                                max.pathway.size = 300)
saveRDS(Pagwas, "Pagwas.rds")

### 1.6.3 Shell commands: Submit scPagwas jobs per trait (with LD file selection)

"

cd ~/pcisatlas/GWAS/scPagwas/run

# DDLLYY traits
for i in ADG100 AGE BF BodyH BodyL BW DFI FCR LEA LMP100 rLH TTN ABW GL HPR LBW LL NBA NHP NMF NMP_repro NSP NWP PSR TNB WSI ABN DEN FSN MOT NMP_sperm SPabn TSN VOL VSL;
do
    LDFILE=DDLLYY
    sed -e 's/TRAITNAME/'$i'/g' ref.aip | sed -e 's/TRAITLD/'$LDFILE'/g' > scPagwas.$i.aip && csub < scPagwas.$i.aip && sleep 1
done

# ZL traits
for i in Ala ALA.C18.3n3 ARA.C20.4n6 Arg Asp C18.1n9 C18.2n6 C20.1 C20.2 Chol Glu Gly His Ile IMF_AA IMP Leu Lys Met Myristic.C14.0 Palmitic.C16.0 Palmitoleic.C16.1 Phe Pro Ser Stearic.C18.0 TAA TFA Thr Tyr Val a-0 a-24 b-0 b-24 DL.24h DL.48h gmIMF IMF_meat_quality L-0 L-24 PH.24h PH.45min;
do
    LDFILE=ZL
    sed -e 's/TRAITNAME/'$i'/g' ref.aip | sed -e 's/TRAITLD/'$LDFILE'/g' > scPagwas.$i.aip && csub < scPagwas.$i.aip && sleep 1
done

# CT traits
for i in abdominal_and_pelvic_cavity_fat_volume abdominal_and_pelvic_cavity_fat_volume_percentage ...; do
    LDFILE=CT
    sed -e 's/TRAITNAME/'$i'/g' ref.aip | sed -e 's/TRAITLD/'$LDFILE'/g' > scPagwas.$i.aip && csub < scPagwas.$i.aip && sleep 1
done

"

### 1.6.4 Shell commands: Merge scPagwas results
"

cd ~/pcisatlas/GWAS/scPagwas/12.sc_data3/02.run

# DDLLYY
head -1 01.run/scPagwas_results.ADG100/ADG100_celltypes_bootstrap_results.csv | sed -e 's/$/,Trait/' > DDLLYY.bs.csv
head -1 01.run/scPagwas_results.ADG100/ADG100_singlecell_scPagwas_score_pvalue.Result.csv | sed -e 's/$/,Trait/' > DDLLYY.trs.csv
for i in ADG100 AGE BF ...; do
    sed '1,2d' ./01.run/scPagwas_results.$i/${i}_celltypes_bootstrap_results.csv | sed -e 's/$/,'$i'/' >> DDLLYY.bs.csv
    sed '1d' ./01.run/scPagwas_results.$i/${i}_singlecell_scPagwas_score_pvalue.Result.csv | sed -e 's/$/,'$i'/' >> DDLLYY.trs.csv
done

# Repeat for ZL and CT ...
"

# SECTION 2: Cell Lineage Level scPagwas Analysis
# ==============================================================================
# scPagwas run at the cell lineage level (5000 cells/lineage) for all traits

## ---- 2.1 Prepare single cell data (subset by CellLineage) ----
setwd("~/pcisatlas/GWAS/scPagwas/celllineage/01.sc_data")
archrPig <- readRDS("~/pcisatlas/data/04.call_peaks/archrPig.Tissue_Date_CellType.peaks.rds")
retained_cells <- subset_cells_archr(archrPig, subset_by = 'CellLineage', n = 5000, seed = 1234)
archrPig <- archrPig[retained_cells, ]
saveRDS(archrPig@cellColData, "archrPig.metadata.rds")
peakMTX <- getMatrixFromProject(ArchRProj = archrPig, useMatrix = "GeneScoreMatrix", binarize = FALSE)
saveRDS(peakMTX, "peakMTX.rds")

library(Seurat)
library(dplyr)
peakMTX <- readRDS("peakMTX.rds")
rownames(peakMTX@assays@data$GeneScoreMatrix) <- peakMTX@elementMetadata@listData$name
peakMTX <- peakMTX@assays@data$GeneScoreMatrix
metadata <- readRDS("archrPig.metadata.rds")
metadata <- metadata[colnames(peakMTX), c('Sample','Tissue','Date','CellType','CellLineage')]
pbmc <- CreateSeuratObject(counts = peakMTX, project = "pbmc3k",
  min.cells = 3, min.features = 200, meta.data = as.data.frame(metadata))
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
Idents(pbmc) <- "CellLineage"
saveRDS(pbmc, "seurat.rds")


## ---- 2.2 Run scPagwas with KEGG pathways ----
setwd("~/pcisatlas/GWAS/scPagwas/celllineage/02.run")
library(scPagwas)
Pagwas <- list()
Single_data <- readRDS('~/pcisatlas/GWAS/scPagwas/celllineage/01.sc_data/seurat.rds')
Genes_by_pathway_kegg <- readRDS("~/pcisatlas/GWAS/scPagwas/kegg/genes.by.pathway_kegg.pig.rds")
Pagwas <- Single_data_input(Pagwas = Pagwas, assay = "RNA",
  Single_data = Single_data, Pathway_list = Genes_by_pathway_kegg)
Pagwas <- Pathway_pcascore_run(Pagwas = Pagwas, Pathway_list = Genes_by_pathway_kegg,
  max.pathway.size = 300)
saveRDS(Pagwas, "Pagwas.rds")


## ---- 2.3 Submit scPagwas jobs per trait ----
# Shell commands: submit jobs for each trait with appropriate LD reference
"
cd ~/pcisatlas/GWAS/scPagwas/celllineage/02.run/01.run.cortGWAS.TRS

# DDLLYY population traits
for i in ADG100 AGE BF BodyH BodyL BW DFI FCR LEA LMP100 rLH TTN ABW GL HPR LBW LL NBA NHP NMF NMP_repro NSP NWP PSR TNB WSI ABN DEN FSN MOT NMP_sperm SPabn TSN VOL VSL;
do
    LDFILE=DDLLYY;
    sed -e \"s/TRAITNAME/\$i/g\" ref.aip | sed -e \"s/TRAITLD/\$LDFILE/g\" > scPagwas.\$i.aip && csub < scPagwas.\$i.aip && sleep 1
done

# ZL population traits
for i in Ala ALA.C18.3n3 ARA.C20.4n6 Arg Asp C18.1n9 C18.2n6 C20.1 C20.2 Chol Glu Gly His Ile IMF_AA IMP Leu Lys Met Myristic.C14.0 Palmitic.C16.0 Palmitoleic.C16.1 Phe Pro Ser Stearic.C18.0 TAA TFA Thr Tyr Val a-0 a-24 b-0 b-24 DL.24h DL.48h gmIMF IMF_meat_quality L-0 L-24 PH.24h PH.45min;
do
    LDFILE=ZL;
    sed -e \"s/TRAITNAME/\$i/g\" ref.aip | sed -e \"s/TRAITLD/\$LDFILE/g\" > scPagwas.\$i.aip && csub < scPagwas.\$i.aip && sleep 1
done

# CT population traits
for i in abdominal_and_pelvic_cavity_fat_volume ... ;  # truncated for brevity
do
    LDFILE=CT;
    sed -e \"s/TRAITNAME/\$i/g\" ref.aip | sed -e \"s/TRAITLD/\$LDFILE/g\" > scPagwas.\$i.aip && csub < scPagwas.\$i.aip && sleep 1
done
"

## ---- 2.4 Merge scPagwas results ----
"
cd ~/pcisatlas/GWAS/scPagwas/celllineage/02.run

# DDLLYY
head -1 01.run.cortGWAS.TRS/scPagwas_results.ADG100/ADG100_celltypes_bootstrap_results.csv | sed -e 's/$/,Trait/' > DDLLYY.bs.csv
head -1 01.run.cortGWAS.TRS/scPagwas_results.ADG100/ADG100_singlecell_scPagwas_score_pvalue.Result.csv | sed -e 's/$/,Trait/' > DDLLYY.trs.csv
for i in ADG100 AGE BF BodyH BodyL BW DFI FCR LEA LMP100 rLH TTN ABW GL HPR LBW LL NBA NHP NMF NMP_repro NSP NWP PSR TNB WSI ABN DEN FSN MOT NMP_sperm SPabn TSN VOL VSL;
do
    sed '1,2d' ./01.run.cortGWAS.TRS/scPagwas_results.\$i/\${i}_celltypes_bootstrap_results.csv | sed -e 's/$/,'\$i'/' >> DDLLYY.bs.csv
    sed '1d' ./01.run.cortGWAS.TRS/scPagwas_results.\$i/\${i}_singlecell_scPagwas_score_pvalue.Result.csv | sed -e 's/$/,'\$i'/' >> DDLLYY.trs.csv
done

# ZL
head -1 01.run.cortGWAS.TRS/scPagwas_results.Ala/Ala_celltypes_bootstrap_results.csv | sed -e 's/$/,Trait/' > ZL.bs.csv
head -1 01.run.cortGWAS.TRS/scPagwas_results.Ala/Ala_singlecell_scPagwas_score_pvalue.Result.csv | sed -e 's/$/,Trait/' > ZL.trs.csv
for i in Ala ALA.C18.3n3 ARA.C20.4n6 Arg Asp C18.1n9 C18.2n6 C20.1 C20.2 Chol Glu Gly His Ile IMF_AA IMP Leu Lys Met Myristic.C14.0 Palmitic.C16.0 Palmitoleic.C16.1 Phe Pro Ser Stearic.C18.0 TAA TFA Thr Tyr Val a-0 a-24 b-0 b-24 DL.24h DL.48h gmIMF IMF_meat_quality L-0 L-24 PH.24h PH.45min;
do
    sed '1,2d' ./01.run.cortGWAS.TRS/scPagwas_results.\$i/\${i}_celltypes_bootstrap_results.csv | sed -e 's/$/,'\$i'/' >> ZL.bs.csv
    sed '1d' ./01.run.cortGWAS.TRS/scPagwas_results.\$i/\${i}_singlecell_scPagwas_score_pvalue.Result.csv | sed -e 's/$/,'\$i'/' >> ZL.trs.csv
done

# CT
head -1 01.run.cortGWAS.TRS/scPagwas_results.abdominal_and_pelvic_cavity_fat_volume/abdominal_and_pelvic_cavity_fat_volume_celltypes_bootstrap_results.csv | sed -e 's/$/,Trait/' > CT.bs.csv
head -1 01.run.cortGWAS.TRS/scPagwas_results.abdominal_and_pelvic_cavity_fat_volume/abdominal_and_pelvic_cavity_fat_volume_singlecell_scPagwas_score_pvalue.Result.csv | sed -e 's/$/,Trait/' > CT.trs.csv
for i in abdominal_and_pelvic_cavity_fat_volume ... ;  # truncated for brevity
do
    sed '1,2d' ./01.run.cortGWAS.TRS/scPagwas_results.\$i/\${i}_celltypes_bootstrap_results.csv | sed -e 's/$/,'\$i'/' >> CT.bs.csv
    sed '1d' ./01.run.cortGWAS.TRS/scPagwas_results.\$i/\${i}_singlecell_scPagwas_score_pvalue.Result.csv | sed -e 's/$/,'\$i'/' >> CT.trs.csv
done
"


# ==============================================================================
# SECTION 3: qgg Tsum GWAS Enrichment (Selection-Overlapped Regions)
# ==============================================================================
# qgg Tsum permutation test for GWAS enrichment of selection signals
# in cell-type-specific cCRE modules

# ---- 2. qgg Tsum GWAS enrichment (selection-overlapped regions) ----

    ##region prepare peak data: intersect selection signals with module peaks
"
cd ~/pcisatlas/selection_sweep/qgg_Tsum/01.peak
for i in M56 M31 M17 M134 M4 M72 M64 M48 M77 M128 M76 M0 M50 M104 M18 M79 M19 M112 M108 M139 M91 M102 M87 M53 M96 M69 M140 M129 M137 M133 M136 M3 M71 M115;
do
  bedtools intersect \\
    -a ~/pcisatlas/selection_sweep/input/ASD.fst_xpehh.01.bed \\
    -b ~/pcisatlas/selection_sweep/qgg_Tsum/module_peaks/\$i.peaks.bed \\
    -wo |awk -v OFS=\"\\t\" '{print \$1,\$2,\$3,\"ASD---\"\$9}' |sort -u > 01.ASD---\$i.bed && \\
  bedtools intersect \\
    -a ~/pcisatlas/selection_sweep/input/WhiteD.fst_xpehh.01.bed \\
    -b ~/pcisatlas/selection_sweep/qgg_Tsum/module_peaks/\$i.peaks.bed \\
    -wo |awk -v OFS=\"\\t\" '{print \$1,\$2,\$3,\"WhiteD---\"\$9}' |sort -u > 01.WhiteD---\$i.bed &
done

cat 01.*bed > 02.all.bed
"
    ##endregion

    ##region get SNPs overlapping these regions
"
cd ~/pcisatlas/selection_sweep/qgg_Tsum/02.snp
for i in ~/pcisatlas/selection_sweep/qgg_Tsum/02.snps/01.gwas/*;
do
  NAME=\`basename \$i\`;
  bedtools intersect -a ../01.peak/02.all.bed -b \$i -wo |awk -v OFS=\"\\t\" '{print \$8,\$4}' > \$NAME.snps &
done
"
    ##endregion

    ##region submit qgg Tsum permutation jobs via csub + ref.aip
"
cd ~/pcisatlas/selection_sweep/qgg_Tsum/03.permutation
for i in ~/pcisatlas/selection_sweep/qgg_Tsum/01.gwas/*;
do
  NAME=\`basename \$i\`;
  sed -e \"s/GROUPNAME/\$NAME/g\" ref.aip > \$NAME.aip && csub < \$NAME.aip && sleep 0.5
done

echo -e \"group\\tm\\tstat\\tp\\tCellType\\tTrait\" > merged.res.txt
for i in *.permu.txt;
do
  NAME=\`basename \$i .permu.txt\`;
  sed '1d' \$i |sed -e \"s/$/\\t\$NAME/\" >> merged.res.txt
done
"
    ##endregion

    ##region load merged results and calculate enrichment ratios
    setwd("~/pcisatlas/selection_sweep/qgg_Tsum/04.plot")
    library(dplyr)
    library(reshape2)
    library(plyr)
    library(pheatmap)
    library(RColorBrewer)
    library(ggplot2)
    library(patchwork)
    library(ggpmisc)
    library(stringr)
    library(ggpubr)
    library(tidyr)

    # load data
    res <- read.table("~/pcisatlas/selection_sweep/qgg_Tsum/03.permutation/merged.res.txt",
                      header = T, stringsAsFactor = F, sep = "\t")
    res$perSNP_Z2 <- res$stat / res$m
    res$Trait <- gsub("-", ".", res$Trait)

    # trait annotations
    chinese_name <- read.table("~/pcisatlas/data/trait_annotation.txt",
                               stringsAsFactor = F, fileEncoding = "UTF-8", sep = "\t")
    chinese_name[which(chinese_name$V2 == "NMP" & chinese_name$V3 == "reproduction"), "V2"] <- "NMP_repro"
    chinese_name[which(chinese_name$V2 == "NMP" & chinese_name$V3 == "sperm"), "V2"] <- "NMP_sperm"
    chinese_name[which(chinese_name$V2 == "IMF" & chinese_name$V3 == "AA"), "V2"] <- "IMF_AA"
    chinese_name[which(chinese_name$V2 == "IMF" & chinese_name$V3 == "meat_quality"), "V2"] <- "IMF_meat_quality"
    res$Trait_Chinese <- mapvalues(res$Trait, chinese_name$V2, chinese_name$V1)
    res$Trait_Category <- mapvalues(res$Trait, chinese_name$V2, chinese_name$V5)

    res$Select_group <- res$group %>% strsplit(., "---") %>% sapply(., function(x){x[[1]]})
    res$Module <- res$group %>% strsplit(., "---") %>% sapply(., function(x){x[[2]]})

    # annotate cell types
    CTinfo2 <- read.table("~/pcisatlas/data/module_celltypes.txt",
                          stringsAsFactor = F, sep = "\t")
    CTinfo2 <- CTinfo2[, 2:3] %>% unique
    res$CellType <- mapvalues(res$Module, CTinfo2$V2, CTinfo2$V3)
    res$negLog10P <- -log10(res$p)
    res$negLog10P[is.infinite(res$negLog10P)] <- 4

    # calculate enrichment ratio
    allback <- read.table("~/pcisatlas/selection_sweep/qgg_Tsum/background.txt",
                          stringsAsFactor = F, sep = "\t")
    allback$V1 <- gsub("-", ".", allback$V1)
    res$background <- mapvalues(res$Trait, allback$V1, allback$V2) %>% as.numeric
    res$Enrichment <- (res$stat / res$m) / res$background
    ##endregion

# clustering of the enrichment matrix. Below is an alternative version using
# ASD clustering order for the Asian panel.

    ##region alternative ordering: clustered by ASD
    res2 <- subset(res,
        (Select_group == "ASD" & Module %in% c("M56","M31","M17","M134","M4","M72","M64","M48","M77","M128","M76","M0","M50","M104","M18","M79","M19","M112","M108","M139","M91","M102")) |
        (Select_group == "WhiteD" & Module %in% c("M112","M87","M53","M91","M96","M69","M140","M129","M137","M133","M136","M64","M3","M71","M115","M18")))
    res2 <- res2 %>% group_by(Trait, Select_group) %>%
        mutate(Enrichment_Zscore = scale(Enrichment) %>% as.numeric()) %>%
        ungroup() %>% as.data.frame

    res3 <- subset(res2, p < 0.05 & Enrichment > 1)
    res3 <- subset(res3, Trait_Chinese %in% c("平均初生重","达100公斤平均日增重","达100公斤日龄",
                                              "100公斤背膘厚","体长","出生重","精子浓度",
                                              "30至120公斤日采食量","30至120公斤校正饲料转化率",
                                              "窝重","100公斤眼肌面积","100公斤瘦肉率","正常精子率","死胎数"))
    traitorder <- c("AGE","DFI","NSP","ABW","ADG100","LEA","FCR","NMP_sperm","BF","LBW","LMP100","BodyL","BW","DEN") %>% rev
    CTorder <- c("Bipolar.cell","Rod.photoreceptor.cell","Di.mesencephalon.inhibitory.neuron","Sncg, Vip",
                 "Goblet.cell","Parietal.cell","Hepatocyte","Leydig.cell","Ciliated.cell",
                 "Ascending.loop.of.Henle.cell","Megakaryocyte","Plasma.cell","T.cell","B.cell",
                 "Kupffer.cell","Dendritic.cell, Macrophage","Somatotrope","Ependymal.cell",
                 "Neural.crest..PNS.glia.","Satellite.cell","Cardiomyocyte","Mesenchymal.cell",
                 "Common.M3","Common.M0")

    res_plot <- res3 %>%
        mutate(Group = if_else(str_detect(Select_group, "ASD"), "Asian", "European")) %>%
        mutate(
            CellType = factor(CellType, levels = CTorder),
            Trait = factor(Trait, levels = traitorder)
        ) %>%
        mutate(
            x = as.numeric(CellType),
            y = as.numeric(Trait)
        ) %>%
        mutate(
            start = if_else(Group == "Asian", 3 * pi / 4, -pi / 4),
            end   = if_else(Group == "Asian", 7 * pi / 4, 3 * pi / 4),
            radius_val = abs(Enrichment_Zscore)
        )

    res_asian <- res_plot %>% filter(Group == "Asian")
    res_euro  <- res_plot %>% filter(Group == "European")
    num_x <- length(CTorder)
    num_y <- length(traitorder)
    global_max_p <- max(res_plot$negLog10P, na.rm = TRUE)

    p_main <- ggplot() +
        coord_fixed(ratio = 1) +
        geom_vline(xintercept = seq(0.5, num_x + 0.5, 1), color = "grey90", size = 0.2) +
        geom_hline(yintercept = seq(0.5, num_y + 0.5, 1), color = "grey90", size = 0.2) +
        new_scale_fill() +
        geom_arc_bar(data = res_asian,
            aes(x0 = x, y0 = y, r0 = 0,
                r = rescale(radius_val, to = c(0.15, 0.45)),
                start = start, end = end, fill = negLog10P),
            color = "white", size = 0.05) +
        scale_fill_gradient(low = "#fee0d2", high = "#de2d26",
            limits = c(0, global_max_p), name = "Asian -log10(p)") +
        new_scale_fill() +
        geom_arc_bar(data = res_euro,
            aes(x0 = x, y0 = y, r0 = 0,
                r = rescale(radius_val, to = c(0.15, 0.45)),
                start = start, end = end, fill = negLog10P),
            color = "white", size = 0.05) +
        scale_fill_gradient(low = "#deebf7", high = "#3182bd",
            limits = c(0, global_max_p), name = "European -log10(p)") +
        geom_point(data = res_plot, aes(x = x, y = y, size = radius_val), alpha = 0) +
        scale_size_continuous(range = c(1, 4), name = "Ratio_Zscore",
            breaks = c(1, 2, 3), labels = c("1", "2", "3")) +
        scale_x_continuous(breaks = 1:num_x, labels = CTorder, expand = c(0, 0)) +
        scale_y_continuous(breaks = 1:num_y, labels = traitorder, expand = c(0, 0)) +
        theme_minimal() +
        theme(
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
            axis.text.y = element_text(size = 7),
            plot.margin = margin(5, 5, 5, 0),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 7)
        ) +
        guides(
            fill = guide_colorbar(barwidth = 0.8, barheight = 4),
            size = guide_legend(
                title = "Ratio_Zscore",
                override.aes = list(alpha = 1, color = "black", size = c(1, 2, 3))
            )
        )

    pdf("splitDotPlot.enrich_Zscore.clustered.pdf", height = 6, width = 6)
    print(p_main)
    dev.off()
    ##endregion


        ##region output trait~celltype pairs for scPagwas-date analysis
        # Filter significant qgg results to find trait-celltype pairs with
        # sufficient cell counts across developmental stages for scPagwas
        tmpdf <- res3[-grep("Common", res3$CellType), c("Trait", "CellType")] %>%
          unique %>%
          separate_rows(CellType, sep = ",\\s*") %>%
          as.data.frame
        metadata <- readRDS("~/pcisatlas/data/cell_metadata.rds")
        CT_date.MTX <- table(metadata$CellType, metadata$Date) %>% unclass %>% as.data.frame
        matrix.trs <- CT_date.MTX[unique(tmpdf$CellType),] %>% unique
        counts.20 <- rowSums(matrix.trs > 20)
        counts.20 <- names(counts.20)[counts.20 > 1]
        tmpdf <- subset(tmpdf, CellType %in% counts.20)
        write.table(tmpdf, "compare.trait_CT.txt", quote=F, sep = "\t", row.name = F, col.name = F)
        write.table(unique(tmpdf[,c('CellType'),drop=F]), "compare.CT.txt", quote=F, sep = "\t", row.name = F, col.name = F)
        ##endregion

    ##endregion

# ---- End of GWAS enrichment analysis ----
# Key outputs:
#   1. scPagwas single-cell data preparation and PCA
#   2. GWAS summary statistics conversion
#   3. KEGG pathway mapping (human -> pig)
#   4. scPagwas job submission per trait
#   5. Cell lineage-level scPagwas (5000 cells/lineage)
