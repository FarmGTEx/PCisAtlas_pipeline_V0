# set packages
library("ArchR")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("viridis")
library("RColorBrewer")
library("gridExtra")
library("Seurat")
library("patchwork")
library("data.table")
library("Signac")
library("pdist")
library("TFBSTools")
library("BSgenome.Sscrofa.UCSC.susScr11")
library("BSgenome.Hsapiens.UCSC.hg38")

# set workers
library(future)
plan("multisession", workers = 10)
options(future.globals.maxSize = 5000 * 1024^2)
future.seed=TRUE
set.seed(1234)

## prepare pig object
# subset object for required cell types
ATAC_project <- readRDS("./Save-ArchR-Project.rds")
sub_ATAC_project <- ATAC_project[rownames(ATAC_project@cellColData[ATAC_project@cellColData$Date == "180d",]),]

ATAC_tissues <- c("Aorta", "LongissimusDorsi", "Lung", "Ovary", "FrontalLobe", "Liver", "Uterus", "RightAtrium", "LeftAtrium", "LeftVentricle", "RightVentricle", "Colon", "Duodenum", "Ileum", "Jejunum", "Esophagus")
sub_ATAC_project <- sub_ATAC_project[rownames(sub_ATAC_project@cellColData[sub_ATAC_project@cellColData$Tissue %in% ATAC_tissues,]),]

ATAC_celltypes <- c("Hepatocyte", "Oligodendrocyte", "Oligodendrocyte.progenitor.cell", "Endocardial.cell", "Microglia", "Satellite.cell", "Type.1.myonuclei", "Type.2.myonuclei", "Alveolar.Type.1.cell", "Alveolar.Type.2.cell", "Enterocyte", "Astrocyte", "Goblet.cell", "Cardiomyocyte", "Sst", "Sst.Chodl", "Lamp5.1", "Lamp5.2", "Vip", "Sncg", "Pvalb", "L2.3.IT", "L5.ET", "L5.IT", "L5.6.NP", "L6.CT", "L6.IT", "L6.IT.Car3", "L6b", "Cholinergic.neuron")
sub_ATAC_project <- sub_ATAC_project[rownames(sub_ATAC_project@cellColData[sub_ATAC_project@cellColData$CellType %in% ATAC_celltypes,]),]
original <- c("Hepatocyte", "Oligodendrocyte", "Oligodendrocyte.progenitor.cell", "Endocardial.cell", "Microglia", "Satellite.cell", "Type.1.myonuclei", "Type.2.myonuclei", "Alveolar.Type.1.cell", "Alveolar.Type.2.cell", "Enterocyte", "Astrocyte", "Goblet.cell", "Cardiomyocyte", "Sst", "Sst.Chodl", "Lamp5.1", "Lamp5.2", "Vip", "Sncg", "Pvalb", "L2.3.IT", "L5.ET", "L5.IT", "L5.6.NP", "L6.CT", "L6.IT", "L6.IT.Car3", "L6b", "Cholinergic.neuron")
change <- c("Hepatocyte", "Oligodendrocyte", "Oligodendrocyte.progenitor.cell", "Endocardial.cell", "Microglia", "Satellite.cell", "Type.1.myonuclei", "Type.2.myonuclei", "Alveolar.Type.1.cell", "Alveolar.Type.2.cell", "Enterocyte", "Astrocyte", "Goblet.cell", "Cardiomyocyte", "IN", "IN", "IN", "IN", "IN", "IN", "IN", "EX", "EX", "EX", "EX", "EX", "EX", "EX", "EX", "EX")
sub_ATAC_project@cellColData$CellType <- mapvalues(sub_ATAC_project@cellColData$CellType, original, change)

# 1000 cells for each cell type
cellnames <- subset_cells_archr(object = sub_ATAC_project, subset_by = 'CellType', n = 1000)
sub_ATAC_project <- sub_ATAC_project[cellnames, ]
ATAC_metadata <- as.data.frame(sub_ATAC_project@cellColData[,c("CellType", "Tissue")])

# generate fragment file
sn_pig_cells <- row.names(ATAC_metadata)
sn_pig_fragments_file <- "./03.pig_all_fragments_sort.gz"
sn_pig_fragments <- CreateFragmentObject(path = sn_pig_fragments_file, cells = sn_pig_cells)

# generate matrix file
pig_peak <- read.table("./00.all.peaks", header = F)
colnames(pig_peak) <- c("pig_peak")
pig_peak_granges <- StringToGRanges(regions = pig_peak$pig_peak)
sn_pig_peak_matrix <- FeatureMatrix(
  fragments = sn_pig_fragments,
  features = pig_peak_granges
)

# pig Seurat object
sn_pig_chrom_assay <- CreateChromatinAssay(
  counts = sn_pig_peak_matrix,
  sep = c("-", "-")
)
sn_pig_pbmc <- CreateSeuratObject(
  counts = sn_pig_chrom_assay,
  assay = "peaks",
  meta.data = ATAC_metadata
)

# add motif matrix
JASPAR2024 <- readJASPARMatrix("./reference_motifs/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.modified.txt")
sn_pig_pbmc <- AddMotifs(
  object = sn_pig_pbmc,
  genome = BSgenome.Sscrofa.UCSC.susScr11,
  pfm = JASPAR2024
)
sn_pig_pbmc <- RunChromVAR(
  object = sn_pig_pbmc,
  genome = BSgenome.Sscrofa.UCSC.susScr11
)


## prepare human object
human_peak <- read.table("./03.human_adult_CRE.bed", header = F)
colnames(human_peak) <- c("human_peak")
human_metadata <- read.table("./06.human_celltype_metadata", header = T, row.names = 1, quote = "", sep = "\t")
original <- c("Hepatocyte", "Oligodendrocyte", "Oligodendrocyte Precursor", "Endocardial Cell", "Microglia", "Satellite Cell", "Type I Skeletal Myocyte", "Type II Skeletal Myocyte", "Alveolar Type 1 (AT1) Cell", "Alveolar Type 2 (AT2) Cell", "Small Intestinal Enterocyte", "Astrocyte 1", "Astrocyte 2", "Small Intestinal Goblet Cell", "Colonic Goblet Cell", "Atrial Cardiomyocyte", "Ventricular Cardiomyocyte", "GABAergic Neuron 1", "GABAergic Neuron 2", "Glutamatergic Neuron 1", "Glutamatergic Neuron 2")
change <- c("Hepatocyte", "Oligodendrocyte", "Oligodendrocyte.progenitor.cell", "Endocardial.cell", "Microglia", "Satellite.cell", "Type.1.myonuclei", "Type.2.myonuclei", "Alveolar.Type.1.cell", "Alveolar.Type.2.cell", "Enterocyte", "Astrocyte", "Astrocyte", "Goblet.cell", "Goblet.cell", "Cardiomyocyte", "Cardiomyocyte", "IN", "IN", "EX", "EX")
human_metadata$cell.type <- mapvalues(human_metadata$cell.type, original, change)
human_metadata$cells <- rownames(human_metadata)
grouped_data <- human_metadata %>% group_by(cell.type)

# 1000 cells fro each cell type
sampled_data <- grouped_data %>%
  group_modify(~ {
    n <- nrow(.x)
    if (n <= 1000) {
      return(.x)
    } else {
      return(sample_n(.x, 1000))
    }
  })
human_metadata <- as.data.frame(sampled_data)
rownames(human_metadata) <- human_metadata$cells

# human matrix
human_cells <- rownames(human_metadata)
human_fragments_file <- "./05.human_all_fragments_sort.gz"
human_fragments <- CreateFragmentObject(path = human_fragments_file, cells = human_cells)
human_peak_granges <- StringToGRanges(regions = human_peak$human_peak)
human_peak_matrix <- FeatureMatrix(
  fragments = human_fragments,
  features = human_peak_granges
)

# human Seurat object
human_chrom_assay <- CreateChromatinAssay(
  counts = human_peak_matrix,
  sep = c("-", "-")
)
human_pbmc <- CreateSeuratObject(
  counts = human_chrom_assay,
  assay = "peaks",
  meta.data = human_metadata
)

# human motif matrix
human_pbmc <- AddMotifs(
  object = human_pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = JASPAR2024
)
human_pbmc <- RunChromVAR(
  object = human_pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38
)


## co-embedding pig and human objects
DefaultAssay(sn_pig_pbmc) <- 'chromvar'
DefaultAssay(human_pbmc) <- 'chromvar'

human_pbmc@meta.data$celltype <- human_pbmc@meta.data$cell.type
sn_pig_pbmc@meta.data$celltype <- sn_pig_pbmc@meta.data$CellType

Idents(object = human_pbmc) <- "celltype"
human_pbmc <- subset(human_pbmc, downsample = 1000)


Idents(object = sn_pig_pbmc) <- "celltype"
sn_pig_pbmc <- subset(sn_pig_pbmc, downsample = 1000)


sn_pig_pbmc@meta.data$set <- "pig"
human_pbmc@meta.data$set <- "human"

sn_pig_pbmc@assays$peaks <- NULL
human_pbmc@assays$peaks <- NULL

# get PCA based on motif z-score
sn_pig_pbmc.feature <- rownames(sn_pig_pbmc[["chromvar"]])
sn_pig_pbmc <- ScaleData(sn_pig_pbmc, features = sn_pig_pbmc.feature)
sn_pig_pbmc <- RunPCA(sn_pig_pbmc, features = sn_pig_pbmc.feature)
sn_pig_pbmc <- RunUMAP(object = sn_pig_pbmc, reduction = 'pca', dims = 1:30)
#ref.se <- FindNeighbors(ref.se, dims = 1:30)

human_pbmc.feature <- rownames(human_pbmc[["chromvar"]])
human_pbmc <- ScaleData(human_pbmc, features = human_pbmc.feature)
human_pbmc <- RunPCA(human_pbmc, features = human_pbmc.feature)
human_pbmc <- RunUMAP(object = human_pbmc, reduction = 'pca', dims = 1:30)
#quy.se <- FindNeighbors(quy.se, dims = 1:30)

se.lst <- list(ref=sn_pig_pbmc, query=human_pbmc)


# find anchors
se.features <- intersect(rownames(sn_pig_pbmc[["chromvar"]]), rownames(human_pbmc[["chromvar"]]))
se.anchors <- FindIntegrationAnchors(object.list = se.lst, anchor.features = se.features)
se.combined <- IntegrateData(anchorset = se.anchors)

# Perform an integrated analysis
DefaultAssay(se.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
se.combined <- ScaleData(se.combined, verbose = FALSE)
se.combined <- RunPCA(se.combined, npcs = 50, verbose = FALSE)
#ElbowPlot(se.combined, ndims = 50)
se.combined <- RunUMAP(se.combined, reduction = "pca", dims = 1:30)
se.combined <- FindNeighbors(se.combined, reduction = "pca", dims = 1:30)
se.combined <- FindClusters(se.combined, resolution = 0.5)
saveRDS(se.combined, "10.pig_human_combined.rds")


## plot UMAPs
library(ggrastr)

# species
integrated_datafrom <- cbind(se.combined@meta.data, se.combined@reductions$umap@cell.embeddings) %>% dplyr::select(UMAP_1, UMAP_2, set)
n <- nrow(integrated_datafrom)
random_indices <- sample(n)
df_shuffled <- integrated_datafrom[random_indices, ]
p1 <- ggplot(data = df_shuffled, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(mapping = aes(color = set), size = 0.5, alpha = 1) +
  theme_bw()+
  scale_color_manual(values = c("#FF7F00", "#B2DF8A")) +
  theme(panel.grid = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 5)))
pdf("11.adult_pig_human_coembedding_species_UMAP.pdf", width = 15, height = 14)
print(p1)
dev.off()

# celltype
integrated_celltype <- cbind(se.combined@meta.data, se.combined@reductions$umap@cell.embeddings) %>% dplyr::select(UMAP_1, UMAP_2, celltype)
mycolors <- c(
  "Endocardial.cell" = "#00CED1",
  "Alveolar.Type.1.cell" = "#2B68CD",
  "Alveolar.Type.2.cell" = "#2F65CD",
  "Goblet.cell" = "#243C8A",
  "Enterocyte" = "#152285",
  "Hepatocyte" = "#000080",
  "Microglia" = "#E57FFF",
  "Satellite.cell" = "#F46BAC",
  "Type.2.myonuclei" = "#D36194",
  "Type.1.myonuclei" = "#B2507D",
  "Cardiomyocyte" = "#8B3A62",
  "Astrocyte" = "#65721D",
  "Oligodendrocyte.progenitor.cell" = "#5F7023",
  "Oligodendrocyte" = "#556B2F",
  "EX" = "#728815",
  "IN" = "#28730D")

p2 <- ggplot(data = integrated_celltype, mapping = aes(x = UMAP_1, y = UMAP_2)) +
  geom_point_rast(mapping = aes(color = celltype), size = 0.5, alpha = 1) +
  theme_bw()+
  scale_color_manual(values = mycolors) +
  theme(panel.grid = element_blank()) + 
  guides(color = guide_legend(override.aes = list(size = 5)))
pdf("12.adult_pig_human_coembedding_celltype_UMAP.pdf", width = 16, height = 14)
print(p2)
dev.off()
