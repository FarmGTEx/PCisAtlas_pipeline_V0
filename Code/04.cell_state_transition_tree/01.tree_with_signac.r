library(ArchR)
addArchRThreads(threads = 22)
options(future.globals.maxSize = 40 * 1024^3)

## step 1: get peaksets for each tissue
archrPig <- readRDS("./archrPig.Tissue_Date_CellType.peaks.rds")

# retain cell types with cells > 10
archrPig@cellColData$Tissue_Date_CellType <- paste(archrPig@cellColData$Tissue, archrPig@cellColData$Date, archrPig@cellColData$CellType, sep = "--")
celltype_counts <- archrPig@cellColData$Tissue_Date_CellType %>% table
good_celltypes <- names(celltype_counts)[which(celltype_counts > 10)]
archrPig <- subset(archrPig@cellColData,Tissue_Date_CellType %in% good_celltypes) %>% rownames %>% archrPig[., ]

# get peaks
tissue = "CerebralCortex"

object <- subset(archrPig@cellColData,GROUP %in% c("Embryo","Head",tissue)) %>% rownames %>% archrPig[.,]
allpeakslist <- c()
for(dates in intersect(c("E14","E18","E22","E29","E35","E43","E55","E93","7d","14d","60d","180d"), unique(object@cellColData$Date))){
    print(paste("##############", dates, "##########"))
    subcelltypes_files <- subset(object@cellColData,Date == dates) %>% .$Tissue_Date_CellType %>% unique %>% paste0(., "-reproduciblePeaks.gr.rds")
    exist_files <- list.files("./PeakCalls") %>% intersect(., subcelltypes_files)
    goodpeaks <- c()
    for(file in exist_files){
        GRpeaks <- readRDS(paste0("./PeakCalls/", file))
        GRpeaks_overlap <- subsetByOverlaps(object@peakSet, GRpeaks)
        if(file == exist_files[1]){
            goodpeaks <- GRpeaks_overlap
        }else{
            goodpeaks <- reduce(c(goodpeaks, GRpeaks_overlap), min.gapwidth = 0)
        }
    }
    allpeakslist <- c(allpeakslist, list(goodpeaks))
}
names(allpeakslist) <- intersect(c("E14","E18","E22","E29","E35","E43","E55","E93","7d","14d","60d","180d"), unique(object@cellColData$Date))
saveRDS(allpeakslist, paste0(tissue, ".allpeaks.stages.rds"))


## step 2: generate co-embedding UMAP
library(future)
plan("multicore", workers = 12)
options(future.globals.maxSize = 50 * 1024 ^ 3)

# add batch group
archrPig@cellColData$BATCH <- paste(archrPig@cellColData$Tissue, archrPig@cellColData$Date)

# run function
plotlist <- c("Sample","Tissue","Date","Tissue2","CellType","CellLineage","BATCH","Fragments","MT_ratio","TSS_enrich","RFiP")
order_dates = c("E14","E18","E22","E29","E35","E43","E55","E93","7d","14d","60d","180d")

object <- subset(archrPig@cellColData,GROUP %in% c("Embryo","Head",tissue)) %>% rownames %>% archrPig[.,]
saveRDS(object@cellColData, paste0(tissue, ".metadata.rds"))
peaklist <- readRDS(paste0("./", tissue, ".allpeaks.stages.rds"))
Subway_map.step_1.UMAP_list.seurat(object, order_date = order_dates, batch_group = "BATCH", saveSeurat = T, plotlist = plotlist, outname = tissue, k_anchor = 5, regionsList = peaklist, useFeatureCutoff = "q0", subset_by = "CellType", subset_count = 500)


## step 3: calculate KNN list
Subway_map.step_2.knn_list(umap_list = paste0("./", tissue, ".adjacent_UMAP_list.rds"),
    metadata = paste0("./", tissue, ".metadata.rds"),
    order_date = order_dates, outname = tissue
)

# 0.2 for cutoff
Subway_map.step_3.knn_res(knn_median_list = paste0(tissue, "_knn_median_list.2.rds"), edge_cutpff = 0.25, outname = tissue)


## step 4: plot trees and heatmaps
Subway_map.step_4.plot_tree(edge = paste0("./",tissue,".edge.txt"), edge_prob = paste0("./",tissue,".edge_prob.txt"),
    celltype_group = NULL, order_date = order_dates, width_para = 16, height_para = 80, outname = tissue, dotsize = 8)
)

Subway_map.step_4.plot_heatmap(edge_all = paste0("./",tissue,".edge_all.rds"),
    order_date = order_dates, margins = c(15, 15), pdf_height = NULL, pdf_width = NULL, orders = T, Key = F, outname = paste0(tissue, ".order")
)
