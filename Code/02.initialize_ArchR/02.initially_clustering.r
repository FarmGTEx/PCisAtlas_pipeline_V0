library(ArchR)
addArchRThreads(threads = 88)
options(future.globals.maxSize = 40 * 1024^3)

archrPig <- readRDS("./Save-ArchR-Project.rds")

## cluster
archrPig <- addIterativeLSI(ArchRProj = archrPig, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2, clusterParams = list(resolution = c(1, 2), sampleCells = 100000, maxClusters = NULL, n.start = 10), varFeatures = 500000, sampleCellsPre = 100000, dimsToUse = 1:30, sampleCellsFinal = 5e5)

# remove batch effect
archrPig@cellColData$Date_Tissue <- paste(archrPig@cellColData$Date, archrPig@cellColData$Tissue, sep = '_')
archrPig <- addHarmony(ArchRProj = archrPig, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Date_Tissue", force = TRUE)
archrPig <- addUMAP(ArchRProj = archrPig, reducedDims = "Harmony",  name = "HarmonyUMAP",  nNeighbors = 30, spread = 3, minDist = 0.3,  metric = "cosine", force = TRUE)
saveRDS(archrPig, 'archrPig.ArchR_object.rds')

# plot
pdf("cluster.pdf", width = 40, height = 10)
for(plotby in colnames(archrPig@cellColData)){
	p <- plotEmbedding(ArchRProj = archrPig, colorBy = "cellColData", name = plotby, embedding = "HarmonyUMAP", colorTitle = plotby) +
	  theme(legend.position = "right", legend.direction = 'vertical', legend.text = element_text(size=8)) +
	  guides(color = guide_legend(override.aes = list(size = 4)))
	print(p)
}
dev.off()
