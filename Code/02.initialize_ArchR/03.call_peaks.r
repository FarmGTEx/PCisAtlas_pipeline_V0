# call peaks and add to object
library(ArchR)
addArchRThreads(threads = 88)
options(future.globals.maxSize = 40 * 1024^3)
library(BSgenome.Sscrofa.UCSC.susScr11)
library(TxDb.Sscrofa.Ensembl.susScr11.V1)
library(org.Spig.eg.db)

# 
archrPig <- readRDS("./archrPig.CT.rds")

# call peaks (CellType + Tissue + Date)
archrPig@cellColData$Tissue_Date_CellType <- paste(archrPig@cellColData$Tissue, archrPig@cellColData$Date, archrPig@cellColData$CellType, sep = "_")
archrPig <- addGroupCoverages(ArchRProj = archrPig,groupBy = "Tissue_Date_CellType",useLabels = TRUE,minCells = 50,maxCells = 500,maxFragments = 30 * 10^6,minReplicates = 2,maxReplicates = 10,sampleRatio = 0.8,kmerLength = 6, force = T)
archrPig <- addReproduciblePeakSet(ArchRProj = archrPig, groupBy = "Tissue_Date_CellType", peakMethod = "Macs2", cutOff = 0.05, reproducibility = "2", peaksPerCell = 1000, maxPeaks = 200000, minCells = 40, excludeChr = c("chrM", "chrY"), genomeSize = 2367939388, shift = -75, extsize = 150, extendSummits = 250, promoterRegion = c(2500, 1000), plot = TRUE, force = T)
archrPig <- addPeakMatrix(archrPig)

# save
saveRDS(archrPig, "archrPig.Tissue_Date_CellType.peaks.rds")
