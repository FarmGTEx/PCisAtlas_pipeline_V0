library(dplyr)
library(plyr)
library(qgg)
library(data.table)

# all gwas Z-scores
gwas_ss <- fread("./FCR_gwas.txt", stringsAsFactor = F, header = F, sep = "\t") %>% as.data.frame
all_z_scores <- gwas_ss[,2]
names(all_z_scores) <- gwas_ss[,1]

# get marker list
markers <- fread("./overlapped_snp.FCR.txt", stringsAsFactor = F, header = F, sep = "\t") %>% as.data.frame
marker_list <- c()
for(celltype in unique(markers$V2)){
  tmpmarkers <- subset(markers,V2 == celltype) %>% .$V1
  marker_list <- c(marker_list, list(tmpmarkers))
}
names(marker_list) <- unique(markers$V2)

# Set test based on sums
mma <- gsea(stat = all_z_scores, sets = marker_list, method = "hyperg", nperm = 10000, ncores = 88, threshold = 0.05)
mma$CellType <- rownames(mma)

# output
write.table(mma, "FCR.permu.txt", row.names = T, col.names = T, quote = F, sep = "\t")
