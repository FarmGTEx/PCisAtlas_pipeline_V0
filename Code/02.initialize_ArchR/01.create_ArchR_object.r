library(ArchR)
addArchRThreads(threads = 88)
options(future.globals.maxSize = 40 * 1024^3)

library(BSgenome.Sscrofa.UCSC.susScr11)
library(TxDb.Sscrofa.Ensembl.susScr11.V1)
library(org.Spig.eg.db)

# annotation
genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Sscrofa.UCSC.susScr11)
geneAnnotation <- createGeneAnnotation(TxDb = TxDb.Sscrofa.Ensembl.susScr11.V1, OrgDb = org.Spig.eg.db)

# get arrow files
ArrowFiles <- read.table('arrowfiles', stringsAsFactor = F)
# check ArrowFiles
for(i in ArrowFiles[,1]){
  if(!file.exists(i))
    print(i)
}
# create archr object
archrPig <- ArchRProject(ArrowFiles = ArrowFiles[,1], geneAnnotation = geneAnnotation, genomeAnnotation = genomeAnnotation, outputDirectory = "./", copyArrows = TRUE)

# save object
saveArchRProject(ArchRProj = archrPig, outputDirectory = "./", load = T)
