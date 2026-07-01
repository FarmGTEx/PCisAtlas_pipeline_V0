###############################################################################
# 01.selection_sweep_enrichment.r
# Selection sweep enrichment analysis for PCisAtlas
# Covers: LOLA enrichment (cCRE modules x selection signatures),
#         Selection prediction heatmaps
#
# Derived from the original Figure5 workflow.
###############################################################################

library(ArchR)
addArchRThreads(threads = 32)
options(future.globals.maxSize = 40 * 1024^3)

set.seed(1234)
source("~/pcisatlas/scripts/rFun.scATAC.r")

# ---- 1. LOLA enrichment analysis (selection signatures in cCRE modules) ----

    ##region prepare selection signatures
"
cd ~/pcisatlas/selection_sweep/00.fst001
cp ~/pcisatlas/selection_sweep/input/ASD.fst_xpehh.01.bed \
    ~/pcisatlas/selection_sweep/input/WhiteD.fst_xpehh.01.bed \
    ~/pcisatlas/selection_sweep/input/new/body_size.chinese_white.fst.01.bed \
    ~/pcisatlas/selection_sweep/input/new/ASD-White.ASD.bed \
    ~/pcisatlas/selection_sweep/input/new/ASD-White.White.bed \
    .
"

"
# get overlapped genes
for i in *bed;
do
    sed 's/^chr//' ~/pcisatlas/reference/genome/Sus_scrofa.Sscrofa11.1.gene.bed |bedtools intersect -a $i -b - -wo |cut -f 5-9 > overlap_genes/$i
done
"
    ##endregion

    ##region get module peak BED files
"
cd ~/pcisatlas/selection_sweep/00.peaks
cp ~/pcisatlas/selection_sweep/input/allpeaks.bed .
cp ~/pcisatlas/selection_sweep/input/modules.bed .
"
    ##endregion

    ##region LOLA database creation

        ##region prepare data
"
cd ~/pcisatlas/selection_sweep/LOLA/DB/regions
for i in ~/pcisatlas/selection_sweep/LOLA/DB/regions/M*.peaks.bed;
do
    NAME=`basename $i .peaks.bed`;
    cut -f 1-3 $i |sed 's/^chr//' > $NAME.bed
done

cd ~/pcisatlas/selection_sweep/LOLA/00.data
sed 's/^chr//' ~/pcisatlas/data/allpeaks.bed > allpeaks.bed
"
        ##endregion

        ##region create LOLA DB files
        setwd("~/pcisatlas/selection_sweep/LOLA/DB")

        db_path <- path.expand("~/pcisatlas/selection_sweep/LOLA/DB")

        # colletion.txt
        col_df <- data.frame(collector = "User", date = Sys.Date(), source = "ATAC-seq", description = "Pig CellType Specific Peaks")
        write.table(col_df, file = file.path(db_path, "collection.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

        # index.txt
        bed_files <- list.files(file.path(db_path, "regions"), pattern = "\\.bed$")
        index_df <- data.frame(filename = bed_files, description = gsub("\\.bed$", "", bed_files))
        write.table(index_df, file = file.path(db_path, "index.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
        ##endregion

    ##endregion

    ##region run LOLA enrichment
    setwd("~/pcisatlas/selection_sweep/LOLA/01.run")
    library(LOLA)
    library(GenomicRanges)

    # load regionDB (指向创建的顶层文件夹)
    regionDB = loadRegionDB("~/pcisatlas/selection_sweep/LOLA/DB")

    # read query BEDs (selection signals)
    asia_query = readBed("~/pcisatlas/selection_sweep/input/ASD.fst_xpehh.01.bed")
    euro_query = readBed("~/pcisatlas/selection_sweep/input/WhiteD.fst_xpehh.01.bed")

    # read universe (all peaks)
    userUniverse = readBed("~/pcisatlas/selection_sweep/LOLA/DB/allpeaks.bed")

    # run LOLA
    res_asia <- runLOLA(userSets = asia_query, userUniverse = userUniverse,
                        regionDB = regionDB, redefineUserSets = TRUE,
                        minOverlap = 1, cores = 8, direction = "enrichment")
    res_euro <- runLOLA(userSets = euro_query, userUniverse = userUniverse,
                        regionDB = regionDB, redefineUserSets = TRUE,
                        minOverlap = 1, cores = 8, direction = "enrichment")

    write.table(as.data.frame(res_asia), "ASD.txt", row.name = F, col.name = T, sep = "\t", quote = F)
    write.table(as.data.frame(res_euro), "WhiteD.txt", row.name = F, col.name = T, sep = "\t", quote = F)
    ##endregion



# NOTE: qgg Tsum GWAS enrichment analysis moved to 02.GWAS_enrichment.r

# ---- 3. Selection prediction heatmaps (brief) ----
# NOTE: requires pre-computed prediction scores from external tool

    ##region prepare data (shell)
"
cd ~/pcisatlas/selection_sweep/prediction/01.prepare_data

# snp
awk -v FS=\"\\t|_\" -v OFS=\"\\t\" '{print \$1,\$2-1,\$2,\$1\"_\"\$2,\$3,\$4,\$5,\$6}' \
  ../00.data/ASW_64_ASD_1276.deltaf_0.5 > ASD.deltaAF.bed
awk -v FS=\"\\t|_\" -v OFS=\"\\t\" '{print \$1,\$2-1,\$2,\$1\"_\"\$2,\$3,\$4,\$5,\$6}' \
  ../00.data/EUW_66_White_798.deltaf_0.5 > EUD.deltaAF.bed

# selection
sed '1d' ~/pcisatlas/selection_sweep/input/ASD_sel.bed \
  |awk -v OFS=\"\\t\" '{if(\$7 <= 0.01 && \$8 == \"xpehh\"){print \$1,\$2-1,\$3}}' \
  |bedtools sort -i - > ASD.fst_xpehh.01.bed
sed '1d' ~/pcisatlas/selection_sweep/input/White_sel.bed \
  |awk -v OFS=\"\\t\" '{if(\$7 <= 0.01 && \$8 == \"xpehh\"){print \$1,\$2-1,\$3}}' \
  |bedtools sort -i - > WhiteD.fst_xpehh.01.bed

sed 's/^>chr/>/' ~/pcisatlas/reference/genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa > Sus_scrofa.Sscrofa11.1.dna.toplevel.fa

# extract regions and SNPs
bedtools intersect -a ASD.deltaAF.bed -b ASD.fst_xpehh.01.bed -wo |cut -f 1-8 |sort -u > 01.ASD.deltaAF.bed
bedtools intersect -a ASD.deltaAF.bed -b ASD.fst_xpehh.01.bed -wo |cut -f 9-11 |sort -u > 01.ASD.fst_xpehh.01.bed
bedtools intersect -a EUD.deltaAF.bed -b WhiteD.fst_xpehh.01.bed -wo |cut -f 1-8 |sort -u > 01.EUD.deltaAF.bed
bedtools intersect -a EUD.deltaAF.bed -b WhiteD.fst_xpehh.01.bed -wo |cut -f 9-11 |sort -u > 01.WhiteD.fst_xpehh.01.bed
"
    ##endregion

    ##region plot prediction heatmaps
    setwd("~/pcisatlas/selection_sweep/prediction/03.plot")
    library(plyr)
    library(dplyr)
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(reshape2)
    library(pheatmap)
    library(RColorBrewer)
    library(gridExtra)

    CTinfo <- read.table("~/pcisatlas/data/CellType.info",
                         stringsAsFactor = F, sep = "\t", header = T)
    CT.anno <- data.frame(row.names = CTinfo$Name.2, CellLineage = CTinfo$summary.2)
    ann_colors <- list(CellLineage = c(`Germ layer cell` = "#CD853F", `Neural cell` = "#2E8B57",
                                       `Endothelial cell` = "#48D1CC", `Epithelial cell` = "#377eb8",
                                       `Immune cell` = "#984ea3", `Stromal cell` = "#FFA500",
                                       `Muscle cell` = "#f781bf", `Endocrine cell` = "#DC143C",
                                       `Germline cell` = "#9C9C9C"))



        ##region combined: ASD vs EUD with shared row order
        res1 <- read.table("~/pcisatlas/selection_sweep/prediction/data/ASD.tsv",
                           header = T, stringsAsFactor = F, sep = "\t", row.names = 1)
        res2 <- read.table("~/pcisatlas/selection_sweep/prediction/data/EUD.tsv",
                           header = T, stringsAsFactor = F, sep = "\t", row.names = 1)

        # cluster by ASD
        p1 <- pheatmap(t(res1),
            color = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)),
            breaks = c(seq(-2, 0, length.out = 50), seq(0.0001, 2, length.out = 50)),
            fontsize_row = 8, fontsize_col = 8,
            cluster_rows = T, cluster_cols = T,
            annotation_row = CT.anno, annotation_colors = ann_colors,
            scale = "column", border_color = "grey60", angle_col = "45",
            silent = TRUE)
        row_order <- p1$tree_row$order
        t_res2_reordered <- t(res2)[row_order, ]
        p2 <- pheatmap(t_res2_reordered,
            color = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100)),
            breaks = c(seq(-2, 0, length.out = 50), seq(0.0001, 2, length.out = 50)),
            fontsize_row = 8, fontsize_col = 8,
            cluster_rows = F, cluster_cols = T,
            annotation_row = CT.anno, annotation_colors = ann_colors,
            scale = "column", border_color = "grey60", angle_col = "45",
            silent = TRUE)

        pdf("Comparison_Heatmaps.ASD_order.pdf", width = 30, height = 20)
        grid.arrange(p1$gtable, p2$gtable, ncol = 2)
        dev.off()

