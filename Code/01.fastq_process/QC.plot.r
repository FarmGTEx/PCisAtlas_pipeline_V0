#!/usr/bin/env Rscript
# _*_ coding: utf-8 _*_
suppressMessages(library('getopt'))
#options(warnings=-1)

command=matrix(c(
  'input', 'i', 1, 'character', 'input file (the cell metadata) [required]',
  'outputname', 'o', 1, 'character', 'prefix of output file [required]',
  'fragments', 'f', 2, 'numeric', 'number of fragments in cells to retain, default "1000"',
  'ratio', 'r', 2, 'numeric', 'peak ratio in cells to retain, default "0.2"',
  'tss', 't', 2, 'numeric', 'TSS enrichment score, default "2"',
  'nsl', 'n', 2, 'numeric', 'nucleosome signal, default "4"',
  'filter', 'w', 0, 'loical', 'whether output filtered metadata',
  'help', 'h', 0,'loical', 'help'),
  byrow=T, ncol=5
)
args=getopt(command)

#help info
if (!is.null(args$help) || is.null(args$input) || is.null(args$outputname)) {
  cat("\n")
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

#set default values
if ( is.null(args$fragments ) ) { args$fragments = 1000 }
if ( is.null(args$ratio ) ) { args$ratio = 0.2 }
if ( is.null(args$tss ) ) { args$tss = 2 }
if ( is.null(args$nsl ) ) { args$nsl = 4 }

##main function
# load packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(ggpmisc)
})

# read in data
metadata <- read.table(args$input, stringsAsFactor = F, sep = '\t')
metadata[is.na(metadata)] <- 0
colnames(metadata) <- c('CellBarcode', 'Fragments', 'FragLen.1-147', 'FragLen.147-294', 'Nucleosome_Signal', 'MT_ratio', 'TSS_enrich', 'FragInPeak', 'RFiP', 'Doublet')[1:ncol(metadata)]
metadata$group = 1

# plot violin
ps <- list()
for(plotgroup in c('Fragments', 'Nucleosome_Signal', 'TSS_enrich', 'RFiP')){
	p <- ggplot(metadata, aes_string(x="group", y=plotgroup)) +
	  geom_violin(fill='#fb8072') +
	  geom_jitter(shape=16, position=position_jitter(0.4), size = 0.2) +
	  xlab(NULL) +
	  theme_bw() +
	  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
	ps <- c(ps, list(p))
}
#pdf("zz.pdf",height=5,width=12)
png(paste0(args$outputname, ".barcode_QC.violin.png"), height=2500, width=1800, res=400)
grid.arrange(grobs=ps,nrow=2)
invisible(dev.off())

# plot histogram
outsummary = data.frame(Cells = 1:6)
ps2 <- list()
for(plotgroup in c('Fragments', 'Nucleosome_Signal', 'TSS_enrich', 'RFiP')){
	binwidths = max(metadata[,plotgroup])/1000
	totalcell = nrow(metadata)
	outsummary[1,1] <- totalcell
	mediancell = median(metadata[,plotgroup])
	if(plotgroup == 'Fragments'){
		goodcell = which(metadata[,plotgroup] > args$fragments & metadata[,plotgroup] < 100000) %>% length
		plotlabels = paste0("\nTotal cells: ", totalcell, "\nGood cells: ", goodcell, "\nMedian fragments: ", mediancell)
		outsummary[2,1] <- goodcell
	}else if(plotgroup == 'Nucleosome_Signal'){
		goodcell = which(metadata[,plotgroup] < args$nsl) %>% length
		plotlabels = paste0("\nTotal cells: ", totalcell, "\nGood cells: ", goodcell)
		outsummary[3,1] <- goodcell
	}else if(plotgroup == 'TSS_enrich'){
		goodcell = which(metadata[,plotgroup] > args$tss) %>% length
		plotlabels = paste0("\nTotal cells: ", totalcell, "\nGood cells: ", goodcell)
		outsummary[4,1] <- goodcell
	}else if(plotgroup == 'RFiP'){
		goodcell = which(metadata[,plotgroup] > args$ratio) %>% length
		plotlabels = paste0("\nTotal cells: ", totalcell, "\nGood cells: ", goodcell)
		outsummary[5,1] <- goodcell
	}
	annotations <- data.frame(xpos = Inf, ypos = Inf, annotateText = plotlabels, hjustvar = 1, vjustvar = 1)
	p <- ggplot(metadata) +
	  aes_string(x=plotgroup) +
	  geom_histogram(binwidth=binwidths, fill="black", color="black", alpha=0) +
#	  geom_vline(aes(xintercept=nxline), colour = "red", linetype = "dashed") +
#	  geom_text_npc(aes(npcx = "right", npcy = "top", label = plotlabels)) + 
#	  annotate(geom = 'text', label = plotlabels, x = -Inf, y = Inf, hjust = 0, vjust = 1) +
	  geom_text(data = annotations, aes(x = xpos, y = ypos, label = annotateText, hjust=hjustvar, vjust=vjustvar)) +
#	  geom_text(data = data.frame(), aes(x = Inf, y = Inf, label = plotlabels), hjust="right", vjust="top") +
	  theme_bw()
	if(plotgroup == 'Fragments'){
		p <- p + geom_vline(aes(xintercept=args$fragments), colour = "red", linetype = "dashed") + geom_vline(aes(xintercept=100000), colour = "red", linetype = "dashed")
	}else if(plotgroup == 'Nucleosome_Signal'){
		p <- p + geom_vline(aes(xintercept=args$nsl), colour = "red", linetype = "dashed")
	}else if(plotgroup == 'TSS_enrich'){
		p <- p + geom_vline(aes(xintercept=args$tss), colour = "red", linetype = "dashed")
	}else{
		p <- p + geom_vline(aes(xintercept=args$ratio), colour = "red", linetype = "dashed")
	}
	ps2 <- c(ps2, list(p))
}
#pdf("zz.pdf",height=6,width=24)
png(paste0(args$outputname, ".barcode_QC.hist.png"), height=2500, width=2500, res=400)
grid.arrange(grobs=ps2,nrow=2)
invisible(dev.off())

# out QC summary
outsummary[6,1] = which(metadata[,'Fragments'] > args$fragments & metadata[,'Fragments'] < 100000 & metadata[,'Nucleosome_Signal'] < args$nsl & metadata[,'TSS_enrich'] > args$tss & metadata[,'RFiP'] > args$ratio) %>% length
rownames(outsummary) <- c('Total_cells', 'Fragments', 'Nucleosome_Signal', 'TSS_enrich', 'RFiP', 'Final_Cells')
write.table(outsummary, paste0(args$outputname, ".QC.summary"), col.names = F, quote = F, sep = '\t')

# plot filtered cells
metadata <- metadata[which(metadata[,'Fragments'] > args$fragments & metadata[,'Fragments'] < 100000 & metadata[,'Nucleosome_Signal'] < args$nsl & metadata[,'TSS_enrich'] > args$tss & metadata[,'RFiP'] > args$ratio),]
if(!is.null(args$filter)){
	write.table(metadata[,1:9], paste0(args$outputname, ".flted.metadata"), row.names = F, col.names = F, quote = F, sep = '\t')
}

# plot violin
ps <- list()
for(plotgroup in c('Fragments', 'Nucleosome_Signal', 'TSS_enrich', 'RFiP')){
	p <- ggplot(metadata, aes_string(x="group", y=plotgroup)) +
	  geom_violin(fill='#fb8072') +
	  geom_jitter(shape=16, position=position_jitter(0.4), size = 0.2) +
	  xlab(NULL) +
	  theme_bw() +
	  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
	ps <- c(ps, list(p))
}
#pdf("zz.pdf",height=5,width=12)
png(paste0(args$outputname, ".barcode_QC.flted.violin.png"), height=2500, width=1800, res=400)
grid.arrange(grobs=ps,nrow=2)
invisible(dev.off())

# plot histogram
ps2 <- list()
for(plotgroup in c('Fragments', 'Nucleosome_Signal', 'TSS_enrich', 'RFiP')){
	binwidths = max(metadata[,plotgroup])/1000
	p <- ggplot(metadata) +
	  aes_string(x=plotgroup) +
	  geom_histogram(binwidth=binwidths, fill="black", color="black", alpha=0) +
	  theme_bw()
	ps2 <- c(ps2, list(p))
}
#pdf("zz.pdf",height=6,width=24)
png(paste0(args$outputname, ".barcode_QC.flted.hist.png"), height=2500, width=2500, res=400)
grid.arrange(grobs=ps2,nrow=2)
invisible(dev.off())

