#!/usr/bin/env Rscript
# _*_ coding: utf-8 _*_
suppressMessages(library('getopt'))
#options(warnings=-1)

command=matrix(c(
  'input', 'i', 1, 'character', 'input file (the cell metadata) [required]',
  'summary', 's', 1, 'character', 'summary file [required]',
  'help', 'h', 0,'loical', 'help'),
  byrow=T, ncol=5
)
args=getopt(command)

#help info
if (!is.null(args$help) || is.null(args$input) || is.null(args$summary)) {
  cat("\n'NonDoublet' in 'Doublet' column is used to determine the which barcodes are used!\n")
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

#set default values

##main function
# load packages

# read in data
metadata <- read.table(args$input, stringsAsFactor = F, sep = '\t')
metadata[is.na(metadata)] <- 0
colnames(metadata) <- c('CellBarcode', 'Fragments', 'FragLen.1-147', 'FragLen.147-294', 'Nucleosome_Signal', 'MT_ratio', 'TSS_enrich', 'FragInPeak', 'RFiP', 'Doublet')[1:ncol(metadata)]
metadata <- subset(metadata,Doublet == "NonDoublet")
fragment_mean <- mean(metadata$Fragments)
fragment_median <- median(metadata$Fragments)

commands1 <- paste0("echo ", paste("Average fragments", fragment_mean, sep = " tab "), " >> ", args$summary)
commands2 <- paste0("echo ", paste("Median fragments", fragment_median, sep = " tab "), " >> ", args$summary)
system(commands1)
system(commands2)
