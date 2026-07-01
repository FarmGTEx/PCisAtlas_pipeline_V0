#!/usr/bin/env Rscript
# _*_ coding: utf-8 _*_
suppressMessages(library('getopt'))
#options(warnings=-1)

command=matrix(c(
  'input', 'i', 1, 'character', 'input file (the bed fragment data) [required]',
  'outputname', 'o', 1, 'character', 'prefix of output file [required]',
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

##main function
# read in data
fragment <- read.table(args$input, stringsAsFactor = F, sep = '\t')
fragment$Fragment_length <- fragment$V3 - fragment$V2
pdf(paste0(args$outputname, ".fragment_len.pdf"), height = 5, width = 5)
hist(fragment$Fragment_length, breaks = 1000, main = args$input, xlab = "Fragment length")
invisible(dev.off())
