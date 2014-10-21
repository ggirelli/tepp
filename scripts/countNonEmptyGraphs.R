#!/usr/bin/env Rscript
#
# ./countNonEmptyGraphs.R graphDirectory

library('igraph')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1) stop('./countNonEmptyGraphs.R graphDirectory')

fs <- list.files(file.path('.', args[1]))
fs <- paste0(args[1], '/', fs)
rm.ids <- c()
for(i in 1:length(fs)) {
	f <- fs[i]
	g <- read.graph(f, format='graphml')
	if(ecount(g) == 0) rm.ids <- append(rm.ids, i)
}

cat(length(fs) - length(rm.ids), ' / ', length(fs), '\n')