#!/usr/bin/env Rscript
#
# ./cleanSSTEP.R graph.dir output.dir
#
# Merges 'name' and 'abe.type' node attributes into new 'name' node attribute.
# The old 'name' is stored as a new 'HUGO' node attribute.

library('igraph')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2) stop('./cleanSSTEP.R graph.dir output.dir')

flist <- list.files(args[1])

for (file in flist) {
	g <- read.graph(file.path(args[1], file), format='graphml')
	V(g)$HUGO <- V(g)$name
	V(g)$name <- paste0(V(g)$HUGO, '~', V(g)$abe.type)
	write.graph(g, file.path(args[2], file), format='graphml')
}

cat('~ END ~\n')
