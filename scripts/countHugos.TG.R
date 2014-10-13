#!/usr/bin/env Rscript
#
# ./countHugos.TG.R graphName outputFile

library('igraph')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2) stop('./countHugos.R graphDirectory outputFile')

cat('> Prepare.\n')
g <- read.graph(file.path('.', args[1]), format='graphml')
print(g)
cat('> Count HUGOs.\n')
res <- c(args[1], length(V(g)), length(unique(V(g)$HUGO)))

cat('> Write output.\n')
write.table(res, args[2], quote=F, row.names=F, col.names=F, sep='\t')

cat('~ END ~\n')
