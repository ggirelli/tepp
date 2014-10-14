#!/usr/bin/env Rscript
#
# ./genalgMaxDistGeneset.R totalGraph1 totalGraph2 label1 label2 nCores
# 01_Graph_Manager.class.R and extendigraph.R are required.

library('GA')
library('igraph')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 5) stop('./genalgMaxDistGeneset.R totalGraph1 totalGraph2 label1 label2 nCores')

nCores <- as.numeric(args[5])
if(nCores == 1 | is.na(nCores)) nCores <- FALSE

g1 <- read.graph(args[1], format='graphml')
g2 <- read.graph(args[2], format='graphml')

genes1 <- unique(V(g1)$HUGO)
genes2 <- unique(V(g2)$HUGO)
genes.tot <- union(genes1, genes2)

fitness = function(chr) {
	source('01_Graph_Manager.class.R')
	
	cat('\n    * Prepare genelist\n')
	genes <- genes.tot[chr == 1]
	genes.rm <- genes.tot[chr == 0]

	cat('\t- Prepare graph #1\n')
	g1.tmp <- delete.vertices(g1, V(g1)[HUGO %in% genes.rm])
	g1.tmp <- delete.vertices(g1.tmp, V(g1)[degree(g1, V(g1)) == 0])
	cat('\t- Graph #1 has', length(unique(V(g1.tmp)$HUGO)), ' nodes\n')

	cat('\t- Prepare graph #2\n')
	g2.tmp <- delete.vertices(g2, V(g2)[HUGO %in% genes.rm])
	g2.tmp <- delete.vertices(g2.tmp, V(g2)[degree(g2, V(g2)) == 0])
	cat('\t- Graph #2 has', length(unique(V(g2.tmp)$HUGO)), ' nodes\n')

	cat('    * Calculate distance\n')
	d <- GraphManager()$calcHIMDist(g1.tmp, g2.tmp, 1)
	cat('\t- d = ', d, '\n')
	return(d)
}

monitor = function(obj) {
	cat(obj@population, '\n')
	cat(obj@fitness, '\n')
}

ga <- ga(
	type = 'binary',
	fitness = fitness,
	nBits = length(genes.tot),
	min = 20,
	max = 30,
	names = genes.tot,
	maxiter = 100,
	popSize = 100,
	parallel = nCores
)
