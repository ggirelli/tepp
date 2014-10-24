#!/usr/bin/env Rscript
#
# ./genalgMaxDistGeneset.R totalGraph1 totalGraph2 label1 label2 nCores
# Graph_Manager.class.R and extendigraph.R are required.

library('GA')
library('igraph')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 6) stop('./genalgMaxDistGeneset.R totalGraph1 totalGraph2 label1 label2 nCores n.iter')

n.iter <- as.numeric(args[6])

nCores <- as.numeric(args[5])
if(nCores == 1 | is.na(nCores)) nCores <- FALSE

g1 <- read.graph(args[1], format='graphml')
g2 <- read.graph(args[2], format='graphml')

genes1 <- unique(V(g1)$HUGO)
genes2 <- unique(V(g2)$HUGO)
genes.tot <- union(genes1, genes2)

fitness = function(chr) {
	source('Graph_Manager.class.R')
	
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

	cat('\t* Calculate distance\n')
	d <- GraphManager()$calcHIMDist(g1.tmp, g2.tmp, 1)
	cat('\t- d = ', d, '\n')
	return(d)
}

cat('> Running GA binary algorithm with:\n')
cat('\t- nBits:', length(genes.tot), '\n')
cat('\t- min:', 0, '\n')
cat('\t- max:', length(genes.tot), '\n')
cat('\t- maxiter:', n.iter, '\n')
cat('\t- popSize:', 100, '\n')
cat('\t- nCores:', nCores, '\n')

system.time({
	ga <- ga(
		type = 'binary',
		fitness = fitness,
		nBits = length(genes.tot),
		min = 0,
		max = length(genes.tot),
		names = genes.tot,
		maxiter = n.iter,
		popSize = 100,
		parallel = nCores
	)
})

write.table(ga@summary, 'ga_summary.dat', row.names=F, quote=F)

svg('ga_convergence.svg')
plot(1:n.iter, ga@summary[,1], type='l', main='Genetic Algorithm', xlab='iteration', ylab='d(ERG+, ERG-)')
lines(1:n.iter, ga@summary[,2], lty=2, col=2)
lines(1:n.iter, ga@summary[,4], lty=1, col=4)
legend(0, 0.25, c('max', 'min', 'mean'), lty=c(1,1,2), col=c(1,4,2))
dev.off()

save('ga', file='ga.drop.Rdata')