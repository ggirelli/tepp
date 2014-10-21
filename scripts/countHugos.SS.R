#!/usr/bin/env Rscript
#
# ./countHugos.SS.R graphDirectory maxHugos outputFile

library('igraph')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3) stop('./countHugos.R graphDirectory maxHugos outputFile')

cat('> Prepare.\n')
flist <- list.files(args[1])
remove.id <- c()
for (i in 1:length(flist)) {
	g <- read.graph(paste0(args[1], '/', flist[i]), format='graphml')
	if (length(E(g)) == 0) remove.id <- append(remove.id, i)
}
remove.id <- sort(remove.id, decreasing=T)
for (id in remove.id) flist <- flist[-id]
n.sample <- length(flist)

res <- c()
hugos <- c()

cat('> Count HUGOs.\n')
for(f in flist) {
	g <- read.graph(file.path('.', paste0(args[1], '/', f)), format='graphml')
	#if(length(V(g)) != length(unique(V(g)$name))) print(f)
	res <- rbind(res, c(f, length(V(g)), length(unique(V(g)$name))))
	hugos <- append(hugos, unique(V(g)$name))
}

cat('> Writing output.\n')
write.table(res, args[3], quote=F, row.names=F, col.names=F, sep='\t')
#write.table(hugos, 'hugo_list.dat', quote=F, row.names=F, col.names=F, sep='\t')

svg('sample_per_gene.svg')
hist(as.numeric(res[,3])/as.numeric(args[2]), main='Samples per %gene', xlab='%genes', ylab="#samples", xlim=c(0,1), breaks=seq(0, 1, by=0.05))
dev.off()

svg('genes_per_sample.svg')
hist(table(hugos)/n.sample, main='Gene per %sample', xlab='%samples', ylab="#genes", xlim=c(0,1), breaks=seq(0, 1, by=0.05))
dev.off()

cat('~ END ~\n')
