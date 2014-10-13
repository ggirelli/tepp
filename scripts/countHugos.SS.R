#!/usr/bin/env Rscript
#
# ./countHugos.SS.R graphDirectory outputFile

library('igraph')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2) stop('./countHugos.R graphDirectory outputFile')

cat('> Prepare.\n')
flist <- list.files(args[1])
remove.id <- c()
for (i in 1:length(flist)) {
	g <- read.graph(paste0(args[1], '/', flist[i]), format='graphml')
	if (length(E(g)) == 0) remove.id <- append(remove.id, i)
}
remove.id <- sort(remove.id, decreasing=T)
for (id in remove.id) flist <- flist[-id]

res <- c()

cat('> Count HUGOs.\n')
for(f in flist) {
	g <- read.graph(file.path('.', paste0(args[1], '/', f)), format='graphml')
	#if(length(V(g)) != length(unique(V(g)$name))) print(f)
	res <- rbind(res, c(f, length(V(g)), length(unique(V(g)$name))))
}

cat('> Writing output.\n')
write.table(res, args[2], quote=F, row.names=F, col.names=F, sep='\t')

svg('genes_per_sample.svg')
hist(as.numeric(res[,3]), main='Gene number across samples', xlab='Number of genes')
dev.off()

cat('~ END ~\n')
