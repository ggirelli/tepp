#!/usr/bin/env Rscript
#
# ./preparePermutations.R paramFile
# 
# SSDNs must be in the folder 'sample-graphs'
# Requires the following files in the same folder:
#     * Graph_Builder.class.R
#     * Graph_Builder.launcher.R
#     * parentDir.Graph_Builder.class.R
#     * parentDir.Graph_Builder.launcher.R
#     * Graph_Manager.class.R
#     * parentDir.Graph_Manager.class.R
#     * extendigraph.R

library(igraph)
library(doParallel)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 1) stop('./preparePermutations.R paramFile')

system.time({

# ----------- #
# READ PARAMS #
# ----------- #

getp = function(table, pname) {
	if(ncol(table) != 2 | nrow(table) == 0) return(NA)
	id <- which(table[,1] == pname)
	if(length(id) != 1) return(NA)
	return(as.character(table[id,2]))
}

cat('> Read parameter file.\n')
ps <- read.table(args[1])

nperm <- getp(ps, 'nperm')
if(is.na(nperm)) {
	stop('Error: nperm required.')
} else {
	nperm <- as.numeric(nperm)
	if(!is.numeric(nperm) | nperm < 1) stop('Error: using wrong nperm.')
	cat(paste0('    * Preparing ', nperm, ' permutation(s).', "\n"))
}

seed <- getp(ps, 'seed')
if(is.na(seed)) {
	seed <- sample(100000000,1)
	cat(paste0("    * Using seed ", seed, '.', "\n"))
} else {
	cat(paste0("    * Using seed ", seed, '.', "\n"))
}

nCores <- getp(ps, 'nCores')
if(is.na(nCores)) {
	nCores <- 1
	cat('    * Using 1 core.\n')
} else {
	nCores <- as.numeric(nCores)
	if(!is.numeric(nCores) | nCores < 1) stop('Error: using wrong nCores.')
	cat(paste0('    * Using ', nCores, ' core(s).', "\n"))
}

graph.dir <- getp(ps, 'graph.dir')
if(is.na(graph.dir)) {
	stop('Error: graph.dir required.')
} else {
	cat(paste0('    * Using graphs in \'', graph.dir, '\'.', "\n"))
}

annotations <- getp(ps, 'annotations')
if(is.na(annotations)) {
	stop('Error: annotations required.')
} else {
	cat(paste0('    * Using annotations from \'', annotations, '\'.', "\n"))
}

sample.col <- getp(ps, 'sample.col')
if(is.na(sample.col)) {
	stop('Error: sample.col required.')
} else {
	cat(paste0('    * Samples from \'', sample.col, '\' column.', "\n"))
	if(sample.col == 'sample') {
		pfcl$Bam.name <- pfcl$sample
		sample.col <- 'Bam.name'
	}
}

gainTab <- getp(ps, 'file-Gain')
lossTab <- getp(ps, 'file-Loss')
pmTab <- getp(ps, 'file-PM')
if(!is.na(pmTab)) {
	if(xor(is.na(gainTab), is.na(lossTab))) {
		if(is.na(gainTab)) {
			cat('    @ Discarding lossTab.\n')
			lossTab <- NA
		}
		if(is.na(lossTab)) {
			cat('    @ Discarding gainTab.\n')
			gainTab <- NA
		}
	}
} else {
	if(xor(is.na(gainTab), is.na(lossTab))) {
		if(is.na(gainTab)) {
			stop('Error: discarding lossTab. No data, abort.')
			lossTab <- NA
		}
		if(is.na(lossTab)) {
			stop('Error: discarding gainTab. No data, abort.')
			gainTab <- NA
		}
	}
}
if(!is.na(gainTab)) {
	cat(paste0('    * Using gains table \'', gainTab, '\'.', "\n"))
}
if(!is.na(lossTab)) {
	cat(paste0('    * Using loss table \'', lossTab, '\'.', "\n"))
}
if(!is.na(pmTab)) {
	cat(paste0('    * Using PM table \'', pmTab, '\'.', "\n"))
}

m.outdir <- getp(ps, 'm.outdir')
if(is.na(m.outdir)) {
	cat('    * ERG- results output directory: \'ERGm\'\n')
	m.outdir <- 'ERGm'
} else {
	cat(paste0('    * ERG- results output directory: \'', m.outdir, '\'', "\n"))
}

p.outdir <- getp(ps, 'p.outdir')
if(is.na(p.outdir)) {
	cat('    * ERG+ results output directory: \'ERGp\'\n')
	p.outdir <- 'ERGp'
} else {
	cat(paste0('    * ERG+ results output directory: \'', p.outdir, '\'', "\n"))
}

gene.label <- getp(ps, 'gene.label')
if(is.na(p.outdir)) stop('Error: gene.label required.')

skip.base <- getp(ps, 'skip.base')
if(!is.na(skip.base)) {
	if(skip.base == 'TRUE') {
		cat("    * Skipping base networks.\n")
		skip.base <- TRUE
	} else {
		skip.base <- FALSE
	}
} else {
	skip.base <- FALSE
}

skip.perm <- getp(ps, 'skip.perm')
if(!is.na(skip.perm)) {
	if(skip.perm == 'TRUE') {
		cat("    * Skipping permuting.\n")
		skip.perm <- TRUE
	} else {
		skip.perm <- FALSE
	}
} else {
	skip.perm <- FALSE
}

skip.build <- getp(ps, 'skip.build')
if(!is.na(skip.build)) {
	if(skip.build == 'TRUE') {
		cat("    * Skipping building permutations.\n")
		skip.build <- TRUE
	} else {
		skip.build <- FALSE
	}
} else {
	skip.build <- FALSE
}

# ----------- #
# PREPARATION #
# ----------- #

cat('> Prepare.\n')
cat('    * Get graph list.\n')
# Retrieve graph list
flist <- list.files(file.path('.', paste0(graph.dir, '/')))
remove.id <- c()
for (i in 1:length(flist)) {
	g <- read.graph(paste0('sample-graphs/', flist[i]), format='graphml')
	if (length(E(g)) == 0) remove.id <- append(remove.id, i)
}
remove.id <- sort(remove.id, decreasing=T)
for (id in remove.id) flist <- flist[-id]

cat('    * Prepare pathology features vector.\n')
# Prepare pathology features vector
pf <- read.table(file.path('.', annotations), header=T)
pfcl <- pf[which(paste0('gra_', eval(parse(text=paste0('pf$', sample.col))), '.graphml') %in% flist),]
pfcl$sample <- eval(parse(text=paste0('pfcl$', sample.col)))
eval(parse(text=paste0('pfcl$', sample.col, ' <- paste0(\'gra_\', eval(parse(text=paste0(\'pfcl$\', sample.col))), \'.graphml\')')))
pfcl$Gleason.score <- pfcl$Gleason_Major + pfcl$Gleason_Minor
pfcl$Gleason.clean <- paste0(pfcl$Gleason_Major,'+',pfcl$Gleason_Minor)
pfcl$ERG.clean <- pfcl$TMPRSS2.ERG_Fusion_Status_FISH
for (i in which(is.na(pfcl$ERG.clean))) {
    pfcl$ERG.clean[i] <- pfcl$ETS.fusion.detected.by.sequencing[i]
}

cat('    * Distinguish classes\n')
# Distinguish classes
erg.m <- eval(parse(text=paste0('pfcl$', sample.col)))[which(pfcl$ERG.clean == 0)]
cat(paste0('> Found ', length(erg.m), ' ERG- samples.', "\n"))
erg.p <- eval(parse(text=paste0('pfcl$', sample.col)))[which(pfcl$ERG.clean > 0)]
cat(paste0('> Found ', length(erg.p), ' ERG+ samples.', "\n"))

if(!skip.base) {
	dir.create(file.path(p.outdir), showWarnings = FALSE)
	dir.create(file.path(m.outdir), showWarnings = FALSE)

	gainTab.out <- NA
	lossTab.out <- NA
	pmTab.out <- NA

	if(!is.na(gainTab)) {
		gainTab.out <- 'erg_clonTab.Gain.txt'
		lossTab.out <- 'erg_clonTab.Loss.txt'

		# Retrieve gain/loss data
		gain <- read.table(file.path('.', gainTab), header=T, sep='\t')
		loss <- read.table(file.path('.', lossTab), header=T, sep='\t')

		# Dividing classes in gain/loss
		cat('> Dividing gain/loss in classes\n')
		erg.p.gain <- gain[which(gain$sample %in% pfcl$sample[which(pfcl$ERG.clean > 0)]),]
		erg.p.loss <- loss[which(loss$sample %in% pfcl$sample[which(pfcl$ERG.clean > 0)]),]
		write.table(erg.p.gain, paste0(p.outdir, '/erg_clonTab.Gain.txt'), quote=F, sep='\t')
		write.table(erg.p.loss, paste0(p.outdir, '/erg_clonTab.Loss.txt'), quote=F, sep='\t')

		erg.m.gain <- gain[which(gain$sample %in% pfcl$sample[which(pfcl$ERG.clean == 0)]),]
		erg.m.loss <- loss[which(loss$sample %in% pfcl$sample[which(pfcl$ERG.clean == 0)]),]
		write.table(erg.m.gain, paste0(m.outdir, '/erg_clonTab.Gain.txt'), quote=F, sep='\t')
		write.table(erg.m.loss, paste0(m.outdir, '/erg_clonTab.Loss.txt'), quote=F, sep='\t')
	}
	if(!is.na(pmTab)) {
		pmTab.out <- 'erg_clonTab.pm.txt'

		# Retrieve gain/loss data
		pm <- read.table(file.path('.', pmTab), header=T, sep='\t')

		# Dividing classes in gain/loss
		cat('> Dividing pm in classes\n')
		erg.p.pm <- pm[which(pm$sample %in% pfcl$sample[which(pfcl$ERG.clean > 0)]),]
		erg.m.pm <- pm[which(pm$sample %in% pfcl$sample[which(pfcl$ERG.clean > 0)]),]
		write.table(erg.p.pm, paste0(p.outdir, '/erg_clonTab.pm.txt'), quote=F, sep='\t')
		write.table(erg.m.pm, paste0(m.outdir, '/erg_clonTab.pm.txt'), quote=F, sep='\t')
	}

	if(!is.na(gainTab.out)) {
		gainTab.out.m <- paste0(m.outdir, '/', gainTab.out)
	} else {
		gainTab.out.m <- ''
	}
	if(!is.na(lossTab.out)) {
		lossTab.out.m <- paste0(m.outdir, '/', lossTab.out)
	} else {
		lossTab.out.m <- ''
	}
	if(!is.na(gainTab.out)) {
		gainTab.out.p <- paste0(p.outdir, '/', gainTab.out)
	} else {
		gainTab.out.p <- ''
	}
	if(!is.na(lossTab.out)) {
		lossTab.out.p <- paste0(p.outdir, '/', lossTab.out)
	} else {
		lossTab.out.p <- ''
	}
	if(!is.na(pmTab.out)) {
		pmTab.out.m <- paste0(m.outdir, '/', pmTab.out)
	} else {
		pmTab.out.m <- ''
	}
	if(!is.na(pmTab.out)) {
		pmTab.out.p <- paste0(p.outdir, '/', pmTab.out)
	} else {
		pmTab.out.p <- ''
	}
	param.names <- c('clusters', 'verbose', 'genes.label', 'output.dir', 'file-Gain', 'file-Loss', 'file-PM')
	m.param.val <- c(nCores, 'TRUE', gene.label, m.outdir, lossTab.out.m, gainTab.out.m, pmTab.out.m)
	p.param.val <- c(nCores, 'TRUE', gene.label, p.outdir, lossTab.out.p, gainTab.out.p, pmTab.out.p)

	write.table(cbind(param.names, m.param.val), paste0('param.ergm.txt'), quote=F, col.names=F, row.names=F, sep=' ')
	write.table(cbind(param.names, p.param.val), paste0('param.ergp.txt'), quote=F, col.names=F, row.names=F, sep=' ')
}

if(!skip.perm) {
	# MAKE PERMUTATIONS
	cat('> Permuting.\n')

	set.seed(seed)

	cores <- makeCluster(nCores)
	registerDoParallel(cores)

	perms <- c()
	res <- foreach (i=1:nperm) %dopar% {
	    library(igraph)
	    samples <- sample(1:sum(c(length(erg.m),length(erg.p))))
	    sample.p <- samples[1:length(erg.p)]
	    sample.m <- samples[(length(erg.p)+1):sum(c(length(erg.m),length(erg.p)))]
	    perms <- append(perms, paste(samples, collapse='_'))

	    dir.create(file.path('.', paste0('s', i)), showWarnings = FALSE)
	    dir.create(file.path(paste0('s', i), p.outdir), showWarnings = FALSE)
	    dir.create(file.path(paste0('s', i), m.outdir), showWarnings = FALSE)

		gainTab.out <- NA
		lossTab.out <- NA
		pmTab.out <- NA

	    if(!is.na(gainTab)) {
	    	gainTab.out <- 'erg_clonTab.Gain.txt'
			lossTab.out <- 'erg_clonTab.Loss.txt'

		    erg.p.gain <- gain[which(gain$sample %in% pfcl$sample[sample.p]),]
		    erg.p.loss <- loss[which(loss$sample %in% pfcl$sample[sample.p]),]
		    write.table(erg.p.gain, paste0('s', i, '/', p.outdir , '/', gainTab.out), quote=F, sep='\t')
		    write.table(erg.p.loss, paste0('s', i, '/', p.outdir , '/', lossTab.out), quote=F, sep='\t')
		    erg.m.gain <- gain[which(gain$sample %in% pfcl$sample[sample.m]),]
		    erg.m.loss <- loss[which(loss$sample %in% pfcl$sample[sample.m]),]
		    write.table(erg.m.gain, paste0('s', i,'/', m.outdir, '/', gainTab.out), quote=F, sep='\t')
		    write.table(erg.m.loss, paste0('s', i,'/', m.outdir, '/', lossTab.out), quote=F, sep='\t')
		}
	    if(!is.na(pmTab)) {
	    	pmTab.out <- 'erg_clonTab.pm.txt'

		    erg.p.pm <- pm[which(pm$sample %in% pfcl$sample[sample.p]),]
		    erg.m.pm <- pm[which(pm$sample %in% pfcl$sample[sample.m]),]
		    write.table(erg.p.pm, paste0('s', i, '/', p.outdir , '/', pmTab.out), quote=F, sep='\t')
		    write.table(erg.m.pm, paste0('s', i, '/', m.outdir , '/', pmTab.out), quote=F, sep='\t')
		}

		if(!is.na(gainTab.out)) {
			gainTab.out.m <- paste0(m.outdir, '/', gainTab.out)
		} else {
			gainTab.out.m <- ''
		}
		if(!is.na(lossTab.out)) {
			lossTab.out.m <- paste0(m.outdir, '/', lossTab.out)
		} else {
			lossTab.out.m <- ''
		}
		if(!is.na(gainTab.out)) {
			gainTab.out.p <- paste0(p.outdir, '/', gainTab.out)
		} else {
			gainTab.out.p <- ''
		}
		if(!is.na(lossTab.out)) {
			lossTab.out.p <- paste0(p.outdir, '/', lossTab.out)
		} else {
			lossTab.out.p <- ''
		}
		if(!is.na(pmTab.out)) {
			pmTab.out.m <- paste0(m.outdir, '/', pmTab.out)
		} else {
			pmTab.out.m <- ''
		}
		if(!is.na(pmTab.out)) {
			pmTab.out.p <- paste0(p.outdir, '/', pmTab.out)
		} else {
			pmTab.out.p <- ''
		}
	    param.names <- c('clusters', 'verbose', 'genes.label', 'output.dir', 'file-Gain', 'file-Loss', 'file-PM')
	    m.param.val <- c(nCores, 'TRUE', gene.label, m.outdir, lossTab.out.m, gainTab.out.m, pmTab.out.m)
	    p.param.val <- c(nCores, 'TRUE', gene.label, p.outdir, lossTab.out.p, gainTab.out.p, pmTab.out.p)

	    write(perms, paste0('s', i, '/perms.dat'))
	    write.table(cbind(param.names, m.param.val), paste0('s', i, '/paramS', i , 'ergm.txt'), quote=F, col.names=F, row.names=F, sep=' ')
	    write.table(cbind(param.names, p.param.val), paste0('s', i, '/paramS', i , 'ergp.txt'), quote=F, col.names=F, row.names=F, sep=' ')
	}
	cat('> Saving permutations log.\n')
	system('cat s*/perms.dat > perms.dat')

	stopCluster(cores)
}

# -------- #
# BUILDING #
# -------- #
cat('> Building dependency networks.\n')

if(!skip.base) {
	cat('	* Building original ERG- dependency graph\n')
	system(paste0('./Graph_Builder.launcher.R -p=param.ergm.txt > log.ergm.dat'))
	cat('	* Building original ERG+ dependency graph\n')
	system(paste0('./Graph_Builder.launcher.R -p=param.ergp.txt > log.ergp.dat'))
}

if(!skip.build) {
	for (i in 1:nperm) {
		system.time({
			cat(paste0('    * Building permutation #', i, "\n"))
			setwd(paste0('s', i))

			cat(paste0("         - Building ERGm dependency graph.\n"))
			system(paste0('../parentDir.Graph_Builder.launcher.R -p=paramS', i, 'ergm.txt > log.s', i, '.ergm.dat'))
			cat(paste0("         - Building ERGp dependency graph.\n"))
			system(paste0('../parentDir.Graph_Builder.launcher.R -p=paramS', i, 'ergp.txt > log.s', i, '.ergp.dat'))

			setwd('..')
		})
	}
}

# ------------------ #
# CALCULATE DISTANCE #
# ------------------ #

cat(paste0('> Calculating distances\n'))

cat('	* Working on original data.\n')
source('./Graph_Manager.class.R')
gm <- read.graph(paste0(p.outdir, '/total_graph.graphml'), format='graphml')
gp <- read.graph(paste0(m.outdir, '/total_graph.graphml'), format='graphml')
ds <- GraphManager()$calcDistances(gm, gp, 1)
write(ds, 'distances.dat')

cat('	* Working on permutations.\n')
cores <- makeCluster(nCores)
registerDoParallel(cores)

res <- foreach(i=1:nperm, .combine='rbind') %dopar% {
	setwd(paste0('s', i))

	source('../parentDir.Graph_Manager.class.R')
	gm <- read.graph(paste0(p.outdir, '/total_graph.graphml'), format='graphml')
	gp <- read.graph(paste0(m.outdir, '/total_graph.graphml'), format='graphml')
	ds <- GraphManager()$calcDistances(gm, gp, 1)
	write.table(ds, 'distances.dat', row.names=F, col.names=F, quote=F, sep='\t')

	setwd('..')

	return(ds)
}

write.table(res, 'dist.perm.dat', row.names=F, col.names=F, quote=F, sep='\t')

stopCluster(cores)

# ------------- #
# PREPARE PLOTS #
# ------------- #

cat(paste0('> Plotting\n'))

t <- res
ori <- ds

hm <- t[,1]
jm <- t[,2]
im <- t[,3]
him <- t[,4]
jim <- t[,5]

if(!file.exists('plots/')) dir.create('plots')

svg(paste0('plots/hamming_hist.', args[3], '.svg'))
hist(hm, xlim=c(0,1), main='Distance distribution', xlab='H', prob=T)
lines(density(hm), col=4, lty=2)
abline(v=ori[1], col=2, lty=3, lwd=1)
dev.off()
svg(paste0('plots/hamming_hist_sr.', args[3], '.svg'))
hist(hm, xlim=c(0,1), breaks=sqrt(length(hm)), main='Distance distribution [SQRT]', xlab='H', prob=T)
lines(density(hm), col=4, lty=2)
abline(v=ori[1], col=2, lty=3, lwd=1)
dev.off()
if(IQR(hm) != 0) {
	svg(paste0('plots/hamming_hist_fd.', args[3], '.svg'))
	hist(hm, xlim=c(0,1), breaks=1/(2*IQR(hm)/length(hm)**(1/3)), main='Distance distribution [FD]', xlab='H', prob=T)
	lines(density(hm), col=4, lty=2)
	abline(v=ori[1], col=2, lty=3, lwd=1)
	dev.off()
}

svg(paste0('plots/jaccard_hist.', args[3], '.svg'))
hist(jm, xlim=c(0,1), main='Distance distribution', xlab='J', prob=T)
lines(density(jm), col=4, lty=2)
abline(v=ori[2], col=2, lty=3, lwd=1)
dev.off()
svg(paste0('plots/jaccard_hist_sr.', args[3], '.svg'))
hist(jm, xlim=c(0,1), breaks=sqrt(length(jm)), main='Distance distribution [SQRT]', xlab='J', prob=T)
lines(density(jm), col=4, lty=2)
abline(v=ori[2], col=2, lty=3, lwd=1)
dev.off()
if(IQR(jm) != 0) {
	svg(paste0('plots/jaccard_hist_fd.', args[3], '.svg'))
	hist(jm, xlim=c(0,1), breaks=1/(2*IQR(jm)/length(jm)**(1/3)), main='Distance distribution [FD]', xlab='J', prob=T)
	lines(density(jm), col=4, lty=2)
	abline(v=ori[2], col=2, lty=3, lwd=1)
	dev.off()
}

svg(paste0('plots/ipsen_hist.', args[3], '.svg'))
hist(im, xlim=c(0,1), main='Distance distribution', xlab='IM', prob=T)
lines(density(im), col=4, lty=2)
abline(v=ori[3], col=2, lty=3, lwd=1)
dev.off()
svg(paste0('plots/ipsen_hist_sr.', args[3], '.svg'))
hist(im, xlim=c(0,1), breaks=sqrt(length(im)), main='Distance distribution [SQRT]', xlab='IM', prob=T)
lines(density(im), col=4, lty=2)
abline(v=ori[3], col=2, lty=3, lwd=1)
dev.off()
if(IQR(im) != 0) {
	svg(paste0('plots/ipsen_hist_fd.', args[3], '.svg'))
	hist(im, xlim=c(0,1), breaks=1/(2*IQR(im)/length(im)**(1/3)), main='Distance distribution [FD]', xlab='IM', prob=T)
	lines(density(im), col=4, lty=2)
	abline(v=ori[3], col=2, lty=3, lwd=1)
	dev.off()
}

svg(paste0('plots/him_hist.', args[3], '.svg'))
hist(him, xlim=c(0,1), main='Distance distribution', xlab='HIM', prob=T)
lines(density(him), col=4, lty=2)
abline(v=ori[4], col=2, lty=3, lwd=1)
dev.off()
svg(paste0('plots/him_hist_sr.', args[3], '.svg'))
hist(him, xlim=c(0,1), breaks=sqrt(length(him)), main='Distance distribution [SQRT]', xlab='HIM', prob=T)
lines(density(him), col=4, lty=2)
abline(v=ori[4], col=2, lty=3, lwd=1)
dev.off()
if(IQR(him) != 0) {
	svg(paste0('plots/him_hist_fd.', args[3], '.svg'))
	hist(him, xlim=c(0,1), breaks=1/(2*IQR(him)/length(him)**(1/3)), main='Distance distribution [FD]', xlab='HIM', prob=T)
	lines(density(him), col=4, lty=2)
	abline(v=ori[4], col=2, lty=3, lwd=1)
	dev.off()
}

svg(paste0('plots/jim_hist.', args[3], '.svg'))
hist(jim, xlim=c(0,1), main='Distance distribution', xlab='H', prob=T)
lines(density(jim), col=4, lty=2)
abline(v=ori[5], col=2, lty=3, lwd=1)
dev.off()
svg(paste0('plots/jim_hist_sr.', args[3], '.svg'))
hist(jim, xlim=c(0,1), breaks=sqrt(length(jim)), main='Distance distribution [SQRT]', xlab='H', prob=T)
lines(density(jim), col=4, lty=2)
abline(v=ori[5], col=2, lty=3, lwd=1)
dev.off()
if(IQR(jim) != 0) {
	svg(paste0('plots/jim_hist_fd.', args[3], '.svg'))
	hist(jim, xlim=c(0,1), breaks=1/(2*IQR(jim)/length(jim)**(1/3)), main='Distance distribution [FD]', xlab='H', prob=T)
	lines(density(jim), col=4, lty=2)
	abline(v=ori[5], col=2, lty=3, lwd=1)
	dev.off()
}

svg(paste0('plots/him_scatterplot.', args[3], '.svg'), width=10, height=10)
plot(hm[lower.tri(hm, diag=F)], im[lower.tri(im, diag=F)], xlab='H', ylab='IM', main='H/IM space', xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,0,.1))
abline(h=0.5, lty=2, col=2)
abline(v=0.5, lty=2, col=2)
dev.off()

svg(paste0('plots/jim_scatterplot.', args[3], '.svg'), width=10, height=10)
plot(jm[lower.tri(jm, diag=F)], im[lower.tri(im, diag=F)], xlab='J', ylab='IM', main='J/IM space', xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,0,.1))
abline(h=0.5, lty=2, col=2)
abline(v=0.5, lty=2, col=2)
dev.off()

# plot

svg('plots/HIMspace.svg', width=10, height=10)
plot(t[,1], t[,3], xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,0,.1), xlab='H', ylab='IM', main='H/IM space after permuting ERG+/- classes')
abline(h=0.5, col=2, lty=2, lwd=2)
abline(v=0.5, col=2, lty=2, lwd=2)
abline(h=ori[3], col=4, lty=3, lwd=2)
abline(v=ori[1], col=4, lty=3, lwd=2)
dev.off()

svg('plots/HIMspace_zoom.svg', width=10, height=10)
plot(t[,1], t[,3], pch=20, col=rgb(0,0,0,.2), xlab='H', ylab='IM', main='H/IM space after permuting ERG+/- classes')
abline(h=ori[3], col=4, lty=3, lwd=2)
abline(v=ori[1], col=4, lty=3, lwd=2)
dev.off()


svg('plots/JIMspace.svg', width=10, height=10)
plot(t[,2], t[,3], xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,0,.1), xlab='J', ylab='IM', main='J/IM space after permuting ERG+/- classes')
abline(h=0.5, col=2, lty=2, lwd=2)
abline(v=0.5, col=2, lty=2, lwd=2)
abline(h=ori[3], col=4, lty=3, lwd=2)
abline(v=ori[2], col=4, lty=3, lwd=2)
dev.off()

svg('plots/JIMspace_zoom.svg', width=10, height=10)
plot(t[,2], t[,3], pch=20, col=rgb(0,0,0,.2), xlab='J', ylab='IM', main='J/IM space after permuting ERG+/- classes')
abline(h=ori[3], col=4, lty=3, lwd=2)
abline(v=ori[2], col=4, lty=3, lwd=2)
dev.off()

cat('~ END ~\n')

})
