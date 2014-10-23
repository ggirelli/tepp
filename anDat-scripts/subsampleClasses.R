#!/usr/bin/env Rscript
#
# ./preparePermutations.R paramFile
# 
# SSDNs must be in the folder 'sample-graphs'
# Requires the following files in the same folder:
#     * Graph_Builder.class.R
#     * Graph_Builder.script.R
#     * parentDir.Graph_Builder.class.R
#     * parentDir.Graph_Builder.script.R
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

skip.prep <- getp(ps, 'skip.prep')
if(!is.na(skip.prep)) {
	if(skip.prep == 'TRUE') {
		cat("    * Skipping preparation.\n")
		skip.prep <- TRUE
	} else {
		skip.prep <- FALSE
	}
} else {
	skip.prep <- FALSE
}

skip.build <- getp(ps, 'skip.build')
if(!is.na(skip.build)) {
	if(skip.build == 'TRUE') {
		cat("    * Skipping building.\n")
		skip.build <- TRUE
	} else {
		skip.build <- FALSE
	}
} else {
	skip.build <- FALSE
}

# ---------- #
# INITIATION #
# ---------- #

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

if(!skip.prep) {

	# ------------------------- #
	# PREPARE ORIGINAL DISTANCE #
	# ------------------------- #
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

	# -------------------------- #
	# PREPARE SUBSAMPLE DISTANC F
	# -------------------------- #
	cat('> Preparing subsamples\n')

	for (i in 1:length(flist)) {
		cat(paste0('	- Subsample #', i, '\n'))
		flist.subsample <- flist[-i]
		pfcl.subsample <-pfcl[-i,]
		sample.rm <- flist[i]

		# Mkdir and cd
		dir.create(paste0('s', i), showWarnings=F)
		setwd(paste0('s', i))

		dir.create(file.path(p.outdir), showWarnings=F)
		dir.create(file.path(m.outdir), showWarnings=F)

		gainTab.out <- NA
		lossTab.out <- NA
		pmTab.out <- NA

		if(!is.na(gainTab)) {
			gainTab.out <- 'erg_clonTab.Gain.txt'
			lossTab.out <- 'erg_clonTab.Loss.txt'

			# Retrieve gain/loss data
			gain <- read.table(file.path('..', gainTab), header=T, sep='\t')
			gain.rm <- which(paste0('gra_', gain$sample, '.graphml') == sample.rm)
			if(length(gain.rm) != 0) gain <- gain[-gain.rm,]
			loss <- read.table(file.path('..', lossTab), header=T, sep='\t')
			loss.rm <- which(paste0('gra_', loss$sample, '.graphml') == sample.rm)
			if(length(loss.rm) != 0) loss <- loss[-loss.rm,]

			# Dividing classes in gain/loss
			cat('		* Dividing gain/loss in classes\n')
			erg.p.gain <- gain[which(gain$sample %in% pfcl.subsample$sample[which(pfcl.subsample$ERG.clean > 0)]),]
			erg.p.loss <- loss[which(loss$sample %in% pfcl.subsample$sample[which(pfcl.subsample$ERG.clean > 0)]),]
			write.table(erg.p.gain, paste0(p.outdir, '/erg_clonTab.Gain.txt'), quote=F, sep='\t')
			write.table(erg.p.loss, paste0(p.outdir, '/erg_clonTab.Loss.txt'), quote=F, sep='\t')

			erg.m.gain <- gain[which(gain$sample %in% pfcl.subsample$sample[which(pfcl.subsample$ERG.clean == 0)]),]
			erg.m.loss <- loss[which(loss$sample %in% pfcl.subsample$sample[which(pfcl.subsample$ERG.clean == 0)]),]
			write.table(erg.m.gain, paste0(m.outdir, '/erg_clonTab.Gain.txt'), quote=F, sep='\t')
			write.table(erg.m.loss, paste0(m.outdir, '/erg_clonTab.Loss.txt'), quote=F, sep='\t')
		}
		if(!is.na(pmTab)) {
			pmTab.out <- 'erg_clonTab.pm.txt'

			# Retrieve gain/loss data
			pm <- read.table(file.path('..', pmTab), header=T, sep='\t')
			pm.rm <- which(paste0('gra_', pm$sample, '.graphml') == sample.rm)
			if(length(pm.rm) != 0) pm <- pm[-pm.rm,]

			# Dividing classes in gain/loss
			cat('		* Dividing pm in classes\n')
			erg.p.pm <- pm[which(pm$sample %in% pfcl.subsample$sample[which(pfcl.subsample$ERG.clean > 0)]),]
			erg.m.pm <- pm[which(pm$sample %in% pfcl.subsample$sample[which(pfcl.subsample$ERG.clean == 0)]),]
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

		setwd('..')
	}
}

# -------- #
# BUILDING #
# -------- #

if(!skip.build) {
	cat('> Building dependency networks and calculating distances.\n')

	cat('	- Working on original data\n')
	cat('		* Building ERG-\n')
	system(paste0('./Graph_Builder.launcher.R -p=param.ergm.txt > log.ergm.dat'))
	cat('		* Building ERG+\n')
	system(paste0('./Graph_Builder.launcher.R -p=param.ergp.txt > log.ergp.dat'))

	cat('	- Working on subsamples\n')
	for (i in 1:length(flist)) {
		system.time({
			cat(paste0('    	* Building subsample #', i, "\n"))
			setwd(paste0('s', i))

			cat(paste0("         	+ Building ERGm dependency graph.\n"))
			system(paste0('../parentDir.Graph_Builder.launcher.R -p=param.ergm.txt > log.s', i, '.ergm.dat'))
			cat(paste0("         	+ Building ERGp dependency graph.\n"))
			system(paste0('../parentDir.Graph_Builder.launcher.R -p=param.ergp.txt > log.s', i, '.ergp.dat'))

			setwd('..')
		})
	}
}

# ------------------ #
# CALCULATE DISTANCE #
# ------------------ #

cat(paste0('    * Measuring distances\n'))

cat('	* Working on original data.\n')
source('./Graph_Manager.class.R')
gm <- read.graph(paste0(p.outdir, '/total_graph.graphml'), format='graphml')
gp <- read.graph(paste0(m.outdir, '/total_graph.graphml'), format='graphml')
ds <- GraphManager()$calcDistances(gm, gp, 1)
write(ds, 'distances.dat')

cat('	* Working on permutations.\n')
cores <- makeCluster(nCores)
registerDoParallel(cores)

res <- foreach(i=1:length(flist), .combine=rbind) %dopar% {
	setwd(paste0('s', i))

	#cat(paste0('    * Measuring distance for permutation #', i, "\n"))
	source('../parentDir.Graph_Manager.class.R')
	gm <- read.graph(paste0(p.outdir, '/total_graph.graphml'), format='graphml')
	gp <- read.graph(paste0(m.outdir, '/total_graph.graphml'), format='graphml')
	ds <- GraphManager()$calcDistances(gm, gp, 1)
	write.table(ds, 'distances.dat', row.names=F, col.names=F, quote=F, sep='\t')

	setwd('..')

	return(ds)
}

write.table(res, 'dist.subs.dat', row.names=F, col.names=F, quote=F, sep='\t')

stopCluster(cores)

cat('~ END ~\n')

})
