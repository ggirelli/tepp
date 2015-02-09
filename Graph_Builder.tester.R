#!/usr/bin/env Rscript

# Read parameters
##############################
args <- commandArgs(trailingOnly=TRUE)

# Read param file
if ( file.exists(args[1]) ) {
	source(args[1])
} else {
	message('No param file detected.')
}

# Load libraries
library(igraph)
library(parallel)

# Write log
##############################
source(file.path(pathToGBclass, 'Graph_Builder.class.R'))

message(paste0('> Testing script (v', GraphBuilder()$version, ') with ', Ncores, ' cores.'))

if(verbose) {
	message('> Verbosity at maximum.')
} else {
	message('> Verbosity at minimum.')
}
if(output.dir != '') {
	message('> Output directory: ', output.dir, '\n')
} else {
	message('\n')
}

message('> PM-data: ', file.list$PM)
message('> Gain-data: ', file.list$Gain)
message('> Loss-data: ', file.list$Loss)
message('> RR-data: ', file.list$RR)

message('> Gene label: \'', genes.label, '\'\n')
if(clean) message('Running clean')

if(length(white.list) != 0 && white.list != '') {
	message('> White list from: ', toString(ff[which(ff[, 1] == 'white.list'), 2]))
	print(as.character(white.list))
} else {
	message('> No white list specified')
}

if(length(black.list) != 0 && black.list != '') {
	message('> Black list from: ', toString(ff[which(ff[, 1] == 'black.list'), 2]))
	print(as.character(black.list))
} else {
	message('> No black list specified')
}


message('\n> Clonal ids:')
print(clonal.val)
message('> Subclonal ids:')
print(subclonal.val)

if(length(attr.table) != 0 && attr.table != '') {
	message('\n> Vertex attributes added to the final graph:')
	print(colnames(attr.table))
}
message('')

# Read Basic Data
######################################
do <- list(
	pm=TRUE,
	gain=TRUE,
	loss=TRUE,
	rr=TRUE
)

#PM
if(is.na(file.list$PM)) {
	do$pm <- FALSE
} else {
	message('> Reading pm.table')
	pm.table <- read.delim(file.list$PM, as.is=T)
}

#Gain
if(is.na(file.list$Gain)) {
	do$gain <- FALSE
} else {
	message('> Reading gain.table')
	gain.table <- read.delim(file.list$Gain, as.is=T)
}

#Loss
if(is.na(file.list$Loss)) {
	do$loss <- FALSE
} else {
	message('> Reading loss.table')
	loss.table <- read.delim(file.list$Loss, as.is=T)
}

#RR
if(is.na(file.list$RR)) {
	do$rr <- FALSE
} else {
	message('> Reading rr.table')
	rr.table <- read.delim(file.list$RR, as.is=T)
}

# Test single-sample graphs
######################################

message('\n# SS CLONAL GRAPHS #')
flist <- list.files(paste0(output.dir, '/sample-graphs/'))
flist <- flist[grepl('^clonal_', flist)]
for (file in flist) {
	message('> Checking graph ', file)

	g <- read.graph(file.path(output.dir, '/sample-graphs/', file), format='graphml')
	
	sample <- unlist(strsplit(file, '.', fixed=T))
	sample <- paste(sample[1:length(sample)-1], collapse='.')
	sample <- unlist(strsplit(sample, '_', fixed=T))
	sample <- paste(sample[2:length(sample)], collapse='_')
	
	table <- c()
	if (do$pm) {
		if (sample %in% pm.table$sample) {
			sub.tab <- cbind(data.frame(pm.table[which(pm.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'pm')
			colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
			if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
		}
	}
	if (do$gain) {
		if (sample %in% gain.table$sample) {
			sub.tab <- cbind(data.frame(gain.table[which(gain.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'gain')
			colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
			if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
		}
	}
	if (do$loss) {
		if (sample %in% loss.table$sample) {
			sub.tab <- cbind(data.frame(loss.table[which(loss.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'loss')
			colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
			if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
		}
	}
	if (do$rr) {
		if (sample %in% rr.table$sample) {
			sub.tab <- cbind(data.frame(rr.table[which(rr.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'rr')
			colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
			if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
		}
	}
	table <- paste(table[,1], table[,2], table[,3], sep='~')

	message(' - Checking vertices')
	checkVertex=function(i, g, data) {
		v <- V(g)[i]

		v.name <- v$name
		v.aberration <- v$abe.type
		v.hugo <- v$HUGO

		# Check name construction
		if (v.name != paste0(v.hugo, '~', v.aberration)) {
			message('Vertex ', v.name, ' [#', i, '] has a wrong name.')
		}

		if (!paste0(v.hugo, '~clonal~', v.aberration) %in% data) message('    * Vertex ', v.hugo, '~', v.aberration, ' [#', i, '] is not present.')
	}
	v.test <- mclapply(1:vcount(g),
		FUN=checkVertex,
		g=g, data=table,
		mc.preschedule=TRUE,
		mc.cores=Ncores
	)
	if(length(which(degree(g, mode='out') != vcount(g))) != 0) message('    * ', length(which(degree(g, mode='out') != vcount(g))), ' vertices have wrong degree.')

	message(' - Checking edges')
	if( length(which(E(g)$weight != 1)) != 0 ) message('    * ', length(which(E(g)$weight != 1)), ' edges have wrong weight.')

}


# Test total graph
######################################

message('\n> Reading graphs.')
gt <- read.graph(file.path(output.dir, 'total_graph.graphml'), format='graphml')
gc <- read.graph(file.path(output.dir, 'clonal_graph.graphml'), format='graphml')
gs <- read.graph(file.path(output.dir, 'subclonal_graph.graphml'), format='graphml')
gn <- read.graph(file.path(output.dir, 'nonclonal_graph.graphml'), format='graphml')

message('\n# CLONAL GRAPH #')
message('#\n# Checking vertices\n#########################\n\n')
checkClonalVertices=function(i, gc) {
	#
	
	v <- V(gc)[i]
	v.name <- v$name
	v.aberration <- v$aberration
	v.hugo <- v$HUGO
	v.clonal <- v$clonal.occ
	v.subclonal <- v$subclonal.occ

	# Check name construction
	if (v.name != paste0(v.hugo, '~', v.aberration)) {
		message('Vertex ', v.name, ' [#', i, '] has a wrong name.')
	}

	# Check if they are clonal it at least one sample
	if ( eval(parse(text=paste0('length(which(', v.aberration, '.table$clonality.status[which(', v.aberration, '.table$', genes.label, '==v.hugo)] %in% clonal.val))'))) == 0 ) {
		message('Vertex ', v.name, ' [#', i, '] is never clonal.')
	}
	
}
vc.test <- mclapply(1:vcount(gc),
	FUN=checkClonalVertices,
	gc=gc,
	mc.preschedule=TRUE,
	mc.cores=Ncores
)

message('\n#\n# Checking edges\n#########################\n\n')
checkClonalEdges=function(i, gc, el) {
	#

	e <- E(gt)[i]

	e.source <- el[i,1]
	e.source.split <- unlist(strsplit(e.source, '~', fixed=T))
	e.source.hugo <- e.source.split[1]
	e.source.aberration <- e.source.split[2]

	e.target <- el[i,2]
	e.target.split <- unlist(strsplit(e.target, '~', fixed=T))
	e.target.hugo <- e.target.split[1]
	e.target.aberration <- e.target.split[2]

	e.weight <- e$weight

	# Check clonal co-occurrence
	source.samples <- eval(parse(text=paste0(e.source.aberration, '.table$sample[intersect(which(', e.source.aberration, '.table$', genes.label, '==e.source.hugo),which(', e.source.aberration, '.table$clonality.status %in% clonal.val))]')))
	target.clonal.samples <- eval(parse(text=paste0(e.target.aberration, '.table$sample[intersect(which(', e.target.aberration, '.table$', genes.label, '==e.target.hugo),which(', e.target.aberration, '.table$clonality.status %in% clonal.val))]')))
	test.clonal <- length(intersect(source.samples, target.clonal.samples))
	if(test.clonal != e.weight) {
		message('Edge #', i, ' from ', e.source, ' to ', e.target, ' has a wrong clonal co-occurrency: ', e.weight, ' instead of ', test.clonal)
	}

}
ec.list <- mclapply(1:ecount(gc),
	FUN=checkClonalEdges,
	gc=gc, el=get.edgelist(gc),
	mc.preschedule=TRUE,
	mc.cores=Ncores
)

message('\n# TOTAL GRAPH #')
message('#\n# Checking vertices\n#########################\n\n')

#(1) Check vertices
checkVertices=function(i, gt) {
	#
	
	v <- V(gt)[i]
	v.name <- v$name
	v.aberration <- v$aberration
	v.hugo <- v$HUGO
	v.clonal <- v$clonal.occ
	v.subclonal <- v$subclonal.occ

	# Check name construction
	if (v.name != paste0(v.hugo, '~', v.aberration)) {
		message('Vertex ', v.name, ' [#', i, '] has a wrong name.')
	}

	# Check presence in Basic data
	if (eval(parse(text=paste0('!do$', v.aberration)))) {
		message('Vertex ', v.name, ' [#', i, '] exists for an absent aberration type.')
	} else if (!v.hugo %in% eval(parse(text=paste0(v.aberration, '.table$', genes.label)))) {
		message('Vertex ', v.name, ' [#', i, '] is not present in the original data.')
	} else {

		#Check clonal occurrency
		cs <- eval(parse(text=paste0(v.aberration, '.table$clonality.status[which(', v.aberration, '.table$', genes.label, ' == v.hugo)]')))
		if(length(which(cs %in% clonal.val)) != v.clonal) {
			message('Vertex ', v.name, ' [#', i, '] clonal occurrence is wrong: ', v.clonal, ' instead of ', length(which(cs %in% clonal.val)))
			#message('Cheking in clonal_graph')
		}
		if(length(which(cs %in% subclonal.val)) != v.subclonal) {
			message('Vertex ', v.name, ' [#', i, '] subclonal occurrence is wrong: ', v.subclonal, ' instead of ', length(which(cs %in% subclonal.val)))
		}
	}
}
v.test <- mclapply(1:vcount(gt),
	FUN=checkVertices,
	gt=gt,
	mc.preschedule=TRUE,
	mc.cores=Ncores
)

message('#\n\n# Checking Edges
	\n#########################\n\n')

#(2) Check edges
checkEdges=function(i, gt, el) {
	#
	
	e <- E(gt)[i]

	e.source <- el[i,1]
	e.source.split <- unlist(strsplit(e.source, '~', fixed=T))
	e.source.hugo <- e.source.split[1]
	e.source.aberration <- e.source.split[2]

	e.target <- el[i,2]
	e.target.split <- unlist(strsplit(e.target, '~', fixed=T))
	e.target.hugo <- e.target.split[1]
	e.target.aberration <- e.target.split[2]

	e.weight <- e$weight
	e.clonal <- e$clonal.cooc
	e.subclonal <- e$subclonal.cooc

	# Check source/target
	source.samples <- eval(parse(text=paste0(e.source.aberration, '.table$sample[intersect(which(', e.source.aberration, '.table$', genes.label, '==e.source.hugo),which(', e.source.aberration, '.table$clonality.status %in% clonal.val))]')))
	target.samples <- eval(parse(text=paste0(e.target.aberration, '.table$sample[intersect(which(', e.target.aberration, '.table$', genes.label, '==e.target.hugo),which(', e.target.aberration, '.table$clonality.status %in% subclonal.val))]')))
	test.weight <- length(intersect(source.samples, target.samples))
	if(test.weight == 0) {
		message('Edge #', i, ' from ', e.source, ' to ', e.target, ' is not present in the Basic Data.')
	}
	if(test.weight != e.weight) {
		message('Edge #', i, ' from ', e.source, ' to ', e.target, ' has a wrong weight.')
	}

	# Check clonal co-occurrence
	target.clonal.samples <- eval(parse(text=paste0(e.target.aberration, '.table$sample[intersect(which(', e.target.aberration, '.table$', genes.label, '==e.target.hugo),which(', e.target.aberration, '.table$clonality.status %in% clonal.val))]')))
	test.clonal <- length(intersect(source.samples, target.clonal.samples))
	if(test.clonal != e.clonal) {
		message('Edge #', i, ' from ', e.source, ' to ', e.target, ' has a wrong clonal co-occurrency: ', e.clonal, ' instead of ', test.clonal)
	}

	# Check subclonal co-occurrence
	source.subclonal.samples <- eval(parse(text=paste0(e.source.aberration, '.table$sample[intersect(which(', e.source.aberration, '.table$', genes.label, '==e.source.hugo),which(', e.source.aberration, '.table$clonality.status %in% subclonal.val))]')))
	test.subclonal <- length(intersect(source.subclonal.samples, target.samples))
	if(test.subclonal != e.subclonal) {
		message('Edge #', i, ' from ', e.source, ' to ', e.target, ' has a wrong subclonal co-occurrency: ', e.subclonal, ' instead of ', test.subclonal)
	}
}
e.test <- mclapply(1:ecount(gt),
	FUN=checkEdges,
	gt=gt, el=get.edgelist(gt),
	mc.preschedule=TRUE,
	mc.cores=Ncores
)