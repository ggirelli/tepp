#!/usr/bin/env Rscript
source('./06_Graph_Manager.class.R')

#-----------------#
# Read parameters #
#-----------------#

args <- commandArgs(trailingOnly=TRUE)

# Check number of parameters
if(length(args)>2) { cat('ERROR: script.R usage #1\n'); quit() }
# Check parameters values and labels
for(arg in args) if(!grepl('^(-v=(TRUE|FALSE|T|F)|-c=[0-9]*|-p=[a-zA-Z0-9()./_~-]*)$', arg)) { cat('ERROR: script.R usage #2\n'); quit() }

# Set default values
clusters <- 4
verbose <- TRUE
file.list <- list()

# Read number of clusters
for(arg in args) if(grepl('^-c=[0-9]*$', arg))  clusters <- as.numeric(strsplit(arg, '-c=')[[1]][2])
if(clusters == 0) clusters <- 1
# Read verbose
for(arg in args) if(grepl('^-v=(TRUE|FALSE|T|F)$', arg)) verbose <- as.logical(strsplit(arg, '-v=')[[1]][2])
# Read param file
for(arg in args) if(grepl('^-p=[a-zA-Z0-9()./_~-]*$', arg)) {
	file.name <- strsplit(arg, '-p=')[[1]][2]
	if(file.exists(file.name)) {
		file <- read.table(file.name, head=FALSE)
		clusters <- as.numeric(as.matrix(file[which(file[, 1] == 'clusters'), 2]))
		verbose <- as.logical(file[which(file[, 1] == 'verbose'), 2])
		file.list$PM <- toString(file[which(file[, 1] == 'file-PM'), 2])
		file.list$Gain <- toString(file[which(file[, 1] == 'file-Gain'), 2])
		file.list$Loss <- toString(file[which(file[, 1] == 'file-Loss'), 2])
		file.list$RR <- toString(file[which(file[, 1] == 'file-RR'), 2])
	}
}

cat('Executing script with', clusters, 'clusters.\n')
if(verbose) {
	cat('Verbosity at maximum.\n')
} else {
	cat('Verbosity at minimum.\n')
}
print(file.list)

#---------#
# Execute #
#---------#

system.time({
	# Declare GraphManager instance
	gm <- GraphManager(clusters=clusters, verbose=verbose)
	
	# PM #
	if(file.list$PM != "") {
		cat('\n######\n# PM #\n######\n')
		gm$builder <- gm$builder$readData(file.path=file.list$PM, sample.column='sample', abe.type='PM')
		gm$builder$build(gm$builder, abe.type="PM")
	}

	# Gain #
	if(file.list$Gain != "") {
		cat('\n########\n# Gain #\n########\n')
		gm$builder <- gm$builder$readData(file.path=file.list$Gain, sample.column='sample', abe.type='Gain')
		gm$builder$build(gm$builder, abe.type="Gain")
	}

	# Loss #
	if(file.list$Loss != "") {
		cat('\n########\n# Loss #\n########\n')
		gm$builder <- gm$builder$readData(file.path=file.list$Loss, sample.column='sample', abe.type='Loss')
		gm$builder$build(gm$builder, abe.type="Loss")
	}

	# RR #
	if(file.list$RR != "") {
		cat('\n######\n# RR #\n######\n')
		gm$builder <- gm$builder$readData(file.path=file.list$RR, sample.column='sample', abe.type='RR')
		gm$builder$build(gm$builder, abe.type="RR")
	}

#	gm$mergeGraphs.noAttr(paste('./sample-graphs/gra_',unique(gm$builder$data[,'sample']),'.graphml',sep=''))
})
