#!/usr/bin/env Rscript
source('./05_Graph_Manager.class.R')

#-----------------#
# Read parameters #
#-----------------#

args <- commandArgs(trailingOnly=TRUE)

# Check number of parameters
if(length(args)>2) { cat('ERROR: script.R usage\n'); quit() }
# Check parameters values and labels
for(arg in args) if(!grepl('^(-v=(TRUE|FALSE|T|F)|-c=[0-9]*)$', arg)) { cat('ERROR: script.R usage\n'); quit() }

# Set default values
clusters <- 4
verbose <- TRUE

# Read number of clusters
for(arg in args) if(grepl('^-c=[0-9]*$', arg))  clusters <- as.numeric(strsplit(arg, '-c=')[[1]][2])
if(clusters == 0) clusters <- 1
# Read verbose
for(arg in args) if(grepl('^-v=(TRUE|FALSE|T|F)$', arg)) verbose <- as.logical(strsplit(arg, '-v=')[[1]][2])

#---------#
# Execute #
#---------#

system.time({
	# Declare GraphManager instance
	gm <- GraphManager(clusters=clusters, verbose=verbose)
	# Read data 
	gm$builder <- gm$builder$readData(file.path='../Tables/20140317.TCGA.246FreezeSamples.PM.txt', sample.column='sample')
	gm$builder$buildGraph(gm$builder,table.out=TRUE)
	gm$mergeGraphs.noAttr(paste('./sample-graphs/gra_',unique(gm$builder$data[,'sample']),'.graphml',sep=''))
})
