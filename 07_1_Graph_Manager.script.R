#!/usr/bin/env Rscript
source('./07_Graph_Manager.class.R')

#-----------------#
# Read parameters #
#-----------------#

args <- commandArgs(trailingOnly=TRUE)

# Check number of parameters
if(length(args) > 2 || length(args) == 0) { cat('ERROR: script.R usage #1\n'); quit() }
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

cat('\nExecuting script with', clusters, 'clusters.\n')
if(verbose) {
  cat('Verbosity at maximum.\n')
} else {
  cat('Verbosity at minimum.\n')
}
cat('\nPM-data: ', file.list$PM, '\n')
cat('Gain-data: ', file.list$Gain, '\n')
cat('Loss-data: ', file.list$Loss, '\n')
cat('RR-data: ', file.list$RR, '\n\n')

#---------#
# Execute #
#---------#

system.time({

  # Declare GraphManager instance
  gm <- GraphManager(clusters=clusters, verbose=verbose)
  cat('\n')

  #-------------------#
  # From MSSA to SSSA #
  #-------------------#

  # PM #
  if(file.list$PM != "") gm <- gm$readData(file.path=file.list$PM, sample.column='sample', abe.type='PM', temp.sample.list=gm$sample.list)
  cat('\n')

  # Gain #
  if(file.list$Gain != "") gm <- gm$readData(file.path=file.list$Gain, sample.column='sample', abe.type='Gain', temp.sample.list=gm$sample.list)
  cat('\n')

  # Loss #
  if(file.list$Loss != "") gm <- gm$readData(file.path=file.list$Loss, sample.column='sample', abe.type='Loss', temp.sample.list=gm$sample.list)
  cat('\n')

  # RR #
  if(file.list$RR != "") gm <- gm$readData(file.path=file.list$RR, sample.column='sample', abe.type='RR', temp.sample.list=gm$sample.list)
  cat('\n')

  #---------------------#
  # Prepare SSMA Graphs #
  #---------------------#

  gm$build(c('PM', 'Gain', 'Loss'), genes.label="Gene.id", clonality.label="clonality.status", clonal.val="clonal", subclonal.val="subclonal", sample.list=gm$sample.list)

})
