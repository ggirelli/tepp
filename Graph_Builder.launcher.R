#!/usr/bin/env Rscript
source('./Graph_Builder.class.R')

gb.version <- 11

#-----------------#
# Read parameters #
#-----------------#

args <- commandArgs(trailingOnly=TRUE)

# Check number of parameters
if( length(args) != 1) stop('ERROR: script.R usage #1\n')

# Set default values
clusters <- 4
verbose <- TRUE
file.list <- list()
genes.label <- "Gene.id"
clonal.val <- c('clonal')
subclonal.val <- c('subclonal')
output.dir <- 'output'
attr.table <- ''
white.list <- list()
black.list <- list()
clean <- FALSE
write.cooc <- FALSE

# Read param file
if ( file.exists(args[1]) ) {
  file.name <- strsplit(args[1], '-p=')[[1]][2]
  if(file.exists(file.name)) {
    ff <- read.table(file.name, head=FALSE)

    clusters <- as.numeric(as.matrix(ff[which(ff[, 1] == 'clusters'), 2]))
    verbose <- as.logical(ff[which(ff[, 1] == 'verbose'), 2])

    output.dir <- toString(ff[which(ff[, 1] == 'output.dir'), 2])

    file.list$PM <- toString(ff[which(ff[, 1] == 'file-PM'), 2])
    if(file.list$PM == '') file.list$PM <- NA
    file.list$Gain <- toString(ff[which(ff[, 1] == 'file-Gain'), 2])
    if(file.list$Gain == '') file.list$Gain <- NA
    file.list$Loss <- toString(ff[which(ff[, 1] == 'file-Loss'), 2])
    if(file.list$Loss == '') file.list$Loss <- NA
    file.list$RR <- toString(ff[which(ff[, 1] == 'file-RR'), 2])
    if(file.list$RR == '') file.list$RR <- NA

    genes.label <- toString(ff[which(ff[, 1] == 'genes.label'), 2])
    white.list <- toString(ff[which(ff[, 1] == 'white.list'), 2])
    if(white.list != '') white.list <- read.table(white.list, header=FALSE, sep='\t')[,1]
    black.list <- toString(ff[which(ff[, 1] == 'black.list'), 2])
    if(black.list != '') black.list <- read.table(black.list, header=FALSE, sep='\t')[,1]
    clean <- as.logical(ff[which(ff[, 1] == 'clean'), 2])
    if(!is.logical(clean) || length(clean) == 0) clean <- FALSE
    write.cooc <- as.logical(ff[which(ff[, 1] == 'write.cooc'), 2])
    if(!is.logical(write.cooc) || length(write.cooc) == 0) write.cooc <- FALSE

    clonal.val <- strsplit(toString(ff[which(ff[, 1] == 'clonal.val'), 2]), ',')[[1]]
    if(length(clonal.val) == 0) clonal.val <- c('clonal')
    subclonal.val <- strsplit(toString(ff[which(ff[, 1] == 'subclonal.val'), 2]), ',')[[1]]
    if(length(subclonal.val) == 0) subclonal.val <- c('subclonal')

    attr.table <- toString(ff[which(ff[, 1] == 'attr.table'), 2])
    if(attr.table != '') attr.table <- read.table(attr.table, header=TRUE, sep='\t')
  }
} else {
  cat('\nNo param file detected.')
}

cat('\nExecuting script (v', gb.version, ') with', clusters, 'clusters.\n')
if(verbose) {
  cat('Verbosity at maximum.\n')
} else {
  cat('Verbosity at minimum.\n')
}
if(output.dir != '') cat('Output directory: ', output.dir, '\n')
cat('\n')

cat('PM-data: ', file.list$PM, '\n')
cat('Gain-data: ', file.list$Gain, '\n')
cat('Loss-data: ', file.list$Loss, '\n')
cat('RR-data: ', file.list$RR, '\n\n')

cat('Gene label: \'', genes.label, '\'\n')
if(clean) {
  cat('Running clean\n')
}
cat('\n')

if(length(white.list) != 0 && white.list != '') {
  cat('White list from: ', toString(ff[which(ff[, 1] == 'white.list'), 2]), '\n')
  print(as.character(white.list))
} else {
  cat('No white list specified\n')
}
if(length(black.list) != 0 && black.list != '') {
  cat('\nBlack list from: ', toString(ff[which(ff[, 1] == 'black.list'), 2]), '\n')
  print(as.character(black.list))
} else {
  cat('\nNo black list specified\n')
}
cat('\n')

cat('Clonal ids:\n')
print(clonal.val)
cat('Subclonal ids:\n')
print(subclonal.val)
cat('\n')

if(length(attr.table) != 0 && attr.table != '') {
  cat('Vertex attributes added to the final graph:\n')
  print(colnames(attr.table))
  cat('\n')
}

#---------#
# Execute #
#---------#

system.time({

  # Declare GraphManager instance
  gb <- GraphBuilder(clusters=clusters, verbose=verbose, genes.label=genes.label, white.list=white.list, black.list=black.list, clonal.val=clonal.val, subclonal.val=subclonal.val, attr.table=attr.table, clean=clean, output.dir=output.dir, write.cooc=write.cooc)
  cat('\n')

  #-------------------#
  # From MSSA to SSSA #
  #-------------------#

  # PM #
  if(length(file.list$PM) != 0 & !is.na(file.list$PM)) gb <- gb$readData(file.path=file.list$PM, sample.column='sample', abe.type='PM', temp.sample.list=gb$sample.list)
  cat('\n')

  # Gain #
  if(length(file.list$Gain) != 0 & !is.na(file.list$Gain)) gb <- gb$readData(file.path=file.list$Gain, sample.column='sample', abe.type='Gain', temp.sample.list=gb$sample.list)
  cat('\n')

  # Loss #
  if(length(file.list$Loss) != 0 & !is.na(file.list$Loss)) gb <- gb$readData(file.path=file.list$Loss, sample.column='sample', abe.type='Loss', temp.sample.list=gb$sample.list)
  cat('\n')

  # RR #
  if(length(file.list$RR) != 0 & !is.na(file.list$RR)) gb <- gb$readData(file.path=file.list$RR, sample.column='sample', abe.type='RR', temp.sample.list=gb$sample.list)
  cat('\n')

  #---------------------#
  # Prepare SSMA Graphs #
  #---------------------#

  gb$build(c('PM', 'Gain', 'Loss'), genes.label="Gene.id", clonality.label="clonality.status", sample.list=gb$sample.list)

})
