# Gabriele Girelli <gabriele@filopoe.it>
# Elaborates clonal/subclonal gene groups into graphs
source("C:/Users/Gire/Desktop/TumorEvolution_project/Rscript/contingency_table_parallel.lib.R")
library('igraph')
library('doParallel')

mkSetWdDir <- function(dir.name) {
  # Prepare the current working directory, i.e.: it creates it or simply sets it
  #
  # Args:
  #   dir.name: absolute path of the working directory
  #
  # Returns:
  #   null
  if (file.exists(dir.name)) {
    setwd(file.path(dir.name))
  } else {
    dir.create(file.path(dir.name))
    setwd(file.path(dir.name))
  }
  return(NULL)
}

# Read data table
splitDataIntoSamples <- function(file.name, clusters) {
  # Splits the data from the table into single-sample files
  #
  # Args:
  #   file.name: the pathway of the data file relative to the wdir
  #   cluster: number of cores to use for parallel computing
  #
  # Returns:
  #   Data contained in the given file.name
  data <- read.table(file.name, header=TRUE, sep="\t", row.names=NULL)
  # Prepare sample list without duplicates
  sample.list <- unique(data$sample)
  cl <- makeCluster(clusters)
  registerDoParallel(cl)
  # Split original data table for each sample in data.frames
  # and save them in temporary directory
  foreach(i = 1:length(sample.list)) %dopar% {
    # On which sample are we working?
    sample.id <- sample.list[i]
    # Get row_ids from original data table for the working_sample
    row.ids <- which(data$sample ==sample.id)
    # Retrieve genes and clonality status for working_sample
    sample.genes <- data$Gene.id[row.ids]
    sample.status <- data$clonality.status[row.ids]
    # Make data.frame and append it to data.frame list
    data.frame <- data.frame(genes = sample.genes, clonality = sample.status)
    write.table(data.frame, file = file.path('.', sample.list[i]))
  }
  stopCluster(cl)
  return(data)
}

mkGraphFromAdjTable <- function(file.name) {
  # Build a graph from an adjacency table contained in a file
  #
  # Args:
  #   file.name: path to the adjacency table file
  #
  # Returns:
  #   null
  data <- read.table(file.name, header = TRUE)
  # Make small contingency matrix
  adjacency.matrix <- matrix(0,
                             nrow = length(data$genes),
                             ncol = length(data$genes)
  )
  rownames(adjacency.matrix)<- data$genes
  colnames(adjacency.matrix)<- data$genes
  # Distinguish clonal and subclonal
  genes.clonal <- which(data$clonality == 'clonal')
  genes.subclonal <- which(data$clonality == 'subclonal')
  # Increment cells in matrices
  for (i in seq(length(genes.clonal))){
    for (j in seq(length(genes.subclonal))){
      adjacency.matrix[genes.clonal[i], genes.subclonal[j]] <- 
        adjacency.matrix[genes.clonal[i], genes.subclonal[j]] + 1
    }
  }
  # Write small matrix
  # write.table(
  #  adjacency.matrix,
  #  file = file.path('.', paste('tab_', sample.id, '.dat', sep = ""))
  #  )
  write.graph(
    graph.adjacency(
      adjacency.matrix,
      mode="directed",
      weighted=TRUE
    ),
    file = file.path('.', paste('gra_', sample.id, '.graphml', sep = "")),
    format = 'graphml'
  )
}

mkSamplesGraph <- function(clusters, samples) {
  # Build graphs of every sample
  #
  # Args:
  #   clusters: number of cores to use for parallel computing
  #
  # Returns:
  #   null
  cl <- makeCluster(clusters)
  registerDoParallel(cl)
  sample.list <- unique(samples)
  # Make graph of every sample
  foreach(i = 1:length(sample.list)) %dopar% {
    # Working sample
    sample.id <- sample.list[i]
    print(paste('- working on sample:', sample.id, sep=" "))
    # Make graph
    mkGraphFromAdjTable(file.path('.', sample.id))
  }
  stopCluster(cl)
  return(NULL)
}

runGraphMaker <- function() {
  # Ask for where to set current working directory
  ans <- 0
  while (ans < 1) {
    print('Specify current directory absolute path:')
    dir.name <- readline('')
    ans <- ifelse(grepl('[^a-z0-9A-Z/\\:._-]', dir.name), 0, 1)
  }
  ans <- 0
  while (ans < 1) {
    print('Specify file name (absolute path or relative to working directory):')
    file.name <- readline('')
    ans <- ifelse(grepl('[^a-z0-9A-Z/\\:._-]', file.name), 0, 1)
  }
  ans <- 0
  while (ans < 1) {
    print('Number of cores to use for parallel computing:')
    clusters <- as.integer(readline(''))
    ans <- ifelse(grepl('[^0-9]', clusters), 0, 1)
  }
  # Act
  print('Preparing environment.')
  mkSetWdDir(dir.name)
  print('Splitting data.')
  data <- splitDataIntoSamples(file.name,clusters)
  print('Make graphs:')
  mkSamplesGraph(clusters, data$sample)
  return(NULL)
}