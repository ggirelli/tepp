GraphMaker <- function() {
  library('igraph')
  library('doParallel')

  # Define empty GraphMaker
  gm <- list(
    file.name = '',
    dir.name = '',
    clusters = 0,
    data = matrix()
  )

  # Define readInput function
  gm$input <- function() {
    # Asks where to set current working directory
    #
    # Returns:
    #   The list of input values
    ans <- 0
    while (ans < 1) {
      print('Specify current directory absolute path:')
      dir.name <- readline('')
      ans <- ifelse(grepl('[^a-z0-9A-Z/\\:._-]', dir.name), 0, 1)
    }
    # Ask the file position
    ans <- 0
    while (ans < 1) {
      print('Specify file name (abs path or relative to working directory):')
      file.name <- readline('')
      ans <- ifelse(grepl('[^a-z0-9A-Z/\\:._-]', file.name), 0, 1)
    }
    # Ask the number of cores
    ans <- 0
    while (ans < 1) {
      print('Number of cores to use for parallel computing:')
      clusters <- readline('')
      ans <- ifelse(grepl('[^0-9]', clusters), 0, 1)
    }
    return(list(
      dir.name = dir.name,
      file.name = file.name,
      clusters = as.integer(clusters)
    ))
  }

  # Define WD function
  gm$WD <- function() {
    # Prepare the current working directory,
    # i.e.: it creates it or simply sets it
    #
    # gm$args:
    #   dir.name
    #
    # Returns:
    #   Message
    if (file.exists(gm$dir.name)) {
      setwd(file.path(gm$dir.name))
    } else {
      dir.create(file.path(gm$dir.name))
      setwd(file.path(gm$dir.name))
    }
    return('SET Working Directory')
  }

  # Define splitData function
  gm$splitData <- function(data = gm$data) {
    # Splits the data from the table into single-sample files
    #
    # gm$args:
    #   file.name
    #   cluster
    #
    # Returns:
    #   Data contained in the given file.name
    print('Splitting Data')
    # Prepare sample list without duplicates
    sample.list <- unique(data$sample)
    # If needed, create output directory
    if (!file.exists('./sample-data/')) {
      dir.create(file.path('./sample-data/'))
    }
    # Declare parallelism
    cl <- makeCluster(gm$clusters)
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
      write.table(
        data.frame,
        file=file.path('./sample-data/', sample.list[i])
      )
    }
    stopCluster(cl)
    # Terminate
    return('Data splitted')
  }

  # Define build function
  gm$build <- function(sample.id, print.table = FALSE) {
    # Build a graph from an adjacency table contained in a file
    #
    # Args:
    #   sample.id: working sample
    #   print.table: whether or not to write also the table file
    #      (default: FALSE)
    #
    # Returns:
    #   Message
    library('igraph')
    data <- read.table(file.path('./sample-data/', sample.id), header = TRUE)
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
    if (print.table) {
      write.table(
        adjacency.matrix,
        file = file.path('./samples-tabs/', paste('tab_', sample.id, '.dat', sep = ""))
      )
    }
    write.graph(
      graph.adjacency(
        adjacency.matrix,
        mode="directed",
        weighted=TRUE
      ),
      file=file.path(
        './sample-graphs/',
        paste('gra_', sample.id, '.graphml', sep = "")
      ),
      format='graphml'
    )
    # Terminate
    return(paste('Built graph for', file.path('.', sample.id), "\n", sep=" "))
  }

  # Define mergeGraphs function
  gm$mergeGraphs <- function(file.list) {
    # Merges all the given .graphml file into a single graph
    #
    # Args:
    #   file.list: an array of file names
    #
    # gm$args:
    #   clusters
    #
    # Return:
    #   Message
    #
    # TODO:
    #   - Issue in merging
    #   - Does not read files
    print("Merging")
    # Declare parallelism
    cl <- makeCluster(gm$clusters)
    registerDoParallel(cl)
    # Get edges
    edges <- foreach(i = 1:length(file.list), .combine=rbind) %dopar% {
      library('igraph')
      # Working sample
      file.name <- file.list[i]
      g <- read.graph(file.path('.', file.name), format='graphml')
      get.edgelist(g)
    }
    stopCluster(cl)
    print("Read edgelist")
    g <- graph.empty()
    print("Made empty graph")
    g <- g + vertices(unique(as.vector(edges)))
    print("Added vertices")
    # Clean edges
    print(length(edges))
    edges <- unique(edges)
    print(length(edges))
    print('Preparing edges')
    cl <- makeCluster(gm$clusters)
    registerDoParallel(cl)
    edges <- foreach(i = 1:length(edges[,1]), .combine=c) %dopar% {
      return(edges[i,])
    }
    stopCluster(cl)
    print(edges)
    g <- g + edges(edges)
    print("Added edges")
    write.graph(g, file=file.path('.', 'total_graph.graphml'), format="graphml")
    print("Written graph")
    return('Merged graphs')
  }

  # Define run function
  gm$run <- function(data = gm$data, merge=TRUE, table=FALSE) {
    # Build graphs of every sample
    #
    # Args:
    #   merge: whether or not to merge all the single-sample graphs in one
    #
    # gm$args:
    #   clusters
    #
    # gm$methods():
    #   splitData()
    #   build()
    #
    # Returns:
    #   Message
    print(gm$splitData())
    # If needed, create output directory
    if (!file.exists('./sample-graphs/')) {
      dir.create(file.path('./sample-graphs/'))
    }
    # Declare parallelism
    cl <- makeCluster(gm$clusters)
    registerDoParallel(cl)
    print('Building graphs')
    sample.list <- unique(gm$data$sample)
    # Make graph of every sample
    f <- foreach(i = 1:length(sample.list)) %dopar% {
      # Working sample
      sample.id <- sample.list[i]
      # Make graph
      gm$build(sample.id=sample.id, print.table=table)
    }
    stopCluster(cl)
    print(f)
    print('All graphs built')
    if (merge) {
      # Merges graphs
      print('Merging graphs')
      gm$mergeGraphs(paste(
        './sample-graphs/gra_',
        unique(gm$data$sample),
        '.graphml',
        sep=''
      ))
    }
    return('All done.')
  }

  # Read input
  tgm <- gm$input()
  gm$dir.name <- tgm$dir.name
  gm$file.name <- tgm$file.name
  gm$clusters <- tgm$clusters
  # Set WD
  print(gm$WD())
  # Read data
  print('Reading data')
  gm$data <- read.table(gm$file.name, header=TRUE, sep="\t", row.names=NULL)
  print('Done')

  # Terminate constructor
  class(gm) <- 'GraphMaker'
  return(gm)
}