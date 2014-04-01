# Class to build graphs
GraphBuilder <- function(clusters=0, verbose=FALSE) {
  library('igraph')
  library('doParallel')
  
  # Define GraphBuilder attributes
  gb <- list(
    # Number of cluster for parallel computing
    clusters = clusters,
    
    # Boolean value to determine verbosity
    verbose = verbose,
    
    readData = function(file.path, header=TRUE, sep='\t', row.names=NULL, sample.column=NULL) {
      # Reads data
      #
      # Args:
      #   file.path: the name of the file which the data are to be
      #		read from or a a readable text-mode connection.
      #	header: a logical value indicating whether the file contains
      #		the names of the variables as its first line.
      #	sep: the field separator character.
      #	row.names: a vector of row names.
      #	sample.column: if file is multi-sample, specify sample
      #		column header.
      #
      # Returns:
      #   Modified GraphBuilder instance
      
      # Read data
      if(gb$verbose) print('Reading whole data')
      data <- read.table(file.path, header=header, sep=sep, row.names=row.names)
      # Save raw data into GraphBuilder instance
      gb$data <- data
      
      # Check for multi-sampling
      if(is.null(sample.column)) {
        # If single-sample file
        return(gb)
      } else {
        # If multi-sample file
        gb <- gb$splitData(data, sample.column)
        return(gb)
      }
    },
    
    splitData = function(data, sample.column, clusters=gb$clusters) {
      # Splits multi-sample data
      #
      # Args:
      #   data: the multi-sample data to split.
      #	sample-column: the name of the sample-id column.
      #	clusters: the number .
      #
      # Returns:
      #	Modified GraphBuilder instance
      
      if(gb$verbose) print('Splitting whole data')
      # Prepare sample list without duplicates
      sample.list <- unique(data[,sample.column])
      # If needed, create output directory
      if (!file.exists('./sample-data/')) {
        dir.create(file.path('./sample-data/'))
      }
      # Declare parallelism
      par <- makeCluster(clusters)
      registerDoParallel(par)
      # Split original data table for each sample in data.frames
      # and save them in temporary directory
      foreach(i=1:length(sample.list)) %dopar% {
        # On which sample are we working?
        sample.id <- sample.list[i]
        # Get row_ids from original data table for the working_sample
        row.ids <- which(data[,sample.column]==sample.id)
        # Write selected rows
        write.table(data[row.ids,], file=file.path('./sample-data/', sample.list[i]))
      }
      stopCluster(par)
      # Terminate
      if(gb$verbose) print('Data splitted')
      gb$split <- TRUE
      gb$isMultiSample <- TRUE
      gb$sample.column <- sample.column
      gb$data <- data
      return(gb)
    },
    
    buildGraph = function(
      gb,
      genes.label="Gene.id",
      clonality.label="clonality.status",
      clonal.val="clonal",
      subclonal.val="subclonal",
      print.table=FALSE
    ) {
      # Builds graphs, can discriminate between single- and multi-sample cases.
      #
      # Args:
      #   gb: current GraphBuilder instance.
      #	  genes.label: label of the gene's id column.
      #	  clonality.label: label of the clonality status column.
      #   clonal.val: value of clonality for clonal genes.
      #   subclonal.val: value of clonality for subclonal genes.
      #   print.table: whether or not to print the adjacency tables.
      #
      # Returns:
      #	Modified GraphBuilder instance

      if(is.null(gb$isMultiSample)) {
        # Single sample data
        # Make empty adjacency matrix
        adjacency.matrix <- matrix(0, nrow=length(gb$data[,genes.label]), ncol=length(gb$data[,genes.label]))
        rownames(adjacency.matrix) <- gb$data[,genes.label]
        colnames(adjacency.matrix) <- gb$data[,genes.label]
        # Distinguish clonal and subclonal
        genes.clonal <- which(gb$data[,clonality.label]==clonal.val)
        genes.subclonal <- which(gb$data[,clonality.label]==subclonal.val)
        # Fill matrix
        for (i in seq(length(genes.clonal))){
          for (j in seq(length(genes.subclonal))){
            adjacency.matrix[genes.clonal[i], genes.subclonal[j]] <-
              adjacency.matrix[genes.clonal[i], genes.subclonal[j]] + 1
          }
        }
        # Write small matrix
        if (print.table) {
          file.name <- paste('tab_', sample.id, '.dat', sep = "")
          write.table(adjacency.matrix, file = file.path('.', file.name))
        }
        g <- graph.adjacency(adjacency.matrix, mode="directed", weighted=TRUE)
        file.name <- paste('gra_', sample.id, '.graphml', sep = "")
        write.graph(g, file = file.path('.', file.name), format='graphml')
        # Terminate
        if(gb$verbose) print('Built graph')
      } else {
        # Multi sample data
		# If needed, create output directory
		if (!file.exists('./sample-graphs/')) dir.create(file.path('./sample-graphs/'))
		if (!file.exists('./sample-data/') && print.table) dir.create(file.path('./sample-tabs/'))
        # Declare parallelism
        par <- makeCluster(gb$clusters)
        registerDoParallel(par)
        if(gb$verbose) print('Building graphs')
        sample.list <- unique(gb$data[,gb$sample.column])
        # Make graph of every sample
        response <- foreach(i=1:length(sample.list)) %dopar% {
          library('igraph')
          # Working sample
          sample.id <- sample.list[i]
          # Read data
          data <- read.table(file.path('./sample-data/', sample.id), header = TRUE)
          # Make empty adjacency matrix
          adjacency.matrix <- matrix(0, nrow=length(data[,genes.label]), ncol=length(data[,genes.label]))
          rownames(adjacency.matrix) <- data[,genes.label]
          colnames(adjacency.matrix) <- data[,genes.label]
          # Distinguish clonal and subclonal
          genes.clonal <- which(data[,clonality.label]==clonal.val)
          genes.subclonal <- which(data[,clonality.label]==subclonal.val)
          # Fill matrix
          for (i in seq(length(genes.clonal))){
            for (j in seq(length(genes.subclonal))){
              adjacency.matrix[genes.clonal[i], genes.subclonal[j]] <-
                adjacency.matrix[genes.clonal[i], genes.subclonal[j]] + 1
            }
          }
          # Write small matrix
          if (print.table) {
            file.name <- paste('tab_', sample.id, '.dat', sep = "")
            write.table(adjacency.matrix, file = file.path('./sample-tabs/', file.name))
          }
          g <- graph.adjacency(adjacency.matrix, mode="directed", weighted=TRUE)
          file.name <- paste('gra_', sample.id, '.graphml', sep = "")
          write.graph(g, file = file.path('./sample-graphs/', file.name), format='graphml')
          # Terminate
          return(paste('Built graph for', sample.id, sep=' '))
        }
        stopCluster(par)
        if(gb$verbose) print(response)
        if(gb$verbose) print('All graphs built')
      }
    }
  )
  
  # Explicitely define GraphBuilder class
  class(gb) <- 'GraphBuilder'
  
  # Return the new GraphBuilder instance
  return(gb)
}

# Class to manage graphs
GraphManager <- function(clusters=0, verbose=FALSE) {
  
  # Define GraphManager attributes
  gm <- list(
    # Number of cluster for parallel computing
    clusters = clusters,
    
    # Boolean value to determine verbosity
    verbose = verbose,

    # GraphBuilder instance
    builder=GraphBuilder(clusters=clusters, verbose=verbose),

    mergeGraphs=function(graph.list) {
		# Merges multiple graphs.
		#
		# Args:
		#   graph.list: list of graphml files.

	    if(gm$verbose) print("Merging")
	    # Declare parallelism
	    cores <- makeCluster(gm$clusters)
	    registerDoParallel(cores)
	    # Get edges
	    edges <- foreach(i=1:length(graph.list), .combine=rbind) %dopar% {
	      library('igraph')
	      file.name <- graph.list[i]
	      get.edgelist(read.graph(file.path('.', file.name), format='graphml'))
	    }
	    stopCluster(cores)
	    # Build new graph
	    if(gm$verbose) print("Building total graph")
	    g <- graph.empty()
	    if(gm$verbose) print("Built empty graph")
	    if(gm$verbose) print("Adding vertices")
	    g <- g + vertices(unique(as.vector(edges)))
	    if(gm$verbose) print("Added vertices")
	    if(gm$verbose) print("Adding edges")
	    edges <- as.vector(t(edges))
	    g <- g + edges(edges)
	    if(gm$verbose) print("Added edges")
	    if(gm$verbose) print("Writing graph")
	    write.graph(g, file=file.path('.', 'total_graph.graphml'), format="graphml")
	    if(gm$verbose) print("Written")
	    if(gm$verbose) print('Graphs merged')
    }
  )
  
  # Explicitely define GraphManager class
  class(gm) <- c('GraphManager','GraphBuilder')
  
  # Return the new GraphManager instance
  return(gm)
}