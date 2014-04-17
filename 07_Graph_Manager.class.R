# Class to build and manage graphs
GraphManager <- function(clusters=0, verbose=FALSE) {
  library('igraph')
  library('doParallel')
  
  # Define GraphBuilder attributes
  gm <- list(
    # Number of cluster for parallel computing
    clusters = clusters,
    
    # Boolean value to determine verbosity
    verbose = verbose,

    # Sample list
    sammple.list = list(),
    
    # Reads data
    readData = function(file.path, header=TRUE, sep='\t', row.names=NULL, sample.column=NULL, abe.type='dummy', temp.sample.list=list()) {
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
      if(gm$verbose) cat('Reading', abe.type, 'data\n')
      data <- read.table(file.path, header=header, sep=sep, row.names=row.names)

      # Save raw data into GraphBuilder instance
      gm$data <- data
      
      # Check for multi-sampling
      if(is.null(sample.column)) {

        # If single-sample file
        return(gm)

      } else {

        # If multi-sample file
        gmnew <- gm$splitData(data, sample.column, abe.type=abe.type, temp.sample.list=temp.sample.list)
        return(gmnew)

      }
    },
    
    # Splits data
    splitData = function(data, sample.column, clusters=gm$clusters, abe.type='dummy', temp.sample.list=list()) {
      # Splits multi-sample data
      #
      # Args:
      #   data: the multi-sample data to split.
      #	sample-column: the name of the sample-id column.
      #	clusters: the number .
      #
      # Returns:
      #	Modified GraphBuilder instance
      
      if(gm$verbose) cat('Splitting', abe.type, 'data\n')

      # Prepare sample list without duplicates
      sample.list <- unique(data[,sample.column])

      # If needed, create output directory
      data.dir <- paste0('./sample-data-', abe.type, '/')
      if (!file.exists(data.dir)) {
        dir.create(file.path(data.dir))
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
        write.table(data[row.ids,], file=file.path(data.dir, sample.list[i]))

      }
      stopCluster(par)

      # Terminate
      if(gm$verbose) cat('Splitted', abe.type, 'data\n')
      gm$split <- TRUE
      gm$isMultiSample <- TRUE
      gm$sample.column <- sample.column
      gm$data <- data
      gm$sample.list <- unique(c(as.character(temp.sample.list), as.character(sample.list)))
      return(gm)
    },

    # Build SSMAs
    buildSSMA = function(sample.id, genes.label="Gene.id", clonality.label="clonality.status", clonal.val="clonal", subclonal.val="subclonal", abe.list=c('PM', 'Gain', 'Loss', 'RR')) {
      # Splits multi-sample data
      #
      # Args:
      #   sample.id: id of the sample to analyze
      #   genes.label: label of the gene's id column.
      #   clonality.label: label of the clonality status column.
      #   clonal.val: value of clonality for clonal genes.
      #   subclonal.val: value of clonality for subclonal genes.
      #   abe.list: list of aberration types
      #
      # Returns:
      #   None, prints graph in ./sample-graphs/

      # Read data
      data <- list()
      for(abe in abe.list) {
        f.name <- eval(parse(text=paste0('"sample-data-', abe, '/', sample.id, '"')))
        if(file.exists(f.name)) eval(parse(text=paste0('data$', abe, ' <- read.table(f.name , header=TRUE)')))
      }

      # Make empty graph
      g <- graph.empty(directed=TRUE)

      # Get clonals
      for(abe in abe.list) {

        genes <- eval(parse(text=paste0('data$', abe, '$', genes.label)))
        clonality <- eval(parse(text=paste0('data$', abe, '$', clonality.label)))
        aberration <- seq(length(genes))
        aberration[] <- abe

        # Add to graph
        g <- add.vertices(g, length(genes), attr=list(name=as.character(genes), clonality.status=as.character(clonality), abe.type=aberration))
      }
      
      # Remove uncertain.clonality
      g <- delete.vertices(g, V(g)[clonality.status == 'uncertain.clonal'])
      g <- delete.vertices(g, V(g)[clonality.status == 'uncertain.subclonal'])
      
      # Prepare edges list
      g <- g + edges(c(t(expand.grid(c(V(g)[clonality.status == 'clonal']), c(V(g)[clonality.status == 'subclonal'])))))
      E(g)$weight <- 1

      # If needed, create output directory
      data.dir <- paste0('./sample-graphs/')
      if (!file.exists(data.dir)) {
        dir.create(file.path(data.dir))
      }

      # Output graph
      write.graph(g, file.path(data.dir, paste0('gra_', sample.id, '.graphml')), format='graphml')
    },

    # Manage building
    build = function(abe.list, genes.label="Gene.id", clonality.label="clonality.status", clonal.val="clonal", subclonal.val="subclonal", sample.list=gm$sample.list) {
      # Splits multi-sample data
      #
      # Args:
      #   abe.list: list of aberration types
      #   genes.label: label of the gene's id column.
      #   clonality.label: label of the clonality status column.
      #   clonal.val: value of clonality for clonal genes.
      #   subclonal.val: value of clonality for subclonal genes.
      #   sample.list: list of samples to analyze
      #
      # Returns:
      #   None

      for(sample in sample.list) {
        cat('\nWorking on', sample)
        gm$buildSSMA(sample, genes.label=genes.label, clonality.label=clonality.label, clonal.val=clonal.val, subclonal.val=subclonal.val, abe.list=abe.list)
        cat('\nPrinted graph:', sample)
      }

      #gb$buildMSMA()
    }
  )
  
  # Explicitely define GraphBuilder class
  class(gm) <- 'GraphManager'
  
  # Return the new GraphBuilder instance
  return(gm)
}