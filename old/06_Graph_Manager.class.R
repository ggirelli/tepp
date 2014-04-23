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
    
    readData = function(file.path, header=TRUE, sep='\t', row.names=NULL, sample.column=NULL, abe.type='dummy') {
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
      if(gb$verbose) cat('Reading whole data\n')
      data <- read.table(file.path, header=header, sep=sep, row.names=row.names)

      # Save raw data into GraphBuilder instance
      gb$data <- data
      
      # Check for multi-sampling
      if(is.null(sample.column)) {

        # If single-sample file
        return(gb)

      } else {

        # If multi-sample file
        gb <- gb$splitData(data, sample.column, abe.type=abe.type)
        return(gb)

      }
    },
    
    splitData = function(data, sample.column, clusters=gb$clusters, abe.type='dummy') {
      # Splits multi-sample data
      #
      # Args:
      #   data: the multi-sample data to split.
      #	sample-column: the name of the sample-id column.
      #	clusters: the number .
      #
      # Returns:
      #	Modified GraphBuilder instance

      if(gb$verbose) cat('Splitting whole data\n')

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
      if(gb$verbose) cat('Data splitted\n')
      gb$split <- TRUE
      gb$isMultiSample <- TRUE
      gb$sample.column <- sample.column
      gb$data <- data
      return(gb)
    },

    build=function(gb, genes.label="Gene.id", clonality.label="clonality.status", clonal.val="clonal", subclonal.val="subclonal", abe.type="dummy") {
      # Builds adjacency matrix, homologous-clonality matrices and graph
      # It can discriminate between single- and multi-sample cases.
      #
      # Args:
      #   gb: current GraphBuilder instance.
      #   genes.label: label of the gene's id column.
      #   clonality.label: label of the clonality status column.
      #   clonal.val: value of clonality for clonal genes.
      #   subclonal.val: value of clonality for subclonal genes.
      #
      # Returns:
      # Modified GraphBuilder instance

      if(is.null(gb$isMultiSample)) {

        #--------------------#
        # Single sample data #
        #--------------------#

        # Make empty adjacency matrix
        adjacency.matrix <- matrix(
          0,
          nrow=length(gb$data[,genes.label]),
          ncol=length(gb$data[,genes.label])
        )
        rownames(adjacency.matrix) <- gb$data[,genes.label]
        colnames(adjacency.matrix) <- gb$data[,genes.label]

        # Distinguish clonal and subclonal
        genes.clonal <- which(gb$data[,clonality.label]==clonal.val)
        genes.subclonal <- which(gb$data[,clonality.label]==subclonal.val)

        # Make empty homologous-clonality matrices
        homo.clonal.matrix <- matrix(0, nrow=length(genes.clonal), ncol=length(genes.clonal))
        rownames(homo.clonal.matrix) <- gb$data[genes.clonal,genes.label]
        colnames(homo.clonal.matrix) <- gb$data[genes.clonal,genes.label]
        homo.subclonal.matrix <- matrix(0, nrow=length(genes.subclonal), ncol=length(genes.subclonal))
        rownames(homo.subclonal.matrix) <- gb$data[genes.subclonal,genes.label]
        colnames(homo.subclonal.matrix) <- gb$data[genes.subclonal,genes.label]

        if(length(genes.clonal) != 0) {
          # Fill matrices
          for(i in seq(length(genes.clonal))){

            if(length(genes.subclonal) != 0) {
              # Fill adjacency matrix
              for(j in seq(length(genes.subclonal))){
                adjacency.matrix[genes.clonal[i], genes.subclonal[j]] <- adjacency.matrix[genes.clonal[i], genes.subclonal[j]] + 1
              }
            }

            # Fill homologous.clonal
            for(j in seq(length(genes.clonal))) {
              # Empty main diagonal
              if(i != j) homo.clonal.matrix[i, j] <- homo.clonal.matrix[i, j] +1
            }
          }
        }

        if(length(genes.subclonal) != 0) {
          # Fill homologous.subclonal
          for(i in seq(length(genes.subclonal))){
            for(j in seq(length(genes.subclonal))) {
              # Empty main diagonal
              if(i != j) homo.subclonal.matrix[i, j] <- homo.subclonal.matrix[i, j] +1
            }
          }
        }

        # Write matrices
        file.name <- paste('adjTab_', sample.id, '-', abe.type, '.dat', sep = "")
        write.table(adjacency.matrix, file = file.path('.', file.name))
        file.name <- paste('homCloTab_', sample.id, '-', abe.type, '.dat', sep = "")
        write.table(homo.subclonal.matrix, file = file.path('.', file.name))
        file.name <- paste('homSubTab_', sample.id, '-', abe.type, '.dat', sep = "")
        write.table(homo.subclonal.matrix, file = file.path('.', file.name))

        # Build and write EvoPath graph
        g <- graph.adjacency(adjacency.matrix, mode="directed", weighted=TRUE)
        file.name <- paste('gra_', sample.id, '-', abe.type, '.graphml', sep = "")
        write.graph(g, file = file.path('.', file.name), format='graphml')
        # Build and write Cooc graphs
        #g <- graph.adjacency(homo.clonal.matrix, mode="undirected", weighted=TRUE)
        #file.name <- paste('graClo_', sample.id, '-', abe.type, '.graphml', sep = "")
        #write.graph(g, file = file.path('.', file.name), format='graphml')
        #g <- graph.adjacency(homo.subclonal.matrix, mode="undirected", weighted=TRUE)
        #file.name <- paste('graSub_', sample.id, '-', abe.type, '.graphml', sep = "")
        #write.graph(g, file = file.path('.', file.name), format='graphml')

        # Terminate
        if(gb$verbose) cat('Built graph\n')

      } else {

        #-------------------#
        # Multi sample data #
        #-------------------#

        # If needed, create output directory
        if (!file.exists(paste0('./sample-graphs-', abe.type, '/'))) dir.create(file.path(paste0('./sample-graphs-', abe.type, '/')))
        if (!file.exists(paste0('./sample-tables-', abe.type, '/'))) dir.create(file.path(paste0('./sample-tables-', abe.type, '/')))

        if(gb$verbose) cat('Building graphs\n')
        sample.list <- unique(gb$data[,gb$sample.column])

        # Declare parallelism
        par <- makeCluster(gb$clusters)
        registerDoParallel(par)
        # Make graph of every sample
        response <- foreach(i=1:length(sample.list)) %dopar% {
          library('igraph')
          
          # Working sample
          sample.id <- sample.list[i]

          # Read data
          data <- read.table(file.path(paste0('./sample-data-', abe.type, '/', sample.id)), header = TRUE)

          # Make empty adjacency matrix
          adjacency.matrix <- matrix(
            0,
            nrow=length(data[,genes.label]),
            ncol=length(data[,genes.label])
          )
          rownames(adjacency.matrix) <- data[,genes.label]
          colnames(adjacency.matrix) <- data[,genes.label]

          # Distinguish clonal and subclonal
          genes.clonal <- which(data[,clonality.label]==clonal.val)
          genes.subclonal <- which(data[,clonality.label]==subclonal.val)

          # Make empty homologous-clonality matrices
          homo.clonal.matrix <- matrix(0, nrow=length(genes.clonal), ncol=length(genes.clonal))
          rownames(homo.clonal.matrix) <- gb$data[genes.clonal,genes.label]
          colnames(homo.clonal.matrix) <- gb$data[genes.clonal,genes.label]
          homo.subclonal.matrix <- matrix(0, nrow=length(genes.subclonal), ncol=length(genes.subclonal))
          rownames(homo.subclonal.matrix) <- gb$data[genes.subclonal,genes.label]
          colnames(homo.subclonal.matrix) <- gb$data[genes.subclonal,genes.label]

          if(length(genes.clonal) != 0) {
            # Fill matrices
            for(i in seq(length(genes.clonal))){

              if(length(genes.subclonal) != 0) {
                # Fill adjacency matrix
                for(j in seq(length(genes.subclonal))) adjacency.matrix[genes.clonal[i], genes.subclonal[j]] <- adjacency.matrix[genes.clonal[i], genes.subclonal[j]] + 1
              }

              # Fill homologous.clonal
              for(j in seq(length(genes.clonal))) {
                # Empty main diagonal
                if(i != j) homo.clonal.matrix[i, j] <- homo.clonal.matrix[i, j] +1
              }
            }

            if(length(genes.subclonal) != 0) {
              # Fill homologous.subclonal
              for(i in seq(length(genes.subclonal))){
                for(j in seq(length(genes.subclonal))) {
                  # Empty main diagonal
                  if(i != j) homo.subclonal.matrix[i, j] <- homo.subclonal.matrix[i, j] +1
                }
              }
            }
          }
          
          # Write matrices
          file.name <- paste('adjTab_', sample.id, '.dat', sep = "")
          write.table(adjacency.matrix, file = file.path(paste0('./sample-tables-', abe.type, '/'), file.name))
          file.name <- paste('homCloTab_', sample.id, '.dat', sep = "")
          write.table(homo.clonal.matrix, file = file.path(paste0('./sample-tables-', abe.type, '/'), file.name))
          file.name <- paste('homSubTab_', sample.id, '.dat', sep = "")
          write.table(homo.subclonal.matrix, file = file.path(paste0('./sample-tables-', abe.type, '/'), file.name))
          # Build and write graph
          g <- graph.adjacency(adjacency.matrix, mode="directed", weighted=TRUE)
          file.name <- paste('gra_', sample.id, '.graphml', sep = "")
          write.graph(g, file = file.path(paste0('./sample-graphs-', abe.type, '/'), file.name), format='graphml')
          # Build and write Cooc graphs
          #g <- graph.adjacency(homo.clonal.matrix, mode="undirected", weighted=TRUE)
          #file.name <- paste('graClo_', sample.id, '.graphml', sep = "")
          #write.graph(g, file = file.path('./sample-graphs', abe.type, '/', file.name), format='graphml')
          #g <- graph.adjacency(homo.subclonal.matrix, mode="undirected", weighted=TRUE)
          #file.name <- paste('graSub_', sample.id, '.graphml', sep = "")
          #write.graph(g, file = file.path('./sample-graphs', abe.type, '/', file.name), format='graphml')

          # Terminate
          return(paste('* Built matrices and graph for', sample.id, sep=' '))

        }
        stopCluster(par)

        if(gb$verbose) for(answ in response) cat(answ,'\n')
        if(gb$verbose) cat('All matrices and graphs built\n')
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

    mergeGraphs.noAttr=function(graph.list) {
      # Merges multiple graphs.
      #
      # Args:
      #   graph.list: list of graphml files.

      if(gm$verbose) cat("Merging\n")

      # Declare parallelism
      cores <- makeCluster(gm$clusters)
      registerDoParallel(cores)
      # Get edges
      edges <- foreach(i=1:length(graph.list), .combine=rbind) %dopar% {
        library('igraph')

        file.name <- graph.list[i]
        g <- read.graph(file.path('.', file.name), format='graphml')
        cbind(get.edgelist(g), get.edge.attribute(g, 'weight'))

      }
      stopCluster(cores)

      # Build new graph
      if(gm$verbose) cat("Building empty graph\n")
      g <- graph.empty()

      # Vertices
      if(gm$verbose) cat("Adding vertices\n")
      g <- g + vertices(unique(as.vector(cbind(edges[,1], edges[,2]))))

      # Edges
      if(gm$verbose) cat("Adding edges\n")
      g <- g + edges(as.vector(t(cbind(edges[,1], edges[,2]))))
      if(gm$verbose) cat("Adding weights\n")
      E(g)$weight <- cbind(as.numeric(edges[,3]))
      g <- simplify(g, remove.multiple=TRUE, remove.loops=FALSE, edge.attr.comb=list(weight="sum", "ignore"))

      

      # Write graph
      if(gm$verbose) cat("Writing graph\n")
      write.graph(g, file=file.path('.', 'total_graph.graphml'), format="graphml")

      # Terminate
      if(gm$verbose) cat('Graphs merged\n')
    },

    graphsDistance.noAttr=function(g.one, g.two, attr.v.id='name') {
      # Calculates the distance between two graphs
      # Edges and vertices attributes are ignored during intersection
      #
      # Args:
      # g.one, g.two: the two graphs
      # attr.v.id: id vertex attribute
      #
      # Returns:
      # The distances as numeric vector

      # Vertices #
      v.merge <- c(eval(parse(text=paste0('V(g.one)$', attr.v.id))), eval(parse(text=paste0('V(g.two)$', attr.v.id))))
      v.merge <- v.merge[which(duplicated(v.merge) | duplicated(v.merge, fromLast=TRUE))]

      # Edges #
      e.merge <- rbind(get.edgelist(g.one),get.edgelist(g.two))
      e.merge.p <- paste(e.merge[,1], e.merge[,2], sep='-->')
      e.merge <- e.merge[which(duplicated(e.merge.p) | duplicated(e.merge.p, fromLast=TRUE)),]

      # Intersection #
      g.merge <- graph.empty()
      g.merge <- g.merge + vertices(v.merge)
      g.merge <- g.merge + edges(c(t(e.merge)))
      g.merge <- delete.vertices(g.merge, which(degree(g.merge, V(g.merge)) == 0))

      # First distance #
      d1 <- 1 - length(V(g.merge)) / max(length(V(g.one)), length(V(g.two)))

      # Second distance #
      d2 <- 1 - sum(degree(g.merge, V(g.merge))/2) / max(sum(degree(g.one, V(g.one))), sum(degree(g.two, V(g.two))))

      return(c(d1,d2))
    }
  )
  
  # Explicitely define GraphManager class
  class(gm) <- c('GraphManager','GraphBuilder')
  
  # Return the new GraphManager instance
  return(gm)
}