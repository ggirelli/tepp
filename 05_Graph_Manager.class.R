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

      # If needed, create output directory
      if (!file.exists('./sample-graphs/')) dir.create(file.path('./sample-graphs/'))
      if (!file.exists('./sample-data/') && print.table) dir.create(file.path('./sample-tabs/'))

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

        # Build and write graph
        g <- graph.adjacency(adjacency.matrix, mode="directed", weighted=TRUE)
        file.name <- paste('gra_', sample.id, '.graphml', sep = "")
        write.graph(g, file = file.path('.', file.name), format='graphml')

        # Terminate
        if(gb$verbose) print('Built graph')

      } else {

        #-------------------#
        # Multi sample data #
        #-------------------#

        if(gb$verbose) print('Building graphs')
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
          data <- read.table(file.path('./sample-data/', sample.id), header = TRUE)

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

          # Build and write graph
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

    mergeGraphs.noattr=function(graph.list) {
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
        g <- read.graph(file.path('.', file.name), format='graphml')
        cbind(get.edgelist(g),get.edge.attribute(g, 'weight'))

	    }
	    stopCluster(cores)

	    # Build new graph
	    if(gm$verbose) print("Building empty graph")
	    g <- graph.empty()
	    if(gm$verbose) print("Adding vertices")
	    g <- g + vertices(unique(as.vector(edges)))

      # Get weights
      if(gm$verbose) print('Preparing edges weights')
      edgelist <- cbind(edges[,1],edges[,2])
      # Duplicated and unique edges
      multi <- which(duplicated(edgelist) | duplicated(edgelist, fromLast=TRUE))
      uni <- which(!(duplicated(edgelist) | duplicated(edgelist, fromLast=TRUE)))
      # Prepare edges array
      weights <- 1:length(edgelist[,1])
      weights[] <- 0
      # Assign weights to unique edges
      weights[uni] <- as.integer(edges[uni,3])
      # Prepare matrix with weights of duplicated edges
      cores <- makeCluster(gm$clusters)
      registerDoParallel(cores)
      t <- foreach(i=1:length(multi), .combine=rbind) %dopar% {

        c(
          multi[i],
          edges[multi[i],1],
          edges[multi[i],2],
          sum(as.integer(edges[which(
            edges[,1]==edges[multi[i],1] & edges[,2]==edges[multi[i],2]
          ),3])))

      }
      stopCluster(cores)

      # Assign weights to duplicated edges
      weights[as.integer(t[which(!duplicated(cbind(t[,2],t[,3]))),1])] <-
        as.integer(t[which(!duplicated(cbind(t[,2],t[,3]))),4])
      # Remove duplicated edges
      weights <- weights[which(weights>=1)]

      # Add edges
      if(gm$verbose) print("Adding edges")
      edgelist <- as.vector(t(unique(cbind(edges[,1],edges[,2]))))
      g <- g + edges(edgelist, weight=weights)

      # Write graph
	    if(gm$verbose) print("Writing graph")
	    write.graph(g, file=file.path('.', 'total_graph.graphml'), format="graphml")

      # Terminate
	    if(gm$verbose) print('Graphs merged')
    },

    mergeGraphCouple.jointAttr=function(
      g.one, g.two,
      attr.v.list=c(), attr.e.list=c(),
      attr.v.id='name',
      attr.v.action=list(), attr.e.action=list()
      ) {
      #
      # Merges two graphs with common edge/vertex attributes
      #
      # Args:
      #   g.one, g.two: the graphs to be merged
      #   attr.v.list: list of vertex attributes to keep in the merged graph
      #       each attribute MUST be present in both g.one and g.two
      #   attr.e.list: list of edge attributes to keep in the merged graph
      #       each attribute MUST be present in both g.one and g.two
      #   attr.v.id: attribute used to identify each vertex (default: 'name')
      #   attr.v.action: how to act with each vertex attribute (e.g.: 'sum')
      #         attr.v.action$ATTR <- ACTION_STRING
      #   attr.e.action: how to act with each edge attribute (e.g.: 'sum')
      #         attr.e.action$ATTR <- ACTION_STRING
      #
      # ACTION_STRING values:
      #   'sum': sum the values after applying as.numeric()
      #
      # Return:
      #    Merge graph
      
      # (1) if attributes in g.one and g.two are different, abort operation
      if(!identical(
        list.vertex.attributes(g.one),
        list.vertex.attributes(g.two)
      )) return(graph.empty()) # vertex
      if(!identical(
        list.edge.attributes(g.one),
        list.edge.attributes(g.two)
      )) return(graph.empty()) #edge

      ###############
      # Start merge #
      ###############
      g.merge <- g.one

      #----------#
      # Vertices #
      #----------#

      # Prepare final graph vertices id by id
      v.merge.id.list <- c()
      for(i in 1:length(V(g.merge))) eval(parse(text=paste(
        'v.merge.id.list <- append(v.merge.id.list,V(g.merge)[',
          i, ']$', attr.v.id, ')', sep=''
      )))

      # Scroll through g.two vertices
      for(i in 1:length(V(g.two))) {

        # Select vertex
        v.temp <- V(g.two)[i]
        eval(parse(text=paste('v.id <- V(g.two)[', i,']$',attr.v.id, sep='')))
        
        # If vertex is already present
        if(v.id %in% v.merge.id.list) {
          
          # Act on the attributes
          for(attr in attr.v.list) {

            # Retrieve the action associated to the attribute
            eval(parse(text=paste(
              'attr.action <- attr.v.action$', attr,
            sep='')))

            # Verify if there is an action associated to the attribute
            if(!is.null(attr.action) && identical(attr.action, 'sum')) {

              # Get g.two attribute value
              eval(parse(text=paste(
                'a.two <- V(g.two)[', i, ']$', attr,
              sep='')))

              # Get g.one attribute
              eval(parse(text=paste(
                'a.one <- V(g.one)[which(v.merge.id.list == v.id)]$',
                attr, sep=''
              )))

              # Assign attributes sum to g.merge vertex
              eval(parse(text=paste(
                'V(g.merge)[which(v.merge.id.list == v.id)]$', attr,
                ' <- as.numeric(a.one) + as.numeric(a.two)'
              , sep='')))
            }
          }
        } else {
          
          attrs <- ''
          for(attr in attr.v.list) {

            # Retrieve attribute value
            eval(parse(text=paste(
              'val <- V(g.two)[', i, ']$', attr, sep=''
            )))

            # Insert required comma
            if(!identical(attrs, '')) attrs <- paste(attrs, ',', sep='')

            # Prepare attributes assignment string
            attrs <- paste(attrs, attr, '=\'', val, '\'', sep=' ')
          }
          
          if(!identical(attrs, '')) attrs <- paste(',', attrs, sep=' ')
          
          # Add all the vertices with their attributes
          
          eval(parse(text=paste(
            'g.merge <- g.merge + vertex(\'', v.id, '\'', attrs, ')',
          sep='')))
        }
      }

      #------#
      # Edge #
      #------#

      # Prepare final graph edges
      e.merge.list <- get.edgelist(g.merge)
      e.merge.list <- paste(e.merge.list[,1], e.merge.list[,2], sep=' <-- ')

      # Scroll through g.two edges
      for(i in 1:length(get.edgelist(g.two)[,1])) {

        # Select edge
        e.temp <- get.edgelist(g.two)[i,]

        # If edge is already present
        if(paste(e.temp[1], e.temp[2], sep=' <-- ') %in% e.merge.list) {

          # Act on the attributes
          for(attr in attr.e.list) {

            # Retrieve the action associated to the attribute
            eval(parse(text=paste(
              'attr.action <- attr.e.action$', attr,
            sep='')))

            # Verify if there is an action associated to the attribute
            if(!is.null(attr.action) && identical(attr.action, 'sum')) {

              # Get g.two attribute value
              eval(parse(text=paste(
                'a.two <- E(g.two)[', i, ']$', attr,
              sep='')))

              # Get g.one attribute
              eval(parse(text=paste(
                'a.one <- E(g.one)[which(e.merge.list == paste(',
                  'e.temp[1], e.temp[2], sep=\' <-- \'))]$', attr,
              sep='')))

              # Assign attributes sum to g.merge vertex
              eval(parse(text=paste(
                'E(g.merge)[which(e.merge.list == paste(',
                  'e.temp[1], e.temp[2], sep=\' <-- \'))]$', attr,
                ' <- as.numeric(a.one) + as.numeric(a.two)'
              , sep='')))
            }
          }
        } else {
          attrs <- ''
          for(attr in attr.e.list) {

            # Retrieve attribute value
            eval(parse(text=paste(
              'val <- E(g.two)[', i, ']$', attr, sep=''
            )))

            # Insert required comma
            if(!identical(attrs, '')) attrs <- paste(attrs, ',', sep='')

            # Prepare attributes assignment string
            attrs <- paste(attrs, attr, '=', val, sep=' ')
          }
          if(!identical(attrs, '')) attrs <- paste(',', attrs, sep=' ')

          # Add all the edges with their attributes
          cmd <- paste(
            'g.merge <- g.merge + edge(c(\'',
            e.temp[1], '\',\'', e.temp[2],
            '\')', attrs, ')',
          sep='')
          eval(parse(text=cmd))

        }
      }

      #-----------#
      # Terminate #
      #-----------#
      return(g.merge)
    }
  )
  
  # Explicitely define GraphManager class
  class(gm) <- c('GraphManager','GraphBuilder')
  
  # Return the new GraphManager instance
  return(gm)
}