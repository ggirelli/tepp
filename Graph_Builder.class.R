# Class to build and manage graphs
GraphBuilder <- function(clusters=0, verbose=FALSE, genes.label="Gene.id", white.list=list(), black.list=list(), clonal.val=c('clonal'), subclonal.val=c('subclonal'), attr.table='', clean=FALSE) {
  library('igraph')
  library('doParallel')
  
  # Define GraphBuilder attributes
  gb <- list(
    
    #------------------#
    # INPUT PARAMETERS #
    #------------------#

    # Number of cluster for parallel computing
    clusters = clusters,
    
    # Boolean value to determine verbosity
    verbose = verbose,

    # Gene.id label
    genes.label = genes.label,

    # Gene white/black-list
    white.list = white.list,
    black.list = black.list,

    # Values of clonality.status to identify clonal/subclonal genes
    clonal.val = clonal.val,
    subclonal.val = subclonal.val,

    # Clean: keeps only whitelisted genes
    clean = clean,

    # Sample list
    sammple.list = list(),

    # Vertex attribute table
    attr.table = attr.table,
    
    #-----------#
    # FUNCTIONS #
    #-----------#

    # Reads data
    readData = function(file.path, header=TRUE, sep='\t', row.names=NULL, sample.column=NULL, genes.label=gb$genes.label, white.list=gb$white.list, black.list=gb$black.list, abe.type='dummy', temp.sample.list=list(), clean=gb$clean) {
      #
      # Args:
      #   file.path: the name of the file which the data are to be
      #   read from or a a readable text-mode connection.
      # header: a logical value indicating whether the file contains
      #   the names of the variables as its first line.
      # sep: the field separator character.
      # row.names: a vector of row names.
      # sample.column: if file is multi-sample, specify sample
      #   column header.
      #
      # Returns:
      #   Modified GraphBuilder instance
      
      # Read data
      if(gb$verbose) cat('Reading', abe.type, 'data\n')
      data <- read.table(file.path, header=header, sep=sep, row.names=row.names)

      # Save raw data into GraphBuilder instance
      gb$data <- data
      
      # Check for multi-sampling
      if( is.null(sample.column)) {

        # Clean single-sample data
        if(clean && length(white.list) != 0) {
          for(i in seq(length(data[,1]))) {
            id <- eval(parse(text=paste0('data$', genes.label, '[i]')))
            if(!(id %in% white.list)) data <- data[-i]
          }
        }
        # Blacklisting
        if(length(blac.list) != 0) {
          for(i in seq(length(data[,1]))) {
            id <- eval(parse(text=paste0('data$', genes.label, '[i]')))
            if(id %in% white.list) data <- data[-i]
          }
        }

        # If single-sample file
        return(gb)

      } else {

        # If multi-sample file
        gbnew <- gb$splitData(data, sample.column, genes.label=genes.label, white.list=white.list, black.list=black.list, abe.type=abe.type, temp.sample.list=temp.sample.list, clean=clean)
        return(gbnew)

      }
    },
    
    # Splits data
    splitData = function(data, sample.column, clusters=gb$clusters, genes.label=gb$genes.label, white.list=gb$white.list, black.list=gb$black.list, abe.type='dummy', temp.sample.list=list(), clean=gb$clean) {
      # Splits multi-sample data
      #
      # Args:
      #   data: the multi-sample data to split.
      # sample-column: the name of the sample-id column.
      # clusters: the number .
      #
      # Returns:
      # Modified GraphBuilder instance
      
      if(gb$verbose) cat('Splitting', abe.type, 'data\n')
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
        row.ids <- which(data[,sample.column] == sample.id)
        # Clean sample
        if(clean && length(white.list) != 0) row.ids <- intersect(row.ids, which(data[,genes.label] %in% white.list))
        # Blacklisting
        if(length(black.list) != 0) row.ids <- intersect(row.ids, which(!(data[,genes.label] %in% black.list)))
        # Write selected rows
        if(length(row.ids) != 0) write.table(data[row.ids,], file=file.path(data.dir, sample.list[i]))

      }
      stopCluster(par)

      # Terminate
      if(gb$verbose) cat('Splitted', abe.type, 'data\n')
      gb$split <- TRUE
      gb$isMultiSample <- TRUE
      gb$sample.column <- sample.column
      gb$data <- data
      gb$sample.list <- unique(c(as.character(temp.sample.list), as.character(sample.list)))
      return(gb)
    },

    buildFinalSSMA = function(abe.list, genes.label=gb$genes.label, clonality.label="clonality.status", clonal.val=gb$clonal.val, subclonal.val=gb$subclonal.val, sample.list=gb$sample.list) {
      # Build SSMAs for final MSMA
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

      # Declare parallelism
      par <- makeCluster(clusters)
      registerDoParallel(par)
      # Execute buildSSMA for each sample
      foreach(sample.id=sample.list) %dopar% {
        library('igraph')
        
        # Read data
        data <- list()
        for(abe in abe.list) {
          f.name <- eval(parse(text=paste0('"sample-data-', abe, '/', sample.id, '"')))
          if(file.exists(f.name)) eval(parse(text=paste0('data$', abe, ' <- read.table(f.name , header=TRUE, sep=" ")')))
        }

        if(length(data) != 0) {
          # Make empty graph
          g <- graph.empty(directed=TRUE)
          g.clonal <- g
          g.subclonal <- g

          # Get clonals
          for(abe in abe.list) {

            genes <- eval(parse(text=paste0('data$', abe, '$', genes.label)))
            clonality <- eval(parse(text=paste0('data$', abe, '$', clonality.label)))
            aberration <- seq(length(genes))
            aberration[] <- tolower(abe)

            # Add to graph
            g <- add.vertices(g, length(genes), attr=list(name=as.character(genes), clonality.status=as.character(clonality), abe.type=aberration))
          }

          # Remove vertices that are neither clonal nor subclonal
          g <- delete.vertices(g, V(g)[intersect(which(!(levels(V(g)$clonality.status) %in% clonal.val)), which(!(levels(V(g)$clonality.status) %in% subclonal.val)))])

          # Prepare edges list
          g <- g + edges(c(t(expand.grid(c(V(g)[clonality.status %in% clonal.val]), c(V(g)[clonality.status %in% subclonal.val])))))
          E(g)$weight <- 1

          # If needed, create output directory
          data.dir <- paste0('./sample-graphs/')
          if (!file.exists(data.dir)) {
            dir.create(file.path(data.dir))
          }

          # Output graph
          write.graph(g, file.path(data.dir, paste0('gra_', sample.id, '.graphml')), format='graphml')
        }
      }
      stopCluster(par)
    },

    buildClonalSSMA = function(abe.list, genes.label=gb$genes.label, clonality.label="clonality.status", clonal.val=gb$clonal.val, subclonal.val=gb$subclonal.val, sample.list=gb$sample.list, v.list=list()) {
      # Build SSMAs for clonality co-occurrency MSMAs
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

      # Declare parallelism
      par <- makeCluster(clusters)
      registerDoParallel(par)
      # Execute buildSSMA for each sample
      foreach(sample.id=sample.list) %dopar% {
        library('igraph')
        
        # Read data
        data <- list()
        for(abe in abe.list) {
          f.name <- eval(parse(text=paste0('"sample-data-', abe, '/', sample.id, '"')))
          if(file.exists(f.name)) {
            temp.data <- read.table(f.name)
            temp.data <- temp.data[which(paste(as.character(eval(parse(text=paste0('temp.data$', genes.label)))), tolower(abe), sep='~') %in% v.list),]
            eval(parse(text=paste0('data$', abe, ' <- temp.data')))
          }
        }

        if(length(data) != 0) {
          # Make empty graph
          g.clonal <- graph.empty(directed=TRUE)
          g.subclonal <- graph.empty(directed=TRUE)
          g.nonclonal <- graph.empty(directed=TRUE)

          # Get clonals
          for(abe in abe.list) {

            genes <- eval(parse(text=paste0('data$', abe, '$', genes.label)))
            clonality <- as.character(eval(parse(text=paste0('data$', abe, '$', clonality.label))))
            aberration <- seq(length(genes))
            aberration[] <- tolower(abe)

            # Add to graph
            g.clonal <- add.vertices(g.clonal, length(genes[which(clonality %in% clonal.val)]), attr=list(name=as.character(genes[which(clonality %in% clonal.val)]), abe.type=aberration[which(clonality %in% clonal.val)]))
            g.subclonal <- add.vertices(g.subclonal, length(genes[which(clonality %in% subclonal.val)]), attr=list(name=as.character(genes[which(clonality %in% subclonal.val)]), abe.type=aberration[which(clonality %in% subclonal.val)]))
            g.nonclonal <- add.vertices(g.nonclonal, length(genes), attr=list(name=as.character(genes), abe.type=aberration, clonality=clonality))
          }

          # Prepare edges list for clonals
          g.clonal <- g.clonal + edges(c(t(expand.grid(V(g.clonal)$name, V(g.clonal)$name))))
          E(g.clonal)$weight <- 1

          # Prepare edges list for subclonals
          g.subclonal <- g.subclonal + edges(c(t(expand.grid(V(g.subclonal)$name, V(g.subclonal)$name))))
          E(g.subclonal)$weight <- 1

          # Prepare edges list for nonclonals
          g.nonclonal <- g.nonclonal + edges(c(t(expand.grid(V(g.nonclonal)[!(clonality %in% append(clonal.val, subclonal.val))]$name, V(g.nonclonal)$name))))
          g.nonclonal <- g.nonclonal + edges(c(t(expand.grid(V(g.nonclonal)$name, V(g.nonclonal)[!(clonality %in% append(clonal.val, subclonal.val))]$name))))
          E(g.nonclonal)$weight <- 1

          # If needed, create output directory
          data.dir <- paste0('./sample-graphs/')
          if (!file.exists(data.dir)) {
            dir.create(file.path(data.dir))
          }

          # Output graph
          if(length(V(g.clonal)) != 0) write.graph(g.clonal, file.path(data.dir, paste0('gra_clonal_', sample.id, '.graphml')), format='graphml')
          if(length(V(g.subclonal)) != 0) write.graph(g.subclonal, file.path(data.dir, paste0('gra_subclonal_', sample.id, '.graphml')), format='graphml')
          if(length(V(g.nonclonal)) != 0) write.graph(g.nonclonal, file.path(data.dir, paste0('gra_nonclonal_', sample.id, '.graphml')), format='graphml')
        }
      }
      stopCluster(par)
    },

    # Build MSMA
    buildMSMA = function(graph.list, directed=TRUE) {
      # Merges SSMAs into a single MSMA summing the edge's weight
      #
      # Args:
      #   graph.list:
      #
      # Return:
      #

      # Declare parallelism
      cores <- makeCluster(gb$clusters)
      registerDoParallel(cores)
      # Get edges
      edges <- foreach(i=1:length(graph.list), .combine=rbind) %dopar% {
        library('igraph')

        file.name <- graph.list[i]
        if(file.exists(file.path('./sample-graphs/', file.name))) {
          g <- read.graph(file.path('./sample-graphs/', file.name), format='graphml')
          edgelist <- get.edgelist(g)
          if(length(edgelist) != 0) {
            edgelist.n <- get.edgelist(g, name=FALSE)
            cbind(edgelist, E(g)$weight, V(g)[edgelist.n[,1]]$abe.type, V(g)[edgelist.n[,2]]$abe.type)
          }
        }
      }
      stopCluster(cores)

      if(length(edges) != 0) {

        # Build new graph
        if(gb$verbose) cat("Building empty graph\n")
        g <- graph.empty(directed=directed)

        # Vertices
        #(2) Then scroll the table (tricky) to get each node and its attributes and add them to the MSMA.
        if(gb$verbose) cat("Preparing sources\n")
        source.table <- paste(edges[,1], edges[,4], sep='~')
        if(gb$verbose) cat("Preparing targets\n")
        target.table <- paste(edges[,2], edges[,5], sep='~')
        if(gb$verbose) cat("Preparing vertices\n")
        vertices.table <- unique(c(source.table, target.table))
        if(gb$verbose) cat("Adding vertices with attributes\n")
        g <- add.vertices(g, nv=length(vertices.table), attr=list(name=vertices.table, aberration=matrix(unlist(strsplit(vertices.table, '~')), nrow=2)[2,], HUGO=matrix(unlist(strsplit(vertices.table, '~')), nrow=2)[1,]))

        # Edges
        if(gb$verbose) cat("Adding edges\n")
        #(3) Then add all the edges identifying the nodes based on their attributes.
        g <- g + edges(c(t(cbind(source.table, target.table))))
        if(gb$verbose) cat("Adding edges' weight\n")
        E(g)$weight <- as.numeric(edges[,3])

        #(4) Finally apply 'simplify' to remove multiple edges and sum their weights.
        if(gb$verbose) cat("Simplifying\n")
        g <- simplify(g, remove.multiple=TRUE, remove.loops=FALSE, edge.attr.comb=list(weight="sum", "ignore"))


        if(gb$verbose) cat('\nMaximum edge weight: ', max(E(g)$weight), '\n')
        if(gb$verbose) cat('Maximum vertex degree: ', max(degree(g, V(g))), '\n')

        # Terminate
        if(gb$verbose) cat('Graphs merged\n')

        return(g)

      } else {

        return(graph.empty())

      }
    },

    setVAttributes = function(graph, attr.table, key.label="HUGO") {
      # Adds vertex attributes to a graph object
      #
      # Args:
      #   graph: the graph
      #   attr.file: the path to the file
      #   key.label: the column/vertex label used to assign the attributes
      #
      # Returns:
      #   The graph with vertex attributes

      # Get key.attribute values
      sources <- eval(parse(text=paste0('attr.table$', key.label)))
      targets <- unique(eval(parse(text=paste0('V(graph)$', key.label))))

      # Iterate through target key.values (HUGO)
      for(key in targets[which(targets %in% sources)]) {
        # Iterate through new attributes to assign
        for(attr in colnames(attr.table)[which(colnames(attr.table) != key.label)]) {
          eval(parse(text=paste0('V(graph)[', key.label, ' == key]$', attr, ' <- attr.table$', attr, '[which(sources == key)]')))
        }
      }

      return(graph)
    },

    # Manage building
    build = function(abe.list, genes.label=gb$genes.label, clonality.label="clonality.status", clonal.val=gb$clonal.val, subclonal.val=gb$subclonal.val, sample.list=gb$sample.list, attr.table=gb$attr.table) {
      # Builds graphs after data read/split
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

      if(gb$verbose) cat('# Preparing SSMAs.\n')
      gb$buildFinalSSMA(abe.list, genes.label=genes.label, clonality.label=clonality.label, clonal.val=clonal.val, subclonal.val=subclonal.val, sample.list=sample.list)
      if(gb$verbose) cat('SSMAs prepared.\n')

      if(gb$verbose) cat("\n# Merging SSMAs into MSMA\n")
      g.total <- gb$buildMSMA(paste0('gra_', sample.list, '.graphml'))

      if(gb$verbose) cat('\n# Preparing Clonality SSMAs.\n')
      gb$buildClonalSSMA(abe.list, genes.label=genes.label, clonality.label=clonality.label, clonal.val=clonal.val, subclonal.val=subclonal.val, sample.list=sample.list, v.list=V(g.total)$name)
      if(gb$verbose) cat('SSMAs prepared.\n')

      if(gb$verbose) cat("\n# Merging SSMAs into MSMA · Clonal co-occurrency\n")
      g.clonal <- gb$buildMSMA(paste0('gra_clonal_', sample.list, '.graphml'))
      if(gb$verbose) cat("\n# Merging SSMAs into MSMA · Subclonal co-occurrency\n")
      g.subclonal <- gb$buildMSMA(paste0('gra_subclonal_', sample.list, '.graphml'))
      if(gb$verbose) cat("\n# Merging SSMAs into MSMA · Uncertain_clonality co-occurrency\n")
      g.nonclonal <- gb$buildMSMA(paste0('gra_nonclonal_', sample.list, '.graphml'))

      if(gb$verbose) cat("\n# Retrieving co-occurrency data\n")
      # Get which edges can be present in the co-occurrency graphs
      el.tot <- get.edgelist(g.total)
      e.in.clo <- (el.tot[,1] %in% V(g.clonal)$name & el.tot[,2] %in% V(g.clonal)$name)
      e.in.sub <- (el.tot[,1] %in% V(g.subclonal)$name & el.tot[,2] %in% V(g.subclonal)$name)
      e.in.non <- (el.tot[,1] %in% V(g.nonclonal)$name & el.tot[,2] %in% V(g.nonclonal)$name)
      # Retrieve correct clonal co-occurrency data
      if(length(V(g.clonal)) != 0) {
        if(gb$verbose) cat(" · Retrieving clonal co-occurrency data\n")
        clonal.cooc <- seq(length(E(g.total)))
        clonal.cooc[] <- 0
        clonal.cooc.ids <- get.edge.ids(g.clonal, t(get.edgelist(g.total)[which(e.in.clo),]), error=FALSE)
        clonal.cooc[which(e.in.clo)[which(clonal.cooc.ids != 0)]] <- E(g.clonal)[clonal.cooc.ids]$weight
        E(g.total)$clonal.cooc <- clonal.cooc
      }
      # Retrieve correct subclonal co-occurrency data
      if(length(V(g.subclonal)) != 0) {
        if(gb$verbose) cat(" · Retrieving subclonal co-occurrency data\n")
        subclonal.cooc <- seq(length(E(g.total)))
        subclonal.cooc[] <- 0
        subclonal.cooc.ids <- get.edge.ids(g.subclonal, t(get.edgelist(g.total)[which(e.in.sub),]), error=FALSE)
        subclonal.cooc[which(e.in.sub)[which(subclonal.cooc.ids != 0)]] <- E(g.subclonal)[subclonal.cooc.ids]$weight
        E(g.total)$subclonal.cooc <- subclonal.cooc
      }
      # Retrieve correct nonclonal co-occurrency data
      if(length(V(g.nonclonal)) != 0) {
        if(gb$verbose) cat(" · Retrieving uncertain_clonality co-occurrency data\n")
        nonclonal.cooc <- seq(length(E(g.total)))
        nonclonal.cooc[] <- 0
        nonclonal.cooc.ids <- get.edge.ids(g.nonclonal, t(get.edgelist(g.total)[which(e.in.non),]), error=FALSE)
        nonclonal.cooc[which(e.in.non)[which(nonclonal.cooc.ids != 0)]] <- E(g.nonclonal)[nonclonal.cooc.ids]$weight
        E(g.total)$nonclonal.cooc <- nonclonal.cooc
      }

      # Retrieve and assign new vertex attributes
      if(length(attr.table) != 0 && attr.table != '') {
        if(gb$verbose) cat("\nAssigning new vertex attributes based on HUGO\n")
        g.total <- gb$setVAttributes(g.total, attr.table)
      }

      if(gb$verbose) cat("\nWriting graph\n")
      write.graph(g.total, file.path('.', 'total_graph.graphml'), format='graphml')
      #write.graph(g.clonal, file.path('.', 'clonal_graph.graphml'), format='graphml')
      #write.graph(g.subclonal, file.path('.', 'subclonal_graph.graphml'), format='graphml')
      #write.graph(g.nonclonal, file.path('.', 'nonclonal_graph.graphml'), format='graphml')

      if(gb$verbose) cat("\nFIN\n\n")
    }
  )
  
  # Explicitely define GraphBuilder class
  class(gb) <- 'GraphBuilder'
  
  # Return the new GraphBuilder instance
  return(gb)
}