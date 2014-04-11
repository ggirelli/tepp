    build=function(gb, genes.label="Gene.id", clonality.label="clonality.status", clonal.val="clonal", subclonal.val="subclonal") {
      # Builds adjacency matrix, homologous-clonality matrices and graph
      # It can discriminate between single- and multi-sample cases.
      #
      # Args:
      #   gb: current GraphBuilder instance.
      #	  genes.label: label of the gene's id column.
      #	  clonality.label: label of the clonality status column.
      #   clonal.val: value of clonality for clonal genes.
      #   subclonal.val: value of clonality for subclonal genes.
      #
      # Returns:
      #	Modified GraphBuilder instance

      # If needed, create output directory
      if (!file.exists('./sample-graphs/')) dir.create(file.path('./sample-graphs/'))
      if (!file.exists('./sample-tables/')) dir.create(file.path('./sample-tables/'))

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

        # Fill matrices
        for(i in seq(length(genes.clonal))){
          # Fill adjacency matrix
          for(j in seq(length(genes.subclonal))){
            adjacency.matrix[genes.clonal[i], genes.subclonal[j]] <- adjacency.matrix[genes.clonal[i], genes.subclonal[j]] + 1
          }
          # Fill homologous.clonal
          for(j in seq(length(genes.clonal))) {
            # Empty main diagonal
            if(i != j) homo.clonal.matrix[genes.clonal[i], genes.clonal[j]] <- homo.clonal.matrix[genes.clonal[i], genes.clonal[j]] +1
          }
        }
        # Fill homologous.subclonal
        for(i in seq(length(genes.subclonal))){
          for(j in seq(length(genes.subclonal))) {
            # Empty main diagonal
            if(i != j) homo.clonal.matrix[genes.subclonal[i], genes.subclonal[j]] <- homo.clonal.matrix[genes.subclonal[i], genes.subclonal[j]] +1
          }
        }

        # Write matrices
        file.name <- paste('adjTab_', sample.id, '.dat', sep = "")
        write.table(adjacency.matrix, file = file.path('.', file.name))
        file.name <- paste('homCloTab_', sample.id, '.dat', sep = "")
        write.table(homo.subclonal.matrix, file = file.path('.', file.name))
        file.name <- paste('homSubTab_', sample.id, '.dat', sep = "")
        write.table(homo.subclonal.matrix, file = file.path('.', file.name))

        # Build and write graph
        g <- graph.adjacency(adjacency.matrix, mode="directed", weighted=TRUE)
        file.name <- paste('gra_', sample.id, '.graphml', sep = "")
        write.graph(g, file = file.path('.', file.name), format='graphml')

        # Terminate
        if(gb$verbose) cat('Built graph\n')

      } else {

        #-------------------#
        # Multi sample data #
        #-------------------#

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

          # Make empty homologous-clonality matrices
          homo.clonal.matrix <- matrix(0, nrow=length(genes.clonal), ncol=length(genes.clonal))
          rownames(homo.clonal.matrix) <- gb$data[genes.clonal,genes.label]
          colnames(homo.clonal.matrix) <- gb$data[genes.clonal,genes.label]
          homo.subclonal.matrix <- matrix(0, nrow=length(genes.subclonal), ncol=length(genes.subclonal))
          rownames(homo.subclonal.matrix) <- gb$data[genes.subclonal,genes.label]
          colnames(homo.subclonal.matrix) <- gb$data[genes.subclonal,genes.label]

          # Fill matrices
          for(i in seq(length(genes.clonal))){
            # Fill adjacency matrix
            for(j in seq(length(genes.subclonal))){
              adjacency.matrix[genes.clonal[i], genes.subclonal[j]] <- adjacency.matrix[genes.clonal[i], genes.subclonal[j]] + 1
            }
            # Fill homologous.clonal
            for(j in seq(length(genes.clonal))) {
              # Empty main diagonal
              if(i != j) homo.clonal.matrix[genes.clonal[i], genes.clonal[j]] <- homo.clonal.matrix[genes.clonal[i], genes.clonal[j]] +1
            }
          }
          # Fill homologous.subclonal
          for(i in seq(length(genes.subclonal))){
            for(j in seq(length(genes.subclonal))) {
              # Empty main diagonal
              if(i != j) homo.clonal.matrix[genes.subclonal[i], genes.subclonal[j]] <- homo.clonal.matrix[genes.subclonal[i], genes.subclonal[j]] +1
            }
          }

          # Write matrices
          file.name <- paste('adjTab_', sample.id, '.dat', sep = "")
          write.table(adjacency.matrix, file = file.path('./sample-tables/', file.name))
          file.name <- paste('homCloTab_', sample.id, '.dat', sep = "")
          write.table(homo.clonal.matrix, file = file.path('./sample-tables/', file.name))
          file.name <- paste('homSubab_', sample.id, '.dat', sep = "")
          write.table(homo.subclonal.matrix, file = file.path('./sample-tables/', file.name))

          # Build and write graph
          g <- graph.adjacency(adjacency.matrix, mode="directed", weighted=TRUE)
          file.name <- paste('gra_', sample.id, '.graphml', sep = "")
          write.graph(g, file = file.path('./sample-graphs/', file.name), format='graphml')

          # Terminate
          return(paste('* Built matrices and graph for', sample.id, sep=' '))

        }
        stopCluster(par)

        if(gb$verbose) for(answ in response) cat(answ,'\n')
        if(gb$verbose) cat('All matrices and graphs built\n')
      }
    }