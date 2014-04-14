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
    get.edgelist(g)

  }
  stopCluster(cores)

  # Build new graph
  if(gm$verbose) cat("Building empty graph\n")
  g <- graph.empty()

  # Vertices
  if(gm$verbose) cat("Adding vertices\n")
  g <- g + vertices(unique(as.vector(edges)))

  # Edges
  if(gm$verbose) cat("Adding edges\n")
  g <- g + edges(as.vector(t(edges)))
  if(gm$verbose) cat("Adding weights\n")
  E(g)$weight <- count.multiple(g, E(g))
  if(gm$verbose) cat("Remove multiple edges\n")
  g <- delete.edges(g, which(is.multiple(g, E(g))))

  # Write graph
  if(gm$verbose) cat("Writing graph\n")
  write.graph(g, file=file.path('.', 'total_graph.graphml'), format="graphml")

  # Terminate
  if(gm$verbose) cat('Graphs merged\n')
}