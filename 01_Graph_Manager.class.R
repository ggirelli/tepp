
source('./extendIgraph.R')

# Class to manage graphml graphs and perform graph operations
GraphManager <- function() {

	# Instantiate Graph Manager
	gm <- list(

		#------------#
		# ATTRIBUTES #
		#------------#



		#-----------#
		# FUNCTIONS #
		#-----------#
		
		# Transform

		undirected.noAttr = function(g) {
			# Transforms an DIRECTED graph into a UNDIRECTED one
			# Disregards edges/vertices attributes
			# 
			# Args:
			# 	g: undirected graph
			# 
			# Returns:
			# 	The UNDIRECTED graph
			
			if(!is.directed(g)) return(F)
			
			# Create undirected empty graph
			gf <- graph.empty(directed=T)

			# Add vertices
			gf <- gf + vertices(paste0(V(g)$name, '~IN'))
			gf <- gf + vertices(paste0(V(g)$name, '~OUT'))

			# Add edges
			el <- get.edgelist(g)
			gf <- gf + edges(c(t(cbind(paste0(el[,1], '~OUT'), paste0(el[,2], '~IN')))))

			# Remove 0-degree vertices
			gf <- delete.vertices(gf, V(gf)[which(degree(gf, V(gf)) == 0)])

			# Return undirected graph
			return(gf)
		},

		undirected = function(g) {
			# Transforms an DIRECTED graph into a UNDIRECTED one
			# Keeps edges/vertices attributes
			# 
			# Args:
			# 	g: undirected graph
			# 
			# Returns:
			# 	The UNDIRECTED graph
			
			if(!is.directed(g)) return(F)
			
			# Create undirected empty graph
			gf <- graph.empty(directed=T)

			# Add vertices
			gf <- gf + vertices(paste0(V(g)$name, '~IN'))
			gf <- gf + vertices(paste0(V(g)$name, '~OUT'))

			# Add vertices attributes
			attr.list <- list.vertex.attributes(g)
			for(attr.name in attr.list[which(attr.list != 'name')]) {
				eval(parse(text=paste0('V(gf)[1:length(V(g))]$', attr.name, ' <- V(g)$', attr.name)))
				eval(parse(text=paste0('V(gf)[length(V(g))+1:length(V(gf))]$', attr.name, ' <- V(g)$', attr.name)))
			}

			# Add edges
			el <- get.edgelist(g)
			gf <- gf + edges(c(t(cbind(paste0(el[,1], '~OUT'), paste0(el[,2], '~IN')))))

			# Add edges attributes
			attr.list <- list.edge.attributes(g)
			for(attr.name in attr.list[which(attr.list != 'name')]) eval(parse(text=paste0('E(gf)$', attr.name, ' <- E(g)$', attr.name)))

			# Remove 0-degree vertices
			gf <- delete.vertices(gf, V(gf)[which(degree(gf, V(gf)) == 0)])

			# Return undirected graph
			return(gf)
		},

		# Measures
		
		clusteringCoefficient = function(v, env, graph) {
			# Calculates clustering coefficient of a certain vertex in a given graph
			# 
			# Args:
			# 	v: vertex
			# 	env: vertex environment
			# 	graph: vertex-containing graph
			# 
			# Return:
			# 	Clustering coefficient

			# Prepare vertex
			class(v) <- 'igraph.vs'; attr(v, 'env') <- env

			# Retrieve neighbors
			neigh <- neighbors(graph, v, mode='all')

			# Return 0 if less than 2 neighbors
			if(length(neigh) < 2) return(0)

			# Retrieve subgraph
			sg <- as.undirected(induced.subgraph(graph, neigh))

			# Calculate clustering coefficient
			cc <- 2 * ecount(sg) / (vcount(sg) * (vcount(sg) - 1))

			# Terminate
			return(cc)
		},

		clusteringCoefficients = function(g) {
			# Calculates the clustering coefficient of the given graph
			# 
			# Args:
			# 	g: graph
			# 
			# Returns
			# 	The clustering coefficient

			# Calculates clustering coefficient for each node
			c.list <- sapply(V(g), FUN=function(v, env, graph) {
				return(gm$clusteringCoefficient(v, env, graph))
			}, env=attr(V(g), 'env'), graph=g)

			# Terminate
			return(mean(c.list))
		},

		# Compare

		calcHammingDist = function(g.one, g.two) {
			# Calculates the Hamming (edit) distance between two UNDIRECTED graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#
			# Returns:
			#	The Hamming distance H(g.one,g.two)

			# Get edges
			el.one <- get.edgelist(g.one)
			el.two <- get.edgelist(g.two)

			# Get number of common edges
			common <- length(intersect(paste0(el.one[,1], '~', el.one[,2]), paste0(el.two[,1], '~', el.two[,2])))
			common <-  common + length(intersect(paste0(el.one[,2], '~', el.one[,1]), paste0(el.two[,1], '~', el.two[,2])))

			# Not normalized distance
			dH.raw <- (length(el.one[,1]) + length(el.two[,1])) - (2 * common)

			# Normalize distance
			max.v <- max(length(V(g.one)), length(V(g.two)))
			K <- (max.v * (max.v - 1))
			K <- K / 2
			dH <- dH.raw / K

			# Return distance
			return(dH)
		},

		calcIpsenDist = function(g.one, g.two, gamma) {
			# Calculates the Ipsen-Mikhailov (spectral) distance between two UNDIRECTED graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#	gamma: parameter corresponding to the HWHM of the calculated Lorentz distributions
			#
			# Returns:
			#	The Ipsen-Mikhailov distance IM(g.one,g.two)
			
			# Read graphs
			gs <- list(g.one, g.two)

			# SpectralDensity function
			specDens = function(omega, omegadef, gamma) {
				k <- 0
				for(i in 2:length(omegadef)) k = k + (gamma / ((omega - omegadef[i])^2 + gamma^2))
				return(k)
			}
			# Normalization constant
			sdK = function(omegadef, gamma) {
				return(integrate(specDens, lower=0, upper=Inf, omegadef=omegadef, gamma=gamma, stop.on.error = FALSE)$value)
			}
			# Normalized spectral density (rho[omega])
			rhoO = function(omega, omegadef, gamma, k) {
				return(specDens(omega, omegadef, gamma) / k)
			}

			# Prepare graphs data
			gs.data <- list()
			for(i in 1:length(gs)) {
				g <- gs[[i]]
				g.data <- list()
				g.data$g <- g

				# Get adjacency matrix
				g.data$adj <- as.matrix(get.adjacency(g))

				# Build laplacian matrix
				g.data$lap <- - g.data$adj
				for(i in seq(length(g.data$adj[,1]))) {
					g.data$lap[i,i] <- g.data$lap[i,i] + degree(g, V(g)[i])
				}

				# Get 'defined' frequencies
				g.data$eigva <- round(sort(eigen(g.data$lap)$values),5)
				g.data$freqs <- sqrt(abs(g.data$eigva))

				# Calculate K
				g.data$k <- sdK(g.data$freqs, gamma)

				# Return graph data
				gs.data <- append(gs.data, list(g.data))
			}

			# IM distance
			sdDiffSq = function(omega, one.omegadef, one.k, two.omegadef, two.k, gamma) {
				return((rhoO(omega, one.omegadef, gamma, one.k) - rhoO(omega, two.omegadef, gamma, two.k))**2)
			}
			dIM <- sqrt(integrate(sdDiffSq, lower=0, upper=Inf, one.omegadef=gs.data[[1]]$freqs, one.k=gs.data[[1]]$k, two.omegadef=gs.data[[2]]$freqs, two.k=gs.data[[2]]$k, gamma=gamma)$value)

			# Return Ipsen-Mikhailov distance
			return(dIM)
		},

		calcHIMDist = function(g.one, g.two, gamma, xi) {
			# Calculates the Ipsen-Mikhailov (spectral) distance between two UNDIRECTED graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#	gamma: parameter corresponding to the HWHM of the calculated Lorentz distributions in IM distance calculation
			#	xi: parameter corresponding to the weight of dIM over dH in the final distance
			#
			# Returns:
			#	The Ipsen-Mikhailov distance IM(g.one,g.two)
			
			dH <- gm$calcHammingDist(g.one, g.two)
			dIM <- gm$calcIpsenDist(g.one, g.two, gamma)
			dHIM <- (1/sqrt(1+xi)) * sqrt(dH**2 + xi * dIM**2)
			return(dHIM)
		},

		# Operations
		
		merge = function(g.one, g.two) {
			# Identify the 'bigger' graph
			if(vcount(g.one) > vcount(g.two)) {

				#-------#
				# NODES #
				#-------#

				# Identify uncommon vertex between the two graphs
				uncommon.vertex <- which(!(V(g.two) %in% V(g.one)))
				if(length(uncommon.vertex) != 0) {
					# Retrieve uncommon vertex attributes list
					attrs <- get.vertex.attributes(V(g.two)[uncommon.vertex])

					# Prepare string for attributes assignment
					s <- ''; for(attr in colnames(attrs)) s <- paste0(s, ', ', attr, '=attrs[,\'', attr, '\']')
					
					# Assign attributes
					eval(parse(text=paste0('g.one <- g.one + vertices((length(V(g.one))+1):(length(V(g.one))+length(attrs[,1]))', s, ')')))
				}

				#-------#
				# EDGES #
				#-------#
				
				# Identify uncommon edges between the two graphs
				uncommon.edge <- which(!(E(g.two) %in% E(g.one)))
				if(length(uncommon.edge) == 0) return(g.one)

				if(length(list.edge.attributes(g.two)) != 0) {
					# Retrieve uncommon edge attributes list
					attrs <- get.edge.attributes(E(g.two)[uncommon.edge])

					# Prepare string for attributes assignment
					s <- ''; for(attr in colnames(attrs)) s <- paste0(s, ', ', attr, '=attrs[,\'', attr, '\']')
				} else {
					s <- ''
				}
				
				# Assign attributes
				edge.list <- c(t(get.edgelist(g.two)[uncommon.edge,]))
				eval(parse(text=paste0('g.one <- g.one + edges(edge.list', s, ')')))

				# Terminate
				return(g.one)
			} else {
				return(gm$merge(g.two, g.one))
			}
		},

		subtract = function(g.one, g.two) {
			# Removes edges and nodes that are common to both graphs from the first graph
			
			common.vertices <- which(V(g.one) %in% V(g.two))
			if(length(common.vertices) != 0 ) g.one <- g.one - vertices(V(g.one)[common.vertices])

			common.edges <- which(E(g.one) %in% E(g.two))
			if(length(common.edges) != 0) g.one <- delete.edges(g.one, E(g.one)[common.edges])

			return(g.one)
		},

		intersect = function(g.one, g.two) {
			# Intersects edges and nodes
			
			uncommon.vertices <- which(!(V(g.one) %in% V(g.two)))
			if(length(uncommon.vertices) != 0 ) g.one <- g.one - vertices(V(g.one)[uncommon.vertices])

			uncommon.edges <- which(!(E(g.one) %in% E(g.two)))
			if(length(uncommon.edges) != 0) g.one <- delete.edges(g.one, E(g.one)[uncommon.edges])

			return(g.one)
		}

	)

	# Assign class attribute
	class(gm) <- 'GraphManager'

	# Return instantiaded Graph Manager
	return(gm)
}
