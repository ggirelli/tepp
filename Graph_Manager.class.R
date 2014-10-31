library('igraph')
source('./extendIgraph.R')

# Class to manage graphml graphs and perform graph operations
GraphManager <- function() {

	# Instantiate Graph Manager
	gm <- list(

		#------------#
		# ATTRIBUTES #
		#------------#

		version = 2.1,

		#-----------#
		# FUNCTIONS #
		#-----------#
		
		# Transform
		#-----------------------------------

		undirected.noAttr = function(g) {
			# Transforms a DIRECTED graph into an UNDIRECTED one
			# Disregards edges/vertices attributes
			# 
			# Args:
			# 	g: undirected graph
			# 
			# Returns:
			# 	The UNDIRECTED graph
			
			if(!is.directed(g)) return(F)
			
			# Create undirected empty graph
			gf <- graph.empty(directed=F)

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
			# Transforms a DIRECTED graph into an UNDIRECTED one
			# Keeps edges/vertices attributes
			# 
			# Args:
			# 	g: undirected graph
			# 
			# Returns:
			# 	The UNDIRECTED graph
			
			if(!is.directed(g)) return(F)
			
			# Create undirected empty graph
			gf <- graph.empty(directed=F)

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

		rename.vertex.attr = function(g, old.name, new.name) {
			# Changes the name of a certain vertex attribute
			# 
			# Args:
			# 	g: graph
			# 	old.name
			# 	new.name
			# 
			# Returns:
			# 	Updated graph
			
			eval(parse(text=paste0('V(g)$', new.name, ' <- V(g)$', old.name)))
			g <- remove.vertex.attribute(g, old.name)
			return(g)
		},

		rename.vertex.attributes = function(g, map) {
			# Changes the name of a given set of vertex attribute
			# 
			# Args:
			# 	g: graph
			# 	old.name
			# 	new.name
			# 
			# Returns:
			# 	Updated graph
			
			for(old in names(map)) g <- gm$rename.vertex.attr(g, old, eval(parse(text=paste0('map$', old))))
			return(g)
		},

		rename.edge.attr = function(g, old.name, new.name) {
			# Changes the name of a certain edge attribute
			# 
			# Args:
			# 	g: graph
			# 	map: list with old.name => new.name
			# 
			# Returns:
			# 	Updated graph
			
			eval(parse(text=paste0('E(g)$', new.name, ' <- E(g)$', old.name)))
			g <- remove.edge.attribute(g, old.name)
			return(g)
		},

		rename.edge.attributes = function(g, map) {
			# Changes the name of a given set edge attribute
			# 
			# Args:
			# 	g: graph
			# 	map: list with old.name => new.name
			# 
			# Returns:
			# 	Updated graph
			
			for(old in names(map)) g <- gm$rename.edge.attr(g, old, eval(parse(text=paste0('map$', old))))
			return(g)
		},

		rename.attributes = function(g, vertex.map=list(), edge.map=list()) {
			# Changes the name of a given set of edge and/or vertex attributes
			# 
			# Args:
			# 	g: graph
			# 	vertex.map: list with old.name => new.name of vertex attributes
			# 	edge.map: list with old.name => new.name of edge attributes
			# 	
			# Returns:
			# 	Updated graph

			if(vertex.map != list()) g <- gm$rename.vertex.attributes(g, vertex.map)
			if(edge.map != list()) g <- gm$rename.edge.attributes(g, edge.map)
			return(g)
		},

		remove.bidirection = function(g, keepLoops=F) {
			# Removes bidirected edges from the graph and cleans it from any zero-degree vertex
			# 
			# Args:
			# 	g: graph
			# 	keepLoops: boolean, whether to keep loops
			# 
			# Returns:
			# 	Updated graph

			# identify bidirectionalities
			bidir.index <- which(duplicated(rbind(el, cbind(el[,2],el[,1]))))-ecount(g)

			# Remove bidirectionalities
			if(keepLoops) {
				g <- delete.edges(g, bidir.index(which(!is.loop(g, bidir.index))))
			} else {
				g <- delete.edges(g, bidir.index)
			}

			# Remove 0-degree nodes
			zero.index <- which(degree(g, V(g)) == 0)
			if(length(zero.index) != 0) {
				g <- delete.vertices(g, zero.index)
			}

			# Terminate
			return(g)
		},

		# Measures
		#-----------------------------------
		
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
		#-----------------------------------

		calcHammingDist = function(g.one, g.two) {
			# Calculates the Hamming (edit) distance between two UNDIRECTED graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#
			# Returns:
			#	The Hamming distance H(g.one,g.two)

			if (is.directed(g.one)) g.one <- gm$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- gm$undirected.noAttr(g.two)

			# Check if size is the same
			add.one <- V(g.two)[which(!(V(g.two)$name %in% V(g.one)$name))]
			add.two <- V(g.one)[which(!(V(g.one)$name %in% V(g.two)$name))]
			g.one <- add.vertices(g.one, length(add.one), attr=list(name=add.one$name))
			g.two <- add.vertices(g.two, length(add.two), attr=list(name=add.two$name))
			n <- length(V(g.one))

			# Get edges
			el.one <- get.edgelist(g.one)
			el.two <- get.edgelist(g.two)

			# Get number of common edges
			common <- length(intersect(paste0(el.one[,1], '~', el.one[,2]), paste0(el.two[,1], '~', el.two[,2])))
			common <- common + length(intersect(paste0(el.one[,2], '~', el.one[,1]), paste0(el.two[,1], '~', el.two[,2])))
			# Not normalized distance
			dH.raw <- (length(el.one[,1]) + length(el.two[,1])) - (2 * common)

			# Normalize distance
			max.v <- max(length(V(g.one)), length(V(g.two)))
			K <- (max.v * (max.v - 1))
			if(common != 0) K <- K / 2
			dH <- dH.raw / K

			# Return distance
			return(dH)
		},

		calcJaccardDist = function(g.one, g.two) {
			# Calculates the Jaccard (edit) distance between two UNDIRECTED graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#
			# Returns:
			#	The Jaccard distance J(g.one,g.two)

			if (is.directed(g.one)) g.one <- gm$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- gm$undirected.noAttr(g.two)

			# Get edges
			el.one <- get.edgelist(g.one)
			el.two <- get.edgelist(g.two)

			# Get number of common edges
			common <- length(intersect(paste0(el.one[,1], '~', el.one[,2]), paste0(el.two[,1], '~', el.two[,2])))
			common <- common + length(intersect(paste0(el.one[,2], '~', el.one[,1]), paste0(el.two[,1], '~', el.two[,2])))

			dJ <- ((length(el.one[,1]) + length(el.two[,1])) - 2 * common) / ((length(el.one[,1]) + length(el.two[,1])) - common)

			# Return distance
			return(dJ)
		},

		calcJaccardSubsetDist = function(g.one, g.two) {
			# Calculates the Jaccard Subset (edit) distance between two UNDIRECTED graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#
			# Returns:
			#	The Jaccard distance J(g.one,g.two)

			if (is.directed(g.one)) g.one <- gm$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- gm$undirected.noAttr(g.two)

			# Get edges
			el.one <- get.edgelist(g.one)
			el.two <- get.edgelist(g.two)

			# Get number of common edges
			common <- length(intersect(paste0(el.one[,1], '~', el.one[,2]), paste0(el.two[,1], '~', el.two[,2])))
			common <- common + length(intersect(paste0(el.one[,2], '~', el.one[,1]), paste0(el.two[,1], '~', el.two[,2])))

			# Calculate similarity and than distance
			JS <- common / min(length(el.one[,1]), length(el.two[,1]))
			dJS <- 1 - JS

			# Return distance
			return(dJS)
		},

		calcIpsenDist = function(g.one, g.two) {
			# Calculates the Ipsen-Mikhailov (spectral) distance between two UNDIRECTED graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#	gamma: parameter corresponding to the HWHM of the calculated Lorentz distributions
			#
			# Returns:
			#	The Ipsen-Mikhailov distance IM(g.one,g.two)
			
			if (is.directed(g.one)) g.one <- gm$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- gm$undirected.noAttr(g.two)

			# Read graphs
			gs <- list(g.one, g.two)

			# Check if size is the same
			add.one <- V(g.two)[which(!(V(g.two)$name %in% V(g.one)$name))]
			add.two <- V(g.one)[which(!(V(g.one)$name %in% V(g.two)$name))]
			g.one <- add.vertices(g.one, length(add.one), attr=list(name=add.one$name))
			g.two <- add.vertices(g.two, length(add.two), attr=list(name=add.two$name))
			n <- length(V(g.one))

			gammaBar = function(gamma, n) {
				return( sqrt(1/(gamma*pi)+(1/(2*gamma*(pi/2+atan(sqrt(n)/gamma))^2))*(pi/2+((gamma*sqrt(n))/((gamma^2)+n))+atan(sqrt(n)/gamma))+((-4*gamma)/((pi/2+atan(sqrt(n)/gamma))*pi*(4*(gamma^2)+n)))*(pi-(gamma/sqrt(n))*log((gamma^2)/((gamma^2)+n))+atan(sqrt(n)/gamma))) - 1)
			}
			gamma <- uniroot(gammaBar, c(1e-16,200), n=n, tol=1e-16)$root

			# SpectralDensity function
			specDens = function(omega, omegadef, gamma) {
				k <- 0
				for(i in 2:length(omegadef)) k <- k + (gamma / ((omega - omegadef[i])^2 + gamma^2))
				return(k)
			}
			# Normalization constant
			sdK = function(omegadef, gamma) {
				return(integrate(specDens, lower=0, upper=Inf, omegadef=omegadef, gamma=gamma, stop.on.error = FALSE, subdivisions=2000)$value)
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
			dIM <- sqrt(integrate(sdDiffSq, lower=0, upper=Inf, one.omegadef=gs.data[[1]]$freqs, one.k=gs.data[[1]]$k, two.omegadef=gs.data[[2]]$freqs, two.k=gs.data[[2]]$k, gamma=gamma, subdivisions=2000)$value)

			# Return Ipsen-Mikhailov distance
			return(dIM)
		},

		calcHIMDist = function(g.one, g.two, xi) {
			# Calculates the HIM between two UNDIRECTED graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#	xi: parameter corresponding to the weight of dIM over dH in the final distance
			#
			# Returns:
			#	The HIM distance (g.one,g.two)
			
			if (is.directed(g.one)) g.one <- gm$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- gm$undirected.noAttr(g.two)

			dH <- gm$calcHammingDist(g.one, g.two)
			dIM <- gm$calcIpsenDist(g.one, g.two)
			dHIM <- (1/sqrt(1+xi)) * sqrt(dH**2 + xi * dIM**2)
			return(dHIM)
		},

		calcJIMDist = function(g.one, g.two, xi) {
			# Calculates the JIM between two UNDIRECTED graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#	xi: parameter corresponding to the weight of dIM over dJ in the final distance
			#
			# Returns:
			#	The JIM distance (g.one,g.two)

			if (is.directed(g.one)) g.one <- gm$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- gm$undirected.noAttr(g.two)

			dJ <- gm$calcJaccardDist(g.one, g.two)
			dIM <- gm$calcIpsenDist(g.one, g.two)
			dJIM <- (1/sqrt(1+xi)) * sqrt(dJ**2 + xi * dIM**2)
			return(dJIM)
		},

		calcJIMSDist = function(g.one, g.two, xi) {
			# Calculates the JIMS between two UNDIRECTED graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#	xi: parameter corresponding to the weight of dIM over dJ in the final distance
			#
			# Returns:
			#	The JIM distance (g.one,g.two)

			if (is.directed(g.one)) g.one <- gm$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- gm$undirected.noAttr(g.two)

			dJS <- gm$calcJaccardSubsetDist(g.one, g.two)
			dIM <- gm$calcIpsenDist(g.one, g.two)
			dJIM <- (1/sqrt(1+xi)) * sqrt(dJS**2 + xi * dIM**2)
			return(dJIM)
		},

		calcDistances = function(g.one, g.two, xi) {
			# Calculates Hamming, Ipsen-Mikhailov and HIM distances between two undirected graphs
			# 
			# Args:
			# 	g.one: first graph
			# 	g.two: second graph
			# 	xi: parameter corresponding to the weight of dIM over dH in the final distance
			# 	
			# Returns:
			# 	A tuble containing respectiveli H, IM and HIM distances
			
			if (is.directed(g.one)) g.one <- gm$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- gm$undirected.noAttr(g.two)
			
			dH <- gm$calcHammingDist(g.one, g.two)
			dJ <- gm$calcJaccardDist(g.one, g.two)
			dJS <- gm$calcJaccardSubsetDist(g.one, g.two)
			dIM <- gm$calcIpsenDist(g.one, g.two)
			dHIM <- (1/sqrt(1+xi)) * sqrt(dH**2 + xi * dIM**2)
			dJIM <- (1/sqrt(1+xi)) * sqrt(dJ**2 + xi * dIM**2)
			dJIMS <- (1/sqrt(1+xi)) * sqrt(dJS**2 + xi * dIM**2)
			return(c(dH, dJ, dJS, dIM, dHIM, dJIM, dJIMS))
		},

		# Operations
		#-----------------------------------

		merge = function(g.one, g.two, vertex.key.label='name', vertex.attr.comb = getIgraphOpt("vertex.attr.comb"), edge.attr.comb = getIgraphOpt("edge.attr.comb")) {
			# Identify the 'bigger' graph
			if(vcount(g.one) > vcount(g.two)) {

				#-------#
				# NODES #
				#-------#

				if(length(list.vertex.attributes(g.two)) != 0) {
					# Retrieve vertex attributes list
					attrs <- get.vertex.attributes(V(g.two))

					# Prepare string for attributes assignment
					attr.names <- colnames(attrs)
					attr.names <- attr.names[which(attr.names != vertex.key.label)]
					attr.list <- list()
					s <- ''
					for(attr in attr.names) {
						attr.list[[attr]] <- get.vertex.attr(attr, V(g.two))
						s <- paste0(s, ', ', attr, '=attr.list$', attr)
					}
				} else {
					s <- ''
				}

				# Assign attributes
				eval(parse(text=paste0('g.one <- g.one + vertices(V(g.two)$', vertex.key.label, s, ')')))

				#-------#
				# EDGES #
				#-------#

				if(length(list.edge.attributes(g.two)) != 0) {
					# Retrieve edge attributes list
					attrs <- get.edge.attributes(E(g.two))

					# Prepare string for attributes assignment
					attr.list <- list()
					s <- ''
					for(attr in colnames(attrs)) {
						attr.list[[attr]] <- get.vertex.attr(attr, V(g.two))
						s <- paste0(s, ', ', attr, '=attr.list$', attr)
					}
				} else {
					s <- ''
				}
				
				# Assign attributes
				edge.list <- c(t(get.edgelist(g.two)))
				eval(parse(text=paste0('g.one <- g.one + edges(edge.list', s, ')')))

				#----------#
				# SIMPLIFY #
				#----------#
				
				vkeys <- eval(parse(text=paste0('V(g.one)$', vertex.key.label)))
				vkeys.uniq <- unique(eval(parse(text=paste0('V(g.one)$', vertex.key.label))))
				vids <- matrix(1:length(vkeys.uniq), nrow=1)
				colnames(vids) <- vkeys.uniq
				
				g.one <- contract.vertices(g.one, vids[,vkeys], vertex.attr.comb=vertex.attr.comb)
				g.one <- simplify(g.one, remove.loops=F, edge.attr.comb=edge.attr.comb)


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
		},

		contains = function(g.one, g.two) {
			# Verifies if g.one contains g.two
			
			if(length(which(!(V(g.two) %in% V(g.one)))) != 0) return(FALSE)
			if(length(which(!(E(g.two) %in% E(g.one)))) != 0) return(FALSE)

			return(TRUE)
		},

		contained = function(g.one, g.two) {
			# Verifies if g.two contains g.one
			return(gm$contains(g.two, g.one))
		}

	)

	# Assign class attribute
	class(gm) <- 'GraphManager'

	# Return instantiaded Graph Manager
	return(gm)
}
