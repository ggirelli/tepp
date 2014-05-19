source('./extendigraph.R')

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
			if(vcount(g.one) > vcount(g.two)) {
				attrs <- get.vertex.attributes(V(g.two)[which(!(V(g.two) %in% V(g.one)))])
				s <- ''
				for(attr in colnames(attrs)) s <- paste0(s, ', ', attr, '=attrs[,\'', attr, '\']')
				print(s)
				eval(parse(text=paste0('g.one <- g.one + vertices((length(V(g.one))+1):(length(V(g.one))+length(attrs[,1]))', s, ')')))
				return(g.one)
			} else {
				attrs <- get.vertex.attributes(V(g.one)[which(!(V(g.one) %in% V(g.two)))])
				s <- ''
				for(attr in colnames(attrs)) s <- paste0(s, ', ', attr, '=attrs[,\'', attr, '\']')
				print(s)
				eval(parse(text=paste0('g.two <- g.two + vertices((length(V(g.two))+1):(length(V(g.two))+length(attrs[,1]))', s, ')')))
				return(g.two)
			}
		}

	)

	# Assign class attribute
	class(gm) <- 'GraphManager'

	# Return instantiaded Graph Manager
	return(gm)
}
