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

		calcHammingDist = function(g.one, g.two) {
			# Calculates the Hamming (edit) distance between two graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#
			# Returns:
			#	The Hamming distance H(g.one,g.two)

			# Get edges
			el.one <- get.edgelist(g1)
			el.two <- get.edgelist(g2)

			# Get number of common edges
			common <- length(intersect(el.one, el.two))

			# Not normalized distance
			d <- (length(el.one) + length(el.two)) - (2 * common)

			# Normalize distance
			max.v <- max(length(V(g1)), length(V(g2)))
			d <- d / (max.v * (max.v - 1))

			# Return distance
			return d
		},

		calcIMDist = function(g.one, g.two) {
			# Calculates the Ipsen-Mikhailov (spectral) distance between two graphs
			#
			# Args:
			#	g.one: first graph
			#	g.two: second graph
			#
			# Returns:
			#	The Ipsen-Mikhailov distance IM(g.one,g.two)
		}

	)

	# Assign class attribute
	class(gm) <- 'GraphManager'

	# Return instantiaded Graph Manager
	return(gm)

}