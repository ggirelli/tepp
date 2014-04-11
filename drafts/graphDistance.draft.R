graphsDistance.noAttr <- function(g.one, g.two, attr.v.id='name') {
	#
	# Calculates the distance between two graphs
	# Edges and vertices attributes are ignored during intersection
	#
	# Args:
	#	g.one, g.two: the two graphs
	#	attr.v.id: id vertex attribute
	#
	# Returns:
	#	The distances as numeric vector
	#

	# Vertices #
	v.merge <- c(eval(parse(text=paste('V(g.one)$', attr.v.id, sep=''))), eval(parse(text=paste('V(g.two)$', attr.v.id, sep=''))))
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

library('igraph')
g <- graph.empty()
g <- g + vertex('a','b','c','e','g','h')
g <- g + edges(c('a','c'))
g <- g + edges(c('b','c'))
g <- g + edges(c('b','g'))
g <- g + edges(c('e','b'))
g <- g + edges(c('e','c'))
g <- g + edges(c('g','c'))
g <- g + edges(c('g','h'))
g <- g + edges(c('h','b'))
g <- g + edges(c('h','c'))
g <- g + edges(c('h','e'))
h <- graph.empty()
h <- h + vertex('a','c','d','e','f','g','i','l')
h <- h + edges(c('a','f'))
h <- h + edges(c('a','g'))
h <- h + edges(c('c','f'))
h <- h + edges(c('d','c'))
h <- h + edges(c('d','f'))
h <- h + edges(c('e','c'))
h <- h + edges(c('e','f'))
h <- h + edges(c('g','c'))
h <- h + edges(c('i','c'))
h <- h + edges(c('i','f'))
h <- h + edges(c('l','a'))
h <- h + edges(c('l','d'))
h <- h + edges(c('l','g'))
cat(graphDistance(g,h), '\n')