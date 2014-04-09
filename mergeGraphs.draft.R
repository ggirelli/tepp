mergeGraphCouple.jointAttr <- function(g.one, g.two, attr.v.list=c(), attr.e.list=c(), attr.v.id='name', attr.v.action=list(), attr.e.action=list(), one.dominant) {
	#
	# Merges two graphs with common edge/vertex attributes
	#
	# Args:
	#	g.one, g.two: the graphs to be merged
	#	attr.v.list: list of vertex attributes to keep in the merged graph
	#				each attribute MUST be present in both g.one and g.two
	#	attr.e.list: list of edge attributes to keep in the merged graph
	#				each attribute MUST be present in both g.one and g.two
	#	attr.v.id: attribute used to identify each vertex (default: 'name')
	#	attr.v.action: how to act with each vertex attribute (e.g.: 'sum')
	#					attr.v.action$ATTR <- ACTION_STRING
	#	attr.e.action: how to act with each edge attribute (e.g.: 'sum')
	#					attr.e.action$ATTR <- ACTION_STRING
	#	one.dominant: TRUE means that in case of conflit without any specified
	#				action, g.one attributes prevail over g.two attributes.
	#				Vice versa if FALSE
	#
	# ACTION_STRING values:
	#	'sum': sum the values after applying as.numeric()
	#
	# Return:
	#	Merge graph
	
	# (1) if attributes in g.one and g.two are different, abort operation
	if(!identical
		(list.vertex.attributes(g.one),
		list.vertex.attributes(g.two)
	)) return(graph.empty()) # vertex
	if(!identical(
		list.edge.attributes(g.one),
		list.edge.attributes(g.two)
	)) return(graph.empty()) #edge

	###############
	# Start merge #
	###############

	# Assign the 'bigger' graph to g.merge
	g.one.length <- length(V(g.one)) + length(get.edgelist(g.one)[,1])
	g.two.length <- length(V(g.two)) + length(get.edgelist(g.two)[,1])
	if(g.one.length >= g.two.length) {
		g.merge <- g.one
	} else {
		g.merge <- g.two
		g.two <- g.one
		g.one <- g.merge
		one.dominant <- !one.dominant
	}

	#----------#
	# Vertices #
	#----------#

	# Prepare final graph vertices id by id
	v.merge.id.list <- c()
	for(i in 1:length(V(g.merge))) eval(parse(text=paste('v.merge.id.list <- append(v.merge.id.list,V(g.merge)[', i, ']$', attr.v.id, ')', sep='')))

	# Scroll through g.two vertices
	for(i in 1:length(V(g.two))) {

		# Select vertex
		v <- V(g.two)[i]
		v.id <- eval(parse(text=paste('V(g.two)[', i,']$',attr.v.id, sep='')))
		
		# If vertex is already present
		if(v.id %in% v.merge.id.list) {
			
			# Act on the attributes
			for(attr in attr.v.list) {

				# Retrieve the action associated to the attribute
				eval(parse(text=paste('attr.action <- attr.v.action$', attr, sep='')))

				# Verify if there is an action associated to the attribute
				if(identical(attr.action, 'sum')) {

					# Get g.two attribute value
					eval(parse(text=paste('a.two <- V(g.two)[', i, ']$', attr, sep='')))

					# Get g.one attribute
					eval(parse(text=paste('a.one <- V(g.one)[which(v.merge.id.list == v.id)]$', attr, sep='')))

					# Assign attributes sum to g.merge vertex
					eval(parse(text=paste('V(g.merge)[which(v.merge.id.list == v.id)]$', attr, ' <- as.numeric(a.one) + as.numeric(a.two)', sep='')))

				} else {

					# Assign attribute (g.two prevails)
					if(!one.dominant) eval(parse(text=paste('V(g.merge)[which(V(g.merge)$', attr.v.id, ' == v.id)]', '$', attr, ' <- V(g.two)[i]$', attr, sep='')))

				}

			}
		} else {
			
			# Prepare attributes assignment string
			attrs <- ''
			for(attr in attr.v.list) attrs <- paste(attrs, ', ', attr, '=\'', eval(parse(text=paste('V(g.two)[i]$', attr, sep=''))), '\'', sep=' ')
			
			# Add all the vertices with their attributes
			
			eval(parse(text=paste('g.merge <- g.merge + vertex(\'', v.id, '\'', attrs, ')', sep='')))
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
				eval(parse(text=paste('attr.action <- attr.e.action$', attr, sep='')))

				# Verify if there is an action associated to the attribute
				if(identical(attr.action, 'sum')) {

					# Get g.two attribute value
					eval(parse(text=paste('a.two <- E(g.two)[i]$', attr, sep='')))

					# Get g.one attribute
					eval(parse(text=paste('a.one <- E(g.one)[which(e.merge.list == paste(e.temp[1], e.temp[2], sep=\' <-- \'))]$', attr, sep='')))

					# Assign attributes sum to g.merge vertex
					eval(parse(text=paste('E(g.merge)[which(e.merge.list == paste(e.temp[1], e.temp[2], sep=\' <-- \'))]$', attr, ' <- as.numeric(a.one) + as.numeric(a.two)', sep='')))
				} else {

					# Assign attribute (g.two prevails)
					if(!one.dominant) eval(parse(text=paste('E(g.merge)[which(e.merge.list == paste(e.temp[1], e.temp[2], sep=\' <-- \'))]$', attr, ' <- as.numeric(paste(\'a.two <- E(g.two)[i]$\', attr, sep=\'\'))', sep='')))

				}
			}
		} else {

			# Prepare attributes assignment string
			attrs <- ''
			for(attr in attr.e.list) attrs <- paste(attrs, ', ', attr, '=', eval(parse(text=paste('E(g.two)[i]$', attr, sep=''))), sep=' ')

			# Add all the edges with their attributes
			eval(parse(text=paste('g.merge <- g.merge + edge(c(\'', e.temp[1], '\',\'', e.temp[2], '\')', attrs, ')', sep='')))

		}
	}

	#-----------#
	# Terminate #
	#-----------#
	return(g.merge)
}

#library('igraph')
#g1 <- graph.ring(10, directed=TRUE)
#g2 <- graph.ring(10, directed=TRUE)
#V(g1)$name <- letters[1:10]
#V(g2)$name <- letters[5:14]
#E(g1)$weight <- 1:10
#E(g2)$weight <- 11:20
#g3 <- mergeGraphs(
#	g1,
#	g2,
#	attr.v.list=c('name'),
#	attr.e.list=c('weight'),
#	attr.e.action=list(weight='sum')
#)

mergeGraphCouple.disjointAttr <- function( g.one, g.two, attr.v.list=c(), attr.e.list=c(), attr.v.id='name', attr.v.action=list(), attr.e.action=list(), one.dominant=TRUE) {
	#
	# Merges two graphs
	#
	# Args:
	#	g.one, g.two: the graphs to be merged
	#	attr.v.list: list of vertex attributes to keep in the merged graph
	#	attr.e.list: list of edge attributes to keep in the merged graph
	#	attr.v.id: attribute used to identify each vertex (default: 'name')
	#	attr.v.action: how to act with some vertex attribute (e.g.: 'sum')
	#					attr.v.action$ATTR <- list(attr=ACTION_STRING,...)
	#	attr.e.action: how to act with some edge attribute (e.g.: 'sum')
	#					attr.e.action$ATTR <- list(attr=ACTION_STRING,...)
	#	one.dominant: TRUE means that in case of conflit without any specified
	#				action, g.one attributes prevail over g.two attributes.
	#				Vice versa if FALSE
	#
	# ACTION_STRING values:
	#	'sum': sum the values after applying as.numeric()
	#
	# Return:
	#	Merge graph

	# Assign the 'bigger' graph to g.merge
	g.one.length <- length(V(g.one)) + length(get.edgelist(g.one)[,1])
	g.two.length <- length(V(g.two)) + length(get.edgelist(g.two)[,1])
	if(g.one.length >= g.two.length) {
		g.merge <- g.one
	} else {
		g.merge <- g.two
		g.two <- g.one
		g.one <- g.merge
		one.dominant <- !one.dominant
	}

	#----------#
	# Vertices #
	#----------#

	# Prepare final graph vertices id by id
	v.merge.id.list <- c()
	for(i in 1:length(V(g.merge))) eval(parse(text=paste('v.merge.id.list <- append(v.merge.id.list,V(g.merge)[i]$', attr.v.id, ')', sep='')))

	# Scroll through g.two vertices
	for(i in 1:length(V.(g.two))) {

		# Selecte vertex (i is the index of both vertex and its attributes)
		v.temp <- V(g.two)[i]
		eval(parse(text=paste('v.id <- V(g.two)[i]$',attr.v.id, sep='')))

		# If vertex is already present
        if(v.id %in% v.merge.id.list) {

        	# Retrieve corresponding vertex from g.merge
        	eval(parse(text=paste('merge.v.temp <- V(g.merge)[which(V(g.merge)$', attr.v.id,' == v.id)]',sep='')))

			# Evaluate the attributes
			for(attr in attr.v.list) {
				if(eval(parse(text=paste('!is.na(merge.v.temp$',
					attr, ')', sep='')))) {

					# Merge attributes
					if(attr %in% attr.v.action) {

						if(eval(parse(text=paste('attr.v.action$', attr, " == 'sum'", sep='')))) paste('V(g.merge)[which(V(g.merge)$', attr.v.id, ' == v.id)]$', attr, ' <- V(g.merge)[which(V(g.merge)$', attr.v.id, ' == v.id)]$', attr, ' + V(g.two)[i]$', attr, sep='')

					} else if(!one.dominant) {

						# Assign attribute (g.two prevails)
						eval(parse(text=paste('V(g.merge)[which(V(g.merge)$', attr.v.id, ' == v.id)]', '$', attr, ' <- V(g.two)[i]$', attr, sep='')))

					}

				} else {

					# Assign attribute (g.two prevails)
					if(!one.dominant) eval(parse(text=paste('V(g.merge)[which(V(g.merge)$', attr.v.id, ' == v.id)]', '$', attr, ' <- V(g.two)[i]$', attr, sep='')))

				}
			}

        } else {

        	# Build attributes assignment string
        	attrs <- ''
        	for(attr in attr.v.list) attrs <- paste(attrs, ', ', attr, '=\'', eval(parse(text=paste('V(g.two)[i]$', attr, sep=''))), '\'', sep='')

        	# Add it with its attributes
        	eval(parse(text=paste('g.merge <- g.merge + vertex(v.id, ', attrs, ')', sep='')))
        }

	}

	#-------#
	# Edges #
	#-------#

	# Prepare final graph edges
	e.merge.list <- get.edgelist(g.merge)
	e.merge.list <- paste(e.merge.list[,1], e.merge.list[,2], sep=' <-- ')

	# Scroll through g.two edges
	for(i in 1:length(get.edgelist(g.two)[,1])) {

		# Select edge (i is the index of both edge and its attributes)
		e.temp <- get.edgelist(g.two)[i,]

		if(paste(e.temp[1], e.temp[2], sep=' <-- ') %in% e.merge.list) {

			# Merge the attributes
			for(attr in attr.e.list) {

			}

		} else {

			# Add it with its attributes

		}

	}

	# Terminate
	return(g.merge)
}
