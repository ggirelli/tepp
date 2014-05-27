"==.igraph.vs" = function(x, y) {
	# Can compare a couple of vertices, or a vertex to a vertex vector.
	# The comparison is carried out on the attributes of each vertex (except for the id).
	
	# If both node lists have length == 1 perform pair comparison
	if(length(x) == 1 && length(y) == 1) {

		# Retrieve attribute name list
		attr.list.x <- list.vertex.attributes(get('graph', attr(x, 'env')))
		attr.list.y <- list.vertex.attributes(get('graph', attr(y, 'env')))

		if(length(attr.list.x) == 0 | length(attr.list.y) == 0) {
			cat('Error: graphs have no node attributes.','\n')
			return(FALSE)
		}

		# Check attribute name list
		if(length(attr.list.x) != length(attr.list.y)) return(FALSE)
		if(length(which(attr.list.x %in% attr.list.y)) != length(attr.list.x)) return(FALSE)

		# Retrieve attribute values
		attr.value.x <- lapply(attr.list.x, x, FUN=function(attr, x) { return(eval(parse(text=paste0('x$', attr)))) })
		attr.value.y <- lapply(attr.list.y, y, FUN=function(attr, y) { return(eval(parse(text=paste0('y$', attr)))) })

		# Verify identity
		if(length(which(vapply(attr.value.x, FUN=function(x,y) { return(x %in% unlist(y)) }, FUN.VALUE=c(logical(1), logical(0)), y=attr.value.y))) == 0) return(FALSE)
		return(TRUE)

	} else { # Other cases:

		if(length(x) == 1) { # x is singular while y is not
			return(vapply(y, FUN=function(x, y, env) {
				if(y == V(get('graph', env))[x]) return(TRUE)
				return(FALSE)
			}, FUN.VALUE=c(logical(1),logical(0)), y=x, env=attr(y, 'env')))
		} else if(length(y) == 1) { # y is singular while x is not
			cat('Error: cannot compare node vector to single node.\n')
			return(NULL)
		} else {
			# Both x and y are note singular
			# They have different length
			if(length(x) != length(y)) {
				cat('Error: can compare only node vectors with the same length.')
				return(NULL)
			}
			# They have the same length
			cat('Compared node vectors in a pair-wise fashion.\n')
			return(vapply(x, FUN=function(x, y, env) {
				if(V(get('graph', env))[x] == y[x]) return(TRUE)
				return(FALSE)
			}, FUN.VALUE=c(logical(1),logical(0)), y=y, env=attr(x, 'env')))
		}
	}
}

"==.igraph.es" = function(x, y) {
	# Can compare a couple of edges, or an edge to an edge vector.
	# The comparison is carried out on the attributes of each edge  and, obiously, also on target/source.
	
	# If both edge lists have length == 1 perform pair comparison
	if(length(x) == 1 && length(y) == 1) {

		# Retrieve attribute name list
		attr.list.x <- list.edge.attributes(get('graph', attr(x, 'env')))
		attr.list.y <- list.edge.attributes(get('graph', attr(y, 'env')))
		
		# Check attribute name list
		if(length(attr.list.x) != length(attr.list.y)) return(FALSE)
		if(length(which(attr.list.x %in% attr.list.y)) != length(attr.list.x)) return(FALSE)

		# Retrieve target/source list
		edge.list.x <- get.edgelist(get('graph', attr(x, 'env')))[as.numeric(x),]
		edge.list.y <- get.edgelist(get('graph', attr(y, 'env')))[as.numeric(y),]

		# Check target/source list
		if(!identical(edge.list.x, edge.list.y)) return(FALSE)

		# If there are edge attributes, check their value
		if(length(attr.list.x) != 0) {
			# Retrieve attribute values
			attr.value.x <- lapply(attr.list.x, x, FUN=function(attr, x) { return(eval(parse(text=paste0('x$', attr)))) })
			attr.value.y <- lapply(attr.list.y, y, FUN=function(attr, y) { return(eval(parse(text=paste0('y$', attr)))) })

			# Verify identity
			if(length(which(vapply(attr.value.x, FUN=function(x,y) { return(x %in% unlist(y)) }, FUN.VALUE=c(logical(1), logical(0)), y=attr.value.y))) == 0) return(FALSE)
		} else {
			return(TRUE)
		}

	} else { # Other cases:

		if(length(x) == 1) { # x is singular while y is not
			return(vapply(y, FUN=function(x, y, env) {
				if(y == E(get('graph', env))[x]) return(TRUE)
				return(FALSE)
			}, FUN.VALUE=c(logical(1),logical(0)), y=x, env=attr(y, 'env')))
		} else if(length(y) == 1) { # y is singular while x is not
			cat('Error: cannot compare edge vector to single edge.\n')
			return(NULL)
		} else {
			# Both x and y are note singular
			# They have different length
			if(length(x) != length(y)) {
				cat('Error: can compare only edge vectors with the same length.')
				return(NULL)
			}
			# They have the same length
			cat('Compared edge vectors in a pair-wise fashion.\n')
			return(vapply(x, FUN=function(x, y, env) {
				if(E(get('graph', env))[x] == y[x]) return(TRUE)
				return(FALSE)
			}, FUN.VALUE=c(logical(1),logical(0)), y=y, env=attr(x, 'env')))
		}
	}
}

"%in%" = function(x, y) {
	# Extended to allow operations with igraph.vs (nodes)
	
	if(class(x) == class(y) && class(x) == 'igraph.vs') {
		if(length(y) == 0 || length(x) == 0) {
			cat('Error: cannot compare NULL\n')
			return(NULL)
		}
		if(length(x) == 1 && length(y) == 1) return(x == y)
		if(length(x) == 1 && length(y) > 1) return(length(which(x == y)) != 0)
		if(length(y) == 1 && length(x) > 1) return(length(which(y == x)) != 0)
		if(length(y) > 1 && length(x) > 1) {
			# Retrieve attribute name list
			attr.list.x <- list.vertex.attributes(get('graph', attr(x, 'env')))
			attr.list.y <- list.vertex.attributes(get('graph', attr(y, 'env')))

			if(length(attr.list.x) == 0 | length(attr.list.y) == 0) {
				cat('Error: graphs have no node attributes.','\n')
				return(FALSE)
			}

			return(vapply(x, FUN=function(x, y, env) {
				class(x) <- 'igraph.vs'; attr(x, 'env') <- env
				class(y) <- 'igraph.vs'
				return(x %in% y)
			}, FUN.VALUE=c(logical(0), logical(1)), y=y, env=attr(x, 'env')))
		}
		cat("Error: wrong lengths\n")
		return(NULL)
	}

	if(class(x) == class(y) && class(x) == 'igraph.es') {
		if(length(y) == 0 || length(x) == 0) {
			cat('Error: cannot compare NULL\n')
			return(NULL)
		}
		if(length(x) == 1 && length(y) == 1) return(x == y)
		if(length(x) == 1 && length(y) > 1) return(length(which(x == y)) != 0)
		if(length(y) == 1 && length(x) > 1) return(length(which(y == x)) != 0)
		if(length(y) > 1 && length(x) > 1) {
			return(vapply(x, FUN=function(x, y, env) {
				class(x) <- 'igraph.es'; attr(x, 'env') <- env
				class(y) <- 'igraph.es'
				return(x %in% y)
			}, FUN.VALUE=c(logical(0), logical(1)), y=y, env=attr(x, 'env')))
		}
		cat("Error: wrong lengths\n")
		return(NULL)
	}

	# Old %in%
	match(x, y, nomatch = 0L) > 0L
}

get.vertex.attr = function(name, v) {
	# Returns a certain attribute of a given vertex
	# 
	# Args:
	# 	name: attribute name
	# 	v: vertex
	# 	
	# Returns:
	# 	Attribute value
	
	return(eval(parse(text=paste0("V(get('graph', attr(v, 'env')))[v]$", name))))
}

get.vertex.attributes = function(v, id=FALSE) {
	# Returns all the attributes of a given vertex
	# 
	# Args:
	# 	v: vertex
	# 	
	# Returns:
	# 	List of attributes values

	# Retrieve vertex attributes list
	vl <- list.vertex.attributes(get('graph', attr(v, 'env')))

	# Check for attributes
	if(length(vl) == 0) {
		cat('No attributes to be retrieved.', '\n')
		return(NULL)
	}

	# Prepare table
	t <- sapply(vl, FUN=function(name, v) { return(get.vertex.attr(name, v)) }, v=v)

	# Add id column name
	if(id) {
		t <- cbind(1:length(t[,1]), t)
		colnames(t)[1] <- 'id'
	}

	# Terminate
	return(t)
}

get.edge.attr = function(name, e) {
	# Returns a certain attribute of a given edge
	# 
	# Args:
	# 	name: attribute name
	# 	e: edge
	# 	
	# Returns:
	# 	Attribute value
	
	return(eval(parse(text=paste0("E(get('graph', attr(e, 'env')))[e]$", name))))
}

get.edge.attributes = function(e, id=FALSE, path=FALSE) {
	# Returns all the attributes of a given edge
	# 
	# Args:
	# 	e: edge
	# 	
	# Returns:
	# 	List of attributes values

	# Retrieve edge attributes list
	el <- list.edge.attributes(get('graph', attr(e, 'env')))

	# Check for attributes
	if(length(el) == 0) {
		cat('No attributes to be retrieved.', '\n')
		return(NULL)
	}

	# Prepare table
	t <- sapply(el, FUN=function(name, e) { return(get.edge.attr(name, e)) }, e=e)

	# Add id column name
	if(id) {
		t <- cbind(1:length(t[,1]), t)
		colnames(t)[1] <- 'id'
	}

	# Add source/target
	if(path) {
		el <- get.edgelist(get('graph', attr(e, 'env')), names=FALSE)
		t <- cbind(t, paste0('n', el[,1]), paste0('n', el[,2]))
		colnames(t)[(length(t[1,])-1):length(t[1,])] <- c('source', 'target')
	}
	
	# Terminate
	return(t)
}

write.graph = function(graph, file, format) {
    if (!is.igraph(graph)) {
        stop("Not a graph object")
    }
    if (!is.character(file) || length(grep("://", file, fixed = TRUE)) > 
        0 || length(grep("~", file, fixed = TRUE)) > 0) {
        tmpfile <- TRUE
        origfile <- file
        file <- tempfile()
    }
    else {
        tmpfile <- FALSE
    }
    res <- switch(format, pajek = write.graph.pajek(graph, file), edgelist = write.graph.edgelist(graph, file), 
        	ncol = write.graph.ncol(graph, file), lgl = write.graph.lgl(graph, 
            file), graphml = write.graph.graphml(graph, 
            file), dimacs = write.graph.dimacs(graph, file), gml = write.graph.gml(graph, file), dot = write.graph.dot(graph, 
            file), leda = write.graph.leda(graph, file), json = write.graph.json(graph, file),
             stop(paste("Unknown file format:", format)))
    if (tmpfile) {
        buffer <- read.graph.toraw(file)
        write.graph.fromraw(buffer, origfile)
    }
    invisible(res)
}

write.graph.json = function(graph, file) {
	library('rjson')

	l <- list(nodes=list(), edges=list())

	# NODES
	l$nodes <- apply(get.vertex.attributes(V(g), id=TRUE), MARGIN=1, FUN=function(x, index) {
		data <- list(id=paste0('n', as.vector(x['id'])))
		for(attr in names(x)[which(names(x) != 'id')]) {
			data <- append(data, eval(parse(text=paste0('x[\'', attr, '\']'))))
		}
		return(list(data=data))
	})

	# EDGES
	l$edges <- apply(get.edge.attributes(E(g), id=TRUE, path=TRUE), MARGIN=1, FUN=function(x, index) {
		data <- list(id=paste0('e', as.vector(x['id'])))
		for(attr in names(x)[which(names(x) != 'id')]) {
			data <- append(data, eval(parse(text=paste0('x[\'', attr, '\']'))))
		}
		return(list(data=data))
	})

	write(toJSON(l), file)
}