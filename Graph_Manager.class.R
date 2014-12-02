library('igraph')

# Class to manage graphml graphs and perform graph operations
GraphManager <- function() {

	# Instantiate Graph Manager
	nm <- list(

		# ---------- #
		# ATTRIBUTES #
		# ---------- #

		version = 2.2,

		# --------- #
		# FUNCTIONS #
		# --------- #
		
		# Conversions
		#-----------------------------------
		
		graph.to.attr.table = function (graph) {
			# Converts a graph into vertex/edge attribute tables
			# 
			# Args:
			# 	graph: network instance
			# 
			# Returns:
			# 	list(nodes=v.attr.table, edges=e.attr.table)
			
			# VERTICES #

			v.attr.list <- list.vertex.attributes(graph)
			v.count <- vcount(graph)

			if ( 0 == v.count ) {

				# Empty network
				v.attr.table <- c()

			} else if ( 1 == v.count ) {
				
				# Single vertex network
				v.attr.table <- c()
				for (attr in v.attr.list) {
					v.attr.table <- append(v.attr.table,
						eval(parse(text=paste0('V(graph)[1]$', attr))))
				}
				names(v.attr.table) <- v.attr.list

			} else {
				
				# 'normal' network
				v.attr.table <- c()
				for (attr in v.attr.list) {
					v.attr.table <- cbind(v.attr.table,
						eval(parse(text=paste0('V(graph)$', attr))))
				}
				colnames(v.attr.table) <- v.attr.list

			}

			# EDGES #
			
			e.attr.list <- list.edge.attributes(graph)
			e.count <- ecount(graph)
			
			if ( 0 == e.count ) {

				# Empty network
				e.attr.table <- c()

			} else if ( 1 == e.count ) {
				
				# Single edge network
				e.attr.table <- c()
				for (attr in e.attr.list) {
					e.attr.table <- append(e.attr.table,
						eval(parse(text=paste0('E(g)[1]$', attr))))
				}
				if ( !is.null(e.attr.table) ) names(e.attr.table) <- e.attr.list

				# Add source/target columns
				e.attr.table <- GraphManager()$add.edges.extremities(e.attr.table, graph, T)

			} else {

				# 'normal' network
				e.attr.table <- c()
				for (attr in e.attr.list) {
					e.attr.table <- cbind(e.attr.table,
						eval(parse(text=paste0('E(g)$', attr))))
				}
				if ( !is.null(e.attr.table) ) colnames(e.attr.table) <- e.attr.list

				# Add source/target columns
				e.attr.table <- GraphManager()$add.edges.extremities(e.attr.table, graph, T)

			}

			# END #
			return(list(nodes=v.attr.table, edges=e.attr.table))
		},

		attr.tables.to.graph = function (v.attr.table, e.attr.table) {
			# Converts the attribute table of  a graph into a graph
			# 
			# Args:
			# 	v.attr.table
			# 	e.attr.table
			# 
			# Returns:
			# 	graph
			
			g <- graph.empty()

			# VERTICES #
			if ( !is.null(nrow(v.attr.table)) ) {

				# Non-empty table
				g <- add.vertices(g, nrow(v.attr.table))
				for (attr in colnames(v.attr.table)) {
					attr.col.id <- which(colnames(v.attr.table) == attr)
					eval(parse(text=paste0('V(g)$', attr, ' <- v.attr.table[, attr.col.id]')))
				}

			} else if ( 0 != length(v.attr.table) ) {

				# Single-row table
				g <- add.vertices(g, 1)
				for (attr in names(v.attr.table)) {
					attr.col.id <- which(names(v.attr.table) == attr)
					eval(parse(text=paste0('V(g)$', attr, ' <- v.attr.table[attr.col.id]')))
				}

			}

			# EDGES #
			if ( !is.null(nrow(e.attr.table)) ) {

				# Non-empty table
				source.col.id <- which(colnames(e.attr.table) == 'source')
				target.col.id <- which(colnames(e.attr.table) == 'target')
				edge.pairwise.list <- c(rbind(
					unlist(lapply(e.attr.table[, source.col.id],
						FUN=function(x,vt) { return(which(vt[, 'id'] == x)) },vt=v.attr.table)),
					unlist(lapply(e.attr.table[, target.col.id],
						FUN=function(x,vt) { return(which(vt[, 'id'] == x)) },vt=v.attr.table))
				))
				g <- add.edges(g, edge.pairwise.list)
				for (attr in colnames(e.attr.table)) {
					attr.col.id <- which(colnames(e.attr.table) == attr)
					eval(parse(text=paste0('E(g)$', attr, ' <- e.attr.table[, attr.col.id]')))
				}

			} else if ( 0 != length(e.attr.table) ) {

				# Single-row table
				source.col.id <- which(names(e.attr.table) == 'source')
				target.col.id <- which(names(e.attr.table) == 'target')
				edge.pairwise.list <- c(
					V(g)[id == e.attr.table[source.col.id]],
					V(g)[id == e.attr.table[target.col.id]]
				)
				g <- add.edges(g, edge.pairwise.list)
				for (attr in names(e.attr.table)) {
					attr.col.id <- which(names(e.attr.table) == attr)
					eval(parse(text=paste0('E(g)$', attr, ' <- e.attr.table[attr.col.id]')))
				}

			}

			# END #
			return(g)
		},

		graph.to.list = function (graph) {
			# Converts a graph into a list, to allow JSON output
			# 
			# Args:
			# 	graph
			# 
			# Returns:
			# 	List graph for JSON output
			
			v.count <- vcount(graph)
			e.count <- ecount(graph)

			graph.list <- list(nodes=list(), edges=list())

			# VERTICES #
			if ( 1 == v.count ) {

				# Single-node graph
				graph.list$nodes <- list(data=
					GraphManager()$get.vertex.attributes(V(graph)[1], graph))

			} else {

				# 'normal' graph
				graph.list$nodes <- lapply(V(graph), FUN=function(v) {
					l <- GraphManager()$get.vertex.attributes(v, graph)
					return(list(data=l))
				})

			}

			# EDGES #
			if ( 1 == e.count ) {

				# Single-edge graph
				graph.list$edges <- list(data=
					GraphManager()$get.edge.attributes(E(graph)[1], graph))

			} else {

				# 'normal' graph
				graph.list$edges <- lapply(E(graph), FUN=function(e) {
					l <- GraphManager()$get.edge.attributes(e, graph)
					return(list(data=l))
				})

			}

			# RESTRUCTURE LIST #
			
			if ( 1 == v.count ) {
				if ( 0 == e.count ) {

					# Single node, no edges
					graph.list$nodes$group <- 'nodes'
					graph.list <- list(graph.list$nodes)

				} else if ( 1 == e.count ) {

					# Single node, single edge
					graph.list$nodes$group <- 'nodes'
					graph.list$edges$group <- 'edges'
					graph.list <- list(graph.list$nodes, graph.list$edges)

				} else {

					# Single node, multiple edges
					graph.list$nodes$group <- 'nodes'
					edges <- graph.list$edges
					graph.list <- list(graph.list$nodes)
					for (edge in edges) {
						edge$group <- 'edges'
						graph.list <- append(graph.list, list(edge))
					}

				}
			} else if ( 1 == e.count ) {

				# Single edge, multiple nodes
				graph.list$edges$group <- 'edges'
				nodes <- graph.list$nodes
				graph.list <- list(graph.list$edges)
				for (node in nodes) {
					node$group <- 'nodes'
					graph.list <- append(graph.list, list(node))
				}

			}

			# END #
			return(graph.list)
		},

		attr.tables.to.list = function (v.attr.table, e.attr.table) {
			# Converts a graph into a list, to allow JSON output
			# Based on attribute tables
			# 
			# Args:
			# 	v.attr.table
			# 	e.attr.table
			# 
			# Returns:
			# 	List graph for JSON output
			
			if ( !is.null(nrow(v.attr.table)) ) {
				v.count <- nrow(v.attr.table)
			} else if ( 0 != length(v.attr.table) ) {
				v.count <- 1
			} else {
				v.count <- 0
			}

			if ( !is.null(nrow(e.attr.table)) ) {
				e.count <- nrow(e.attr.table)
			} else if ( 0 != length(e.attr.table) ) {
				e.count <- 1
			} else {
				e.count <- 0
			}


			graph.list <- list(nodes=list(), edges=list())
			v.attr.table <- GraphManager()$expand.attr.table(v.attr.table, c('x','y'))

			# VERTICES #
			if ( 1 == v.count ) {

				# Single-node graph
				l <- list()
				for (col in names(v.attr.table)) {
					col.id <- which(col == names(v.attr.table))
					l[col] <- v.attr.table[col.id]
				}
				graph.list$nodes <- list(data=l, position=list(x=as.numeric(l['x']),y=as.numeric(l['y'])))

			} else if ( 0 != v.count ) {

				# 'normal' graph
				graph.list$nodes <- lapply(1:nrow(v.attr.table), FUN=function(x, v.attr.table) {
					l <- list()
					for (col in colnames(v.attr.table)) {
						col.id <- which(col == colnames(v.attr.table))
						l[col] <- v.attr.table[x, col.id]
					}

					return(list(data=l, position=list(x=as.numeric(l['x']),y=as.numeric(l['y']))))
				}, v.attr.table=v.attr.table)

			}
			# EDGES #
			if ( 1 == e.count ) {

				# Single-edge graph
				l <- list()
				for (col in names(e.attr.table)) {
					col.id <- which(col == names(e.attr.table))
					l[col] <- e.attr.table[col.id]
				}
				graph.list$edges <- list(data=l)

			} else if ( 0 != e.count ) {

				# 'normal' graph
				graph.list$edges <- lapply(1:nrow(e.attr.table), FUN=function(x, e.attr.table) {
					l <- list()
					for (col in colnames(e.attr.table)) {
						col.id <- which(col == colnames(e.attr.table))
						l[col] <- e.attr.table[x, col.id]
					}
					return(list(data=l))
				}, e.attr.table=e.attr.table)

			}

			# RESTRUCTURE LIST #
			
			if ( 1 == v.count ) {
				if ( 0 == e.count ) {

					# Single node, no edges
					graph.list$nodes$group <- 'nodes'
					graph.list <- list(graph.list$nodes)

				} else if ( 1 == e.count ) {

					# Single node, single edge
					graph.list$nodes$group <- 'nodes'
					graph.list$edges$group <- 'edges'
					graph.list <- list(graph.list$nodes, graph.list$edges)

				} else {

					# Single node, multiple edges
					graph.list$nodes$group <- 'nodes'
					edges <- graph.list$edges
					graph.list <- list(graph.list$nodes)
					for (edge in edges) {
						edge$group <- 'edges'
						graph.list <- append(graph.list, list(edge))
					}

				}
			} else if ( 1 == e.count ) {

				# Single edge, multiple nodes
				graph.list$edges$group <- 'edges'
				nodes <- graph.list$nodes
				graph.list <- list(graph.list$edges)
				for (node in nodes) {
					node$group <- 'nodes'
					graph.list <- append(graph.list, list(node))
				}

			}

			# END #
			return(graph.list)
		},

		graph.list.to.attr.tables = function (graph.list) {
			# Converts a graph from JSON to attr.tables
			# 
			# Args:
			# 	graph.list in JSON format
			# 
			# Returns:
			# 	list(nodes=v.attr.table, edges=e.attr.table)

			v.attr.table <- c()
			e.attr.table <- c()

			if ( !is.null(graph.list$nodes) ) {

				# List divided into 'nodes' and 'edges'
				for (node in graph.list$nodes) {
					single.row <- c()
					for (v in node$data) {
						single.row <- append(single.row, v)
					}
					names(single.row) <- names(node$data)
					v.attr.table <- rbind(v.attr.table, single.row)
				}
				for (edge in graph.list$edges) {
					single.row <- c()
					for (v in edge$data) {
						single.row <- append(single.row, v)
					}
					names(single.row) <- names(edge$data)
					e.attr.table <- rbind(e.attr.table, single.row)
				}

			} else {

				# Mixed-list
				for (el in graph.list) {
					if ( 'nodes' == el$group ) {
						node <- el
						single.row <- c()
						for (v in node$data) {
							single.row <- append(single.row, v)
						}
						names(single.row) <- names(node$data)
						v.attr.table <- rbind(v.attr.table, single.row)
					}
					if ( 'edges' == el$group ) {
						edge <- el
						single.row <- c()
						for (v in edge$data) {
							single.row <- append(single.row, v)
						}
						names(single.row) <- names(edge$data)
						e.attr.table <- rbind(e.attr.table, single.row)
					}
				}

			}

			# END #
			return(list(nodes=v.attr.table, edges=e.attr.table))
		},

		graph.list.to.graph = function (graph.list) {
			# ***STILL WORKING ON THIS ONE***
			# To be used when a graph.list contains list attributes
			# 
			# Args:
			# 	graph.list
			# 
			# Returns:
			# 	graph

			g <- graph.empty()
			
			if ( !is.null(graph.list$nodes) ) {

				for (node in graph.list$nodes) {
					v <- vertex()
					for (attr in names(node$data)) {
						if ( 'list' == class(node$data[attr]) ) {
							eval(parse(text=paste0('v$', attr, ' <- as.character(node$data[attr])')))
						} else {
							eval(parse(text=paste0('v$', attr, ' <- node$data[attr]')))
						}
					}
					g <- g + v
				}

				for (edge in graph.list$edges) {
					e <- edge(V(g)[id == edge$data$source], V(g)[id == edge$data$target])
					for (attr in names(edge$data)[which(!names(edge$data) %in% c('source','target'))]) {
						if ( 'list' == class(edge$data[attr]) ) {
							eval(parse(text=paste0('e$', attr, ' <- as.character(edge$data[attr])')))
						} else {
							eval(parse(text=paste0('e$', attr, ' <- edge$data[attr]')))
						}
					}
					g <- g + e
				}

			} else {

				for (el in graph.list) {

					if ( 'nodes' == el$group ) {
						node <- el
						v <- vertex()
						for (attr in names(node$data)) {
							if ( 'list' == class(node$data[attr]) ) {
								eval(parse(text=paste0('v$', attr, ' <- as.character(node$data[attr])')))
							} else {
								eval(parse(text=paste0('v$', attr, ' <- node$data[attr]')))
							}
						}
						g <- g + v
					}

					if ( 'edges' == el$group ) {
						edge <- el
						e <- edge(V(g)[id == edge$data$source], V(g)[id == edge$data$target])
						for (attr in names(edge$data)[which(!names(edge$data) %in% c('source','target'))]) {
							if ( 'list' == class(edge$data[attr]) ) {
								eval(parse(text=paste0('e$', attr, ' <- as.character(edge$data[attr])')))
							} else {
								eval(parse(text=paste0('e$', attr, ' <- edge$data[attr]')))
							}
						}
						g <- g + e
					}

				}

			}

			# END #
			return(g)
		},

		# Table management
		#-----------------------------------

		get.vertex.attributes = function (v, graph) {
			# Retrieves the (attribute,value) couples of a given vertex
			# 
			# Args:
			# 	v: the vertex
			# 	graph
			# 
			# Returns:
			# 	NULL if no vertex attributes were found
			# 	A list with k=attr.name and v=attr.val

			l <- list()

			# Vertex attributes list
			v.attr.list <- list.vertex.attributes(graph)

			if ( 0 == length(v.attr.list) ) {

				# No vertex attributes
				l <- NULL

			} else {

				# Save attribute values in l
				for (attr in v.attr.list) {
					l[attr] <- eval(parse(text=paste0('V(graph)[v]$', attr)))
				}

			}

			# END #
			return(l)
		},

		get.edge.attributes = function (e, graph) {
			# Retrieves the (attribute,value) couples of a given edge
			# 
			# Args:
			# 	e: the edge
			# 	graph
			# 
			# Returns:
			# 	NULL if no edge attributes were found
			# 	A list with k=attr.name and v=attr.val

			l <- list()

			# Edge attributes list
			e.attr.list <- list.edge.attributes(graph)

			if ( 0 == length(e.attr.list) ) {

				# No edge attributes
				l <- NULL

			} else {

				# Save attribute values in l
				for (attr in e.attr.list) {
					l[attr] <- eval(parse(text=paste0('E(graph)[e]$', attr)))
				}

			}

			# END #
			return(l)
		},

		expand.attr.table = function (table, attr.list) {
			# Expands the column of an attribute table
			# 
			# Args:
			# 	table
			# 	attr.list
			# 
			# Returns:
			# 	The expanded table
			
			if ( !is.null(nrow(table)) ) {
				
				# Non-empty table
				for (attr in attr.list) {
					if ( is.null(attr) ) next
					if ( !attr %in% colnames(table) ) {
						table <- cbind(table,NA)
						colnames(table)[ncol(table)] <- attr
					}
				}
				
			} else if ( 0 != length(table) ) {
				
				# Single-row table
				for (attr in attr.list) {
					if ( is.null(attr) ) next
					if ( !attr %in% names(table) ) {
						table <- append(table,NA)
						names(table)[length(table)] <- attr
					}
				}

			}

			# END #
			return(table)
		},

		add.collapsed.col = function (table, attr.list, col.name, sep) {
			# Collapses columns of a table
			# 
			# Args:
			# 	table
			# 	attr.list
			# 	sep: character for sep
			# 
			# Returns:
			# 	The updated table
			
			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				collapsed <- c()
				for (attr in attr.list) {
					collapsed <- paste0(collapsed, sep,
						table[, which(attr == colnames(table)) ])
				}
				table <- cbind(table, collapsed)
				colnames(table)[ncol(table)] <- col.name

			} else if ( 0 != length(table) ) {

				# Single-row table
				collapsed <- ''
				for (attr in attr.list) {
					collapsed <- paste0(collapsed, sep,
						table[ which(attr == names(table)) ])
				}
				table <- append(table, collapsed)
				names(table)[length(table)] <- col.name

			}

			# END #
			return(table)
		},

		add.prefix.to.col = function (table, col, prefix) {
			# Adds a prefix to every cell of a given column
			# 
			# Args:
			# 	table
			# 	col: column label
			# 	prefix
			# 
			# Returns:
			# 	The updated table

			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				if ( col %in% colnames(table) ) {
					col.id <- which(colnames(table) == col)
					table[, col.id] <- paste0(prefix, table[, col.id])
				}

			} else if ( 0 != length(table) ) {

				# Single-row table
				if ( col %in% names(table) ) {
					col.id <- which(names(table) == col)
					table[col.id] <- paste0(prefix, table[col.id])
				}

			}

			# END #
			return(table)
		},

		count.rows = function (table) {
			#

			c.count <- 0

			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				c.count <- nrow(table)

			} else if ( 0 != length(table) ) {

				# Single-row table
				c.count <- 1

			}

			# END #
			return(c.count)
		},

		get.col = function (table, col) {
			# Extracts a column from a table
			# 
			# Args:
			# 	table
			# 	col
			# 
			# Returns:
			# 	The column
			
			data <- NULL

			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				if ( col %in% colnames(table) ) {
					data <- table[, which(col == colnames(table))]
				}

			} else if ( 0 != length(table) ) {

				# Single-row table
				if ( col %in% names(table) ) {
					data <- table[which(col == names(table))]
				}

			}

			# END #
			return(data)
		},

		get.col.names = function (table) {
			# Retrieves the list of colnames
			# 
			# Args:
			# 	table
			# 
			# Returns:
			# 	The list of colnames

			c.list <- NULL

			if ( !is.null(nrow(table)) ) {
				c.list <- colnames(table)
			} else if ( 0 != length(table) ) {
				c.list <- names(table)
			}

			# END #
			return(c.list)
		},

		add.col.names = function (table, col.names) {
			# Adds colnames to a table
			# ncol(table) must be equal to length(col.names)
			# 
			# Args:
			# 	table
			# 	col.names
			# 	
			# Returns:
			# 	The updated table

			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				if ( length(col.names) == ncol(table) ) {
					colnames(table) <- col.names
				}

			} else if ( 0 != length(table) ) {

				# Single-row table
				if ( length(col.names) == length(table) ) {
					names(table) <- col.names
				}

			}

			# END #
			return(table)
		},

		rename.col = function (table, old.name, new.name) {

			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				colnames <- colnames(table)
				if ( old.name %in% colnames && !new.name %in% colnames ) {
					col.id <- which(old.name == colnames)
					colnames(table)[col.id] <- new.name
				}

			} else if ( 0 != length(table) ) {

				# Single-row table
				colnames <- names(table)
				if ( old.name %in% colnames && !new.name %in% colnames ) {
					col.id <- which(old.name == colnames)
					names(table)[col.id] <- new.name
				}

			}

			# END #
			return(table)
		},

		rm.cols = function (table, col.list) {
			# Removes columns from table if present
			# 
			# Args:
			# 	table
			# 	col.list
			# 	
			# Returns:
			# 	The updated table
			
			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				for (col in col.list) {
					if ( col %in% colnames(table) ) {
						rm <- which(colnames(table) == col)
						if ( 0 != length(rm) ) table <- table[, -rm]
					}
				}

			} else if ( 0 != length(table) ) {

				# Single-row table
				for (col in col.list) {
					if ( col %in% names(table) ) {
						rm <- which(names(table) == col)
						if ( 0 != length(rm) ) table <- table[-rm]
					}
				}

			}

			# END #
			return(table)
		},

		rm.rows.based.on.identity = function (table, identity.col, identity.to.rm) {
			# Removes rows from table based on identity
			# 
			# Args:
			# 	table
			# 	identity.col name
			# 	identity.to.rm {vector}
			# 
			# Returns:
			# 	The updated table or NULL

			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				if ( !identity.col %in% colnames(table) ) return(NULL)

				identity.col.id <- which(identity.col == colnames(table))
				rows.to.rm <- which(table[, identity.col.id] %in% identity.to.rm)
				if ( 0 != length(rows.to.rm) ) {
					table <- table[-rows.to.rm, ]
				}

			} else if ( 0 != length(table) ) {

				# Single-row table
				if ( !identity.col %in% names(table) ) return(NULL)

				identity.col.id <- which(identity.col == names(table))
				if ( table[identity.col.id] %in% identity.to.rm ) return(NULL)

			}

			# END #
			return(table)
		},

		rm.duplicated.identity = function (table, identity.col) {
			# Removes rows with duplicated identity
			# 
			# Args:
			# 	table
			# 	identity.col
			# 
			# Returns:
			# 	The updated table

			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				identity.col.id <- which(identity.col == colnames(table))
				identity.dups <- which(duplicated(table[, identity.col.id]))
				if ( 0 != length(identity.dups) ) table <- table[-identity.dups, ]

			}

			# END #
			return(table)
		},

		sort.table.cols = function (table) {
			# Orders a table columns
			# 
			# Args:
			# 	table
			# 
			# Returns:
			# 	The updated table
			
			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				table <- table[, order(colnames(table))]

			} else if ( 0 != length(table) ) {

				# Single-row table
				table <- table[order(names(table))]

			}

			# END #
			return(table)
		},

		# Edge extremities
		#-----------------------------------

		add.edges.extremities = function (e.attr.table, graph, names) {
			# Adds source/target columns to an e.attr.table
			# 
			# Args:
			# 	e.attr.table
			# 	graph
			# 	names: boolean for names or numerical ids
			# 
			# Returns:
			# 	The updated table
			
			e.count <- ecount(graph)

			if ( 0 == e.count ) {

				# Empty network
				e.attr.table <- c()

			} else if ( 1 == e.count ) {

				# Single edge network
				e.attr.table.tmp <- GraphManager()$rm.cols(e.attr.table, c('source', 'target'))

				add.col.name <- FALSE
				if ( is.null(nrow(e.attr.table.tmp)) ) {
					add.col.name <- TRUE
					col.name <- names(e.attr.table)[1]
				}

				e.attr.table <- append(e.attr.table.tmp, get.edgelist(graph, names=names))

				if ( add.col.name) {
					if ( 0 != length(col.name) ) names(e.attr.table)[1] <- col.name
				}
				if ( is.null(names(e.attr.table)) ) {
					names(e.attr.table) <- c('source', 'target')
				} else {
					names(e.attr.table)[length(e.attr.table)-1] <- 'source'
					names(e.attr.table)[length(e.attr.table)] <- 'target'
				}

			} else {
				
				# 'normal' network
				e.attr.table.tmp <- GraphManager()$rm.cols(e.attr.table, c('source', 'target'))

				add.col.name <- FALSE
				if ( is.null(nrow(e.attr.table.tmp)) ) {
					add.col.name <- TRUE
					col.name <- colnames(e.attr.table)[1]
				}
				e.attr.table <- cbind(e.attr.table.tmp, get.edgelist(graph, names=names))
				
				if ( add.col.name) {
					if ( 0 != length(col.name) ) colnames(e.attr.table)[1] <- col.name
				}
				if ( is.null(colnames(e.attr.table)) ) {
					colnames <- c('source', 'target')
				} else {
					colnames(e.attr.table)[ncol(e.attr.table)-1] <- 'source'
					colnames(e.attr.table)[ncol(e.attr.table)] <- 'target'
				}

			}

			# END #
			return(e.attr.table)
		},

		update.row.ids = function (table) {
			# ***TO APPLY ONLY ON TABLES SHRINKED BASED ON IDENTITY***
			# Updates the id of the rows of the given table
			# Id starts from '0'
			# 
			# Args:
			# 	table
			# 
			# Returns:
			# 	The updated table with numerical ids column 'id'

			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				table <- GraphManager()$expand.attr.table(table, c('id'))
				table[, which('id' == colnames(table))] <- 0:(nrow(table) - 1)

			} else if ( 0 != length(table) ) {

				# Single-row table
				table <- GraphManager()$expand.attr.table(table, c('id'))
				table[which('id' == names(table))] <- 0

			}

			# END #
			return(table)
		},

		update.row.ids.and.extremities = function (e.attr.table, e.prefix,
			v.attr.table, v.prefix, v.identity.col) {
			# Updates ids and extremities of attribute tables
			# v.attr.table must contain 'id'
			# e.attr.table must contain 'id', 'source' and 'target'
			# Extremities must be identity based
			# 
			# Args:
			# 	e.attr.table
			# 	e.prefix
			# 	v.attr.table
			# 	v.prefix
			# 	v.identity.col
			# 
			# Returns:
			# 	The updated tables: list(nodes=v.attr.table, edges=e.attr.table)

			# VERTICES #

			v.attr.table <- GraphManager()$update.row.ids(v.attr.table)
			v.attr.table <- GraphManager()$add.prefix.to.col(v.attr.table, 'id', 'v')

			# EDGES #
			
			e.attr.table <- GraphManager()$check.extremities(e.attr.table, v.attr.table, v.identity.col)
			e.attr.table <- GraphManager()$convert.extremities.to.v.id.based.on.table(e.attr.table,
				v.attr.table, v.identity.col)
			e.attr.table <- GraphManager()$update.row.ids(e.attr.table)
			e.attr.table <- GraphManager()$add.prefix.to.col(e.attr.table, 'id', 'e')


			# END #
			return(list(nodes=v.attr.table, edges=e.attr.table))
		},

		convert.extremities.to.v.identity = function (e.attr.table, v.attr.table,
			v.identity.col, graph) {
			# Converts edge extremities in vertex identities
			# The actual extremities are first replaced with numerical v.id
			# 
			# Args:
			# 	e.attr.table
			# 	v.attr.table
			# 	v.identity.col
			# 	graph
			# 
			# Returns:
			# 	The updated e.attr.table

			if ( is.null(v.attr.table) ) {
				return(NULL)
			} else if ( is.null(nrow(v.attr.table)) ) {
				if ( 0 == length(v.attr.table) ) return(NULL)
			}

			v.attr.table <- GraphManager()$update.row.ids(v.attr.table)
			e.attr.table <- GraphManager()$convert.extremities.to.v.id(e.attr.table,
				v.attr.table, v.identity.col, graph)

			if ( !is.null(nrow(e.attr.table)) ) {

				# Non-empty table

				# Replace numerical v.id with v.identity
				source.col.id <- which('source' == colnames(e.attr.table))
				e.attr.table[, source.col.id] <- GraphManager()$get.col(v.attr.table,
					v.identity.col)[as.numeric(e.attr.table[, source.col.id])]
				target.col.id <- which('target' == colnames(e.attr.table))
				e.attr.table[, target.col.id] <- GraphManager()$get.col(
						v.attr.table, v.identity.col)[as.numeric(e.attr.table[, target.col.id])]

			} else if ( 0!= length(e.attr.table) ) {

				# Single-row table

				# Replace numerical v.id with v.identity
				source.col.id <- which('source' == names(e.attr.table))
				e.attr.table[source.col.id] <- GraphManager()$get.col(v.attr.table,
					v.identity.col)[as.numeric(e.attr.table[source.col.id])]
				target.col.id <- which('target' == names(e.attr.table))
				e.attr.table[target.col.id] <- GraphManager()$get.col(v.attr.table,
					v.identity.col)[as.numeric(e.attr.table[target.col.id])]

			}

			# END #
			return(e.attr.table)
		},

		convert.extremities.to.v.id = function (e.attr.table, v.attr.table,
			v.identity.col, graph) {
			# Converts edge extremities in numerical v.id based on the graph
			# 
			# Args:
			# 	e.attr.table
			# 	v.attr.table
			# 	v.identity.col
			# 	graph
			# 
			# Returns:
			# 	The updated e.attr.table

			if ( is.null(v.attr.table) ) {
				return(NULL)
			} else if ( is.null(nrow(v.attr.table)) ) {
				if ( 0 == length(v.attr.table) ) return(NULL)
			}

			v.attr.table <- GraphManager()$update.row.ids(v.attr.table)

			if ( !is.null(nrow(e.attr.table)) ) {

				# Non-empty table
				e.attr.table.tmp <- GraphManager()$rm.cols(e.attr.table, c('source', 'target'))

				add.col.name <- FALSE
				if ( is.null(nrow(e.attr.table.tmp)) ) {
					add.col.name <- TRUE
					col.name <- colnames(e.attr.table)[1]
				}

				# Add numerical v.id as extremities
				e.attr.table <- GraphManager()$add.edges.extremities(e.attr.table.tmp, graph, F)

				if ( add.col.name) {
					if ( 0 != length(col.name) ) colnames(e.attr.table)[1] <- col.name
				}

			} else if ( 0 != length(e.attr.table) ) {

				# Single-row table
				e.attr.table.tmp <- GraphManager()$rm.cols(e.attr.table, c('source,', 'target'))

				add.col.name <- FALSE
				if ( is.null(nrow(e.attr.table.tmp)) ) {
					add.col.name <- TRUE
					col.name <- names(e.attr.table)[1]
				}

				# Add numerical v.id as extremities
				e.attr.table <- GraphManager()$add.edges.extremities(e.attr.table.tmp, graph, F)

				if ( add.col.name) {
					if ( 0 != length(col.name) ) names(e.attr.table)[1] <- col.name
				}

			}

			# END #
			return(e.attr.table)
		},

		convert.extremities.to.v.id.based.on.table = function (e.attr.table, v.attr.table,
			v.identity.col) {
			# Converts edge extremities in numerical v.id based on the v.attr.table
			# The current source/target must be v.identity based
			# 
			# Args:
			# 	e.attr.table
			# 	v.attr.table
			# 	v.identity.col
			# 	graph
			# 
			# Returns:
			# 	The updated e.attr.table

			if ( is.null(v.attr.table) ) {
				return(NULL)
			} else if ( is.null(nrow(v.attr.table)) ) {
				if ( 0 == length(v.attr.table) ) return(NULL)
			}
			
			if ( !is.null(nrow(e.attr.table)) ) {
				
				# Non-empty table
				e.attr.table <- GraphManager()$expand.attr.table(e.attr.table,
					c('source', 'target'))
				source.col.id <- which('source' == colnames(e.attr.table))
				target.col.id <- which('target' == colnames(e.attr.table))
				
				if ( !is.null(nrow(v.attr.table)) ) {
					
					# Non-empty table
					v.identity.col.id <- which(v.identity.col == colnames(v.attr.table))
					v.id.col.id <- which('id' == colnames(v.attr.table))
					
					e.attr.table[, source.col.id] <- unlist(lapply(e.attr.table[, source.col.id],
						FUN=function(x, v.id.col.id) {
							return(v.attr.table[which(
								v.attr.table[, v.identity.col.id] == x), v.id.col.id])
						}, v.id.col.id=v.id.col.id))
					e.attr.table[, target.col.id] <- unlist(lapply(e.attr.table[, target.col.id],
						FUN=function(x, v.id.col.id) {
							return(v.attr.table[which(
								v.attr.table[, v.identity.col.id] == x), v.id.col.id])
						}, v.id.col.id=v.id.col.id))

				} else if ( 0 != length(v.attr.table) ) {

					# Single.row table
					e.attr.table[, source.col.id] <- v.attr.table[v.id.col.id]
					e.attr.table[, target.col.id] <- v.attr.table[v.id.col.id]

				}
				

			} else if ( 0 != length(e.attr.table) ) {

				# Single-row table
				e.attr.table <- GraphManager()$expand.attr.table(e.attr.table,
					c('source', 'target'))
				source.col.id <- which('source' == names(e.attr.table))
				target.col.id <- which('target' == names(e.attr.table))

				if ( !is.null(nrow(v.attr.table)) ) {

					# Non-empty table
					v.identity.col.id <- which(v.identity.col == colnames(v.attr.table))
					v.id.col.id <- which('id' == colnames(v.attr.table))

					e.attr.table[source.col.id] <- v.attr.table[which(
						v.attr.table[, v.identity.col.id] == e.attr.table[source.col.id]
						), v.id.col.id]
					e.attr.table[target.col.id] <- v.attr.table[which(
						v.attr.table[, v.identity.col.id] == e.attr.table[target.col.id]
						), v.id.col.id]

				} else if ( 0 != length(v.attr.table) ) {

					# Single-row table
					e.attr.table[source.col.id] <- v.attr.table[v.id.col.id]
					e.attr.table[target.col.id] <- v.attr.table[v.id.col.id]

				}

			}

			# END #
			return(e.attr.table)
		},

		check.extremities = function (e.attr.table, v.attr.table, v.identity.col) {
			# Removes rows from e.attr.table that lost one or both extremities
			# 'source' and 'target' columns must be present in e.attr.table
			# Extremities must be identity-based
			# 
			# Args:
			# 	e.attr.table
			# 	v.attr.table
			# 	v.identity.table
			# 
			# Returns:
			# 	The updated table

			# RETRIEVE V.IDENTITIES #

			if ( !is.null(nrow(v.attr.table)) ) {

				# Non-empty table
				if ( !v.identity.col %in% colnames(v.attr.table) ) return(NULL)

				v.identity.col.id <- which(v.identity.col == colnames(v.attr.table))
				v.identity <- v.attr.table[, v.identity.col.id]

			} else if ( 0 != length(v.attr.table) ) {

				# Single-row table
				if ( !v.identity.col %in% names(v.attr.table) ) return(NULL)

				v.identity.col.id <- which(v.identity.col == names(v.attr.table))
				v.identity <- v.attr.table[v.identity.col.id]

			} else {

				# No vertices, no edges
				return(NULL)

			}
			
			# CHECK EDGES #
			
			if ( !is.null(nrow(e.attr.table)) ) {

				# Non-empty table
				if ( 0 != length(which(!c('source', 'target') %in% colnames(e.attr.table) ))) {
					return(NULL)
				}

				source.col.id <- which('source' == colnames(e.attr.table))
				rows.to.rm <- which(!e.attr.table[, source.col.id] %in% v.identity)
				if ( 0!= length(rows.to.rm) ) {
					e.attr.table <- e.attr.table[-rows.to.rm, ]
				}
				target.col.id <- which('target' == colnames(e.attr.table))
				rows.to.rm <- which(!e.attr.table[, target.col.id] %in% v.identity)
				if ( 0!= length(rows.to.rm) ) {
					e.attr.table <- e.attr.table[-rows.to.rm, ]
				}

			} else if ( 0!= length(e.attr.table) ) {

				# Single-row table
				if ( 0 != length(which(!c('source', 'target') %in% names(e.attr.table))) ) {
					return(NULL)
				}

				source.col.id <- which('source' == names(e.attr.table))
				if ( !e.attr.table[source.col.id] %in% v.identity ) return(NULL)
				target.col.id <- which('target' == names(e.attr.table))
				if ( !e.attr.table[target.col.id] %in% v.identity ) return(NULL)

			}

			# END #
			return(e.attr.table)
		},

		# Manage multiple networks
		#-----------------------------------

		append.to.table.list = function (table.list, table) {
			# Append a table to a table.list
			# 
			# Args:
			# 	table.list
			# 	table
			# 
			# Returns:
			# 	The updated table.list
			
			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				return(append(table.list, list(table)))

			} else if ( 0 != length(table) ) {

				# Single-row table
				return(append(table.list, list(table)))
			}

			# END #
			return(table.list)
		},

		get.colnames.from.table.list = function (table.list) {
			# Returns the colnames of the tables in a table.list
			# 
			# Args:
			# 	table.list
			# 
			# Returns:
			# 	The colnames of the tables in a table.list	

			colnames <- c()
			for (table in table.list) {

				if ( !is.null(nrow(table)) ) {

					# Non-empty table
					colnames <- append(colnames, colnames(table))

				} else if ( 0 != length(table) ) {

					# Single-row table
					colnames <- append(colnames, names(table))

				}

			}

			# END #
			return(unique(colnames))
		},

		merge.tables.from.table.list = function (table.list) {
			# Executes the rbind of the tables in a table.list
			# 
			# Args:
			# 	table.list
			# 	
			# Returns:
			# 	The rbound table.list

			end.table <- c()
			colnames <- GraphManager()$get.colnames.from.table.list(table.list)

			# RBIND TABLES #
			for (table in table.list) {

				if ( !is.null(nrow(table)) ) {

					# Non-empty table
					if ( 0 != length(which(!colnames(table) %in% colnames)) ) {
						table <- GraphManager()$expand.attr.table(table, colnames)
						table <- sort.table.cols(table)
					}
					end.table <- rbind(end.table, table)

				} else if ( 0 != length(table) ) {

					# Single-row table
					if ( 0 != length(which(!names(table) %in% colnames)) ) {
						table <- GraphManager()$expand.attr.table(table, colnames)
						table <- sort.table.cols(table)
					}
					end.table <- rbind(end.table, table)

				}

			}

			# ADD COLNAMES #
			if ( !is.null(nrow(end.table)) ) {

				# Non-empty table
				colnames(end.table) <- colnames

			} else if ( 0!= length(end.table) ) {

				# Single-row table
				names(end.table) <- colnames

			}

			# END #
			return(end.table)
		},

		extract.subtable.based.on.identity = function (table, identity.col, identity.val) {
			# Extract a subtable based on identity
			# 
			# Args:
			# 	table
			# 	identity.col name
			# 	identity.val value
			# 
			# Returns:
			# 	NULL if the identity.col or the identity.val were not found
			# 	The extracted subtable

			if( !is.null(nrow(table)) ) {

				# Non-empty table
				if ( !identity.col %in% colnames(table) ) return(NULL)

				identity.list <- table[, which(colnames(table) == identity.col)]
				if ( !identity.val %in% identity.list ) return(NULL)

				subtable <- table[which(identity.list == identity.val),]

			} else if ( 0 != length(table) ) {

				# Single-row table
				if ( !identity.col %in% names(table) ) return(NULL)
				if ( !identity.val == table[which(names(table)) == identity.col] ) return(NULL)
				subtable <- table

			} else {

				# Empty table
				subtable <- NULL

			}

			# END #
			return(subtable)
		},

		apply.fun.based.on.identity = function (table, identity.col, behaviors,
			add.count, add.count.label) {
			# Applies preset function to table columns based on identity
			# 
			# Args:
			# 	table
			# 	identity.col
			# 	behaviors: the preset functions, a list with k=col and v=preset
			# 	add.count: whether to add a row.count column {Boolean}
			# 	add.count.label
			# 
			# Returns:
			# 	The updated/shrinked table
			
			end.table <- c()

			if ( !is.null(nrow(table)) ) {

				# Non-empty table
				if ( !identity.col %in% colnames(table) ) return(NULL)

				identity.list <- table[, which(colnames(table) == identity.col)]
				identity.list.unique <- unique(identity.list)

				for (id in identity.list.unique) {
					subtable <- GraphManager()$extract.subtable.based.on.identity(
						table, identity.col, id)

					if ( !is.null(subtable) ) {
						if ( !is.null(nrow(subtable)) ) {

							# Non-empty subtable
							single.row <- c()

							# Remove cols from behaviors
							if ( 'id' %in% names(behaviors) ) {
								behaviors['id'] <- NULL
							}
							if ( 'source' %in% names(behaviors) ) {
								behaviors['source'] <- NULL
							}
							if ( 'target' %in% names(behaviors) ) {
								behaviors['target'] <- NULL
							}

							# Apply behaviors
							for (col in names(behaviors)) {
								if ( col %in% colnames(subtable) ) {
									if ( 'sum' == behaviors[col] ) {
										single.row <- append(single.row,
											sum(as.numeric(subtable[,
												which(colnames(subtable) == col)])))
									} else if ( 'prod' == behaviors[col] ) {
										single.row <- append(single.row,
											prod(as.numeric(subtable[,
												which(colnames(subtable) == col)])))
									} else if ( 'min' == behaviors[col] ) {
										single.row <- append(single.row,
											min(as.numeric(subtable[,
												which(colnames(subtable) == col)])))
									} else if ( 'max' == behaviors[col] ) {
										single.row <- append(single.row,
											max(as.numeric(subtable[,
												which(colnames(subtable) == col)])))
									} else if ( 'random' == behaviors[col] ) {
										single.row <- append(single.row,
											subtable[round(runif(1, 0, 1) * nrow(subtable)),
											which(colnames(subtable) == col)])
									} else if ( 'mean' == behaviors[col] ) {
										single.row <- append(single.row,
											mean(as.numeric(subtable[,
												which(colnames(subtable) == col)])))
									} else if ( 'median' == behaviors[col] ) {
										single.row <- append(single.row,
											median(as.numeric(subtable[,
												which(colnames(subtable) == col)])))
									} else if ( 'concat' == behaviors[col] ) {
										single.row <- append(single.row,
											paste(subtable[,
												which(colnames(subtable) == col)], collapse='_'))
									} else if ( 'first' == behaviors[col] ) {
										single.row <- append(single.row,
											subtable[1, which(colnames(subtable) == col)])
									} else if ( 'last' == behaviors[col] ) {
										single.row <- append(single.row,
											subtable[nrow(subtable),
											which(colnames(subtable) == col)])
									}
									if ( 'ignore' != behaviors[col] ) {
										names(single.row)[length(single.row)] <- col
									}
								}
							}

							# Add row count column
							if ( add.count ) {
								single.row <- append(single.row, nrow(subtable))
								names(single.row)[length(single.row)] <- add.count.label
							}

							# Add missing attributes
							for (col in colnames(subtable)) {
								if ( !col %in% names(behaviors) ) {
									col.id <- which(col == colnames(subtable))
									single.row <- append(single.row, subtable[1, col.id])
									names(single.row)[length(single.row)] <- col
								}
							}

							# Re-order columns, just in case
							single.row <- GraphManager()$sort.table.cols(single.row)

						} else if ( 0 != length(subtable) ) {

							# Single-row subtable
							single.row <- subtable

							# Remove cols from behaviors
							if ( 'id' %in% names(behaviors) ) {
								behaviors['id'] <- NULL
							}
							if ( 'source' %in% names(behaviors) ) {
								behaviors['source'] <- NULL
							}
							if ( 'target' %in% names(behaviors) ) {
								behaviors['target'] <- NULL
							}

							# Remove ignored cols
							cols.to.rm <- c()
							for (col in names(single.row)) {
								if ( col %in% names(behaviors) ) {
									if ( 'ignore' == behaviors[col] ) {
										cols.to.rm <- append(cols.to.rm, col)
									}
								}
							}
							single.row <- GraphManager()$rm.cols(single.row, cols.to.rm)

							# Add row count column
							if ( add.count ) {
								single.row <- append(single.row, 1)
								names(single.row)[length(single.row)] <- add.count.label
							}

							# Re-order columns, just in case
							single.row <- GraphManager()$sort.table.cols(single.row)
						}
						end.table <- rbind(end.table, single.row)

					}

				}

			} else if ( 0 != length(table) ) {

				# Single-row table
				single.row <- table

				# Remove ignored cols
				cols.to.rm <- c()
				for (col in names(single.row)) {
					if ( col %in% names(behaviors) ) {
						if ( 'ignore' == behaviors[col] ) {
							cols.to.rm <- append(cols.to.rm, col)
						}
					}
				}
				single.row <- GraphManager()$rm.cols(single.row, cols.to.rm)

				# Add row count column
				if ( add.count ) {
					single.row <- append(single.row, 1)
					names(single.row)[length(single.row)] <- add.count.label
				}

				# Re-order columns, just in case
				single.row <- GraphManager()$sort.table.cols(single.row)

				end.table <- single.row

			}

			# END #
			return(end.table)
		},

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
			
			for(old in names(map)) g <- GraphManager()$rename.vertex.attr(g, old, eval(parse(text=paste0('map$', old))))
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
			
			for(old in names(map)) g <- GraphManager()$rename.edge.attr(g, old, eval(parse(text=paste0('map$', old))))
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

			if(vertex.map != list()) g <- GraphManager()$rename.vertex.attributes(g, vertex.map)
			if(edge.map != list()) g <- GraphManager()$rename.edge.attributes(g, edge.map)
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
				return(GraphManager()$clusteringCoefficient(v, env, graph))
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

			if (is.directed(g.one)) g.one <- GraphManager()$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- GraphManager()$undirected.noAttr(g.two)

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
			K <- (max.v * (max.v - 1))/2
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

			if (is.directed(g.one)) g.one <- GraphManager()$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- GraphManager()$undirected.noAttr(g.two)

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

			if (is.directed(g.one)) g.one <- GraphManager()$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- GraphManager()$undirected.noAttr(g.two)

			# Get edges
			el.one <- get.edgelist(g.one)
			el.two <- get.edgelist(g.two)

			# Get number of common edges
			common <- length(intersect(paste0(el.one[,1], '~', el.one[,2]), paste0(el.two[,1], '~', el.two[,2])))
			common <- common + length(intersect(paste0(el.one[,2], '~', el.one[,1]), paste0(el.two[,1], '~', el.two[,2])))

			# Calculate similarity and than distance
			num <- common
			den <- min(length(el.one[,1]), length(el.two[,1]))
			if ( den == 0 ) {
				if ( length(el.one[,1]) == length(el.two[,1]) ) {
					JS <- 1
				} else {
					JS <- 0
				}
			} else {
				JS <- num / den
			}
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
			
			if (is.directed(g.one)) g.one <- GraphManager()$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- GraphManager()$undirected.noAttr(g.two)

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
			
			if (is.directed(g.one)) g.one <- GraphManager()$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- GraphManager()$undirected.noAttr(g.two)

			dH <- GraphManager()$calcHammingDist(g.one, g.two)
			dIM <- GraphManager()$calcIpsenDist(g.one, g.two)
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

			if (is.directed(g.one)) g.one <- GraphManager()$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- GraphManager()$undirected.noAttr(g.two)

			dJ <- GraphManager()$calcJaccardDist(g.one, g.two)
			dIM <- GraphManager()$calcIpsenDist(g.one, g.two)
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

			if (is.directed(g.one)) g.one <- GraphManager()$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- GraphManager()$undirected.noAttr(g.two)

			dJS <- GraphManager()$calcJaccardSubsetDist(g.one, g.two)
			dIM <- GraphManager()$calcIpsenDist(g.one, g.two)
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
			
			if (is.directed(g.one)) g.one <- GraphManager()$undirected.noAttr(g.one)
			if (is.directed(g.two)) g.two <- GraphManager()$undirected.noAttr(g.two)
			
			dH <- GraphManager()$calcHammingDist(g.one, g.two)
			dJ <- GraphManager()$calcJaccardDist(g.one, g.two)
			dJS <- GraphManager()$calcJaccardSubsetDist(g.one, g.two)
			dIM <- GraphManager()$calcIpsenDist(g.one, g.two)
			dHIM <- (1/sqrt(1+xi)) * sqrt(dH**2 + xi * dIM**2)
			dJIM <- (1/sqrt(1+xi)) * sqrt(dJ**2 + xi * dIM**2)
			dJIMS <- (1/sqrt(1+xi)) * sqrt(dJS**2 + xi * dIM**2)
			return(c(dH, dJ, dJS, dIM, dHIM, dJIM, dJIMS))
		}

	)

	# Assign class attribute
	class(nm) <- 'GraphManager'

	# Return instantiated Graph Manager
	return(nm)
}
