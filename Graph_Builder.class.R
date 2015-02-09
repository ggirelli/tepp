# Class to build and manage graphs
GraphBuilder <- function(
	clusters=0,
	verbose=FALSE,
	genes.label="Gene.id",
	white.list=list(),
	black.list=list(),
	clonal.val=c('clonal'),
	subclonal.val=c('subclonal'),
	attr.table='',
	clean=FALSE,
	output.dir='.',
	write.cooc=FALSE,
	pathToClass='.'
) {
	library('igraph')
	library('parallel')

	if(output.dir != '.')
		if(!file.exists(output.dir))
			dir.create(file.path(output.dir))
	
	# Define GraphBuilder attributes
	gb <- list(

		# ---------- #
		# ATTRIBUTES #
		# ---------- #
		
		version = 12,

		#------------------#
		# INPUT PARAMETERS #
		#------------------#

		# output directory
		output.dir = output.dir,

		# Number of cluster for parallel computing
		clusters = clusters,
		
		# Boolean value to determine verbosity
		verbose = verbose,

		# Gene.id label
		genes.label = genes.label,

		# Gene white/black-list
		white.list = white.list,
		black.list = black.list,

		# Values of clonality.status to identify clonal/subclonal genes
		clonal.val = clonal.val,
		subclonal.val = subclonal.val,

		# Clean: keeps only whitelisted genes
		clean = clean,

		# Wether to write the co-occurrency graphs in graphml format
		write.cooc = write.cooc,

		# Sample list
		sample.list = list(),

		# Vertex attribute table
		attr.table = attr.table,

		# Absolute path to GB class
		pathToClass = pathToClass,

		#
		nodeInNameSep = '~',
		
		#-----------#
		# FUNCTIONS #
		#-----------#

		# Reads data
		readData = function(
			file.path,
			header=TRUE,
			sep='\t',
			row.names=NULL,
			sample.column=NULL,
			genes.label=gb$genes.label,
			white.list=gb$white.list,
			black.list=gb$black.list,
			abe.type='dummy',
			temp.sample.list=list(),
			clean=gb$clean
		) {
			#
			# Args:
			#   file.path: the name of the file which the data are to be
			#   read from or a a readable text-mode connection.
			# header: a logical value indicating whether the file contains
			#   the names of the variables as its first line.
			# sep: the field separator character.
			# row.names: a vector of row names.
			# sample.column: if file is multi-sample, specify sample
			#   column header.
			#
			# Returns:
			#   Modified GraphBuilder instance
			
			# Read data
			if(gb$verbose) cat('Reading', abe.type, 'data\n')
			data <- read.table(file.path, header=header, sep=sep, row.names=row.names)

			# Save raw data into GraphBuilder instance
			gb$data <- data
			
			# Check for multi-sampling
			if( is.null(sample.column)) {

				# Clean single-sample data
				if(clean && length(white.list) != 0) {
					for(i in seq(length(data[,1]))) {
						id <- eval(parse(text=paste0('data$', genes.label, '[i]')))
						if(!(id %in% white.list)) data <- data[-i]
					}
				}
				# Blacklisting
				if(length(black.list) != 0) {
					for(i in seq(length(data[,1]))) {
						id <- eval(parse(text=paste0('data$', genes.label, '[i]')))
						if(id %in% white.list) data <- data[-i]
					}
				}

				# If single-sample file
				return(gb)

			} else {

				# If multi-sample file
				gbnew <- gb$splitData(
					data, sample.column,
					genes.label=genes.label,
					white.list=white.list,
					black.list=black.list,
					abe.type=abe.type,
					temp.sample.list=temp.sample.list,
					clean=clean
				)
				return(gbnew)

			}
		},

		splitSingleSampleData=function(
			i, data, sample.list,
			sample.column, data.dir,
			clean=gb$clean,
			white.list=gb$white.list,
			black.list=gb$black.list,
			genes.label=gb$genes.label
		) {
			# Splits original data table for each sample in data.frames
			# and save them in temporary directory
			# 
			# Args:
			# 	i: sample index
			# 	data:
			# 	sample.list:
			# 	sample.column:
			# 	clean:
			# 	white.list:
			# 	black.list:
			# 	genes.label:
			# 	data.dir:
			# 	
			# Returns:
			# 	NULL
			# 
			
			# On which sample are we working?
			sample.id <- sample.list[i]
			# Get row_ids from original data table for the working_sample
			row.ids <- which(data[,sample.column] == sample.id)
			# Clean sample
			if(clean && length(white.list) != 0)
                row.ids <- intersect(row.ids, which(data[,genes.label] %in% white.list))
			# Blacklisting
			if(length(black.list) != 0)
                row.ids <- intersect(row.ids, which(!(data[,genes.label] %in% black.list)))
			# Write selected rows
			if(length(row.ids) != 0)
                write.table(data[row.ids,], file=file.path(data.dir, sample.id))
		},
		
		# Splits data
		splitData = function(
			data, sample.column,
			clusters=gb$clusters,
			genes.label=gb$genes.label,
			white.list=gb$white.list,
			black.list=gb$black.list,
			abe.type='dummy',
			temp.sample.list=list(),
			clean=gb$clean,
			output.dir=gb$output.dir
		) {
			# Splits multi-sample data
			#
			# Args:
			#   data: the multi-sample data to split.
			# sample-column: the name of the sample-id column.
			# clusters: the number .
			#
			# Returns:
			# Modified GraphBuilder instance
			if(gb$verbose) cat('Splitting', abe.type, 'data\n')

			# Prepare sample list without duplicates
			sample.list <- unique(data[,sample.column])

			# If needed, create $output directory
			data.dir <- paste0(output.dir, '/sample-data-', abe.type, '/')
			if (!file.exists(data.dir)) {
				dir.create(file.path(data.dir))
			}

			mclapply(1:length(sample.list),
				FUN=gb$splitSingleSampleData,
				data, sample.list, sample.column, data.dir,
				clean, white.list, black.list,
				genes.label,
				mc.cores=clusters,
				mc.preschedule=TRUE
			)

			# Terminate
			if(gb$verbose) cat('Splitted', abe.type, 'data\n')
			gb$split <- TRUE
			gb$isMultiSample <- TRUE
			gb$sample.column <- sample.column
			gb$data <- data
			gb$sample.list <- unique(c(as.character(temp.sample.list), as.character(sample.list)))
			return(gb)
		},

		buildFinalSSMA = function(
			abe.list,
			genes.label=gb$genes.label,
			clonality.label="clonality.status",
			clonal.val=gb$clonal.val,
			subclonal.val=gb$subclonal.val,
			sample.list=gb$sample.list,
			output.dir=gb$output.dir,
			pathToClass=gb$pathToClass
		) {
			# Build SSMAs for final MSMA
			#
			# Args:
			#   abe.list: list of aberration types
			#   genes.label: label of the gene's id column.
			#   clonality.label: label of the clonality status column.
			#   clonal.val: value of clonality for clonal genes.
			#   subclonal.val: value of clonality for subclonal genes.
			#   sample.list: list of samples to analyze
			#   output.dir:
			#   pathToClass:
			#
			# Returns:
			#   None

			buildSSMA=function(sample.id,
				pathToClass, output.dir, abe.list,
				genes.label, clonality.label,
				clonal.val, subclonal.val
			) {
				# Execute buildSSMA for each sample
				# 
				# Args:
				# 	sample.id:
				# 	
				# Returns:
				# 	NULL
				# 	
				
				# Load library and classes
				library('igraph')
				source(file.path(pathToClass, 'Graph_Builder.class.R'))

				# Read data
				data <- list()
				for (i in 1:length(abe.list)) {
					abe <- abe.list[i]
					f.name <- eval(parse(text=paste0(
                        '"', output.dir, '/sample-data-', abe, '/', sample.id, '"'
                    )))
					if(file.exists(f.name)) eval(parse(text=paste0(
                        'data$', abe, ' <- read.table(f.name , header=TRUE, sep=" ")'
                    )))
					
					# Clean data from duplicates
					data <- GraphBuilder()$rmDuplicatedGenes(
                        data, abe, clonality.label,
                        genes.label=genes.label,
                        clonal.val=clonal.val,
                        subclonal.val=subclonal.val
                    )
				}

				if(length(data) != 0) {
					# Make empty graph
					g <- graph.empty(directed=TRUE)
					g.clonal <- g
					g.subclonal <- g

					# Get clonals
					for(abe in abe.list) {

						genes <- eval(parse(text=paste0('data$', abe, '$', genes.label)))
						clonality <- eval(parse(text=paste0('data$', abe, '$', clonality.label)))
						aberration <- seq(length(genes))
						aberration[] <- tolower(abe)

						# Add to graph
                        attr.list <- list(
                        	name=paste0(as.character(genes), '~', aberration),
                            HUGO=as.character(genes),
                            clonality.status=as.character(clonality),
                            abe.type=aberration
                        )
						g <- add.vertices(g, length(genes), attr=attr.list)
					}

					# Remove vertices that are neither clonal nor subclonal
					g <- delete.vertices(g, V(g)[intersect(
                        which(!(levels(V(g)$clonality.status) %in% clonal.val)),
                        which(!(levels(V(g)$clonality.status) %in% subclonal.val))
                    )])

					# Prepare edges list
					g <- g + edges(c(t(expand.grid(
                        c(V(g)[clonality.status %in% clonal.val]),
                        c(V(g)[clonality.status %in% subclonal.val])
                    ))))
					E(g)$weight <- 1

					# If needed, create output directory
					data.dir <- paste0(output.dir, '/sample-graphs/')
					if (!file.exists(data.dir)) {
						dir.create(file.path(data.dir))
					}

					# Output graph
					write.graph(g, file.path(data.dir, paste0(sample.id, '.graphml')),
                        format='graphml')
				}
			}
			mclapply(sample.list,
				FUN=buildSSMA,
				pathToClass, output.dir, abe.list,
				genes.label, clonality.label,
				clonal.val, subclonal.val,
				mc.cores=clusters,
				mc.preschedule=TRUE
			)

		},

		buildClonalSSMA = function(
			abe.list,
			genes.label=gb$genes.label,
			clonality.label="clonality.status",
			clonal.val=gb$clonal.val,
			subclonal.val=gb$subclonal.val,
			sample.list=gb$sample.list,
			v.list=list(),
			output.dir=gb$output.dir,
			pathToClass=gb$pathToClass
		) {
			# Build SSMAs for clonality co-occurrency MSMAs
			#
			# Args:
			#   abe.list: list of aberration types
			#   genes.label: label of the gene's id column.
			#   clonality.label: label of the clonality status column.
			#   clonal.val: value of clonality for clonal genes.
			#   subclonal.val: value of clonality for subclonal genes.
			#   sample.list: list of samples to analyze
			#   v.list:
			#   output.dir:
			#   pathToClass:
			#
			# Returns:
			#   None
			
			buildSSMA=function(sample.id,
				pathToClass, output.dir,
				abe.list, v.list,
				genes.label, clonality.label,
				clonal.val, subclonal.val
			) {
				# Execute buildSSMA for each sample
				# 
				# Args:
				# 	sample.id:
				# 
				# Returns:
				# 	NULL
				# 
				
				library('igraph')
				source(file.path(pathToClass, 'Graph_Builder.class.R'))
				
				# Read data
				data <- list()
				for (i in 1:length(abe.list)) {
					abe <- abe.list[i]
					# For eacha berration type
					f.name <- paste0(output.dir, '/sample-data-', abe, '/', sample.id)
					if(file.exists(f.name)) {
						temp.data <- read.table(f.name)
						temp.data <- temp.data[which(
                            paste(as.character(eval(parse(text=paste0(
                                'temp.data$', genes.label
                            )))), tolower(abe), sep='~') %in% v.list
                        ),]
						eval(parse(text=paste0('data$', abe, ' <- temp.data')))
					}

					# Clean data from duplicates
					data <- GraphBuilder()$rmDuplicatedGenes(
                        data, abe, clonality.label,
                        genes.label=genes.label,
                        clonal.val=clonal.val,
                        subclonal.val=subclonal.val
                    )
				}

				if(length(data) != 0) {
					# Make empty graph
					g.clonal <- graph.empty(directed=TRUE)
					g.subclonal <- graph.empty(directed=TRUE)
					g.nonclonal <- graph.empty(directed=TRUE)

					# Get clonals
					for(abe in abe.list) {

						genes <- eval(parse(text=paste0('data$', abe, '$', genes.label)))
						clonality <- as.character(eval(parse(text=paste0(
                            'data$', abe, '$', clonality.label
                        ))))
						aberration <- seq(length(genes))
						aberration[] <- tolower(abe)

						# Add to graph
                        attr.list <- list(
                            name=paste0(as.character(genes[which(clonality %in% clonal.val)]),
                            	'~', aberration[which(clonality %in% clonal.val)]),
                            HUGO=as.character(genes[which(clonality %in% clonal.val)]),
                            abe.type=aberration[which(clonality %in% clonal.val)]
                        )
						g.clonal <- add.vertices(
                            g.clonal,
                            length(genes[which(clonality %in% clonal.val)]),
                            attr=attr.list
                        )
                        attr.list <- list(
                            name=paste0(as.character(genes[which(clonality %in% subclonal.val)]),
                            	'~', aberration[which(clonality %in% subclonal.val)]),
                            HUGO=as.character(genes[which(clonality %in% subclonal.val)]),
                            abe.type=aberration[which(clonality %in% subclonal.val)]
                        )
						g.subclonal <- add.vertices(
                            g.subclonal,
                            length(genes[which(clonality %in% subclonal.val)]),
                            attr=attr.list
                        )
                        attr.list <- list(
                            name=paste0(as.character(genes), '~', aberration),
                            HUGO=as.character(genes),
                            abe.type=aberration,
                            clonality=clonality
                        )
						g.nonclonal <- add.vertices(g.nonclonal, length(genes), attr=attr.list)
					}

					# Prepare edges list for clonals
					g.clonal <- g.clonal + edges(c(t(expand.grid(
                        V(g.clonal)$name, V(g.clonal)$name
                    ))))
					E(g.clonal)$weight <- 1

					# Prepare edges list for subclonals
					g.subclonal <- g.subclonal + edges(c(t(expand.grid(
                        V(g.subclonal)$name, V(g.subclonal)$name
                    ))))
					E(g.subclonal)$weight <- 1

					# Prepare edges list for nonclonals
					g.nonclonal <- g.nonclonal + edges(c(t(expand.grid(
                        V(g.nonclonal)[!(clonality %in% append(clonal.val, subclonal.val))]$name,
                        V(g.nonclonal)$name
                    ))))
					g.nonclonal <- g.nonclonal + edges(c(t(expand.grid(
                        V(g.nonclonal)$name,
                        V(g.nonclonal)[!(clonality %in% append(clonal.val, subclonal.val))]$name
                    ))))
					E(g.nonclonal)$weight <- 1

					# If needed, create $output directory
					data.dir <- paste0(output.dir, '/sample-graphs/')
					if (!file.exists(data.dir)) {
						dir.create(file.path(data.dir))
					}

					# Output graph
					if(length(V(g.clonal)) != 0) write.graph(
                        g.clonal,
                        file.path(data.dir, paste0('clonal_', sample.id, '.graphml')),
                        format='graphml'
                    )
					if(length(V(g.subclonal)) != 0) write.graph(
                        g.subclonal,
                        file.path(data.dir, paste0('subclonal_', sample.id, '.graphml')),
                        format='graphml'
                    )
					if(length(V(g.nonclonal)) != 0) write.graph(
                        g.nonclonal,
                        file.path(data.dir, paste0('nonclonal_', sample.id, '.graphml')),
                        format='graphml'
                    )
				}
			}
			mclapply(sample.list,
				FUN=buildSSMA,
				pathToClass, output.dir,
				abe.list, v.list,
				genes.label, clonality.label,
				clonal.val, subclonal.val,
				mc.cores=clusters,
				mc.preschedule=TRUE
			)

		},

		# Build MSMA
		buildMSMA = function(
			graph.list,
			directed=TRUE,
			output.dir=gb$output.dir,
			doOcc=FALSE,
			remove.loops=FALSE,
			clusters=gb$clusters
		) {
			# Merges SSMAs into a single MSMA summing the edge's weight
			#
			# Args:
			#   graph.list:
			#
			# Return:
			#
			
			getEdges=function(i,
				graph.list, output.dir
			) {
				# Gets the edges
				# 
				# Args:
				# 	i:
				# 	graph.list:
				# 	output.dir:
				# 
				# Returns:
				# 	NULL
				
				library('igraph')

				file.name <- graph.list[i]
				if(file.exists(file.path(output.dir, '/sample-graphs/', file.name))) {
					g <- read.graph(
                        file.path(output.dir, '/sample-graphs/', file.name),
                        format='graphml'
                    )
					edgelist <- get.edgelist(g)
					if(length(edgelist) != 0) {
						edgelist.n <- get.edgelist(g, name=FALSE)
						cbind(
                            V(g)[edgelist.n[,1]]$HUGO,
                            V(g)[edgelist.n[,2]]$HUGO,
                            E(g)$weight,
                            V(g)[edgelist.n[,1]]$abe.type,
                            V(g)[edgelist.n[,2]]$abe.type,
                            file.name
                        )
					}
				}
			}
			e.list <- mclapply(1:length(graph.list),
				FUN=getEdges,
				graph.list, output.dir,
				mc.cores=clusters,
				mc.preschedule=TRUE
			)
			edges <- do.call(rbind, e.list)
			colnames(edges) <- c(
                'source', 'target', 'weight',
                'source.abe.type', 'target.abe.type', 'sample'
            )

			if(length(edges) != 0) {

				# Build new graph
				if(gb$verbose) cat("Building empty graph\n")
				g <- graph.empty(directed=directed)

				# Vertices
				#(2) Then scroll the table (tricky) to get each node
                # and its attributes and add them to the MSMA.
				if(gb$verbose) cat("Preparing sources\n")
				source.table <- paste0(edges[,1], '~', edges[,4])

				if(gb$verbose) cat("Preparing targets\n")
				target.table <- paste0(edges[,2], '~', edges[,5])

				if(gb$verbose) cat("Preparing vertices\n")
				vertices.table <- unique(c(source.table, target.table))

				if(gb$verbose) cat("Adding vertices with attributes\n")
				tmpMat <- matrix(unlist(strsplit(vertices.table, '~')), nrow=2)
                attr.list <- list(
                    name=vertices.table,
                    aberration=tmpMat[2,],
                    HUGO=tmpMat[1,]
                )
				g <- add.vertices(g, nv=length(vertices.table), attr=attr.list)

				if(doOcc) {
					if(gb$verbose) cat("Retrieving vertices occurrencies\n")

					countNodePerSample=function(n.list,
						edges, sep, Ncores=clusters
					) {
						# Assembles the edges table and the list of source/target
						# 
						# Args:
						# 	n.list:
						# 	edges:
						# 	sep:
						# 	Ncores:
						# 	
						# Returns:
						# 	
						# 
						
						nlist <- mclapply(unique(paste0(n.list, sep, edges[,6])),
							FUN=function(x) {
								paste(unlist(strsplit(x, '~', fixed=T))[1:2], collapse='~')
							}, mc.preschedule=T,
							mc.cores=Ncores)
						return(table(unlist(nlist)))
					}
					# Retrieve clonal occurrencies
					source.counts <- countNodePerSample(source.table, edges, gb$nodeInNameSep)
					# Retrieve subclonal occurrencies
					target.counts <- countNodePerSample(target.table, edges, gb$nodeInNameSep)

					if(gb$verbose) cat("Adding vertices occurrencies\n")
					orderOccurrencies=function(name,
						source.counts, target.counts
					) {
						# Orders the occurrency data previously retrieved
						# 
						# Args:
						# 	name:
						# 	source.counts:
						# 	target.counts:
						# 
						# Returns:
						# 	
						# 
						
						if (!name %in% names(source.counts)) {
							clonal.occ <- 0
						} else {
							clonal.occ <- source.counts[name]
						}
						if (!name %in% names(target.counts)) {
							subclonal.occ <- 0
						} else {
							subclonal.occ <- target.counts[name]
						}
						return(c(clonal.occ, subclonal.occ))
					}
					occList <- mclapply(V(g)$name,
						FUN=orderOccurrencies, source.counts, target.counts,
						mc.preschedule=TRUE,
						mc.cores=clusters
					)
					occTab <- do.call(rbind, occList)
					V(g)$clonal.occ <- occTab[,1]
					V(g)$subclonal.occ <- occTab[,2]
				}

				# Edges
				if(gb$verbose) cat("Adding edges\n")
				#(3) Then add all the edges identifying the nodes based on their attributes.
				g <- g + edges(c(t(cbind(source.table, target.table))))
				if(gb$verbose) cat("Adding edges' weight\n")
				E(g)$weight <- as.numeric(edges[,3])

				#(4) Finally apply 'simplify' to remove multiple edges and sum their weights.
				if(gb$verbose) cat("Simplifying\n")
				g <- simplify(
                    g, remove.multiple=TRUE,
                    remove.loops=remove.loops,
                    edge.attr.comb=list(weight="sum", "ignore")
                )

				if(gb$verbose) cat('\nMaximum edge weight: ', max(E(g)$weight), '\n')
				if(gb$verbose) cat('Maximum vertex degree: ', max(degree(g, V(g))), '\n')

				# Terminate
				if(gb$verbose) cat('Graphs merged\n')

				return(g)

			} else {

				return(graph.empty())

			}
		},

		rmDuplicatedGenes.SCNA = function(
			data, abe, clonality.label,
			genes.label=gb$genes.label,
			clonal.val=gb$clonal.val,
			subclonal.val=gb$subclonal.val
		) {
			# Removes duplicated genes from the SCNA table
			#
            # Args:
            #   data: multi-table
            #   abe: aberration type [SCNA]
            #   clonality.label: label of clonality status column
            #   genes.label: label of gene ID column
            #   clonal.val: clonality status value(s) for 'clonal'
            #   subclonal.val: clonality status value(s) for 'subclonal'
            # 
            # Returns:
            #   The multi-table w/o duplicated genes in the SCNA tables
            #
			
			dups <- which(duplicated(eval(parse(text=paste0('data$', abe, '$', genes.label)))))
			while(length(dups) != 0) {
				i <- dups[1]
				dup.ids <- which(eval(parse(text=paste0(
                    'data$', abe, '$', genes.label, ' == data$', abe, '$', genes.label, '[i]'
                ))))
				sub.ids <- dup.ids[which(eval(parse(text=paste0(
                    'data$', abe, '$', clonality.label
                )))[dup.ids] %in% subclonal.val)]
				clo.ids <- dup.ids[which(eval(parse(text=paste0(
                    'data$', abe, '$', clonality.label
                )))[dup.ids] %in% clonal.val)]
				non.ids <- dup.ids[which(!eval(parse(text=paste0(
                    'data$', abe, '$', clonality.label
                )))[dup.ids] %in% union(clonal.val, subclonal.val))]
				if (length(sub.ids) != 0 & length(clo.ids) == 0) {
					# Keep one of the duplicates
					rm.ids <- dup.ids 
					rm.ids <- rm.ids[-which.max(eval(parse(text=paste0(
                        'data$', abe, '$perc.overlap[rm.ids]'
                    ))))]
					eval(parse(text=paste0(
                        'data$', abe, ' <- data$', abe, '[-rm.ids,]'
                    )))
				}
				if (length(sub.ids) == 0 & length(clo.ids) != 0) {
					# Keep one of the duplicates
					rm.ids <- dup.ids
					rm.ids <- rm.ids[-which.max(eval(parse(text=paste0(
                        'data$', abe, '$perc.overlap[rm.ids]'
                    ))))]
					eval(parse(text=paste0('data$', abe, ' <- data$', abe, '[-rm.ids,]')))
				}
				if (length(sub.ids) != 0 & length(clo.ids) != 0) {
					# Keep one of the subclonals
					sub.ids <- sub.ids[-which.max(eval(parse(text=paste0(
                        'data$', abe, '$perc.overlap[sub.ids]'
                    ))))]
					eval(parse(text=paste0('data$', abe, ' <- data$', abe, '[-clo.ids,]')))
					if(length(sub.ids) != 0) eval(parse(text=paste0(
                        'data$', abe, ' <- data$', abe, '[-sub.ids,]'
                    )))
				}
				if(length(sub.ids) == 0 & length(clo.ids) == 0) {
					# Remove non-clonals
					eval(parse(text=paste0('data$', abe, ' <- data$', abe, '[-non.ids,]')))
				}
				# Change counter
				dups <- which(duplicated(eval(parse(text=paste0('data$', abe, '$', genes.label)))))
			}
			return(data)
		},

		rmDuplicatedGenes.PM = function(
			data, abe, clonality.label,
			genes.label=gb$genes.label,
			clonal.val=gb$clonal.val,
			subclonal.val=gb$subclonal.val
		) {
            # Removes duplicated genes from the PM table
            #
            # Args:
            #   data: multi-table
            #   abe: aberration type [PM]
            #   clonality.label: label of clonality status column
            #   genes.label: label of gene ID column
            #   clonal.val: clonality status value(s) for 'clonal'
            #   subclonal.val: clonality status value(s) for 'subclonal'
            # 
            # Returns:
            #   The multi-table w/o duplicated genes in the PM table
            #
			
			dups <- which(duplicated(as.character(eval(parse(text=paste0(
                'data$', abe, '$', genes.label
            ))))))
			while(length(dups) != 0) {
				i <- dups[1]
				dup.ids <- which(eval(parse(text=paste0(
                    'data$', abe, '$', genes.label, ' == data$', abe, '$', genes.label, '[i]'
                ))))
				sub.ids <- dup.ids[which(eval(parse(text=paste0(
                    'data$', abe, '$', clonality.label
                )))[dup.ids] %in% subclonal.val)]
				clo.ids <- dup.ids[which(eval(parse(text=paste0(
                    'data$', abe, '$', clonality.label
                )))[dup.ids] %in% clonal.val)]
				non.ids <- dup.ids[which(!eval(parse(text=paste0(
                    'data$', abe, '$', clonality.label
                )))[dup.ids] %in% union(clonal.val, subclonal.val))]
				if (length(sub.ids) != 0 & length(clo.ids) == 0) {
					# Keep one of the duplicates
					rm.ids <- dup.ids 
					rm.ids <- rm.ids[-1]
					eval(parse(text=paste0('data$', abe, ' <- data$', abe, '[-rm.ids,]')))
				}
				if (length(sub.ids) == 0 & length(clo.ids) != 0) {
					# Keep one of the duplicates
					rm.ids <- dup.ids
					rm.ids <- rm.ids[-1]
					eval(parse(text=paste0('data$', abe, ' <- data$', abe, '[-rm.ids,]')))
				}
				if (length(sub.ids) != 0 & length(clo.ids) != 0) {
					# Keep one of the subclonals
					sub.ids <- sub.ids[-1]
					eval(parse(text=paste0('data$', abe, ' <- data$', abe, '[-clo.ids,]')))
					if(length(sub.ids) != 0) eval(parse(text=paste0(
                        'data$', abe, ' <- data$', abe, '[-sub.ids,]'
                    )))
				}
				if(length(sub.ids) == 0 & length(clo.ids) == 0) {
					# Remove non-clonals
					eval(parse(text=paste0('data$', abe, ' <- data$', abe, '[-non.ids,]')))
				}
				# Change counter
				dups <- which(duplicated(eval(parse(text=paste0('data$', abe, '$', genes.label)))))
			}
			return(data)
		},

		rmDuplicatedGenes = function(
			data, abe, clonality.label,
			genes.label=gb$genes.label,
			clonal.val=gb$clonal.val,
			subclonal.val=gb$subclonal.val
		) {
            # Removes duplicated genes from the [SCNA|PM] table
            #
            # Args:
            #   data: multi-table
            #   abe: aberration type [SCNA|PM]
            #   clonality.label: label of clonality status column
            #   genes.label: label of gene ID column
            #   clonal.val: clonality status value(s) for 'clonal'
            #   subclonal.val: clonality status value(s) for 'subclonal'
            # 
            # Returns:
            #   The multi-table w/o duplicated genes
            #   
            # @TODO:
            #   Implement rmDuplicated for RR
            # 
			
			gb <- GraphBuilder(
                genes.label=gb$genes.label,
                clonal.val=gb$clonal.val,
                subclonal.val=gb$subclonal.val
            )
			if(abe == 'PM') {
				return(gb$rmDuplicatedGenes.PM(data, abe, clonality.label))
			} else if(abe == 'Loss' | abe == 'Gain') {
				return(gb$rmDuplicatedGenes.SCNA(data, abe, clonality.label))
			} else {
				return(data)
			}
		},

		setVAttributes = function(
			graph, attr.table,
			key.label="HUGO"
		) {
			# Adds vertex attributes to a graph object
			#
			# Args:
			#   graph: the graph
			#   attr.file: the path to the file
			#   key.label: the column/vertex label used to assign the attributes
			#
			# Returns:
			#   The graph with vertex attributes

			# Get key.attribute values
			sources <- eval(parse(text=paste0('attr.table$', key.label)))
			targets <- unique(eval(parse(text=paste0('V(graph)$', key.label))))

			# Iterate through target key.values (HUGO)
			for(key in targets[which(targets %in% sources)]) {
				# Iterate through new attributes to assign
				for(attr in colnames(attr.table)[which(colnames(attr.table) != key.label)]) {
					eval(parse(text=paste0(
                        'V(graph)[', key.label, ' == key]$', attr,
                        ' <- attr.table$', attr, '[which(sources == key)]'
                    )))
				}
			}

			return(graph)
		},

		# Manage building
		build = function(
			abe.list,
			genes.label=gb$genes.label,
			clonality.label="clonality.status",
			clonal.val=gb$clonal.val,
			subclonal.val=gb$subclonal.val,
			sample.list=gb$sample.list,
			attr.table=gb$attr.table,
			output.dir=gb$output.dir,
			write.cooc=gb$write.cooc
		) {
			# Builds graphs after data read/split
			#
			# Args:
			#   abe.list: list of aberration types
			#   genes.label: label of the gene's id column.
			#   clonality.label: label of the clonality status column.
			#   clonal.val: value of clonality for clonal genes.
			#   subclonal.val: value of clonality for subclonal genes.
			#   sample.list: list of samples to analyze
			#
			# Returns:
			#   None
			
			if(gb$verbose) cat('# Preparing SSMAs.\n')
			gb$buildFinalSSMA(
				abe.list,
				genes.label=genes.label,
				clonality.label=clonality.label,
				clonal.val=clonal.val,
				subclonal.val=subclonal.val,
				sample.list=sample.list
			)
			if(gb$verbose) cat('SSMAs prepared.\n')

			if(gb$verbose) cat("\n# Merging SSMAs into MSMA\n")
			g.total <- gb$buildMSMA(paste0(sample.list, '.graphml'), doOcc=TRUE)

			if(gb$verbose) cat('\n# Preparing Clonality SSMAs.\n')
			gb$buildClonalSSMA(
				abe.list,
				genes.label=genes.label,
				clonality.label=clonality.label,
				clonal.val=clonal.val,
				subclonal.val=subclonal.val,
				sample.list=sample.list,
				v.list=V(g.total)$name
			)
			if(gb$verbose) cat('SSMAs prepared.\n')

			if(gb$verbose) cat("\n# Merging SSMAs into MSMA · Clonal co-occurrency\n")
			g.clonal <- gb$buildMSMA(paste0('clonal_', sample.list, '.graphml'), remove.loops=TRUE)
			if(gb$verbose) cat("\n# Merging SSMAs into MSMA · Subclonal co-occurrency\n")
			g.subclonal <- gb$buildMSMA(paste0('subclonal_', sample.list, '.graphml'), remove.loops=TRUE)
			if(gb$verbose) cat("\n# Merging SSMAs into MSMA · Uncertain_clonality co-occurrency\n")
			g.nonclonal <- gb$buildMSMA(paste0('nonclonal_', sample.list, '.graphml'), remove.loops=TRUE)

			if(gb$verbose) cat("\n# Retrieving co-occurrency data\n")
			# Get which edges can be present in the co-occurrency graphs
			el.tot <- get.edgelist(g.total)
			e.in.clo <- (el.tot[,1] %in% V(g.clonal)$name & el.tot[,2] %in% V(g.clonal)$name)
			e.in.sub <- (el.tot[,1] %in% V(g.subclonal)$name & el.tot[,2] %in% V(g.subclonal)$name)
			e.in.non <- (el.tot[,1] %in% V(g.nonclonal)$name & el.tot[,2] %in% V(g.nonclonal)$name)
			# Retrieve correct clonal co-occurrency data
			if(length(V(g.clonal)) != 0) {
				if(gb$verbose) cat(" · Retrieving clonal co-occurrency data\n")
				clonal.cooc <- rep(0, length(E(g.total)))
				clonal.cooc.ids <- get.edge.ids(
                    g.clonal,
                    t(get.edgelist(g.total)[which(e.in.clo),]),
                    error=FALSE
                )
				clonal.cooc[
                    which(e.in.clo)[which(clonal.cooc.ids != 0)]
                ] <- E(g.clonal)[clonal.cooc.ids]$weight
				E(g.total)$clonal.cooc <- clonal.cooc
			}
			# Retrieve correct subclonal co-occurrency data
			if(length(V(g.subclonal)) != 0) {
				if(gb$verbose) cat(" · Retrieving subclonal co-occurrency data\n")
				subclonal.cooc <- rep(0, length(E(g.total)))
				subclonal.cooc.ids <- get.edge.ids(
                    g.subclonal,
                    t(get.edgelist(g.total)[which(e.in.sub),]),
                    error=FALSE
                )
				subclonal.cooc[
                    which(e.in.sub)[which(subclonal.cooc.ids != 0)]
                ] <- E(g.subclonal)[subclonal.cooc.ids]$weight
				E(g.total)$subclonal.cooc <- subclonal.cooc
			}
			# Retrieve correct nonclonal co-occurrency data
			if(length(V(g.nonclonal)) != 0) {
				if(gb$verbose) cat(" · Retrieving uncertain_clonality co-occurrency data\n")
				nonclonal.cooc <- rep(0, length(E(g.total)))
				nonclonal.cooc.ids <- get.edge.ids(
                    g.nonclonal,
                    t(get.edgelist(g.total)[which(e.in.non),]),
                    error=FALSE
                )
				nonclonal.cooc[
                    which(e.in.non)[which(nonclonal.cooc.ids != 0)]
                ] <- E(g.nonclonal)[nonclonal.cooc.ids]$weight
				E(g.total)$nonclonal.cooc <- nonclonal.cooc
			}

			# Retrieve and assign new vertex attributes
			#if(length(attr.table) != 0 && attr.table != '') {
			#  if(gb$verbose) cat("\nAssigning new vertex attributes based on HUGO\n")
			#  g.total <- gb$setVAttributes(g.total, attr.table)
			#}

			if(gb$verbose) cat("\nWriting graph\n")
			write.graph(g.total, file.path(output.dir, 'total_graph.graphml'), format='graphml')

			if(!write.cooc) {
				E(g.clonal)$weight <- as.character(E(g.clonal)$weight)
				write.graph(
                    g.clonal,
                    file.path(output.dir, 'clonal_graph.lgl'),
                    format='lgl',
                    names='name',
                    weights='weight'
                )
				E(g.subclonal)$weight <- as.character(E(g.subclonal)$weight)
				write.graph(
                    g.subclonal,
                    file.path(output.dir, 'subclonal_graph.lgl'),
                    format='lgl', names='name', weights='weight')
				E(g.nonclonal)$weight <- as.character(E(g.nonclonal)$weight)
				write.graph(
                    g.nonclonal,
                    file.path(output.dir, 'nonclonal_graph.lgl'),
                    format='lgl',
                    names='name',
                    weights='weight'
                )
			} else {
				write.graph(
                    g.clonal,
                    file.path(output.dir, 'clonal_graph.graphml'),
                    format='graphml'
                )
				write.graph(
                    g.subclonal,
                    file.path(output.dir, 'subclonal_graph.graphml'),
                    format='graphml'
                )
				write.graph(
                    g.nonclonal,
                    file.path(output.dir, 'nonclonal_graph.graphml'),
                    format='graphml'
                )
			}

			if(gb$verbose) cat("\nFIN\n\n")
		}
	)
	
	# Explicitely define GraphBuilder class
	class(gb) <- 'GraphBuilder'
	
	# Return the new GraphBuilder instance
	return(gb)
}