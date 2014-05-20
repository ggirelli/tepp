library('igraph')
library('doParallel')

SyntheticGraphData = function() {

	# Instantiate Graph Manager
	sgd <- list(

		#------------#
		# ATTRIBUTES #
		#------------#



		#-----------#
		# FUNCTIONS #
		#-----------#
		
		generateFullData = function(n, directed=TRUE, nsample=1, limit=100) {
			cat('Preparing full graph','\n')
			g <- graph.full(n, directed=directed)
			V(g)$name <- paste0('dummy', c(V(g)))

			cat('Preparing empty sample table','\n')
			sample.list <- list(list(clonal = c(), subclonal = c()))
			edge.list <- get.edgelist(g)
			edge.list.length <- length(E(g))

#			for(i in 1:edge.list.length) {
#				sample.target <- sample(1:nsample, 1)
#				counter <- 0
#				while((edge.list[i,1] %in% sample.list[[sample.target]]$subclonal || edge.list[i,2] %in% sample.list[[sample.target]]$clonal) && counter <= 100) {
#					sample.target <- sample(1:nsample, 1)
#					counter <- counter + 1
#				}
#				if(counter == limit) {
#					cat('Error: infinite loop, limit reached.', '\n')
#					return(NULL)
#				} else {
#					cat('Assigning edge ', i, ' of ', edge.list.length, ' to sample ', sample.target, '\n')
#					if(!(edge.list[i,1] %in% sample.list[[sample.target]]$clonal)) {
#						sample.list[[sample.target]]$clonal <- append(sample.list[[sample.target]]$clonal, edge.list[i,1])
#					}
#					if(!(edge.list[i,2] %in% sample.list[[sample.target]]$subclonal)) {
#						sample.list[[sample.target]]$subclonal <- append(sample.list[[sample.target]]$subclonal, edge.list[i,2])
#					}
#				}
#			}

			for(i in sample(1:edge.list.length)) {
				sample.target <- 1
				while((edge.list[i,1] %in% sample.list[[sample.target]]$subclonal || edge.list[i,2] %in% sample.list[[sample.target]]$clonal) && counter <= 100) {
					sample.target <- sample.target + 1
					if(sample.target > length(sample.list)) sample.list <- append(sample.list, list(list(clonal = c(), subclonal = c())))
				}
				cat('Assigning edge ', i, ' of ', edge.list.length, ' to sample ', sample.target, '\n')
				if(!(edge.list[i,1] %in% sample.list[[sample.target]]$clonal)) {
					sample.list[[sample.target]]$clonal <- append(sample.list[[sample.target]]$clonal, edge.list[i,1])
				}
				if(!(edge.list[i,2] %in% sample.list[[sample.target]]$subclonal)) {
					sample.list[[sample.target]]$subclonal <- append(sample.list[[sample.target]]$subclonal, edge.list[i,2])
				}
			}

			retro.table <- c()
			for(i in 1:length(sample.list)) {
				genes <- sample.list[[1]]$clonal
				status <- genes; status[] <- 'clonal'
				sample <- genes; sample[] <- paste0('sample', i)
				retro.table <- rbind(retro.table, cbind(genes,status,sample))

				genes <- sample.list[[1]]$subclonal
				status <- genes; status[] <- 'subclonal'
				sample <- genes; sample[] <- paste0('sample', i)
				retro.table <- rbind(retro.table, cbind(genes,status,sample))
			}

			rownames(retro.table) <- 1:length(retro.table[,1])
			colnames(retro.table) <- c('Gene.id', 'clonality.status', 'sample')

			return(retro.table)
		}

	)

	# Assign class attribute
	class(sgd) <- 'SyntheticGraphData'

	# Return instantiaded Graph Manager
	return(sgd)
}
