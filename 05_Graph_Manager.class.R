# Class to build graphs
GraphBuilder <- function(clusters=0, verbose=FALSE) {
	library('igraph')
	library('doParallel')

	# Define GraphBuilder attributes
	gb <- list(
		clusters = clusters,
		verbose = verbose,
		readData = function(file.path, header=TRUE, sep='\t', row.names=NULL, sample.column=NULL) {
			# Reads data
			#
			# Args:
			#   file.path: the name of the file which the data are to be
			#		read from or a a readable text-mode connection.
			#	header: a logical value indicating whether the file contains
			#		the names of the variables as its first line.
			#	sep: the field separator character.
			#	row.names: a vector of row names.
			#	sample.column: if file is multi-sample, specify sample
			#		column header.
			#
			# Returns:
			#   Modified GraphBuilder instance

			# Read data
			if(gb$verbose) print('Reading whole data')
			data <- read.table(file.path, header=header, sep=sep, row.names=row.names)
			# Save raw data into GraphBuilder instance
			gb$data <- data

			# Check for multi-sampling
			if(is.null(sample.column)) {
				# If single-sample file
				return(gb)
			} else {
				# If multi-sample file
				gb <- gb$splitData(data, sample.column)
				return(gb)
			}
		},
		splitData = function(data, sample.column, clusters=gb$clusters) {
			# Splits multi-sample data
			#
			# Args:
			#   data: the multi-sample data to split.
			#	sample-column: the name of the sample-id column.
			#	clusters: the number .
			#
			# Returns:
			#	Modified GraphBuilder instance

			if(gb$verbose) print('Splitting whole data')
			# Prepare sample list without duplicates
			sample.list <- unique(data[,sample.column])
			# If needed, create output directory
			if (!file.exists('./sample-data/')) {
				dir.create(file.path('./sample-data/'))
			}
			# Declare parallelism
			cl <- makeCluster(clusters)
			registerDoParallel(cl)
			# Split original data table for each sample in data.frames
			# and save them in temporary directory
			foreach(i = 1:length(sample.list)) %dopar% {
				# On which sample are we working?
				sample.id <- sample.list[i]
				# Get row_ids from original data table for the working_sample
				row.ids <- which(data[,sample.column]==sample.id)
				# Write selected rows
				write.table(data[row.ids,], file=file.path('./sample-data/', sample.list[i]))
			}
			stopCluster(cl)
			# Terminate
			if(gb$verbose) print('Data splitted')
			gb$split <- TRUE
			gb$isMultiSample <- TRUE
			gb$data <- data
			return(gb)
		}
	)

	# Explicitely define GraphBuilder class
	class(gb) <- 'GraphBuilder'

	# Return the new GraphBuilder instance
	return(gb)
}

# Class to manage graphs
GraphManager <- function(clusters=0, verbose=FALSE) {

	# Define GraphManager attributes
	gm <- list(
		builder=GraphBuilder(clusters=clusters, verbose=verbose), # GraphBuilder instance
		try=function() print(1)
	)

	# Explicitely define GraphManager class
	class(gm) <- c('GraphManager','GraphBuilder')

	# Return the new GraphManager instance
	return(gm)
}