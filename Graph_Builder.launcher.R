#!/usr/bin/env Rscript


# Read parameters
##############################
args <- commandArgs(trailingOnly=TRUE)

# Read param file
if ( file.exists(args[1]) ) {
	source(args[1])
} else {
	cat('\nNo param file detected.')
}


# Write log
##############################
source(file.path(pathToGBclass, 'Graph_Builder.class.R'))

cat('\nExecuting script (v', GraphBuilder()$version, ') with', Ncores, ' cores.\n')

if(verbose) cat('Verbosity at maximum.\n')
else cat('Verbosity at minimum.\n')

if(output.dir != '') cat('Output directory: ', output.dir, '\n\n')
else cat('\n')

cat('PM-data: ', file.list$PM, '\n')
cat('Gain-data: ', file.list$Gain, '\n')
cat('Loss-data: ', file.list$Loss, '\n')
cat('RR-data: ', file.list$RR, '\n\n')

cat('Gene label: \'', genes.label, '\'\n\n')
if(clean) cat('Running clean\n\n')

if(length(white.list) != 0 && white.list != '') {
	cat('White list from: ', toString(ff[which(ff[, 1] == 'white.list'), 2]), '\n')
	print(as.character(white.list))
} else cat('No white list specified\n')

if(length(black.list) != 0 && black.list != '') {
	cat('\nBlack list from: ', toString(ff[which(ff[, 1] == 'black.list'), 2]), '\n')
	print(as.character(black.list))
	cat('\n')
} else cat('\nNo black list specified\n\n')


cat('Clonal ids:\n')
print(clonal.val)
cat('Subclonal ids:\n')
print(subclonal.val)
cat('\n')

if(length(attr.table) != 0 && attr.table != '') {
	cat('Vertex attributes added to the final graph:\n')
	print(colnames(attr.table))
	cat('\n')
}






# Execute
##############################
system.time({

	# GRAPH BUILDER instance
	######################################
	
	gb <- GraphBuilder(
		clusters=Ncores,
		verbose=verbose,
		genes.label=genes.label,
		white.list=white.list,
		black.list=black.list,
		clonal.val=clonal.val,
		subclonal.val=subclonal.val,
		attr.table=attr.table,
		clean=clean,
		output.dir=output.dir,
		write.cooc=write.cooc
	)
	cat('\n')



	# MSSA -> SSSA
	######################################

	# PM #
	if(length(file.list$PM) != 0 & !is.na(file.list$PM)) {
		gb <- gb$readData(
			file.path=file.list$PM,
			sample.column='sample',
			abe.type='PM',
			temp.sample.list=gb$sample.list
		)
		cat('\n')
	}

	# Gain #
	if(length(file.list$Gain) != 0 & !is.na(file.list$Gain)) {
		gb <- gb$readData(
			file.path=file.list$Gain,
			sample.column='sample',
			abe.type='Gain',
			temp.sample.list=gb$sample.list
		)
		cat('\n')
	}

	# Loss #
	if(length(file.list$Loss) != 0 & !is.na(file.list$Loss)) {
		gb <- gb$readData(
			file.path=file.list$Loss,
			sample.column='sample',
			abe.type='Loss',
			temp.sample.list=gb$sample.list
		)
		cat('\n')
	}

	# RR #
	if(length(file.list$RR) != 0 & !is.na(file.list$RR)) {
		gb <- gb$readData(
			file.path=file.list$RR,
			sample.column='sample',
			abe.type='RR',
			temp.sample.list=gb$sample.list
		)
		cat('\n')
	}



	# SSSA -> SSMA -> MSMA
	######################################

	gb$build(
		c('PM', 'Gain', 'Loss'),
		genes.label="Gene.id",
		clonality.label="clonality.status",
		sample.list=gb$sample.list
	)

})
