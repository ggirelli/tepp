#!/usr/bin/env Rscript

# Read parameters
##############################
args <- commandArgs(trailingOnly=TRUE)

# Read param file
if ( file.exists(args[1]) ) {
	source(args[1])
} else {
	message('No param file detected.')
}

# Load libraries
library(igraph)
library(parallel)
library(data.table)



# Write log
##############################
source(file.path(pathToGBclass, 'Graph_Builder.class.R'))

message(paste0('> Testing script (v', GraphBuilder()$version, ') with ', Ncores, ' cores.'))

if(verbose) {
	message('> Verbosity at maximum.')
} else {
	message('> Verbosity at minimum.')
}
if(output.dir != '') {
	message('> Output directory: ', output.dir, '\n')
} else {
	message('\n')
}

if ( test.cleaned == TRUE ) {
	file.list=clean.list
	message('> Using *CLEAN* version of the tables.')
}
message('> PM-data: ', file.list$PM)
message('> Gain-data: ', file.list$Gain)
message('> Loss-data: ', file.list$Loss)
message('> RR-data: ', file.list$RR)

message('> Gene label: \'', genes.label, '\'\n')
if(clean) message('Running clean')

if(length(white.list) != 0 && white.list != '') {
	message('> White list from: ', toString(ff[which(ff[, 1] == 'white.list'), 2]))
	print(as.character(white.list))
} else {
	message('> No white list specified')
}

if(length(black.list) != 0 && black.list != '') {
	message('> Black list from: ', toString(ff[which(ff[, 1] == 'black.list'), 2]))
	print(as.character(black.list))
} else {
	message('> No black list specified')
}


message('\n> Clonal ids:')
print(clonal.val)
message('> Subclonal ids:')
print(subclonal.val)

if(length(attr.table) != 0 && attr.table != '') {
	message('\n> Vertex attributes added to the final graph:')
	print(colnames(attr.table))
}
message('')

# Clean Basic Data
######################################

rmDups.SCNA = function (table,
	sample.column, clonality.label, genes.label,
	clonal.val, subclonal.val, Ncores
) {

	samples <- eval(parse(text=paste0('table$', sample.column)))
	tables <- mclapply(unique(samples),
		FUN=function(sample, table, samples) {
			table.sample <- table[which(samples == sample),]

			gene.list <- as.character(eval(parse(text=paste0('table.sample$', genes.label))))
			gene.dups <- unique(gene.list[which(duplicated(gene.list))])

			findDupIDs.SCNA = function(gene,
				data=table.sample,
				clonality.label=gb$clonality.label,
				genes.label=gb$genes.label,
				clonal.val=gb$clonal.val,
				subclonal.val=gb$subclonal.val,
				perc.overlap.label='perc.overlap'
			) {
				# 
				# Args:
				# 	gene: a gene to check for duplicates
				# 	data: the data table
				# 
				# Return:
				# 	A list of IDs to be removed
				# 	
				
				all.val <- union(clonal.val, subclonal.val)

				genes <- as.character(eval(parse(text=paste0('data$', genes.label))))
				clonality.status <- eval(parse(text=paste0('data$', clonality.label)))
				perc.overlap <- eval(parse(text=paste0('data$', perc.overlap.label)))

				dup.ids <- which(genes == gene)
				dup.sub.ids <- dup.ids[which(clonality.status[dup.ids] %in% subclonal.val)]
				dup.clo.ids <- dup.ids[which(clonality.status[dup.ids] %in% clonal.val)]
				dup.non.ids <- dup.ids[which(!clonality.status[dup.ids] %in% all.val)]

				if (0 == length(dup.sub.ids) ) {
					if ( 0 == length(dup.clo.ids) ) {
						if ( 0 == length(dup.non.ids) ) {
							# ERROR
						} else {
							max.id <- dup.non.ids[which.max(perc.overlap[dup.non.ids])]

							rm.ids <- dup.ids
							rm.ids <- rm.ids[-which(rm.ids == max.id)]
							return(rm.ids)
						}
					} else {
						max.id <- dup.clo.ids[which.max(perc.overlap[dup.clo.ids])]

						rm.ids <- dup.ids
						rm.ids <- rm.ids[-which(rm.ids == max.id)]
						return(rm.ids)
					}
				} else {
					max.id <- dup.sub.ids[which.max(perc.overlap[dup.sub.ids])]

					rm.ids <- dup.ids
					rm.ids <- rm.ids[-which(rm.ids == max.id)]
					return(rm.ids)
				}
			}

			rm.ids <- unlist(lapply(gene.dups,
				FUN=findDupIDs.SCNA,
				clonality.label=clonality.label,
				genes.label=genes.label,
				clonal.val=clonal.val,
				subclonal.val=subclonal.val
			))

			if ( 0 != length(rm.ids) ) table.sample <- table.sample[-rm.ids,]
			return(table.sample)
		},
		table=table,
		samples=samples,
		mc.cores=Ncores
	)

	message('> Assembling...')
	table2 <- rbindlist(tables)
	return(as.data.frame(table2, stringsAsFactors=F))
}

rmDups.PM = function (table,
	sample.column, clonality.label, genes.label,
	clonal.val, subclonal.val, Ncores
) {

	samples <- eval(parse(text=paste0('table$', sample.column)))
	tables <- mclapply(unique(samples),
		FUN=function(sample, table, samples) {
			table.sample <- table[which(samples == sample),]

			gene.list <- as.character(eval(parse(text=paste0('table.sample$', genes.label))))
			gene.dups <- unique(gene.list[which(duplicated(gene.list))])

			findDupIDs.SCNA = function(gene,
				data=table.sample,
				clonality.label=gb$clonality.label,
				genes.label=gb$genes.label,
				clonal.val=gb$clonal.val,
				subclonal.val=gb$subclonal.val,
				perc.overlap.label='perc.overlap'
			) {
				# 
				# Args:
				# 	gene: a gene to check for duplicates
				# 	data: the data table
				# 
				# Return:
				# 	A list of IDs to be removed
				# 	
				
				all.val <- union(clonal.val, subclonal.val)

				genes <- as.character(eval(parse(text=paste0('data$', genes.label))))
				clonality.status <- eval(parse(text=paste0('data$', clonality.label)))
				perc.overlap <- eval(parse(text=paste0('data$', perc.overlap.label)))

				dup.ids <- which(genes == gene)
				dup.sub.ids <- dup.ids[which(clonality.status[dup.ids] %in% subclonal.val)]
				dup.clo.ids <- dup.ids[which(clonality.status[dup.ids] %in% clonal.val)]
				dup.non.ids <- dup.ids[which(!clonality.status[dup.ids] %in% all.val)]

				if (0 == length(dup.sub.ids) ) {
					if ( 0 == length(dup.clo.ids) ) {
						if ( 0 == length(dup.non.ids) ) {
							# ERROR
						} else {
							rm.ids <- dup.ids
							rm.ids <- rm.ids[-which(rm.ids == dup.non.ids[1])]
							return(rm.ids)
						}
					} else {
						rm.ids <- dup.ids
						rm.ids <- rm.ids[-which(rm.ids == dup.clo.ids[1])]
						return(rm.ids)
					}
				} else {
					rm.ids <- dup.ids
					rm.ids <- rm.ids[-which(rm.ids == dup.sub.ids[1])]
					return(rm.ids)
				}
			}

			rm.ids <- unlist(lapply(gene.dups,
				FUN=findDupIDs.SCNA,
				data=table.sample,
				clonality.label=clonality.label,
				genes.label=genes.label,
				clonal.val=clonal.val,
				subclonal.val=subclonal.val
			))

			if ( 0 != length(rm.ids) ) table.sample <- table.sample[-rm.ids,]
			return(table.sample)
		},
		table=table,
		samples=samples,
		mc.cores=Ncores
	)

	message('> Assembling...')
	table2 <- rbindlist(tables)
	return(as.data.frame(table2, stringsAsFactors=F))
}

# Read Basic Data
######################################

do <- list(
	pm=TRUE,
	gain=TRUE,
	loss=TRUE,
	rr=TRUE
)

#PM
if(is.na(file.list$PM)) {
	do$pm <- FALSE
} else {
	message('> Reading pm.table')
	pm.table <- read.delim(file.list$PM, as.is=T)
	if ( test.clean == TRUE & !test.cleaned == TRUE ) {
		message('> Cleaning pm.table')
		pm.table <- rmDups.PM(pm.table,
			sample.column, clonality.label, genes.label,
			clonal.val, subclonal.val, Ncores
		)
		message('> Writing...')
		write.table(pm.table, paste0(file.list$PM, '.clean'), row.names=F, quote=F, sep='\t')
	}
}

#Gain
if(is.na(file.list$Gain)) {
	do$gain <- FALSE
} else {
	message('> Reading gain.table')
	gain.table <- read.delim(file.list$Gain, as.is=T)
	if ( test.clean == TRUE & !test.cleaned == TRUE ) {
		message('> Cleaning gain.table')
		gain.table <- rmDups.SCNA(gain.table,
			sample.column, clonality.label, genes.label,
			clonal.val, subclonal.val, Ncores
		)
		message('> Writing...')
		write.table(gain.table, paste0(file.list$Gain, '.clean'), row.names=F, quote=F, sep='\t')
	}
}

#Loss
if(is.na(file.list$Loss)) {
	do$loss <- FALSE
} else {
	message('> Reading loss.table')
	loss.table <- read.delim(file.list$Loss, as.is=T)
	if ( test.clean == TRUE & !test.cleaned == TRUE ) {
		message('> Cleaning loss.table')
		loss.table <- rmDups.SCNA(loss.table,
			sample.column, clonality.label, genes.label,
			clonal.val, subclonal.val, Ncores
		)
		message('> Writing...')
		write.table(loss.table, paste0(file.list$Loss, '.clean'), row.names=F, quote=F, sep='\t')
	}
}

#RR
if(is.na(file.list$RR)) {
	do$rr <- FALSE
} else {
	message('> Reading rr.table')
	rr.table <- read.delim(file.list$RR, as.is=T)
}

# Test single-sample graphs
######################################

if(doSingle[1] != FALSE) {

	if(1 %in% doSingle) {
		message('\n# SS CLONAL GRAPHS #')
		flist <- list.files(paste0(output.dir, '/sample-graphs/'))
		flist <- flist[grepl('^clonal_', flist)]
		for (file in flist) {
			message('> Checking graph ', file)

			g <- read.graph(file.path(output.dir, '/sample-graphs/', file), format='graphml')
			
			sample <- unlist(strsplit(file, '.', fixed=T))
			sample <- paste(sample[1:length(sample)-1], collapse='.')
			sample <- unlist(strsplit(sample, '_', fixed=T))
			sample <- paste(sample[2:length(sample)], collapse='_')
			
			table <- c()
			if (do$pm) {
				if (sample %in% pm.table$sample) {
					sub.tab <- cbind(data.frame(pm.table[which(pm.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'pm')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$gain) {
				if (sample %in% gain.table$sample) {
					sub.tab <- cbind(data.frame(gain.table[which(gain.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'gain')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$loss) {
				if (sample %in% loss.table$sample) {
					sub.tab <- cbind(data.frame(loss.table[which(loss.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'loss')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$rr) {
				if (sample %in% rr.table$sample) {
					sub.tab <- cbind(data.frame(rr.table[which(rr.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'rr')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			table <- paste(table[,1], table[,2], table[,3], sep='~')

			message(' - Checking vertices')
			checkVertex=function(i, g, data) {
				v <- V(g)[i]

				v.name <- v$name
				v.aberration <- v$abe.type
				v.hugo <- v$HUGO

				# Check name construction
				if (v.name != paste0(v.hugo, '~', v.aberration)) {
					warning('Vertex ', v.name, ' [#', i, '] has a wrong name.')
					message('Vertex ', v.name, ' [#', i, '] has a wrong name.')
				}

				if (!paste0(v.hugo, '~clonal~', v.aberration) %in% data) {
					warning('    * Vertex ', v.hugo, '~', v.aberration, ' [#', i, '] is not present.')
					message('    * Vertex ', v.hugo, '~', v.aberration, ' [#', i, '] is not present.')
				}
			}
			v.test <- mclapply(1:vcount(g),
				FUN=checkVertex,
				g=g, data=table,
				mc.preschedule=TRUE,
				mc.cores=Ncores
			)
			if(length(which(degree(g, mode='out') != vcount(g))) != 0) {
				warning('    * ', length(which(degree(g, mode='out') != vcount(g))), ' vertices have wrong degree.')
				message('    * ', length(which(degree(g, mode='out') != vcount(g))), ' vertices have wrong degree.')
			}

			message(' - Checking edges')
			if( length(which(E(g)$weight != 1)) != 0 ) {
				warning('    * ', length(which(E(g)$weight != 1)), ' edges have wrong weight.')
				message('    * ', length(which(E(g)$weight != 1)), ' edges have wrong weight.')
			}

		}
	}

	if(2 %in% doSingle) {
		message('\n# SS SUBCLONAL GRAPHS #')
		flist <- list.files(paste0(output.dir, '/sample-graphs/'))
		flist <- flist[grepl('^subclonal_', flist)]
		for (file in flist) {
			message('> Checking graph ', file)

			g <- read.graph(file.path(output.dir, '/sample-graphs/', file), format='graphml')
			
			sample <- unlist(strsplit(file, '.', fixed=T))
			sample <- paste(sample[1:length(sample)-1], collapse='.')
			sample <- unlist(strsplit(sample, '_', fixed=T))
			sample <- paste(sample[2:length(sample)], collapse='_')
			
			table <- c()
			if (do$pm) {
				if (sample %in% pm.table$sample) {
					sub.tab <- cbind(data.frame(pm.table[which(pm.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'pm')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$gain) {
				if (sample %in% gain.table$sample) {
					sub.tab <- cbind(data.frame(gain.table[which(gain.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'gain')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$loss) {
				if (sample %in% loss.table$sample) {
					sub.tab <- cbind(data.frame(loss.table[which(loss.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'loss')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$rr) {
				if (sample %in% rr.table$sample) {
					sub.tab <- cbind(data.frame(rr.table[which(rr.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'rr')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			table <- paste(table[,1], table[,2], table[,3], sep='~')

			message(' - Checking vertices')
			checkVertex=function(i, g, data) {
				v <- V(g)[i]

				v.name <- v$name
				v.aberration <- v$abe.type
				v.hugo <- v$HUGO

				# Check name construction
				if (v.name != paste0(v.hugo, '~', v.aberration)) {
					warning('Vertex ', v.name, ' [#', i, '] has a wrong name.')
					message('Vertex ', v.name, ' [#', i, '] has a wrong name.')
				}

				if (!paste0(v.hugo, '~subclonal~', v.aberration) %in% data) {
					warning('    * Vertex ', v.hugo, '~', v.aberration, ' [#', i, '] is not present.')
					message('    * Vertex ', v.hugo, '~', v.aberration, ' [#', i, '] is not present.')
				}
			}
			v.test <- mclapply(1:vcount(g),
				FUN=checkVertex,
				g=g, data=table,
				mc.preschedule=TRUE,
				mc.cores=Ncores
			)
			if(length(which(degree(g, mode='out') != vcount(g))) != 0) {
				warning('    * ', length(which(degree(g, mode='out') != vcount(g))), ' vertices have wrong degree.')
				message('    * ', length(which(degree(g, mode='out') != vcount(g))), ' vertices have wrong degree.')
			}

			message(' - Checking edges')
			if( length(which(E(g)$weight != 1)) != 0 ) {
				warning('    * ', length(which(E(g)$weight != 1)), ' edges have wrong weight.')
				message('    * ', length(which(E(g)$weight != 1)), ' edges have wrong weight.')
			}

		}
	}

	if(3 %in% doSingle) {
		message('\n# SS NONCLONAL GRAPHS #')
		flist <- list.files(paste0(output.dir, '/sample-graphs/'))
		flist <- flist[grepl('^nonclonal_', flist)]
		for (file in flist) {
			message('> Checking graph ', file)

			g <- read.graph(file.path(output.dir, '/sample-graphs/', file), format='graphml')
			
			sample <- unlist(strsplit(file, '.', fixed=T))
			sample <- paste(sample[1:length(sample)-1], collapse='.')
			sample <- unlist(strsplit(sample, '_', fixed=T))
			sample <- paste(sample[2:length(sample)], collapse='_')
			
			table <- c()
			if (do$pm) {
				if (sample %in% pm.table$sample) {
					sub.tab <- cbind(data.frame(pm.table[which(pm.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'pm')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$gain) {
				if (sample %in% gain.table$sample) {
					sub.tab <- cbind(data.frame(gain.table[which(gain.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'gain')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$loss) {
				if (sample %in% loss.table$sample) {
					sub.tab <- cbind(data.frame(loss.table[which(loss.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'loss')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$rr) {
				if (sample %in% rr.table$sample) {
					sub.tab <- cbind(data.frame(rr.table[which(rr.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'rr')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			table <- paste(table[,1], table[,2], table[,3], sep='~')

			message(' - Checking vertices')
			checkVertex=function(i, g, data) {
				v <- V(g)[i]

				v.name <- v$name
				v.aberration <- v$abe.type
				v.hugo <- v$HUGO
				v.clonality <- v$clonality

				# Check name construction
				if (v.name != paste0(v.hugo, '~', v.aberration)) {
					warning('Vertex ', v.name, ' [#', i, '] has a wrong name.')
					message('Vertex ', v.name, ' [#', i, '] has a wrong name.')
				}

				if (!paste0(v.hugo, '~', v.clonality, '~', v.aberration) %in% data) {
					warning('    * Vertex ', v.hugo, '~', v.aberration, ' [#', i, '] is not present.')
					message('    * Vertex ', v.hugo, '~', v.aberration, ' [#', i, '] is not present.')
				}
			}
			v.test <- mclapply(1:vcount(g),
				FUN=checkVertex,
				g=g, data=table,
				mc.preschedule=TRUE,
				mc.cores=Ncores
			)
			v.nonclonal <- V(g)[!V(g)$clonality %in% append(clonal.val, subclonal.val)]
			if(length(which(degree(g, v.nonclonal, mode='out') != vcount(g))) != 0) {
				warning('    * ', length(which(degree(g, mode='out') != vcount(g))), ' vertices have wrong degree.')
				message('    * ', length(which(degree(g, mode='out') != vcount(g))), ' vertices have wrong degree.')
				print(degree(g, v.nonclonal, mode='out'))
				print(g)
				print(table(V(g)$clonality))
				quit()
			}

			message(' - Checking edges')
			if( length(which(E(g)$weight != 1)) != 0 ) {
				warning('    * ', length(which(E(g)$weight != 1)), ' edges have wrong weight.')
				message('    * ', length(which(E(g)$weight != 1)), ' edges have wrong weight.')
			}

		}
	}

	if(4 %in% doSingle) {
		message('\n# SS DEPENDENCY GRAPHS #')
		flist <- list.files(paste0(output.dir, '/sample-graphs/'))
		flist <- flist[!grepl('^clonal_', flist)]
		flist <- flist[!grepl('^subclonal_', flist)]
		flist <- flist[!grepl('^nonclonal_', flist)]
		for (file in flist) {
			message('> Checking graph ', file)

			g <- read.graph(file.path(output.dir, '/sample-graphs/', file), format='graphml')
			
			sample <- unlist(strsplit(file, '.', fixed=T))
			sample <- paste(sample[1:length(sample)-1], collapse='.')
			
			table <- c()
			if (do$pm) {
				if (sample %in% pm.table$sample) {
					sub.tab <- cbind(data.frame(pm.table[which(pm.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'pm')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$gain) {
				if (sample %in% gain.table$sample) {
					sub.tab <- cbind(data.frame(gain.table[which(gain.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'gain')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$loss) {
				if (sample %in% loss.table$sample) {
					sub.tab <- cbind(data.frame(loss.table[which(loss.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'loss')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			if (do$rr) {
				if (sample %in% rr.table$sample) {
					sub.tab <- cbind(data.frame(rr.table[which(rr.table$sample == sample),c(genes.label, 'clonality.status')], stringsAsFactors=FALSE), 'rr')
					colnames(sub.tab) <- c(genes.label, 'clonality.status', 'abe.type')
					if(nrow(sub.tab) != 0) table <- rbind(table, sub.tab)
				}
			}
			table <- paste(table[,1], table[,2], table[,3], sep='~')

			message(' - Checking vertices')
			checkVertex=function(i, g, data) {
				v <- V(g)[i]

				v.name <- v$name
				v.aberration <- v$abe.type
				v.hugo <- v$HUGO
				v.clonality <- v$clonality.status

				# Check name construction
				if (v.name != paste0(v.hugo, '~', v.aberration)) {
					warning('Vertex ', v.name, ' [#', i, '] has a wrong name.')
					message('Vertex ', v.name, ' [#', i, '] has a wrong name.')
				}

				if (!paste0(v.hugo, '~', v.clonality, '~', v.aberration) %in% data) {
					warning('    * Vertex ', v.hugo, '~', v.aberration, ' [#', i, '] is not present.')
					message('    * Vertex ', v.hugo, '~', v.aberration, ' [#', i, '] is not present.')
				}
			}
			v.test <- mclapply(1:vcount(g),
				FUN=checkVertex,
				g=g, data=table,
				mc.preschedule=TRUE,
				mc.cores=Ncores
			)
			v.clonal <- V(g)[V(g)$clonality.status %in% clonal.val]
			v.subclonal <- V(g)[V(g)$clonality.status %in% subclonal.val]
			if(length(which(degree(g,v.clonal) != length(v.subclonal))) != 0) {
				warning('    * ', length(which(degree(g, mode='out') != vcount(g))), ' vertices have wrong degree.')
				message('    * ', length(which(degree(g, mode='out') != vcount(g))), ' vertices have wrong degree.')
			}
			
			message(' - Checking edges')
			if( length(which(E(g)$weight != 1)) != 0 ) {
				warning('    * ', length(which(E(g)$weight != 1)), ' edges have wrong weight.')
				message('    * ', length(which(E(g)$weight != 1)), ' edges have wrong weight.')
			}

		}
	}

}



# Test total graph
######################################

#message('\n> Reading graphs.')
#gs <- read.graph(file.path(output.dir, 'subclonal_graph.graphml'), format='graphml')
#gn <- read.graph(file.path(output.dir, 'nonclonal_graph.graphml'), format='graphml')

if(doClonal) {

	message('\n# CLONAL GRAPH #')

	message('\n> Reading graph.')
	gc <- read.graph(file.path(output.dir, 'clonal_graph.graphml'), format='graphml')

	message('#\n# Checking vertices\n#########################\n\n')
	checkClonalVertices=function(i, gc) {
		#
		
		v <- V(gc)[i]
		v.name <- v$name
		v.aberration <- v$aberration
		v.hugo <- v$HUGO
		v.clonal <- v$clonal.occ
		v.subclonal <- v$subclonal.occ

		# Check name construction
		if (v.name != paste0(v.hugo, '~', v.aberration)) {
			message('Vertex ', v.name, ' [#', i, '] has a wrong name.')
		}

		# Check if they are clonal it at least one sample
		if ( eval(parse(text=paste0('length(which(', v.aberration, '.table$clonality.status[which(', v.aberration, '.table$', genes.label, '==v.hugo)] %in% clonal.val))'))) == 0 ) {
			message('Vertex ', v.name, ' [#', i, '] is never clonal.')
		}
		
	}
	vc.test <- mclapply(1:vcount(gc),
		FUN=checkClonalVertices,
		gc=gc,
		mc.preschedule=TRUE,
		mc.cores=Ncores
	)

	message('\n#\n# Checking edges\n#########################\n\n')
	checkClonalEdges=function(i, gc, el) {
		#

		e <- E(gc)[i]

		e.source <- el[i,1]
		e.source.split <- unlist(strsplit(e.source, '~', fixed=T))
		e.source.hugo <- e.source.split[1]
		e.source.aberration <- e.source.split[2]

		e.target <- el[i,2]
		e.target.split <- unlist(strsplit(e.target, '~', fixed=T))
		e.target.hugo <- e.target.split[1]
		e.target.aberration <- e.target.split[2]

		e.weight <- e$weight

		# Check clonal co-occurrence
		source.samples <- eval(parse(text=paste0(e.source.aberration, '.table$sample[intersect(which(', e.source.aberration, '.table$', genes.label, '==e.source.hugo),which(', e.source.aberration, '.table$clonality.status %in% clonal.val))]')))
		target.clonal.samples <- eval(parse(text=paste0(e.target.aberration, '.table$sample[intersect(which(', e.target.aberration, '.table$', genes.label, '==e.target.hugo),which(', e.target.aberration, '.table$clonality.status %in% clonal.val))]')))
		test.clonal <- length(intersect(source.samples, target.clonal.samples))
		if(test.clonal != e.weight) {
			message('Edge #', i, ' from ', e.source, ' to ', e.target, ' has a wrong clonal co-occurrency: ', e.weight, ' instead of ', test.clonal)
		}

	}
	ec.list <- mclapply(1:ecount(gc),
		FUN=checkClonalEdges,
		gc=gc, el=get.edgelist(gc),
		mc.preschedule=TRUE,
		mc.cores=Ncores
	)

}


if(doTotal) {

	message('\n# TOTAL GRAPH #')

	message('> Reading graph.')
	gt <- read.graph(file.path(output.dir, 'total_graph.graphml'), format='graphml')

	message('#\n# Checking vertices\n#########################\n\n')

	#(1) Check vertices
	checkVertices=function(i, gt) {
		#
		
		v <- V(gt)[i]
		v.name <- v$name
		v.aberration <- v$aberration
		v.hugo <- v$HUGO
		v.clonal <- v$clonal.occ
		v.subclonal <- v$subclonal.occ

		# Check name construction
		if (v.name != paste0(v.hugo, '~', v.aberration)) {
			warning('Vertex ', v.name, ' [#', i, '] has a wrong name.')
			message('Vertex ', v.name, ' [#', i, '] has a wrong name.')
		}

		# Check presence in Basic data
		if (eval(parse(text=paste0('!do$', v.aberration)))) {
			warning('Vertex ', v.name, ' [#', i, '] exists for an absent aberration type.')
			message('Vertex ', v.name, ' [#', i, '] exists for an absent aberration type.')
		} else if (!v.hugo %in% eval(parse(text=paste0(v.aberration, '.table$', genes.label)))) {
			warning('Vertex ', v.name, ' [#', i, '] is not present in the original data.')
			message('Vertex ', v.name, ' [#', i, '] is not present in the original data.')
		} else {

			#Check clonal occurrency
			cs <- eval(parse(text=paste0(v.aberration, '.table$clonality.status[which(', v.aberration, '.table$', genes.label, ' == v.hugo)]')))
			if(length(which(cs %in% clonal.val)) != v.clonal) {
				warning('Vertex ', v.name, ' [#', i, '] clonal occurrence is wrong: ', v.clonal, ' instead of ', length(which(cs %in% clonal.val)))
				message('Vertex ', v.name, ' [#', i, '] clonal occurrence is wrong: ', v.clonal, ' instead of ', length(which(cs %in% clonal.val)))
				#message('Cheking in clonal_graph')
			}
			if(length(which(cs %in% subclonal.val)) != v.subclonal) {
				warning('Vertex ', v.name, ' [#', i, '] subclonal occurrence is wrong: ', v.subclonal, ' instead of ', length(which(cs %in% subclonal.val)))
				message('Vertex ', v.name, ' [#', i, '] subclonal occurrence is wrong: ', v.subclonal, ' instead of ', length(which(cs %in% subclonal.val)))
			}
		}
	}
	v.test <- mclapply(1:vcount(gt),
		FUN=checkVertices,
		gt=gt,
		mc.preschedule=TRUE,
		mc.cores=Ncores
	)

	message('#\n# Checking Edges\n#########################\n\n')

	#(2) Check edges
	checkEdges=function(i, gt, el) {
		#
		
		e <- E(gt)[i]

		e.source <- el[i,1]
		e.source.split <- unlist(strsplit(e.source, '~', fixed=T))
		e.source.hugo <- e.source.split[1]
		e.source.aberration <- e.source.split[2]

		e.target <- el[i,2]
		e.target.split <- unlist(strsplit(e.target, '~', fixed=T))
		e.target.hugo <- e.target.split[1]
		e.target.aberration <- e.target.split[2]

		e.weight <- e$weight
		e.clonal <- e$clonal.cooc
		e.subclonal <- e$subclonal.cooc

		# Check source/target
		source.samples <- eval(parse(text=paste0(e.source.aberration, '.table$sample[intersect(which(', e.source.aberration, '.table$', genes.label, '==e.source.hugo),which(', e.source.aberration, '.table$clonality.status %in% clonal.val))]')))
		target.samples <- eval(parse(text=paste0(e.target.aberration, '.table$sample[intersect(which(', e.target.aberration, '.table$', genes.label, '==e.target.hugo),which(', e.target.aberration, '.table$clonality.status %in% subclonal.val))]')))
		test.weight <- length(intersect(source.samples, target.samples))
		if(test.weight == 0) {
			warning('Edge #', i, ' from ', e.source, ' to ', e.target, ' is not present in the Basic Data.')
			message('Edge #', i, ' from ', e.source, ' to ', e.target, ' is not present in the Basic Data.')
		}
		if(test.weight != e.weight) {
			warning('Edge #', i, ' from ', e.source, ' to ', e.target, ' has a wrong weight.')
			message('Edge #', i, ' from ', e.source, ' to ', e.target, ' has a wrong weight.')
		}

		# Check clonal co-occurrence
		target.clonal.samples <- eval(parse(text=paste0(e.target.aberration, '.table$sample[intersect(which(', e.target.aberration, '.table$', genes.label, '==e.target.hugo),which(', e.target.aberration, '.table$clonality.status %in% clonal.val))]')))
		test.clonal <- length(intersect(source.samples, target.clonal.samples))
		if(test.clonal != e.clonal) {
			warning('Edge #', i, ' from ', e.source, ' to ', e.target, ' has a wrong clonal co-occurrency: ', e.clonal, ' instead of ', test.clonal)
			message('Edge #', i, ' from ', e.source, ' to ', e.target, ' has a wrong clonal co-occurrency: ', e.clonal, ' instead of ', test.clonal)
		}

		# Check subclonal co-occurrence
		source.subclonal.samples <- eval(parse(text=paste0(e.source.aberration, '.table$sample[intersect(which(', e.source.aberration, '.table$', genes.label, '==e.source.hugo),which(', e.source.aberration, '.table$clonality.status %in% subclonal.val))]')))
		test.subclonal <- length(intersect(source.subclonal.samples, target.samples))
		if(test.subclonal != e.subclonal) {
			warning('Edge #', i, ' from ', e.source, ' to ', e.target, ' has a wrong subclonal co-occurrency: ', e.subclonal, ' instead of ', test.subclonal)
			message('Edge #', i, ' from ', e.source, ' to ', e.target, ' has a wrong subclonal co-occurrency: ', e.subclonal, ' instead of ', test.subclonal)
		}
	}
	e.test <- mclapply(1:ecount(gt),
		FUN=checkEdges,
		gt=gt, el=get.edgelist(gt),
		mc.preschedule=TRUE,
		mc.cores=Ncores
	)

}

warnings()
