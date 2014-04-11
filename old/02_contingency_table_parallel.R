# Gabriele Girelli <gabriele@filopoe.it>
# Elaborates clonal/subclonal gene groups into graphs

# Ask for where to set current working directory
print('Specify current directory absolute path:')
dir.name <- readline('')
print('Specify file name (absolute path 
  or path relative to working directory):')
file.name <- readline('')

# Make directory if it does not exist
if (file.exists(dir.name)) {
  setwd(file.path(dir.name))
} else {
  dir.create(file.path(dir.name))
  setwd(file.path(dir.name))
}

# Read data table
data <- read.table(file.name, header=TRUE, sep="\t", row.names=NULL)
# Prepare sample list without duplicates
sample.list <- unique(data$sample)
cl <- makeCluster(3)
registerDoParallel(cl)
# Split original data table for each sample in data.frames
# and save them in temporary directory
foreach(i = 1:length(sample.list)) %dopar% {
  # On which sample are we working?
  sample.id <- sample.list[i]
  # Get row_ids from original data table for the working_sample
  row.ids <- which(data$sample ==sample.id)
  # Retrieve genes and clonality status for working_sample
  sample.genes <- data$Gene.id[row.ids]
  sample.status <- data$clonality.status[row.ids]
  # Make data.frame and append it to data.frame list
  data.frame <- data.frame(genes = sample.genes, clonality = sample.status)
  write.table(data.frame, file = file.path('.', sample.list[i]))
}
stopCluster(cl)

cl <- makeCluster(3)
registerDoParallel(cl)
# Prepare gene list without duplicates
genelist <- unique(Gene.id)
# Empty data
# Start filling giant matrix and make small matrices
foreach(i = 1:length(sample.list)) %dopar% {
  library('igraph')
  # On which sample are we working?
  sample.id <- sample.list[i]
  # Read sample data
  data <- read.table(file.path('.', sample.id), header = TRUE)
  # Make small contingency matrix
  adjacency.matrix <- matrix(0,
    nrow = length(data$genes),
    ncol = length(data$genes)
  )
  rownames(adjacency.matrix)<- data$genes
  colnames(adjacency.matrix)<- data$genes
  # Distinguish clonal and subclonal
  genes.clonal <- which(data$clonality == 'clonal')
  genes.subclonal <- which(data$clonality == 'subclonal')
  # Increment cells in matrices
  for (i in seq(length(genes.clonal))){
    for (j in seq(length(genes.subclonal))){
      adjacency.matrix[genes.clonal[i], genes.subclonal[j]] <- 
        adjacency.matrix[genes.clonal[i], genes.subclonal[j]] + 1
    }
  }
  # Write small matrix
  # write.table(adjacency.matrix,
  #  file = file.path('.', paste('tab_', sample.id, '.dat', sep = "")))
  write.graph(
    graph.adjacency(adjacency.matrix, mode="directed", weighted=TRUE),
    file = file.path('.', paste('gra_', sample.id, '.graphml', sep = "")),
    format = 'graphml')
}
stopCluster(cl)