# List of paths to the different tables
# PM = Point Mutations
# Gain = CN > 2
# Loss = CN < 2
# RR = Rearrangement
file.list <- list(
	PM="",
	Gain="",
	Loss="",
	RR=""
)

# Label of the column containing the gene IDs
genes.label <- "Gene.id"

# String vector containing the clonality values corresponding to 'clonal' and 'subclonal'
clonal.val <- c('clonal')
subclonal.val <- c('subclonal')

# Path to the output directory
output.dir <- 'output'

# NOT YET IMPLEMENTED
attr.table <- ''

# List of genes to keep (white) and to remove (black)
white.list <- list()
black.list <- list()

# Clean mode removes all the genes that are not in the white.list
clean <- FALSE

# Whether to write as output also the co-occurrency multi-sample networks in LGL format
write.cooc <- FALSE

# Verbose mode
verbose <- TRUE

# Number of cores
clusters <- 4