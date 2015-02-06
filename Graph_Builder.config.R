# Path to directory containing Graph_Builder.class.R
pathToGBclass <- "."

# List of paths to the different tables
# PM = Point Mutations
# Gain = CN > 2
# Loss = CN < 2
# RR = Rearrangement
# Give value NA when excluding a certain aberration type from the analysis
file.list <- list(
	PM="Tables/20140317.TCGA.246FreezeSamples.PM.txt",
	Gain="Tables/20140317.TCGA.246FreezeSamples.Loss.txt",
	Loss="Tables/20140317.TCGA.246FreezeSamples.Gain.txt",
	RR=NA
)

# Label of the column containing the gene IDs
genes.label <- "Gene.id"

# String vector containing the clonality values corresponding to 'clonal' and 'subclonal'
clonal.val <- c('clonal')
subclonal.val <- c('subclonal')

# Path to the output directory
output.dir <- 'output'

# Tables with attribute to be added to the nodes [NOT YET IMPLEMENTED]
attr.table <- NULL

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
Ncores <- 4