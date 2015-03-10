# Path to directory containing Graph_Builder.class.R
pathToGBclass="/TEP/"

# List of paths to the different tables
# PM = Point Mutations
# Gain = CN > 2
# Loss = CN < 2
# RR = Rearrangement
# Give value NA when excluding a certain aberration type from the analysis
file.list=list(
	PM="/TEP_Runs/Example/BasicData/PM_clonalityTable.csv",
	Gain="/TEP_Runs/Example/BasicData/genesSCNA_clonalityTable.Gain.csv",
	Loss="/TEP_Runs/Example/BasicData/genesSCNA_clonalityTable.Loss.csv",
	RR=NA
)

# Label of the column containing the gene IDs
genes.label="Gene.id"

# Label of the column containing the sample IDs
sample.column='sample'

# Label of the column containing the clonality.status
clonality.label='clonality.status'

# String vector containing the clonality values corresponding to 'clonal' and 'subclonal'
clonal.val=c('clonal')
subclonal.val=c('subclonal')

# Path to the output directory
output.dir='/TEP_Runs/Example/Results/'

# Tables with attribute to be added to the nodes [NOT YET IMPLEMENTED]
attr.table=NULL

# List of genes to keep (white) and to remove (black)
white.list=list()
black.list=list()

# Clean mode removes all the genes that are not in the white.list
clean=F

# Whether to write as output also the co-occurrency multi-sample networks in LGL format
write.cooc=T

# Verbose mode
verbose=T

# Number of cores
Ncores=40


# PARAMS FOR TESTER #

# doSingle = whether to check single-sample graphs
# FALSE: don't check
# array:
# 	1: clonal
# 	2: subclonal
# 	3: nonclonal
# 	4: dependency
doSingle=F
#doSingle=1:4

# Whether to check total clonal co-occurrency graph
doClonal=T

# Whether to check total dependency graph
doTotal=T

# Whether to remove duplicated aberrations
test.clean=T

# Whether to use cleaned version of the files
test.cleaned=T
clean.list=list(PM=NA, Gain=NA, Loss=NA, RR=NA)
if ( !is.na(file.list$PM) ) clean.list$PM=paste0(file.list$PM, '.clean')
if ( !is.na(file.list$Gain) ) clean.list$Gain=paste0(file.list$Gain, '.clean')
if ( !is.na(file.list$Loss) ) clean.list$Loss=paste0(file.list$Loss, '.clean')
if ( !is.na(file.list$RR) ) clean.list$RR=paste0(file.list$RR, '.clean')
