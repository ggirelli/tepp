# Ask for where to set current working directory
print( 'Specify current directory absolute path:' )
dirname <- readline( '' )
# Make directory if it does not exist
if ( file.exists( dirname ) ) {
  setwd( file.path( dirname  ) )
} else {
  dir.create( file.path( dirname ) )
  setwd( file.path( dirname  ) )
}

# Read data table
print( 'Specify file name (absolute path or path relative to working directory):' )
filename <- readline( '' )
data <- read.table( filename, header=TRUE, sep="\t", row.names=NULL )
attach(data)
# Prepare sample list without duplicates
samplelist <- unique( sample )

# Split original data table for each sample in data.frames and save them in temporary directory
for ( i in seq( length( samplelist ) ) ) {
  # On which sample are we working?
  idsample <- samplelist[i]
  # Get row_ids from original data table for the working_sample
  idrows <- which( sample == idsample )
  # Retrieve genes and clonality status for working_sample
  samplegenes <- Gene.id[idrows]
  samplestatus <- clonality.status[idrows]
  # Make data.frame and append it to data.frame list
  dataframe <- data.frame( genes = samplegenes, clonality = samplestatus )
  write.table( dataframe, file = file.path( '.', samplelist[i] ) )
}

# Prepare gene list without duplicates
genelist <- unique( Gene.id )
# Empty data
detach( data )
# Make giant matrix
giantmatrix <- matrix( 0, nrow = length( genelist ), ncol = length( genelist ) )
rownames( giantmatrix ) <- genelist
colnames( giantmatrix ) <- genelist
# Start filling giant matrix and make small matrices
for ( i in seq( length( samplelist ) ) ) {
  # On which sample are we working?
  idsample <- samplelist[i]
  # Read sample data
  data <- read.table( file.path( '.', idsample ), header = TRUE )
  attach( data )  
  # Make small contingency matrix
  smallmatrix <- matrix( 0, nrow = length( genes ), ncol = length( genes ) )
  rownames( smallmatrix ) <- genes
  colnames( smallmatrix ) <- genes
  # Distinguish clonal and subclonal
  clonal <- which( clonality == 'clonal' )
  subclonal <- which( clonality == 'subclonal' )
  # Increment cells in matrices
  for( i in seq( length( clonal ) ) ) {
    for( j in seq( length( subclonal ) ) ) {
      giantmatrix[clonal[i],subclonal[j]] <- giantmatrix[clonal[i],subclonal[j]] + 1
      smallmatrix[clonal[i],subclonal[j]] <- smallmatrix[clonal[i],subclonal[j]] + 1
    }
  }
  # Write small matrix
  # write.table( smallmatrix, file = file.path( '.', paste( 'tab_', idsample, '.dat', sep = "" ) ) )
  write.graph( graph.adjacency( smallmatrix, mode="directed", weighted=TRUE ), file = file.path( '.', paste( 'gra_', idsample, '.graphml', sep = "" ) ), format = 'graphml' )
  detach( data )
}

# Write giant matrix
# write.table( giantmatrix, file = file.path( '.', 'total_tab.dat' ) )
write.graph( graph.adjacency( giantmatrix, mode="directed", weighted=TRUE ), file = file.path( '.', 'total_graph.graphml' ), format = 'graphml' )
