source('./05_Graph_Manager.class.R')

system.time({
  gm <- GraphManager(clusters=4, verbose=TRUE)
  gm$builder <- gm$builder$readData(file.path='../Tables/20140317.TCGA.246FreezeSamples.PM.txt', sample.column='sample')
  gm$builder$buildGraph(gm$builder,table.out=TRUE)
  gm$mergeGraphs.noAttr(paste('./sample-graphs/gra_',unique(gm$builder$data[,'sample']),'.graphml',sep=''))
})
