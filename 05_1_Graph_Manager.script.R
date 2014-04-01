gm <- GraphManager(clusters=4, verbose=TRUE)
gm$builder <- gm$builder$readData(file.path='C:\\Users\\Gire\\Desktop\\TumorEvolution_project\\Tables\\20140317.TCGA.246FreezeSamples.PM.txt', sample.column='sample')
gm$builder$buildGraph(gm$builder,print.table=TRUE)
gm$mergeGraphs(paste('./sample-graphs/gra_',unique(gm$builder$data[,'sample']),'.graphml',sep=''))