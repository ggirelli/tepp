#---------------#
# mix functions #
#---------------#

#mergeNumericLists <- function(l1, l2) {
#  for(i in 1:length(l1)) {
#    if (names(l1)[i] != "" && !is.null(l2[[names(l1)[i]]])) {
#      l1[names(l1)[i]] <- l1[[names(l1)[i]]] + l2[[names(l1)[i]]]
#    }
#  }
#  for(i in 1:length(l2)) {
#    if (names(l2)[i] != "" && is.null(l1[[names(l2)[i]]])) {
#      l1[names(l2)[i]] <- l2[[names(l2)[i]]]
#    }
#  }
#  return(l1)
#}
#
#print(mergeNumericLists(a,b))

#-------------#
# MSC PROBLEM #
#-------------#
 
# library('igraph')
# 
# g <- graph.empty()
# g <- g + vertex(letters[1:6])
# g <- g + edges(c('a','e','a','c','b','d','b','e','b','f','c','d','c','e','c','f','d','e'))
# 
# h <- graph.empty()
# h <- h + vertex(letters[1:5])
# h <- h + vertex('f')
# h <- h + edges(c('a','d','b','c','b','d','b','e','c','d','c','e','d','e'))
# h <- h + edges(c('b','f','c','f'))
# 
# gh <- graph.intersection(g,h)
# i <- graph.empty()
# for(v in V(gh)$name) if(v %in% c(get.edgelist(gh)[,1],get.edgelist(gh)[,2])) i <- i + vertex(v)
# i <- i + edges(t(get.edgelist(gh)))
# 
# d <- 1 - (length(V(i)) / max(length(V(g)),length(V(h))))
