library('igraph')

g <- graph.ring(5)
V(g)$name <- letters[1:5]
V(g)$weight <- 1:5
g <- g + vertex('a', weight=1)

h <- graph.ring(5, directed=T)
V(h)$name <- letters[1:5]
V(h)$weight <- 1:5
