library('igraph')

g <- graph.empty(directed=F)
g <- g + vertices(seq(8))
g <- g + edges(c(1, 2, 1, 5, 1, 8, 2, 7, 2, 8, 3, 6, 3, 7, 4, 5, 4, 6, 7, 8))

adj <- get.adjacency(g)
lap <- - adj
for(i in seq(length(adj[,1]))) lap[i,i] <- lap[i,i] + degree(g, V(g)[i])
eig <- eigen(lap)$values
fre <- abs(sqrt(abs(eig)))

gamma <- 0.4450034

K.part <- c()
errs <- c()
subs <- c()
for(i in seq(length(fre))) {
	igr <- integrate(
		f=function(x, xdef, ga) {
			return(1 / ((x - xdef)^2 + ga^2))
		},
		lower=0,
		upper=Inf,
		xdef=fre[i],
		ga=gamma
	)
	K.part <- append(K.part, igr$value)
	errs <- append(errs, igr$abs.error)
	subs <- append(subs, igr$subdivisions)
}
K <- 1 / (gamma * sum(K.part))

y <- c()
for(x in seq(0,8,by=0.0001)) {
	y <- append(y, K * sum( gamma / ( (x - fre)^2 + gamma^2 ) ))
}
plot(seq(0,8,by=0.0001), y, type='l')
abline(v=eig, col=2)