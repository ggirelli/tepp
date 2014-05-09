library('igraph')

# Make graph
g <- graph.empty(directed=F)
g <- g + vertices(seq(10))
g <- g + edges(c(1,2,1,3,1,4,1,5,1,6,1,7,2,3,3,4,4,5,5,6,6,7,3,8,7,8,8,9,9,10))
#tkplot(g)

# Get adjacency matrix and build lap matrix
adj <- as.matrix(get.adjacency(g))
lap <- - adj
for(i in seq(length(adj[,1]))) lap[i,i] <- lap[i,i] + degree(g, V(g)[i])
# Get 'defined' frequencies
eigva <- round(sort(eigen(lap)$values),5)
freq <- sqrt(abs(eigva))

# Make graph
g2 <- graph.empty(directed=F)
g2 <- g2 + vertices(seq(10))
g2 <- g2 + edges(c(1,2,1,3,1,4,1,5,1,6,1,7,2,3,3,4,4,5,5,6,6,7,4,8,5,8,8,9,9,10))
#tkplot(g2)

# Get adjacency matrix and build lap matrix
adj2 <- as.matrix(get.adjacency(g2))
lap2 <- - adj2
for(i in seq(length(adj2[,1]))) lap2[i,i] <- lap2[i,i] + degree(g2, V(g2)[i])
# Get 'defined' frequencies
eigva2 <- round(sort(eigen(lap2)$values),5)
freq2 <- sqrt(abs(eigva2))

# HWHM parameter
gamma <- 0.4450034

# SpectralDensity function
specDens = function(omega, omegadef, gamma) {
  k <- 0
  for(i in 2:length(omegadef)) k = k + (gamma / ((omega - omegadef[i])^2 + gamma^2))
  return(k)
}
# Normalization constant
sdK = function(omegadef, gamma) {
  return(integrate(specDens, lower=0, upper=Inf, omegadef=omegadef, gamma=gamma, stop.on.error = FALSE)$value)
}
# Normalized spectral density (rho[omega])
rhoO = function(omega, omegadef, gamma, k) {
  return(specDens(omega, omegadef, gamma) / k)
}

# Spectral plot
k <- sdK(freq, gamma)
y <- c()
for(x in seq(0,4,by=0.001)) y <- append(y, rhoO(x, freq, gamma, k))
plot(seq(0,4,by=0.001), y, type='l'); abline(v=eigva, col=2)
k2 <- sdK(freq2, gamma)
y2 <- c()
for(x in seq(0,4,by=0.001)) y2 <- append(y2, rhoO(x, freq2, gamma, k2))
lines(seq(0,4,by=0.001), y2, type='l', lty=2); abline(v=eigva2, col=2, lty=2)

#IM distance
sdDiffSq = function(omega, omegadef1, k1, omegadef2, k2, gamma) {
  return((rhoO(omega, omegadef1, gamma, k1) - rhoO(omega, omegadef2, gamma, k2))**2)
}
dIM <- sqrt(integrate(sdDiffSq, lower=0, upper=Inf, omegadef1=freq, k1=k, omegadef2=freq2, k2=k2, gamma=gamma)$value)

calcHammingDist = function(g.one, g.two) {
  # Calculates the Hamming (edit) distance between two graphs
  # The two graphs MUST HAVE the same nodes
  #
  # Args:
  # g.one: first graph
  # g.two: second graph
  #
  # Returns:
  # The Hamming distance H(g.one,g.two)

  # Get edges
  adj.one <- get.adjacency(g.one)
  adj.two <- get.adjacency(g.two)

  d <- sum(abs(adj.one - adj.two))
  # Normalize distance
  max.v <- max(length(adj.one[,1]), length(adj.two[,1]))
  d <- d / (max.v * (max.v - 1))

  # Return distance
  return(d)
}
dH <- calcHammingDist(g, g2)

eta <- 1

dHIM <- (1/sqrt(1+eta)) * sqrt(dH**2 + eta * dIM**2)