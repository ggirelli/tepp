library('igraph')

#---------#
# prepare #
#---------#

# Make graph
g <- graph.empty(directed=F)
g <- g + vertices(seq(10))
g <- g + edges(c(1,2,1,3,1,4,1,5,1,6,1,7,2,3,3,4,4,5,5,6,6,7,4,8,5,8,8,9,9,10))
#tkplot(g)

# Get adjacency matrix and build lap matrix
adj <- as.matrix(get.adjacency(g))
lap <- - adj
for(i in seq(length(adj[,1]))) lap[i,i] <- lap[i,i] + degree(g, V(g)[i])
# Get 'defined' frequencies
eigva <- round(sort(eigen(lap)$values),5)
freq <- sqrt(abs(eigva))

# HWHM parameter
gamma <- 0.08




#---------#
# my code #
#---------#

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
y <- c()
for(x in seq(0,4,by=0.001)) y <- append(y, rhoO(x, freq, gamma, sdK(freq, gamma)))
plot(seq(0,4,by=0.001), y, type='l'); abline(v=eigva, col=2)




#----------#
# nettools #
#----------#

## Useful function for computing Ipsen distance
##--------------------------------------------------
spec <- function(mm){
  sort(eigen(mm)$values)
}

## D2 - ipsen02evolutionary
lorentz <- function(omega,mygamma,given_omega){
  l <-0
  for(i in 2:length(given_omega)){
    l = l + mygamma/( (omega-given_omega[i])**2+mygamma**2)                          }
  return(l)
}

K <- function(mygamma,given_omega){
  return(1/integrate(lorentz,lower=0,upper=Inf,mygamma=mygamma,given_omega=given_omega, stop.on.error = FALSE)$value)
}

rho <- function(omega, mygamma, ll){
  ll[[2]]*lorentz(omega,mygamma,ll[[1]])
}

## Code
##-----------------------------------------------
myomega <- sqrt(abs(round(spec(lap),5)))
mygamma <- gamma
myk <- K(mygamma,myomega)
y <- c()
for(x in seq(0,4,by=0.001)) y <- append(y, rho(x, mygamma, list(myomega, myk)))
plot(seq(0,4,by=0.001), y, type='l'); abline(v=eigva, col=2)