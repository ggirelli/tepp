#!/usr/bin/env Rscript
#
# ./plotOccurrency.R total_graph suffix output_dir

library('igraph')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3) stop('./plotOccurrency.R total_graph suffix output_dir')

# Mkdir
if(!file.exists(args[3])) dir.create(args[3])

# Read
g <- read.graph(args[1], format='graphml')

co <- as.numeric(V(g)$clonal.occ)
so <- as.numeric(V(g)$subclonal.occ)

x.max <- max(max(co), max(so))
y.max <- max(max(density(co)$y), max(density(so)$y))
breaks.clo <- seq(-0.5, max(co)+1)
breaks.sub <- seq(-0.5, max(co)+1)
h.clo <- hist(co, main='Clonal status distribution', xlab='number of samples', xlim=c(0, x.max), breaks=breaks.clo, prob=T, plot=F)
h.sub <- hist(so, main='Subclonal status distribution', xlab='number of samples', xlim=c(0, x.max), breaks=breaks.sub, prob=T, plot=F)

svg(paste0(args[3], '/clonal_occ_hist.', args[2], '.svg'))
hist(co, main='Clonal status distribution', xlab='number of samples', xlim=c(0, x.max), ylim=c(0, y.max), breaks=breaks.clo, prob=T)
lines(density(co), col=4, lty=2)
dev.off()

svg(paste0(args[3], '/subclonal_occ_hist.', args[2], '.svg'))
hist(so, main='Subclonal status distribution', xlab='number of samples', xlim=c(0, x.max), ylim=c(0, y.max), breaks=breaks.sub, prob=T)
lines(density(so), col=4, lty=2)
dev.off()

t <- cbind(co, so)
n.start <- length(co)
t <- t[-which(co == 0),]
t <- t[-which(so == 0),]
n.stop <- length(t[,1])
print(paste0('Removed ', n.start-n.stop, ' aberrations from ratio plot.'))

ratio <- t[,1]/t[,2]
lratio <- log(ratio)

x.lim <- c(min(lratio), max(lratio))
x.den <- density(lratio)$y
breaks <- seq(round(min(lratio))-0.5, round(max(lratio))+0.5, by=0.2)
h.ratio <- hist(lratio, xlim=x.lim, breaks=breaks, plot=F)
y.lim <- c(min(x.den), max(x.den, max(h.ratio$density)))

svg(paste0(args[3], '/lratio_hist.', args[2], '.svg'))
hist(lratio, xlim=x.lim, ylim=y.lim, breaks=breaks, main='Clonality status distribution', xlab='log(clonal_samples/subclonal_samples)', prob=T)
lines(density(lratio), col=4, lty=2)
dev.off()

