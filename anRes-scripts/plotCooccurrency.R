#!/usr/bin/env Rscript
#
# ./plotCooccurrency.R total_graph n.sample suffix output_dir

library('igraph')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) stop('./plotCooccurrency.R total_graph n.sample suffix output_dir')

# Mkdir
if(!file.exists(args[4])) dir.create(args[4])

# Read
g <- read.graph(args[1], format='graphml')

cs <- as.numeric(E(g)$weight)
cc <- as.numeric(E(g)$clonal.cooc)
ss <- as.numeric(E(g)$subclonal.cooc)
nn <- as.numeric(E(g)$nonclonal.cooc)

x.max <- max(c(max(cs/as.numeric(args[2])), max(cc/as.numeric(args[2])), max(ss/as.numeric(args[2])), max(nn/as.numeric(args[2]))))
y.max <- max(c(max(density(cs/as.numeric(args[2]))$y), max(density(cc/as.numeric(args[2]))$y), max(density(ss/as.numeric(args[2]))$y), max(density(nn/as.numeric(args[2]))$y)))
breaks.cs <- seq(0, 1, by=0.01)
breaks.cc <- seq(0, 1, by=0.01)
breaks.ss <- seq(0, 1, by=0.01)
breaks.nn <- seq(0, 1, by=0.01)
h.cs <- hist(cs/as.numeric(args[2]), breaks=breaks.cs, plot=F)
h.cc <- hist(cc/as.numeric(args[2]), breaks=breaks.cc, plot=F)
h.ss <- hist(ss/as.numeric(args[2]), breaks=breaks.ss, plot=F)
h.nn <- hist(nn/as.numeric(args[2]), breaks=breaks.nn, plot=F)

svg(paste0(args[4], '/edge_hist.', args[3], '.svg'))
hist(cs/as.numeric(args[2]), main='Edge distribution', xlab='%samples', ylim=c(0, y.max), xlim=c(0, 1), breaks=breaks.cs, prob=T)
lines(density(cs/as.numeric(args[2])), col=4, lty=2)
dev.off()

svg(paste0(args[4], '/clonal_cocc_hist.', args[3], '.svg'))
hist(cc/as.numeric(args[2]), main='Clonal cooccurrency distribution', xlab='%samples', ylim=c(0, y.max), xlim=c(0, 1), breaks=breaks.cc, prob=T)
lines(density(cc/as.numeric(args[2])), col=4, lty=2)
dev.off()

svg(paste0(args[4], '/subclonal_cocc_hist.', args[3], '.svg'))
hist(ss/as.numeric(args[2]), main='Subclonal cooccurrency distribution', xlab='%samples', ylim=c(0, y.max), xlim=c(0, 1), breaks=breaks.ss, prob=T)
lines(density(ss/as.numeric(args[2])), col=4, lty=2)
dev.off()

svg(paste0(args[4], '/nonclonal_cocc_hist.', args[3], '.svg'))
hist(nn/as.numeric(args[2]), main='Nonclonal cooccurrency distribution', xlab='%samples', ylim=c(0, y.max), xlim=c(0, 1), breaks=breaks.nn, prob=T)
lines(density(nn/as.numeric(args[2])), col=4, lty=2)
dev.off()

t <- cbind(cs, cs+cc+ss+nn)
n.start <- length(cs)
if(length(which(cs == 0)) != 0) t <- t[-which(cs == 0),]
if(length(which((cs+cc+ss+nn) == 0)) != 0) t <- t[-which((cs+cc+ss+nn) == 0),]
n.stop <- length(t[,1])
ratio <- t[,1]/t[,2]
lratio <- log2(ratio)

x.lim <- c(min(lratio), max(lratio))
x.den <- density(lratio)$y
breaks <- seq(round(min(lratio))-0.5, round(max(lratio))+0.5, by=0.2)
h.ratio <- hist(lratio, breaks=breaks, plot=F)
y.lim <- c(min(x.den), max(x.den, max(h.ratio$density)))

svg(paste0(args[4], '/lratio_hist.', args[3], '.svg'))
hist(lratio, xlim=c(-max(abs(lratio)), 0.5), ylim=y.lim, breaks=breaks, main='Cooccurrency distribution', xlab='log2(dependency/co_occurrency)', prob=T)
lines(density(lratio), col=4, lty=2)
dev.off()

svg(paste0(args[4], '/lratio_scatter.', args[3], '.svg'))
plot(lratio, (t[,1] + t[,2])/as.numeric(args[2]), xlim=c(-max(abs(lratio)), 0.5), ylim=c(0,1), main='Cooccurrency distribution in samples', xlab='log2(dependency/co_occurrency)', ylab='%samples', pch=20, col=rgb(0,0,0,.15))
dev.off()

t <- cbind(cs, cs+cc+ss+nn)
diff <- t[,1] - t[,2]
x.lim <- c(-max(abs(diff)), max(abs(diff)))
x.den <- density(diff)$y
breaks <- seq(round(min(diff))-0.5, round(max(diff))+0.5, by=0.2)
h.diff <- hist(diff, breaks=breaks, plot=F)
y.lim <- c(min(x.den), max(x.den, max(h.ratio$density)))

svg(paste0(args[4], '/diff_hist.', args[3], '.svg'))
hist(diff, xlim=c(-max(abs(lratio)), 0.5), ylim=y.lim, breaks=breaks, main='Cooccurrency distribution', xlab='dependency - co_occurrency', prob=T)
lines(density(diff), col=4, lty=2)
dev.off()

svg(paste0(args[4], '/diff_scatter.', args[3], '.svg'))
plot(diff, (t[,1] + t[,2])/as.numeric(args[2]), xlim=c(-max(abs(lratio)), 0.5), ylim=c(0,1), main='Cooccurrency distribution in samples', xlab='dependency - co_occurrency', ylab='%samples', pch=20, col=rgb(0,0,0,.15))
dev.off()