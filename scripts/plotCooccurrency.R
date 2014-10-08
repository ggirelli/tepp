#!/usr/bin/env Rscript
#
# ./plotCooccurrency.R total_graph suffix output_dir

library('igraph')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3) stop('./plotCooccurrency.R total_graph suffix output_dir')

# Mkdir
if(!file.exists(args[3])) dir.create(args[3])

# Read
g <- read.graph(args[1], format='graphml')

cs <- as.numeric(E(g)$weight)
cc <- as.numeric(E(g)$clonal.cooc)
ss <- as.numeric(E(g)$subclonal.cooc)
nn <- as.numeric(E(g)$nonclonal.cooc)

x.max <- max(c(max(cs), max(cc), max(ss), max(nn)))
y.max <- max(c(max(density(cs)$y), max(density(cc)$y), max(density(ss)$y), max(density(nn)$y)))
breaks.cs <- seq(-0.5, max(cs)+1)
breaks.cc <- seq(-0.5, max(cc)+1)
breaks.ss <- seq(-0.5, max(ss)+1)
breaks.nn <- seq(-0.5, max(nn)+1)
h.cs <- hist(cs, main='Edge distribution', xlab='number of samples', xlim=c(0, x.max), breaks=breaks.cs, prob=T, plot=F)
h.cc <- hist(cc, main='Clonal cooccurrency distribution', xlab='number of samples', xlim=c(0, x.max), breaks=breaks.cc, prob=T, plot=F)
h.ss <- hist(ss, main='Subclonal cooccurrency distribution', xlab='number of samples', xlim=c(0, x.max), breaks=breaks.ss, prob=T, plot=F)
h.nn <- hist(nn, main='Nonclonal cooccurrency distribution', xlab='number of samples', xlim=c(0, x.max), breaks=breaks.nn, prob=T, plot=F)

svg(paste0(args[3], '/edge_hist.', args[2], '.svg'))
hist(cs, main='Edge distribution', xlab='number of samples', ylim=c(0, y.max), xlim=c(0, x.max), breaks=breaks.cs, prob=T)
lines(density(cs), col=4, lty=2)
dev.off()

svg(paste0(args[3], '/clonal_cocc_hist.', args[2], '.svg'))
hist(cc, main='Clonal cooccurrency distribution', xlab='number of samples', ylim=c(0, y.max), xlim=c(0, x.max), breaks=breaks.cc, prob=T)
lines(density(cc), col=4, lty=2)
dev.off()

svg(paste0(args[3], '/subclonal_cocc_hist.', args[2], '.svg'))
hist(ss, main='Subclonal cooccurrency distribution', xlab='number of samples', ylim=c(0, y.max), xlim=c(0, x.max), breaks=breaks.ss, prob=T)
lines(density(ss), col=4, lty=2)
dev.off()

svg(paste0(args[3], '/nonclonal_cocc_hist.', args[2], '.svg'))
hist(nn, main='Nonclonal cooccurrency distribution', xlab='number of samples', ylim=c(0, y.max), xlim=c(0, x.max), breaks=breaks.nn, prob=T)
lines(density(nn), col=4, lty=2)
dev.off()

t <- cbind(cs, cs+cc+ss+nn)
n.start <- length(cs)
if(length(which(cs == 0)) != 0) t <- t[-which(cs == 0),]
if(length(which((cs+cc+ss+nn) == 0)) != 0) t <- t[-which((cs+cc+ss+nn) == 0),]
n.stop <- length(t[,1])
ratio <- t[,1]/t[,2]
lratio <- log(ratio)

x.lim <- c(min(lratio), max(lratio))
x.den <- density(lratio)$y
breaks <- seq(round(min(lratio))-0.5, round(max(lratio))+0.5, by=0.2)
h.ratio <- hist(lratio, xlim=x.lim, breaks=breaks, plot=F)
y.lim <- c(min(x.den), max(x.den, max(h.ratio$density)))

svg(paste0(args[3], '/lratio_hist.', args[2], '.svg'))
hist(lratio, xlim=x.lim, ylim=y.lim, breaks=breaks, main='Cooccurrency distribution', xlab='log(dependency/co_occurrency)', prob=T)
lines(density(lratio), col=4, lty=2)
dev.off()

