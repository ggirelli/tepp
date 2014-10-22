#!/usr/bin/env Rscript
#
# ./plotOccurrency.R total_graph n.samples suffix output_dir

library('igraph')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) stop('./plotOccurrency.R total_graph n.samples suffix output_dir')

# Mkdir
if(!file.exists(args[4])) dir.create(args[4])

# Read
g <- read.graph(args[1], format='graphml')

co <- as.numeric(V(g)$clonal.occ)
so <- as.numeric(V(g)$subclonal.occ)

x.max <- max(max(co/as.numeric(args[2])), max(so/as.numeric(args[2])))
breaks.clo <- seq(0, 1, by=0.01)
breaks.sub <- seq(0, 1, by=0.01)
h.clo <- hist(co/as.numeric(args[2]), breaks=breaks.clo, plot=F)
h.sub <- hist(so/as.numeric(args[2]), breaks=breaks.sub, plot=F)
y.max <- max(max(density(co/as.numeric(args[2]))$y), max(density(so/as.numeric(args[2]))$y), max(h.clo$density), max(h.sub$density))

svg(paste0(args[4], '/clonal_occ_hist.', args[3], '.svg'))
hist(co/as.numeric(args[2]), main='Clonal status distribution', xlab='%samples', xlim=c(0, 1), ylim=c(0, y.max), breaks=breaks.clo, prob=T)
lines(density(co/as.numeric(args[2])), col=4, lty=2)
dev.off()

svg(paste0(args[4], '/subclonal_occ_hist.', args[3], '.svg'))
hist(so/as.numeric(args[2]), main='Subclonal status distribution', xlab='%samples', xlim=c(0, 1), ylim=c(0, y.max), breaks=breaks.sub, prob=T)
lines(density(so/as.numeric(args[2])), col=4, lty=2)
dev.off()

t <- cbind(co, so)
n.start <- length(co)
if(length(which(co == 0)) != 0) t <- t[-which(co == 0),]
if(length(which(so == 0)) != 0) t <- t[-which(so == 0),]
n.stop <- length(t[,1])
print(paste0('Removed ', n.start-n.stop, ' aberrations from ratio plot.'))

ratio <- t[,1]/t[,2]
lratio <- log2(ratio)

x.lim <- c(min(lratio), max(lratio))
x.den <- density(lratio)$y
breaks <- seq(round(min(lratio))-0.5, round(max(lratio))+0.5, by=0.2)
h.ratio <- hist(lratio, breaks=breaks, plot=F)
y.lim <- c(min(x.den), max(x.den, max(h.ratio$density)))

svg(paste0(args[4], '/lratio_hist.', args[3], '.svg'))
hist(lratio, xlim=c(-max(abs(lratio)), max(abs(lratio))), ylim=y.lim, breaks=breaks, main='Clonality status distribution', xlab='log2(clonal_samples/subclonal_samples)', prob=T)
lines(density(lratio), col=4, lty=2)
dev.off()

svg(paste0(args[4], '/lratio_scatter.', args[3], '.svg'))
plot(lratio, (t[,1] + t[,2])/as.numeric(args[2]), xlim=c(-max(abs(lratio)), max(abs(lratio))), main='Clonality status distribution in samples', ylim=c(0,1), ylab='%samples', xlab='log2(clonal_samples/subclonal_samples)', pch=20, col=rgb(0,0,0,.15))
dev.off()

t <- cbind(co, so)
diff <- t[,1] - t[,2]
x.lim <- c(-max(abs(diff)), max(abs(diff)))
x.den <- density(diff)$y
breaks <- seq(round(min(diff))-0.5, round(max(diff))+0.5, by=0.2)
h.diff <- hist(diff, breaks=breaks, plot=F)
y.lim <- c(min(x.den), max(x.den, max(h.diff$density)))

svg(paste0(args[4], '/diff_hist.', args[3], '.svg'))
hist(diff, xlim=c(-max(abs(lratio)), max(abs(lratio))), ylim=y.lim, breaks=breaks, main='Clonality status distribution', xlab='clonal_samples - subclonal_samples', prob=T)
lines(density(diff), col=4, lty=2)
dev.off()

svg(paste0(args[4], '/diff_scatter.', args[3], '.svg'))
plot(diff, (t[,1] + t[,2])/as.numeric(args[2]), xlim=c(-max(abs(lratio)), max(abs(lratio))), ylim=c(0,1), main='Clonality status distribution in samples', xlab='clonal_samples - subclonal_samples', ylab='%samples', pch=20, col=rgb(0,0,0,.15))
dev.off()