#!/usr/bin/env Rscript
#
# ./calcSampleCouplesDistances.R numberOfCores graphDirectory suffix annotationFile
# 
# Graph_Manager.class.R and extendigraph.R are required.

library('doParallel')
library('igraph')
source('Graph_Manager.class.R')

# Read params
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) stop('./calcSampleCouplesDistances.R numberOfCores graphDirectory suffix annotationFile')

# Prepare file list (no empty DNs)
print('> Prepare.')
flist <- list.files(args[2])
remove.id <- c()
for (i in 1:length(flist)) {
	g <- read.graph(paste0('sample-graphs/', flist[i]), format='graphml')
	if (length(E(g)) == 0) remove.id <- append(remove.id, i)
}
remove.id <- sort(remove.id, decreasing=T)
for (id in remove.id) flist <- flist[-id]
n <- length(flist)

# Prepare couples
print('> Build couples.')
couples <- cbind(rep(1:n,n), sort(rep(1:n,n)))

# Start parallel distance calculation
print('> Start parallel')
nclusters <- as.numeric(args[1])
par <- makeCluster(nclusters)
registerDoParallel(par)

l <- 1:length(couples[,1])
distances <- foreach (i=l, .combine=rbind) %dopar% {
	library('igraph')

	t <- system.time({
		gm <- GraphManager()
		xi <- 1

		f.one <- flist[couples[i,1]]
		f.two <- flist[couples[i,2]]

		g.one <- read.graph(paste0('sample-graphs/', f.one), format='graphml')
		g.two <- read.graph(paste0('sample-graphs/', f.two), format='graphml')

		V(g.one)$HUGO <- V(g.one)$name
		V(g.two)$HUGO <- V(g.two)$name
		V(g.one)$name <- paste0(V(g.one)$name, '~', V(g.one)$abe.type)
		V(g.two)$name <- paste0(V(g.two)$name, '~', V(g.two)$abe.type)

		if (is.directed(g.one)) g.one <- gm$undirected.noAttr(g.one)
		if (is.directed(g.two)) g.two <- gm$undirected.noAttr(g.two)

		ds <- gm$calcDistances(g.one, g.two, xi)
		ds <- append(c(couples[i,1], couples[i,2]), ds)
	})
	
	return(ds)
}

stopCluster(par)
print('> End parallel.')

# Output distance

distances <- data.frame(distances)
colnames(distances) <- c('g.one', 'g.two', 'dh', 'dj', 'dim', 'dhim', 'djim')
write.table(distances, paste0('dist.', args[3], '.dat'), quote=F, row.names=F)

# Print plots

print('> Prepare plots.')
hm <- matrix(as.numeric(distances$dh), n, n)
jm <- matrix(as.numeric(distances$dj), n, n)
im <- matrix(as.numeric(distances$dim), n, n)
him <- matrix(as.numeric(distances$dhim), n, n)
jim <- matrix(as.numeric(distances$djim), n, n)
hm <- hm[lower.tri(hm, diag=F)]
jm <- jm[lower.tri(jm, diag=F)]
im <- im[lower.tri(im, diag=F)]
him <- him[lower.tri(him, diag=F)]
jim <- jim[lower.tri(jim, diag=F)]

if(!file.exists('plots/')) dir.create('plots')

y.max <- max(c(hist(hm, plot=F)$density, hist(im, plot=F)$density, hist(him, plot=F)$density, hist(jim, plot=F)$density, density(hm)$y, density(im)$y, density(him)$y, density(jim)$y))
svg(paste0('plots/hamming_hist.', args[3], '.svg'))
hist(hm, xlim=c(0,1), ylim=c(0, y.max), main='Distance distribution', xlab='H', prob=T)
lines(density(hm), col=4, lty=2)
dev.off()
svg(paste0('plots/jaccard_hist.', args[3], '.svg'))
hist(jm, xlim=c(0,1), ylim=c(0, y.max), main='Distance distribution', xlab='J', prob=T)
lines(density(jm), col=4, lty=2)
dev.off()
svg(paste0('plots/ipsen_hist.', args[3], '.svg'))
hist(im, xlim=c(0,1), ylim=c(0, y.max), main='Distance distribution', xlab='IM', prob=T)
lines(density(im), col=4, lty=2)
dev.off()
svg(paste0('plots/him_hist.', args[3], '.svg'))
hist(him, xlim=c(0,1), ylim=c(0, y.max), main='Distance distribution', xlab='HIM', prob=T)
lines(density(him), col=4, lty=2)
dev.off()
svg(paste0('plots/jim_hist.', args[3], '.svg'))
hist(jim, xlim=c(0,1), ylim=c(0, y.max), main='Distance distribution', xlab='JIM', prob=T)
lines(density(jim), col=4, lty=2)
dev.off()

y.max <- max(c(hist(hm, breaks=sqrt(length(hm)), plot=F)$density, hist(im, breaks=sqrt(length(im)), plot=F)$density, hist(him, breaks=sqrt(length(him)), plot=F)$density, hist(jim, breaks=sqrt(length(jim)), plot=F)$density, density(hm)$y, density(im)$y, density(him)$y, density(jim)$y))
svg(paste0('plots/hamming_hist_sr.', args[3], '.svg'))
hist(hm, xlim=c(0,1), ylim=c(0,y.max), breaks=sqrt(length(hm)), main='Distance distribution [SQRT]', xlab='H', prob=T)
lines(density(hm), col=4, lty=2)
dev.off()
svg(paste0('plots/jaccard_hist_sr.', args[3], '.svg'))
hist(jm, xlim=c(0,1), ylim=c(0,y.max), breaks=sqrt(length(jm)), main='Distance distribution [SQRT]', xlab='J', prob=T)
lines(density(jm), col=4, lty=2)
dev.off()
svg(paste0('plots/ipsen_hist_sr.', args[3], '.svg'))
hist(im, xlim=c(0,1), ylim=c(0,y.max), breaks=sqrt(length(im)), main='Distance distribution [SQRT]', xlab='IM', prob=T)
lines(density(im), col=4, lty=2)
dev.off()
svg(paste0('plots/him_hist_sr.', args[3], '.svg'))
hist(him, xlim=c(0,1), ylim=c(0,y.max), breaks=sqrt(length(him)), main='Distance distribution [SQRT]', xlab='HIM', prob=T)
lines(density(him), col=4, lty=2)
dev.off()
svg(paste0('plots/jim_hist_sr.', args[3], '.svg'))
hist(jim, xlim=c(0,1), ylim=c(0,y.max), breaks=sqrt(length(jim)), main='Distance distribution [SQRT]', xlab='JIM', prob=T)
lines(density(jim), col=4, lty=2)
dev.off()

if(IQR(hm) != 0 & IQR(jm) != 0 & IQR(im) != 0 & IQR(him) != 0 & IQR(jim) != 0) {

	y.max <- max(c(hist(hm, breaks=1/(2*IQR(hm)/length(hm)**(1/3)), plot=F)$density, hist(im, breaks=1/(2*IQR(im)/length(im)**(1/3)), plot=F)$density, hist(him, breaks=1/(2*IQR(him)/length(him)**(1/3)), plot=F)$density, hist(jim, breaks=1/(2*IQR(jim)/length(jim)**(1/3)), plot=F)$density, density(hm)$y, density(im)$y, density(him)$y, density(jim)$y))
	svg(paste0('plots/hamming_hist_fd.', args[3], '.svg'))
	hist(hm, xlim=c(0,1), ylim=c(0,y.max), breaks=1/(2*IQR(hm)/length(hm)**(1/3)), main='Distance distribution [FD]', xlab='H', prob=T)
	lines(density(hm), col=4, lty=2)
	dev.off()
	svg(paste0('plots/jaccard_hist_fd.', args[3], '.svg'))
	hist(jm, xlim=c(0,1), ylim=c(0,y.max), breaks=1/(2*IQR(jm)/length(jm)**(1/3)), main='Distance distribution [FD]', xlab='J', prob=T)
	lines(density(jm), col=4, lty=2)
	dev.off()
	svg(paste0('plots/ipsen_hist_fd.', args[3], '.svg'))
	hist(im, xlim=c(0,1), ylim=c(0,y.max), breaks=1/(2*IQR(im)/length(im)**(1/3)), main='Distance distribution [FD]', xlab='IM', prob=T)
	lines(density(im), col=4, lty=2)
	dev.off()
	svg(paste0('plots/him_hist_fd.', args[3], '.svg'))
	hist(him, xlim=c(0,1), ylim=c(0,y.max), breaks=1/(2*IQR(him)/length(him)**(1/3)), main='Distance distribution [FD]', xlab='HIM', prob=T)
	lines(density(him), col=4, lty=2)
	dev.off()
	svg(paste0('plots/jim_hist_fd.', args[3], '.svg'))
	hist(jim, xlim=c(0,1), ylim=c(0,y.max), breaks=1/(2*IQR(jim)/length(jim)**(1/3)), main='Distance distribution [FD]', xlab='JIM', prob=T)
	lines(density(jim), col=4, lty=2)ù
	dev.off()
}

svg(paste0('plots/him_scatterplot.', args[3], '.svg'), width=10, height=10)
plot(hm[lower.tri(hm, diag=F)], im[lower.tri(im, diag=F)], xlab='H', ylab='IM', main='H/IM space', xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,0,.1))
abline(h=0.5, lty=2, col=2)
abline(v=0.5, lty=2, col=2)
dev.off()

svg(paste0('plots/jim_scatterplot.', args[3], '.svg'), width=10, height=10)
plot(jm[lower.tri(jm, diag=F)], im[lower.tri(im, diag=F)], xlab='J', ylab='IM', main='J/IM space', xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,0,.1))
abline(h=0.5, lty=2, col=2)
abline(v=0.5, lty=2, col=2)
dev.off()

print('> Prepare heatmaps.')
if(!file.exists('heatmaps/')) dir.create('heatmaps')

pv <- read.table(args[4], header=T)
pvcl <- pv[which(paste0('gra_', pv$Bam.name, '.graphml') %in% flist),]
pvcl$Bam.name <- paste0('gra_', pvcl$Bam.name, '.grapml')
pvcl$Gleason.score <- pvcl$Gleason_Major + pvcl$Gleason_Minor
pvcl$Gleason.clean <- paste0(pvcl$Gleason_Major,'+',pvcl$Gleason_Minor)
pvcl$ERG.clean <- pvcl$TMPRSS2.ERG_Fusion_Status_FISH
for (i in which(is.na(pvcl$ERG.clean))) pvcl$ERG.clean[i] <- pvcl$ETS.fusion.detected.by.sequencing[i]

library(heatmap.plus)

# Heatmap notations

color.map.col1 <- sapply(pvcl$Age, function(x) {
	if(is.na(x)) return('white')
	if (30<=x && x<40) {
		return('green')
	} else if (40<=x && x<50) {
		return('deepskyblue')
	} else if (50<=x && x<60) {
		return('darkgreen')
	} else if (60<=x && x<70) {
		return('darkorange')
	} else if (70<=x && x<80) {
		return('gold')
	}
})
color.map.col2 <- sapply(pvcl$Serum_PSA_at_diagnosis, function(x) {
	if(is.na(x)) return('white')
	if(x <= 4) {
		return('cyan3')
	} else if (x > 10) {
		return('brown1')
	} else {
		return('blueviolet')
	}
})
color.map.col3 <- sapply(1:n, function(x, major, minor) {
	if(is.na(major[x])) return('white')
	if(is.na(minor[x])) return('white')
	if(major[x]+minor[x] <= 6) {
		return('gold')
	} else if (major[x]+minor[x] == 7) {
		if(major[x] == 3 && minor[x] == 4) {
			return('firebrick1')
		} else if (major[x] == 4 & minor[x] == 3) {
			return('darkmagenta')
		} else {
			return('gray14')
		}
	} else if (major[x]+minor[x] >= 8) {
		return('forestgreen')
	}
}, major=pvcl$Gleason_Major, minor=pvcl$Gleason_Minor)
color.map.col4 <- sapply(pvcl$ERG.clean, function(x) {
	if(is.na(x)) return('white')
	if(x == 0) return('deepskyblue')
	if(x > 0) return('darkorange')
})

hm <- matrix(as.numeric(distances$dh), n, n)
jm <- matrix(as.numeric(distances$dj), n, n)
im <- matrix(as.numeric(distances$dim), n, n)
him <- matrix(as.numeric(distances$dhim), n, n)
jim <- matrix(as.numeric(distances$djim), n, n)

# HAMMING
hc <- hclust(as.dist(hm))
svg(paste0('heatmaps/hamming_dendrogram.', args[3], '.svg'))
plot(hc, hang=-1, xlab='sample', main='Hamming-based Dendrogram', cex=0.3)
dev.off()

svg(paste0('heatmaps/hamming_heat.', args[3], '.svg'))
color.map.col <- matrix(unlist(c(color.map.col1, color.map.col2, color.map.col3, color.map.col4)), nrow=length(flist), ncol=4)
colnames(color.map.col) <- c('Age', 'PSA', 'GS', 'ERG')
rownames(hm) <- colnames(hm)
heatmap.plus::heatmap.plus(hm, Rowv=as.dendrogram(hc), Colv=rev(as.dendrogram(hc)), na.rm=F, symm=T, margins=c(5,12), ColSideColors=color.map.col, cexRow=0.3, cexCol=0.3, main='Hamming')
legend(0.86,0.3,legend=c('30-40', '40-50', '50-60', '60-70', '70-80'),pch=15,col=c('green', 'deepskyblue', 'darkgreen', 'darkorange', 'gold'),cex=0.7,title="Age")
legend(0.86,0.5,legend=c("6-","3+4","4+3","8+"),pch=15,col=c("gold", "firebrick1", "darkmagenta", "forestgreen"),cex=0.7,title="GS")
legend(0.86,0.7,legend=c("< 4","4~10","> 10"),pch=15,col=c("cyan3", "blueviolet", "brown1"),cex=0.7,title="PSA")
legend(0.86,0.9,legend=c("-","+"),pch=15,col=c("deepskyblue", "darkorange"),cex=0.7,title="ERG")
dev.off()

# JACCARD
jc <- hclust(as.dist(jm))
svg(paste0('heatmaps/jaccard_dendrogram.', args[3], '.svg'))
plot(jc, hang=-1, xlab='sample', main='Jaccard-based Dendrogram', cex=0.3)
dev.off()

svg(paste0('heatmaps/jaccard_heat.', args[3], '.svg'))
color.map.col <- matrix(unlist(c(color.map.col1, color.map.col2, color.map.col3, color.map.col4)), nrow=length(flist), ncol=4)
colnames(color.map.col) <- c('Age', 'PSA', 'GS', 'ERG')
rownames(jm) <- colnames(jm)
heatmap.plus::heatmap.plus(jm, Rowv=as.dendrogram(jc), Colv=rev(as.dendrogram(jc)), na.rm=F, symm=T, margins=c(5,12), ColSideColors=color.map.col, cexRow=0.3, cexCol=0.3, main='Jaccard')
legend(0.86,0.3,legend=c('30-40', '40-50', '50-60', '60-70', '70-80'),pch=15,col=c('green', 'deepskyblue', 'darkgreen', 'darkorange', 'gold'),cex=0.7,title="Age")
legend(0.86,0.5,legend=c("6-","3+4","4+3","8+"),pch=15,col=c("gold", "firebrick1", "darkmagenta", "forestgreen"),cex=0.7,title="GS")
legend(0.86,0.7,legend=c("< 4","4~10","> 10"),pch=15,col=c("cyan3", "blueviolet", "brown1"),cex=0.7,title="PSA")
legend(0.86,0.9,legend=c("-","+"),pch=15,col=c("deepskyblue", "darkorange"),cex=0.7,title="ERG")
dev.off()

# IPSEN
ic <- hclust(as.dist(im))
svg(paste0('heatmaps/ipsen_dendrogram.', args[3], '.svg'))
plot(ic, hang=-1, xlab='sample', main='Ipsen-based Dendrogram', cex=0.3)
dev.off()

svg(paste0('heatmaps/ipsen_heat.', args[3], '.svg'))
color.map.col <- matrix(unlist(c(color.map.col1, color.map.col2, color.map.col3, color.map.col4)), nrow=length(flist), ncol=4)
colnames(color.map.col) <- c('Age', 'PSA', 'GS', 'ERG')
rownames(im) <- colnames(im)
heatmap.plus::heatmap.plus(im, Rowv=as.dendrogram(ic), Colv=rev(as.dendrogram(ic)), na.rm=F, symm=T, margins=c(5,12), ColSideColors=color.map.col, cexRow=0.3, cexCol=0.3, main='Ipsen')
legend(0.86,0.3,legend=c('30-40', '40-50', '50-60', '60-70', '70-80'),pch=15,col=c('green', 'deepskyblue', 'darkgreen', 'darkorange', 'gold'),cex=0.7,title="Age")
legend(0.86,0.5,legend=c("6-","3+4","4+3","8+"),pch=15,col=c("gold", "firebrick1", "darkmagenta", "forestgreen"),cex=0.7,title="GS")
legend(0.86,0.7,legend=c("< 4","4~10","> 10"),pch=15,col=c("cyan3", "blueviolet", "brown1"),cex=0.7,title="PSA")
legend(0.86,0.9,legend=c("-","+"),pch=15,col=c("deepskyblue", "darkorange"),cex=0.7,title="ERG")
dev.off()

# HIM
hic <- hclust(as.dist(him))
svg(paste0('heatmaps/him_dendrogram.', args[3], '.svg'))
plot(hic, hang=-1, xlab='sample', main='HIM-based Dendrogram', cex=0.3)
dev.off()

svg(paste0('heatmaps/him_heat.', args[3], '.svg'))
color.map.col <- matrix(unlist(c(color.map.col1, color.map.col2, color.map.col3, color.map.col4)), nrow=length(flist), ncol=4)
colnames(color.map.col) <- c('Age', 'PSA', 'GS', 'ERG')
rownames(him) <- colnames(him)
heatmap.plus::heatmap.plus(him, Rowv=as.dendrogram(hic), Colv=rev(as.dendrogram(hic)), na.rm=F, symm=T, margins=c(5,12), ColSideColors=color.map.col, cexRow=0.3, cexCol=0.3, main='HIM')
legend(0.86,0.3,legend=c('30-40', '40-50', '50-60', '60-70', '70-80'),pch=15,col=c('green', 'deepskyblue', 'darkgreen', 'darkorange', 'gold'),cex=0.7,title="Age")
legend(0.86,0.5,legend=c("6-","3+4","4+3","8+"),pch=15,col=c("gold", "firebrick1", "darkmagenta", "forestgreen"),cex=0.7,title="GS")
legend(0.86,0.7,legend=c("< 4","4~10","> 10"),pch=15,col=c("cyan3", "blueviolet", "brown1"),cex=0.7,title="PSA")
legend(0.86,0.9,legend=c("-","+"),pch=15,col=c("deepskyblue", "darkorange"),cex=0.7,title="ERG")
dev.off()

# JIM
jic <- hclust(as.dist(jim))
svg(paste0('heatmaps/jim_dendrogram.', args[3], '.svg'))
plot(jic, hang=-1, xlab='sample', main='JIM-based Dendrogram', cex=0.3)
dev.off()

svg(paste0('heatmaps/jim_heat.', args[3], '.svg'))
color.map.col <- matrix(unlist(c(color.map.col1, color.map.col2, color.map.col3, color.map.col4)), nrow=length(flist), ncol=4)
colnames(color.map.col) <- c('Age', 'PSA', 'GS', 'ERG')
rownames(jim) <- colnames(jim)
heatmap.plus::heatmap.plus(jim, Rowv=as.dendrogram(jic), Colv=rev(as.dendrogram(jic)), na.rm=F, symm=T, margins=c(5,12), ColSideColors=color.map.col, cexRow=0.3, cexCol=0.3, main='JIM')
legend(0.86,0.3,legend=c('30-40', '40-50', '50-60', '60-70', '70-80'),pch=15,col=c('green', 'deepskyblue', 'darkgreen', 'darkorange', 'gold'),cex=0.7,title="Age")
legend(0.86,0.5,legend=c("6-","3+4","4+3","8+"),pch=15,col=c("gold", "firebrick1", "darkmagenta", "forestgreen"),cex=0.7,title="GS")
legend(0.86,0.7,legend=c("< 4","4~10","> 10"),pch=15,col=c("cyan3", "blueviolet", "brown1"),cex=0.7,title="PSA")
legend(0.86,0.9,legend=c("-","+"),pch=15,col=c("deepskyblue", "darkorange"),cex=0.7,title="ERG")
dev.off()

print('~ END ~')
