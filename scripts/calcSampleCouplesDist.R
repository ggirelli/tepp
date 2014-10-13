#!/usr/bin/env Rscript
#
# ./calcSampleCouplesDistances.R numberOfCores graphDirectory suffix annotationFile
# 
# 01_Graph_Manager.class.R and extendigraph.R are required.

library('doParallel')
library('igraph')
source('01_Graph_Manager.class.R')

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) stop('./calcSampleCouplesDistances.R numberOfCores graphDirectory suffix annotationFile')

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
distances <- matrix(0, n, n)

print('> Build couples.')
couples <- cbind(rep(1:n,n), sort(rep(1:n,n)))

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

distances <- data.frame(distances)
colnames(distances) <- c('g.one', 'g.two', 'dh', 'dj', 'dim', 'dhim', 'djim')
write.table(distances, paste0('dist.', args[3], '.dat'), quote=F, row.names=F)

print('> Prepare plots.')
hm <- matrix(as.numeric(distances$dh), n, n)
jm <- matrix(as.numeric(distances$dj), n, n)
im <- matrix(as.numeric(distances$dim), n, n)
him <- matrix(as.numeric(distances$dhim), n, n)
jim <- matrix(as.numeric(distances$djim), n, n)

if(!file.exists('plots/')) dir.create('plots')

svg(paste0('plots/hamming_hist.', args[3], '.svg'))
hist(hm, xlim=c(0,1), main='Distance distribution', xlab='H', prob=T, ylim=c(0,80))
lines(density(hm), col=4, lty=2)
dev.off()
svg(paste0('plots/hamming_hist_sr.', args[3], '.svg'))
hist(hm, xlim=c(0,1), breaks=sqrt(length(hm)), main='Distance distribution [SQRT]', xlab='H', prob=T, ylim=c(0,80))
lines(density(hm), col=4, lty=2)
dev.off()
if(IQR(hm) != 0) {
	svg(paste0('plots/hamming_hist_fd.', args[3], '.svg'))
	hist(hm, xlim=c(0,1), breaks=1/(2*IQR(hm)/length(hm)**(1/3)), main='Distance distribution [FD]', xlab='H', prob=T, ylim=c(0,80))
	lines(density(hm), col=4, lty=2)
	dev.off()
}

svg(paste0('plots/jaccard_hist.', args[3], '.svg'))
hist(jm, xlim=c(0,1), main='Distance distribution', xlab='J', prob=T, ylim=c(0,80))
lines(density(jm), col=4, lty=2)
dev.off()
svg(paste0('plots/jaccard_hist_sr.', args[3], '.svg'))
hist(jm, xlim=c(0,1), breaks=sqrt(length(jm)), main='Distance distribution [SQRT]', xlab='J', prob=T, ylim=c(0,80))
lines(density(jm), col=4, lty=2)
dev.off()
if(IQR(jm) != 0) {
	svg(paste0('plots/jaccard_hist_fd.', args[3], '.svg'))
	hist(jm, xlim=c(0,1), breaks=1/(2*IQR(jm)/length(jm)**(1/3)), main='Distance distribution [FD]', xlab='J', prob=T, ylim=c(0,80))
	lines(density(jm), col=4, lty=2)
	dev.off()
}

svg(paste0('plots/ipsen_hist.', args[3], '.svg'))
hist(im, xlim=c(0,1), main='Distance distribution', xlab='IM', prob=T, ylim=c(0,80))
lines(density(im), col=4, lty=2)
dev.off()
svg(paste0('plots/ipsen_hist_sr.', args[3], '.svg'))
hist(im, xlim=c(0,1), breaks=sqrt(length(im)), main='Distance distribution [SQRT]', xlab='IM', prob=T, ylim=c(0,80))
lines(density(im), col=4, lty=2)
dev.off()
if(IQR(im) != 0) {
	svg(paste0('plots/ipsen_hist_fd.', args[3], '.svg'))
	hist(im, xlim=c(0,1), breaks=1/(2*IQR(im)/length(im)**(1/3)), main='Distance distribution [FD]', xlab='IM', prob=T, ylim=c(0,80))
	lines(density(im), col=4, lty=2)
	dev.off()
}

svg(paste0('plots/him_hist.', args[3], '.svg'))
hist(him, xlim=c(0,1), main='Distance distribution', xlab='HIM', prob=T, ylim=c(0,80))
lines(density(him), col=4, lty=2)
dev.off()
svg(paste0('plots/him_hist_sr.', args[3], '.svg'))
hist(him, xlim=c(0,1), breaks=sqrt(length(him)), main='Distance distribution [SQRT]', xlab='HIM', prob=T, ylim=c(0,80))
lines(density(him), col=4, lty=2)
dev.off()
if(IQR(him) != 0) {
	svg(paste0('plots/him_hist_fd.', args[3], '.svg'))
	hist(him, xlim=c(0,1), breaks=1/(2*IQR(him)/length(him)**(1/3)), main='Distance distribution [FD]', xlab='HIM', prob=T, ylim=c(0,80))
	lines(density(him), col=4, lty=2)
	dev.off()
}

svg(paste0('plots/jim_hist.', args[3], '.svg'))
hist(jim, xlim=c(0,1), main='Distance distribution', xlab='H', prob=T, ylim=c(0,80))
lines(density(jim), col=4, lty=2)
dev.off()
svg(paste0('plots/jim_hist_sr.', args[3], '.svg'))
hist(jim, xlim=c(0,1), breaks=sqrt(length(jim)), main='Distance distribution [SQRT]', xlab='H', prob=T, ylim=c(0,80))
lines(density(jim), col=4, lty=2)
dev.off()
if(IQR(jim) != 0) {
	svg(paste0('plots/jim_hist_fd.', args[3], '.svg'))
	hist(jim, xlim=c(0,1), breaks=1/(2*IQR(jim)/length(jim)**(1/3)), main='Distance distribution [FD]', xlab='H', prob=T, ylim=c(0,80))
	lines(density(jim), col=4, lty=2)
	dev.off()
}

svg(paste0('plots/him_scatterplot.', args[3], '.svg'))
plot(hm, im, xlab='H', ylab='IM', main='H/IM space', xlim=c(0,1), ylim=c(0,1), pch=20, col='gray30')
abline(h=0.5, lty=2, col=2)
abline(v=0.5, lty=2, col=2)
dev.off()

svg(paste0('plots/jim_scatterplot.', args[3], '.svg'))
plot(jm, im, xlab='J', ylab='IM', main='J/IM space', xlim=c(0,1), ylim=c(0,1), pch=20, col='gray30')
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
	if (45<=x && x<55) {
		return('deepskyblue')
	} else if (55<=x && x<65) {
		return('darkgreen')
	} else if (65<=x && x<75) {
		return('darkorange')
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

# HAMMING
hc <- hclust(as.dist(hm))
svg(paste0('heatmaps/hamming_dendrogram.', args[3], '.svg'))
plot(hc, hang=-1, xlab='sample', main='Hamming-based Dendrogram', cex=0.3)
dev.off()

svg(paste0('heatmaps/hamming_heat.', args[3], '.svg'))
color.map.col <- cbind(color.map.col1, color.map.col2, color.map.col3, color.map.col4)
color.map.row <- cbind(color.map.col4, color.map.col3, color.map.col2, color.map.col1)
colnames(color.map.col) <- c('Age', 'PSA', 'GS', 'ERG')
colnames(color.map.row) <- c('ERG', 'GS', 'PSA', 'Age')
rownames(hm) <- colnames(hm)
heatmap.plus::heatmap.plus(hm, Rowv=as.dendrogram(hc), Colv=rev(as.dendrogram(hc)), na.rm=F, symm=T, margins=c(5,12), ColSideColors=color.map.col, cexRow=0.3, cexCol=0.3, main='Hamming')
legend(0.86,0.9,legend=c("45-55","55-65","65-75"),pch=15,col=c('deepskyblue', 'darkgreen', 'darkorange'),cex=0.7,title="Age")
legend(0.86,0.7,legend=c("6-","3+4","4+3","8+"),pch=15,col=c("gold", "firebrick1", "darkmagenta", "forestgreen"),cex=0.7,title="GS")
legend(0.86,0.5,legend=c("< 4","4~10","> 10"),pch=15,col=c("cyan3", "blueviolet", "brown1"),cex=0.7,title="PSA")
legend(0.86,0.3,legend=c("-","+"),pch=15,col=c("deepskyblue", "darkorange"),cex=0.7,title="ERG")
dev.off()

# JACCARD
jc <- hclust(as.dist(jm))
svg(paste0('heatmaps/jaccard_dendrogram.', args[3], '.svg'))
plot(jc, hang=-1, xlab='sample', main='Jaccard-based Dendrogram', cex=0.3)
dev.off()

svg(paste0('heatmaps/jaccard_heat.', args[3], '.svg'))
color.map.col <- cbind(color.map.col1, color.map.col2, color.map.col3, color.map.col4)
color.map.row <- cbind(color.map.col4, color.map.col3, color.map.col2, color.map.col1)
colnames(color.map.col) <- c('Age', 'PSA', 'GS', 'ERG')
colnames(color.map.row) <- c('ERG', 'GS', 'PSA', 'Age')
rownames(jm) <- colnames(jm)
heatmap.plus::heatmap.plus(jm, Rowv=as.dendrogram(jc), Colv=rev(as.dendrogram(jc)), na.rm=F, symm=T, margins=c(5,12), ColSideColors=color.map.col, cexRow=0.3, cexCol=0.3, main='Hamming')
legend(0.86,0.9,legend=c("45-55","55-65","65-75"),pch=15,col=c('deepskyblue', 'darkgreen', 'darkorange'),cex=0.7,title="Age")
legend(0.86,0.7,legend=c("6-","3+4","4+3","8+"),pch=15,col=c("gold", "firebrick1", "darkmagenta", "forestgreen"),cex=0.7,title="GS")
legend(0.86,0.5,legend=c("< 4","4~10","> 10"),pch=15,col=c("cyan3", "blueviolet", "brown1"),cex=0.7,title="PSA")
legend(0.86,0.3,legend=c("-","+"),pch=15,col=c("deepskyblue", "darkorange"),cex=0.7,title="ERG")
dev.off()

# IPSEN
ic <- hclust(as.dist(im))
svg(paste0('heatmaps/ipsen_dendrogram.', args[3], '.svg'))
plot(ic, hang=-1, xlab='sample', main='Ipsen-based Dendrogram', cex=0.3)
dev.off()

svg(paste0('heatmaps/ipsen_heat.', args[3], '.svg'))
color.map.col <- cbind(color.map.col1, color.map.col2, color.map.col3, color.map.col4)
color.map.row <- cbind(color.map.col4, color.map.col3, color.map.col2, color.map.col1)
colnames(color.map.col) <- c('Age', 'PSA', 'GS', 'ERG')
colnames(color.map.row) <- c('ERG', 'GS', 'PSA', 'Age')
rownames(im) <- colnames(im)
heatmap.plus::heatmap.plus(im, Rowv=as.dendrogram(ic), Colv=rev(as.dendrogram(ic)), na.rm=F, symm=T, margins=c(5,12), ColSideColors=color.map.col, cexRow=0.3, cexCol=0.3, main='Ipsen')
legend(0.86,0.9,legend=c("45-55","55-65","65-75"),pch=15,col=c('deepskyblue', 'darkgreen', 'darkorange'),cex=0.7,title="Age")
legend(0.86,0.7,legend=c("6-","3+4","4+3","8+"),pch=15,col=c("gold", "firebrick1", "darkmagenta", "forestgreen"),cex=0.7,title="GS")
legend(0.86,0.5,legend=c("< 4","4~10","> 10"),pch=15,col=c("cyan3", "blueviolet", "brown1"),cex=0.7,title="PSA")
legend(0.86,0.3,legend=c("-","+"),pch=15,col=c("deepskyblue", "darkorange"),cex=0.7,title="ERG")
dev.off()

# HIM
hic <- hclust(as.dist(him))
svg(paste0('heatmaps/him_dendrogram.', args[3], '.svg'))
plot(hic, hang=-1, xlab='sample', main='HIM-based Dendrogram', cex=0.3)
dev.off()

svg(paste0('heatmaps/him_heat.', args[3], '.svg'))
color.map.col <- cbind(color.map.col1, color.map.col2, color.map.col3, color.map.col4)
color.map.row <- cbind(color.map.col4, color.map.col3, color.map.col2, color.map.col1)
colnames(color.map.col) <- c('Age', 'PSA', 'GS', 'ERG')
colnames(color.map.row) <- c('ERG', 'GS', 'PSA', 'Age')
rownames(him) <- colnames(him)
heatmap.plus::heatmap.plus(him, Rowv=as.dendrogram(hic), Colv=rev(as.dendrogram(hic)), na.rm=F, symm=T, margins=c(5,12), ColSideColors=color.map.col, cexRow=0.3, cexCol=0.3, main='HIM')
legend(0.86,0.9,legend=c("45-55","55-65","65-75"),pch=15,col=c('deepskyblue', 'darkgreen', 'darkorange'),cex=0.7,title="Age")
legend(0.86,0.7,legend=c("6-","3+4","4+3","8+"),pch=15,col=c("gold", "firebrick1", "darkmagenta", "forestgreen"),cex=0.7,title="GS")
legend(0.86,0.5,legend=c("< 4","4~10","> 10"),pch=15,col=c("cyan3", "blueviolet", "brown1"),cex=0.7,title="PSA")
legend(0.86,0.3,legend=c("-","+"),pch=15,col=c("deepskyblue", "darkorange"),cex=0.7,title="ERG")
dev.off()

# JIM
jic <- hclust(as.dist(jim))
svg(paste0('heatmaps/jim_dendrogram.', args[3], '.svg'))
plot(jic, hang=-1, xlab='sample', main='JIM-based Dendrogram', cex=0.3)
dev.off()

svg(paste0('heatmaps/jim_heat.', args[3], '.svg'))
color.map.col <- cbind(color.map.col1, color.map.col2, color.map.col3, color.map.col4)
color.map.row <- cbind(color.map.col4, color.map.col3, color.map.col2, color.map.col1)
colnames(color.map.col) <- c('Age', 'PSA', 'GS', 'ERG')
colnames(color.map.row) <- c('ERG', 'GS', 'PSA', 'Age')
rownames(jim) <- colnames(jim)
heatmap.plus::heatmap.plus(jim, Rowv=as.dendrogram(jic), Colv=rev(as.dendrogram(jic)), na.rm=F, symm=T, margins=c(5,12), ColSideColors=color.map.col, cexRow=0.3, cexCol=0.3, main='JIM')
legend(0.86,0.9,legend=c("45-55","55-65","65-75"),pch=15,col=c('deepskyblue', 'darkgreen', 'darkorange'),cex=0.7,title="Age")
legend(0.86,0.7,legend=c("6-","3+4","4+3","8+"),pch=15,col=c("gold", "firebrick1", "darkmagenta", "forestgreen"),cex=0.7,title="GS")
legend(0.86,0.5,legend=c("< 4","4~10","> 10"),pch=15,col=c("cyan3", "blueviolet", "brown1"),cex=0.7,title="PSA")
legend(0.86,0.3,legend=c("-","+"),pch=15,col=c("deepskyblue", "darkorange"),cex=0.7,title="ERG")
dev.off()

print('~ END ~')
