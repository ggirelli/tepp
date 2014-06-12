#!/usr/bin/env Rscript

cat('Reading Gleason scores.\n')
gt <- read.table('TCGAsamples_Gleason.csv', header=T)
temp <- sapply(as.character(gt$Final.Gleason), FUN=function(x) { return(eval(parse(text=x))); })

cat('Splitting Gleason scores.\n')
l <- list()
l$six <- gt$sample.id[which(temp == 6)]
l$seven <- gt$sample.id[which(temp == 7)]
l$eight <- gt$sample.id[which(temp == 8)]
l$nine <- gt$sample.id[which(temp == 9)]
l$ten <- gt$sample.id[which(temp == 10)]

cat('Reading data.\n')
st <- list()
st$PM <- read.table('20140317.TCGA.246FreezeSamples.PM.txt', header=T, row.names=NULL)
st$Loss <- read.table('20140317.TCGA.246FreezeSamples.Loss.correct.txt', header=T, row.names=NULL)
st$Gain <- read.table('20140317.TCGA.246FreezeSamples.Gain.correct.txt', header=T, row.names=NULL)

cat('Formatting data sample.id.\n')
PM.sample <- sapply(strsplit(as.character(st$PM$sample), '.', fixed=T), function(x) { return(substr(x[2], 1, 16)) })
Loss.sample <- sapply(strsplit(as.character(st$Loss$sample), '.', fixed=T), function(x) { return(substr(x[2], 1, 16)) })
Gain.sample <- sapply(strsplit(as.character(st$Gain$sample), '.', fixed=T), function(x) { return(substr(x[2], 1, 16)) })

cat('Creating output directory.\n')
dirname <- 'GleasonSplitted'
if(!file.exists(dirname))
	dir.create(dirname)
setwd(dirname)

cat('Output for GS6.\n')
write.table(st$PM[which(PM.sample %in% l$six),], '20140317.TCGA.246FreezeSamplesPM.GS6.txt')
write.table(st$Loss[which(Loss.sample %in% l$six),], '20140317.TCGA.246FreezeSamplesLoss.GS6.txt')
write.table(st$Gain[which(Gain.sample %in% l$six),], '20140317.TCGA.246FreezeSamplesGain.GS6.txt')

cat('Output for GS7.\n')
write.table(st$PM[which(PM.sample %in% l$seven),], '20140317.TCGA.246FreezeSamplesPM.GS7.txt')
write.table(st$Loss[which(Loss.sample %in% l$seven),], '20140317.TCGA.246FreezeSamplesLoss.GS7.txt')
write.table(st$Gain[which(Gain.sample %in% l$seven),], '20140317.TCGA.246FreezeSamplesGain.GS7.txt')

cat('Output for GS8+.\n')
write.table(rbind(st$PM[which(PM.sample %in% l$eight),], st$PM[which(PM.sample %in% l$nine),], st$PM[which(PM.sample %in% l$ten),]), '20140317.TCGA.246FreezeSamplesPM.GS8+.txt')
write.table(rbind(st$Loss[which(Loss.sample %in% l$eight),], st$Loss[which(Loss.sample %in% l$nine),], st$Loss[which(Loss.sample %in% l$ten),]), '20140317.TCGA.246FreezeSamplesLoss.GS8+.txt')
write.table(rbind(st$Gain[which(Gain.sample %in% l$eight),], st$Gain[which(Gain.sample %in% l$nine),], st$Gain[which(Gain.sample %in% l$ten),]), '20140317.TCGA.246FreezeSamplesGain.GS8+.txt')