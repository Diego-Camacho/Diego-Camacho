load("sconeFinished.RDS")  

rleGC_med <- function(counts, gcContent, binSize=4e3){
  logCounts <- log1p(counts)
  meds <- rowMedians(as.matrix(logCounts))
  rle <- sweep(as.matrix(logCounts),1,meds,FUN="-") # log scale
  gcBins <- Hmisc::cut2(gcContent, g=round(nrow(logCounts)/binSize))
  rle_med <- vector(length=ncol(logCounts))
  for(kk in 1:nlevels(gcBins)){
    rleBin <- rle[gcBins == levels(gcBins)[kk],]
    binMeds <- matrixStats::colMedians(rleBin)
    rle_med[kk] <- mean(binMeds^2)
  }
  rleMedsGC <- var(rle_med)
}

rleGC_iqr <- function(counts, gcContent, binSize=4e3){
  logCounts <- log1p(counts)
  meds <- rowMedians(as.matrix(logCounts))
  rle <- sweep(as.matrix(logCounts),1,meds,FUN="-") # log scale
  gcBins <- Hmisc::cut2(gcContent, g=round(nrow(logCounts)/binSize))
  rle_iqr <- vector(length=ncol(logCounts))
  for(kk in 1:nlevels(gcBins)){
    rleBin <- rle[gcBins == levels(gcBins)[kk],]
    binIQRs <- matrixStats::colIQRs(rleBin)
    rle_iqr[kk] <- var(binIQRs)
  }
  rleMedsGC <- var(rle_iqr)
}

rleGCh5 <- function(h5File, gcContent, type="med", ...){
  h5List <- h5ls(h5File)
  rleGCVals <- sapply(2:(nrow(h5List)-1), function(normId){
    normCounts <- h5read(h5File, h5List$name[normId])
    if(type == "med"){
      rleValue <- rleGC_med(counts=normCounts, gcContent=gcContent, ...)
    } else if(type == "iqr"){
      rleValue <- rleGC_iqr(counts=normCounts, gcContent=gcContent, ...)
    }
    return(rleValue)
  })
  names(rleGCVals) <- h5List$name[2:(nrow(h5List)-1)]
  return(rleGCVals)
}


h5File="/u/home/d/dr2camac/project-geschwind/Diego_Internship_Summary/ATAC-seq/data/scoreSconefit.h5"


## add GC-bias through RLE to scores
rleGC_med <- rleGCh5(h5File = h5File, 
                  gcContent = gcContentScone,
                  type="med")
rleGC_iqr <- rleGCh5(h5File = h5File, 
                  gcContent = gcContentScone,
                  type="iqr")
scores <- get_scores(sconeFit)
scores <- cbind(scores, 
                rleGC_med = -rleGC_med[rownames(scores)],
                rleGC_iqr = -rleGC_iqr[rownames(scores)])
## recalculate score ranks with new metrics
scoreT <- t(scores)
ranked_scores <- apply(scoreT,1,rank)
mean_score_rank <- rowMeans(ranked_scores)
mean_score_rank <- mean_score_rank[order(mean_score_rank, decreasing=TRUE)]
mean_score_rank

pc_obj <- prcomp(apply(scoreT,1,rank),
                center = TRUE,scale = FALSE)
methodLabelsTmp <- strsplit(rownames(pc_obj$x), split=",")
methodNames <- unlist(lapply(methodLabelsTmp, function(x){
  paste0(x[2:3], collapse=",")
}))
methods <- unlist(lapply(strsplit(methodNames, split=","), "[[", 1))
ruvs <- unlist(lapply(strsplit(methodNames, split=","), function(x) strsplit(x[2],split="_")))[seq(2,length(methodNames)*2, by=2)]
methodNames <- paste(methods,ruvs, sep = "_")

pdf("../../results/04_NormalizationBenchmark/sconeNormTrayectories.pdf")
bp_obj <- biplot_color(pc_obj,y = -mean_score_rank,expand = .6, cex=2)
text(x=bp_obj[,1], y=bp_obj[,2], labels=methodNames, cex=2/3)

library(ggplot2)
library(tidyverse)
getMeanScoreRank <- function(scores){
  scoreT <- t(scores)
  ranked_scores <- apply(scoreT,1,rank)
  mean_score_rank <- rowMeans(ranked_scores)
  mean_score_rank <- mean_score_rank[order(mean_score_rank, decreasing=TRUE)]
  return(mean_score_rank)
}
gcNormMethods <- c("FQ-FQ", "GC-FQ", "cqn_length", "cqn", "GC-FQ_smooth")
rmMethods <- c("gcqn_mean", "gcqn_median_10", "gcqn_median_20",  "gcqn_median_100",
               "gcqn_median_permuted")

load("sconeTrials")

scores <- lapply(scores, function(scoreMat){
  scoreNames <- rownames(scoreMat)
  rmID <- unlist(sapply(rmMethods, function(method){
    grep(x=scoreNames, pattern=method)
  }))
  return(scoreMat[-rmID,])
})

scoreRanks <- lapply(scores, getMeanScoreRank)
rankVector <- do.call(c,scoreRanks)
scoreDfLong <- data.frame(method=names(rankVector),
                          dataset=unlist(mapply(rep,datasets,unlist(lapply(scoreRanks,length)))),
                          rank=rankVector)


pal <- wesanderson::wes_palette("Zissou1", n=40, type="continuous")
id0 <- grep(x=scoreDfLong$method, pattern="no_uv")
scoreDf0 <- scoreDfLong[id0,]
scoreDf0$methodShort <- unlist(lapply(strsplit(as.character(scoreDf0$method), split=","), "[[", 2))
scoreDf0$methodShort[scoreDf0$methodShort == "gcqn_median_50"] <- "GC-FQ"
scoreDf0$methodShort[scoreDf0$methodShort == "edaseq"] <- "FQ-FQ"
scoreDf0$methodShort[scoreDf0$methodShort == "gcqn_smooth"] <- "GC-FQ_smooth"
scoreDf0$methodShort[scoreDf0$methodShort == "fq"] <- "FQ"
scoreDf0$methodShort[scoreDf0$methodShort == "tmm"] <- "TMM"
scoreDf0$methodShort[scoreDf0$methodShort == "deseq"] <- "DESeq2"
scoreDf0$methodShort[scoreDf0$methodShort == "uq"] <- "UQ"
scoreDf0$methodShort[scoreDf0$methodShort == "none"] <- "None"
scoreDf0$methodShort[scoreDf0$methodShort == "sum"] <- "Sum"
scoreDf0$methodShort <- factor(scoreDf0$methodShort)
## scale ranks
scoreDf0 <- scoreDf0 %>% group_by(dataset) %>% mutate(scaledRank=scale(rank))
## take average of scaled ranks
avRank0 <- scoreDf0 %>% group_by(methodShort) %>% summarize(avRank=mean(scaledRank))
oo0 <- order(avRank0$avRank, decreasing=TRUE)
orderedMethodShort0 <- avRank0$methodShort[oo0]
scoreDfOrdered0 <- scoreDf0
scoreDfOrdered0$methodShort <- factor(scoreDfOrdered0$methodShort, levels=orderedMethodShort0)
#pdf("../../results/04_NormalizationBenchmark/heatMapFullEval.pdf")

pFull <- ggplot(scoreDfOrdered0, aes(y=methodShort, x=dataset, fill=scaledRank)) + 
  geom_tile() + 
  scale_fill_gradientn(colours = pal, name = 'Scaled Scone\nrank', limits=c(-3,2.5)) +
  ggtitle("Full evaluation") +
  ylab("") + xlab("Dataset") +
  theme(axis.text.x = element_text(angle = 30, vjust= 0.7, size=9),
        axis.text.y = element_text(size=9)) + 
  theme(panel.background = element_rect(fill="white")) + 
  theme(axis.text.y=element_text(colour=ifelse(orderedMethodShort0 %in% gcNormMethods, "dodgerblue", "black")))


