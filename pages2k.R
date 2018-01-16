#load data
setwd("/Users/sp/Desktop/")
aSector <- read.csv("pages2k_marine_bacon_atlanticsector.csv", stringsAsFactors=FALSE)
pSector <- read.csv("pages2k_marine_bacon_pacificsector.csv", stringsAsFactors=FALSE)
SSTreplaced <- read.csv("pages2k_marine_bacon_SSTreplaced.csv")
O18replaced <- read.csv("pages2k_marine_bacon_O18replaced.csv")

#pull out depth and age
aSector <- data.frame(aSector[, c("CoreNo", "Depth..cm.", "New.Age..CE.")])
pSector <- data.frame(pSector[, c("CoreNo", "Depth..cm.", "New.Age..CE.")])

#combine datasets
depthDat <- rbind(aSector, pSector)

for (i in 1:51) {
  num <- which(depthDat$CoreNo == 41)
  coreDat <- data.frame(matrix(data = NA, length(num), 3))
  colnames(coreDat) <- c("coreNo", "depth", "age")
  coreDat[,1] <- depthDat[num, 1]
  coreDat[,2] <- depthDat[num, 2]
  coreDat[,3] <- depthDat[num, 3]
  #add linear regression
  fit <- lm(depth ~ age, data = coreDat)
  sedRate <- -100 * fit[["age"]]
}

