#Sediment rate calculation for pacific and atlantic pages2k data using linear regression 

#load data
setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/contour_data/")
aSector <- read.csv("pages2k_marine_bacon_atlanticsector.csv", stringsAsFactors=FALSE)
pSector <- read.csv("pages2k_marine_bacon_pacificsector.csv", stringsAsFactors=FALSE)
SSTreplaced <- read.csv("pages2k_marine_bacon_SSTreplaced.csv")
O18replaced <- read.csv("pages2k_marine_bacon_O18replaced.csv")

#pull out depth and age
aSector <- data.frame(aSector[, c("CoreNo", "Depth..cm.", "New.Age..CE.", "CoreId")])
pSector <- data.frame(pSector[, c("CoreNo", "Depth..cm.", "New.Age..CE.", "CoreId")])

#pull out core ID
SSTreplaced <- data.frame(SSTreplaced[, c("dataSetName", "lat", "lon")])
O18replaced <- data.frame(O18replaced[, c("dataSetName", "lat", "lon")])

#combine datasets
depthDat <- rbind(aSector, pSector)
sstO18 <- rbind(SSTreplaced, O18replaced)

#-----atlantic sector-----
aSecRateDat <- matrix(data = NA, 24, 5)
colnames(aSecRateDat) <- c("coreNo", "sedRate","coreID", "lat", "lon")
for (i in 1:24) {
    num <- which(aSector$CoreNo == i)
    aCoreDat <- data.frame(matrix(data = NA, length(num), 4))
    colnames(aCoreDat) <- c("coreNo", "depth", "age","coreID")
    aCoreDat[,1] <- aSector[num, 1]
    aCoreDat[,2] <- aSector[num, 2]
    aCoreDat[,3] <- aSector[num, 3]
    aCoreDat[,4] <- aSector[num, 4]
    
    if (i == 20 | i == 21 | i == 22 | i == 23 | i == 24) {
      fit <- 0
      sedRate <- 0
      } else {
    fit <- lm(depth ~ age, data = aCoreDat) #linear regression
    sedRate <- round(-1000 * fit[["coefficients"]][["age"]], digits = 1)
      }
    
    coreID <- aCoreDat[1,4]
    
    aSecRateDat[i,1] <- i
    aSecRateDat[i,2] <- sedRate
    aSecRateDat[i,3] <- aCoreDat[1,4]
    
    if (i == 8) {
      aSecRateDat[i,4] <- 2.5
      aSecRateDat[i,5] <- 9.38
    } else {
      row <- grep(coreID, sstO18[,1])
      aSecRateDat[i,4] <- sstO18[row[1], 2]
      aSecRateDat[i,5] <- sstO18[row[1], 3]
    }
}

#-----pacific sector-----
pSecRateDat <- matrix(data = NA, 51, 5)
colnames(pSecRateDat) <- c("coreNo", "sedRate","coreID", "lat", "lon")
for (j in 41:51) {
  num <- which(pSector$CoreNo == j)
  pCoreDat <- data.frame(matrix(data = NA, length(num), 4))
  colnames(pCoreDat) <- c("coreNo", "depth", "age","coreID")
  pCoreDat[,1] <- pSector[num, 1]
  pCoreDat[,2] <- pSector[num, 2]
  pCoreDat[,3] <- pSector[num, 3]
  pCoreDat[,4] <- pSector[num, 4]
  
  if (j == 43 | j == 46 | j == 47 | j == 48 | j == 49) {
    fit <- 0
    sedRate <- 0
  } else {
    fit <- lm(depth ~ age, data = pCoreDat) #linear regression
    sedRate <- round(-1000 * fit[["coefficients"]][["age"]], digits = 1)
  }
  
  coreID <- pCoreDat[1,4]
  
  pSecRateDat[j,1] <- j
  pSecRateDat[j,2] <- sedRate
  pSecRateDat[j,3] <- pCoreDat[1,4]
  if (j == 48 | j == 49) {
    pSecRateDat[j,4] <- 4.67
    pSecRateDat[j,5] <- -77.96
  } else {
    row <- grep(coreID, sstO18[,1])
    pSecRateDat[j,4] <- sstO18[row[1], 2]
    pSecRateDat[j,5] <- sstO18[row[1], 3]
  }
}

pSecRateDat <- na.omit(pSecRateDat)

#combine data
setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/contour_data/")
rateDat <- rbind(aSecRateDat, pSecRateDat)
write.csv(rateDat,"pages2kRateDat.csv")
write.csv(aSecRateDat, "aSecRateDat.csv")
write.csv(pSecRateDat, "pSecRateDat.csv")

#contour plot
aSecRateDat <- data.frame(aSecRateDat[, c("sedRate", "lat", "lon")], stringsAsFactors = FALSE)
sres <- as.numeric(as.character(aSecRateDat$sedRate))
lat <- as.numeric(as.character(aSecRateDat$lat))
lon <- as.numeric(as.character(aSecRateDat$lon))
lati<-seq(round(min(lat)),round(max(lat)),1)		#interpolated latitude to 1º
loni<-seq(round(min(lon)),round(max(lon)),1)		#interpolated longitude to 1º

int<-interp(lon,lat,sres,xo=loni,yo=lati, duplicate = "mean")

hised<-which(int$z>=20, arr.ind = TRUE)								#below gives coordinates of indicies with sedrates greater than 12. 	
hilat<-floor(hised/length(loni)) +1				#this seemed like a trick to get the correct column (e.g. latidude)
hilon<-hised-((hilat-1)*length(loni)) +1		#a similar trick to get longitude

contour(int$x, int$y,int$z,las=1,xlab="longitude",ylab="latitude", 
               main = "Sedimentation Rate (cm/ka) - pages2k", 
               xlim = c(-100, 160), ylim = c(-70,70))
map(add = TRUE, xlim = c(-100, 160), ylim = c(-70,70))
points(loni[hilon],lati[hilat],pch=8, 
       xlim = c(-100, 160), ylim = c(-70,70))

output<-cbind(loni[hilon],lati[hilat])
write.csv(output,"hiSedRateSites.csv")