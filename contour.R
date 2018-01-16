
library(akima)

#load sediment rate datasets
setwd("/Users/sp/Desktop/PSR_paleo/core_data/")
barash <- read.csv("Barash_1977.csv")
barashProxy <- data.frame(barash[, c("Longitude", "Latitude", "Sed.rate..cm.ka.")])
colnames(barashProxy) <- c("lon", "lat", "SampleRes")
barashProxy$SampleRes <- barashProxy$SampleRes / 10 #cm/ka to samples/century

setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/proxy_data/")
O18 <- read.csv("O18ProxyDat.csv")
SST <- read.csv("SSTProxyDat.csv")
O18proxy <- na.omit(data.frame(O18[, c("lon", "lat", "SampleRes")]))
SSTproxy <- na.omit(data.frame(SST[, c("lon", "lat", "SampleRes")]))

#combine datasets
proxyComp <- rbind(O18proxy, SSTproxy, barashProxy)

sres <- proxyComp$SampleRes
lon <- proxyComp$lon
lat <- proxyComp$lat
lati<-seq(round(min(lat)),round(max(lat)),1)		#interpolated latitude to 1º
loni<-seq(round(min(lon)),round(max(lon)),1)		#interpolated longitude to 1º

int<-interp(lon,lat,sres,xo=loni,yo=lati, duplicate = "mean")

hised<-which(int$z>=20, arr.ind = TRUE)								#below gives coordinates of indicies with sedrates greater than 12. 	
hilat<-floor(hised/length(loni)) +1				#this seemed like a trick to get the correct column (e.g. latidude)
hilon<-hised-((hilat-1)*length(loni)) +1		#a similar trick to get longitude

filled.contour(int$x, int$y,int$z,las=1,xlab="longitude",ylab="latitude", 
               col = topo.colors(19), xlim = c(-100, 160), ylim = c(-70,70))
points(loni[hilon],lati[hilat],pch=8, 
       xlim = c(-100, 160), ylim = c(-70,70))

output<-cbind(loni[hilon],lati[hilat])
write.csv(output,"hiSedRateSites.csv")