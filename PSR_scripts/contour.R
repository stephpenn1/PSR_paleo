library(akima)
library(maps)
#-----load sediment rate datasets-----

#load Barash data (in cm/ka)
setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/contour_data/")
barash <- read.csv("Barash_1977.csv")
barashProxy <- data.frame(barash[, c("Longitude", "Latitude", "Sed.rate..cm.ka.")])
colnames(barashProxy) <- c("lon", "lat", "sedRate")

#load O18 and SST proxy data (in samples/century)
setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/proxy_data/")
O18 <- read.csv("O18ProxyDat.csv")
SST <- read.csv("SSTProxyDat.csv")
O18proxy <- na.omit(data.frame(O18[, c("lon", "lat", "SampleRes")]))
SSTproxy <- na.omit(data.frame(SST[, c("lon", "lat", "SampleRes")]))
O18proxy$SampleRes <- O18proxy$SampleRes / 10  #samples/century to cm/ka
SSTproxy$SampleRes <- SSTproxy$SampleRes / 10
colnames(O18proxy) <- c("lon", "lat", "sedRate")
colnames(SSTproxy) <- c("lon", "lat", "sedRate")

#combine datasets
proxyComp <- rbind(O18proxy, SSTproxy)

#load sector data (in cm/ka)
setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/contour_data/")
aSecRateDat <- read.csv("aSecRateDat.csv")
pSecRateDat <- read.csv("pSecRateDat.csv")
rateDat <- read.csv("pages2kRateDat.csv")

aSecRateDat <- data.frame(aSecRateDat[, c("lat", "lon", "sedRate")], stringsAsFactors = FALSE)

proxyComp <- data.frame(proxyComp[, c("lat", "lon", "SampleRes")])
colnames(proxyComp) <- c("lat", "lon","sedRate")
rateDat <- data.frame(rateDat[, c("lat", "lon", "sedRate")])

#atlantic contour
aComp <- rbind(barashProxy, aSecRateDat)

comp <- rbind(aComp, proxyComp)

#contour plot
sres <- as.numeric(as.character(aComp$sedRate))
lat <- as.numeric(as.character(aComp$lat))
lon <- as.numeric(as.character(aComp$lon))
lati<-seq(round(min(lat)),round(max(lat)),1)		#interpolated latitude to 1º
loni<-seq(round(min(lon)),round(max(lon)),1)		#interpolated longitude to 1º

int<-interp(lon,lat,sres,xo=loni,yo=lati, duplicate = "mean")

hised<-which(int$z>=10, arr.ind = TRUE)		#below gives coordinates of indicies with sedrates greater than 12. 	
hilat<-floor(hised/length(loni)) +1				#this seemed like a trick to get the correct column (e.g. latidude)
hilon<-hised-((hilat-1)*length(loni)) -7	#a similar trick to get longitude

levels <- c(0, 5, 10, 20)
png(filename="/Users/sp/Desktop/PSR_paleo/PSR_data/contour_data/test.png")
contour(int$x, int$y,int$z,las=1,xlab="longitude",ylab="latitude",
        main = "Sedimentation Rate (cm/ka) - Atlantic",
        levels = levels, col = "royalblue3", xlim = c(-100, 160), ylim = c(-70,70))
map(add = TRUE, xlim = c(-100, 160), ylim = c(-70,70))
points(loni[hised[,1]],lati[hised[,2]],pch=16, xlim = c(-100, 160), ylim = c(-70,70))
dev.off()

output<-cbind(loni[hised[,1]],lati[hised[,2]])
colnames(output) <- c("lon", "lat")
write.csv(output,"hiSedRateSites.csv")
