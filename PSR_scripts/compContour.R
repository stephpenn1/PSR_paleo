library(ggplot2)
library(reshape2)
library(ggmap)
library(akima)

setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/contour_data/")
proxyComp <- read.csv("proxyRateDat.csv")
rateDat <- read.csv("pages2kRateDat.csv")

proxyComp <- data.frame(proxyComp[, c("lat", "lon", "SampleRes")])
colnames(proxyComp) <- c("lat", "lon","sedRate")
rateDat <- data.frame(rateDat[, c("lat", "lon", "sedRate")])
comp <- rbind(proxyComp, rateDat)

#contour plot
sres <- as.numeric(as.character(comp$sedRate))
lat <- as.numeric(as.character(comp$lat))
lon <- as.numeric(as.character(comp$lon))
lati<-seq(round(min(lat)),round(max(lat)),1)		#interpolated latitude to 1º
loni<-seq(round(min(lon)),round(max(lon)),1)		#interpolated longitude to 1º

int<-interp(lon,lat,sres,xo=loni,yo=lati, duplicate = "mean")

hised<-which(int$z>=20, arr.ind = TRUE)		#below gives coordinates of indicies with sedrates greater than 12. 	
hilat<-floor(hised/length(loni)) +1				#this seemed like a trick to get the correct column (e.g. latidude)
hilon<-hised-((hilat-1)*length(loni)) -7	#a similar trick to get longitude

levels <- c(0, 1.25, 2.5, 5, 10, 15, 20, 30, 90, 120, 300)
png(filename="/Users/sp/Desktop/PSR_paleo/PSR_data/contour_data/test.png")
contour(int$x, int$y,int$z,las=1,xlab="longitude",ylab="latitude",
               main = "Sedimentation Rate (cm/ka) - composite", nlevels =10, labcex = 0.7,
               levels = levels, col = "firebrick1", xlim = c(-100, 160), ylim = c(-70,70))
map(add = TRUE, xlim = c(-100, 160), ylim = c(-70,70))
points(loni[hised[,1]],lati[hised[,2]],pch=16, xlim = c(-100, 160), ylim = c(-70,70))
dev.off()

par(new=TRUE)
plot(lon, lat, xlim = c(-100, 160), ylim = c(-70,70), col = "red")

image(int$x, int$y,int$z)
map(add = TRUE, xlim = c(-100, 160), ylim = c(-70,70), wrap = FALSE)

#try now with ggplot
fld <- with(compTest, interp(x = lon, y = lat, z = sres, duplicate = "mean"))

compTest <- data.frame(cbind(lon,lat,sres))

#prepare data
df <- melt(fld$z, na.rm = TRUE)
names(df) <- c("x", "y", "sedRate")
df$Longitude <- fld$x[df$x]
df$Latitude <- fld$y[df$y]
map <- get_map(location = c(-80,-35,10,65), source = "stamen", maptype = "watercolor")

map <- get_map(location = c(-100,-90,80,80), source = "google", maptype = "watercolor")
long<-c(-90,8)
lat<-c(-60,65)
c<- ggmap(map)  
c<- c + geom_point(data=comp, aes(x=lon, y=lat,  colour = sres))
c<- c + scale_colour_gradient(low = "white",high = "red")
c<- c + labs(colour = "Sedimentation Rate cm/ka")
