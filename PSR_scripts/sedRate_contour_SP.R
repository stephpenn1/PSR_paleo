#Contour plot for late-holocene marine sedimentation rate data in the Atlantic Ocean
#Created 11-9-17

library(maps)
library(ggplot2)
library(ggmap)
setwd("/Users/sp/Desktop/PSR_paleo/core_data/")
data <- read.csv("Barash_1977.csv")


data <- data.frame(data[, c("Longitude", "Latitude", "Sed.rate..cm.ka.")])
#sedRate <- (matrix(data$Sed.rate..cm.ka.))
#Latitude <- as.vector(matrix(data$Latitude)) + 180
#Longitude <- as.vector(matrix(data$Longitude)) + 180

#test 1
map <- get_map(location = c(-100,-90,80,80), source = "google", maptype = "watercolor")
long<-c(-90,8)
lat<-c(-60,65)
c<- ggmap(map)  
c<- c + geom_point(data=data, aes(x=Longitude, y=Latitude,  colour = Sed.rate..cm.ka.))
c<- c + scale_colour_gradient(low = "pink",high = "red")
c<- c + labs(colour = "Sedimentation Rate cm/ka")
c<- c + contour(x=Longitude, y=Latitude, z = sedRate)

#test 2
#c<- ggmap(map)  
c<- ggplot(data, aes(x=Longitude, y=Latitude,  z = sedRate))
c<- c + geom_point(data=data, aes(x=Longitude, y=Latitude,  colour = Sed.rate..cm.ka.))
c<- c + scale_colour_gradient(low = "pink",high = "red")
c<- c + labs(colour = "Sedimentation Rate cm/ka")
c<- c + stat_contour()

#test 3
library(akima)
fld <- with(df, interp(x = Longitude, y = Latitude, z = Sed.rate..cm.ka.))
ggmap(map)
map('world', xlim = c(-85,15),
      ylim = c(-40,65))
par(new = T)
contour(x = fld$x, y = fld$y, z = fld$z,
               #color.palette =
              #   colorRampPalette(c("white", "red")),
               col = "blue",
               xlab = "Longitude",
               ylab = "Latitude",
               xlim = c(-85,15),
               ylim = c(-40,65) , 
               main = "Sedimentation Rate \n Barash data")
               #key.title = title(main = "Rate (cm/ka)", cex.main = 1))

#test 4 (Barash data)
library(ggplot2)
library(reshape2)
library(ggmap)
library(akima)

setwd("/Users/sp/Desktop/PSR_paleo/core_data/")
data <- read.csv("Barash_1977.csv")
df <- data.frame(data[, c("Longitude", "Latitude", "Sed.rate..cm.ka.")])
fld <- with(df, interp(x = Longitude, y = Latitude, z = Sed.rate..cm.ka.))

#prepare data
#df$Longitude <- df$Longitude + 180
#df$Latitude <- df$Latitude + 180
df <- melt(fld$z, na.rm = TRUE)
names(df) <- c("x", "y", "sedRate")
df$Longitude <- fld$x[df$x]
df$Latitude <- fld$y[df$y]

map <- get_map(location = c(-80,-35,10,65), source = "stamen", maptype = "watercolor")

b <- ggplot(data = df, aes(x = Longitude, y = Latitude, z = sedRate)) +
  xlim(-80,10) +
  ylim(-35,65) +
  geom_tile(aes(fill = sedRate)) +
  stat_contour(col = "black") +
  ggtitle("Sedimentation Rates \n Barash data") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_continuous(name = "Rate (cm/ka)",
                        low = "blue", high = "red") +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10, vjust = -0.5),
        axis.title.y = element_text(size = 10, vjust = 0.2),
        legend.text = element_text(size = 10))


#test with O18 and SST proxy data
setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/proxy_data/")
O18proxyDat <- read.csv("O18ProxyDat.csv")
SSTproxyDat <- read.csv("SSTProxyDat.csv")

#create data frame of sediment rates and coordinates
O18proxy <- na.omit(data.frame(O18proxyDat[, c("lon", "lat", "SampleRes")]))
SSTproxy <- (data.frame(MgCa[, c("lon", "lat", "SampleRes")]))
fldO18 <- with(O18proxy, interp(x = lon, y = lat, z = SampleRes, duplicate = "mean"))
fldSST <- with(SSTproxy, interp(x = lon, y = lat, z = SampleRes, duplicate = "mean"))

dfO18 <- melt(fldO18$z, na.rm = TRUE)
names(dfO18) <- c("x", "y", "SampleRes")
dfO18$lon <- fldO18$x[dfO18$x]
dfO18$lat <- fldO18$y[dfO18$y]

dfSST <- melt(fldSST$z, na.rm = TRUE)
names(dfSST) <- c("x", "y", "SampleRes")
dfSST$lon <- fldSST$x[dfSST$x]
dfSST$lat <- fldSST$y[dfSST$y]

s <- ggplot(data = dfSST, aes(x = lon, y = lat, z = SampleRes)) +
  geom_tile(aes(fill = SampleRes)) +
  stat_contour(col = "black") +
  ggtitle("Sedimentation Rates \n SST proxy") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_continuous(name = "Rate (cm/ka)",
                        low = "blue", high = "red") +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10, vjust = -0.5),
        axis.title.y = element_text(size = 10, vjust = 0.2),
        legend.text = element_text(size = 10))

p <- ggplot(data = dfO18, aes(x = lon, y = lat, z = SampleRes)) +
  geom_tile(aes(fill = SampleRes)) +
  stat_contour(col = "black") +
  ggtitle("Sedimentation Rates \n O18 proxy") +
  xlab("Longitude") +
  ylab("Latitude") +
  scale_fill_continuous(name = "Rate (cm/ka)",
                        low = "blue", high = "red") +
  theme(plot.title = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 10),
        axis.title.x = element_text(size = 10, vjust = -0.5),
        axis.title.y = element_text(size = 10, vjust = 0.2),
        legend.text = element_text(size = 10))


map <- get_map(location = c(-150,-50,60,90), source = "stamen", maptype = "watercolor")
c<- ggmap(map) 
c<- c + geom_point(data=O18proxy, aes(x=lon, y=lat,  colour = SampleRes))
c<- c + scale_colour_gradient(low = "pink",high = "red")
c<- c + labs(colour = "Sedimentation Rate cm/ka")

