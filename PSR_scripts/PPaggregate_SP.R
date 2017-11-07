#Stephanie Pennington | PSR summer research
#add bivalve and sclero sponge data back into O18 pseudoproxy timeseries after smoothing marine sediment data (before binning)
#Created 7-17-17
#edited 8-7-17 to fix index issue and remove MgCa

library(miscTools)

#read in pseudoproxy timeseries data
setwd("/home/spenn1/PSR_paleo/PSR_data/pseudoproxy/model_data/")
GISSgTckLM_O18<-read.csv("GISSgTckLM_O18.csv")
GISSgTKckLM_O18<-read.csv("GISSgTKckLM_O18.csv")
GISSgTcsLM_O18<-read.csv("GISSgTcsLM_O18.csv")

#read in smoothed marine sediment pseudoproxy data
setwd("/home/spenn1/PSR_paleo/PSR_data/pseudoproxy/smooth")
GISSgTckLM_O18_PPsmooth<-read.csv("GISSgTckLM_O18_PPsmooth.csv")
GISSgTKckLM_O18_PPsmooth<-read.csv("GISSgTKckLM_O18_PPsmooth.csv")
GISSgTcsLM_O18_PPsmooth<-read.csv("GISSgTcsLM_O18_PPsmooth.csv")

#change from data frame to matrix 
GISSgTckLM_O18_PPsmooth_m<-data.matrix(GISSgTckLM_O18_PPsmooth)
GISSgTKckLM_O18_PPsmooth_m<-data.matrix(GISSgTKckLM_O18_PPsmooth)
GISSgTcsLM_O18_PPsmooth_m<-data.matrix(GISSgTcsLM_O18_PPsmooth)

#find non marine sediment indices from metadata
otherIndex<- which(metadataO18$archiveType != "marine sediment" & metadataO18$archiveType != "marine sediments")

#pull out bivalve and sclerosponge archive type data from original PPs
GISSgTckLM_O18_PPother<-GISSgTckLM_O18[otherIndex,]
GISSgTKckLM_O18_PPother<-GISSgTKckLM_O18[otherIndex,]
GISSgTcsLM_O18_PPother<-GISSgTcsLM_O18[otherIndex,]


#GISSgTckLM_O18
data<-GISSgTckLM_O18_PPsmooth_m
for (i in 1:length(otherIndex)) {
  new<-insertRow(data, otherIndex[i],GISSgTckLM_O18_PPother[i,]) #insert rows into smoothed matrix
  data<-new
}
GISSgTckLM_O18_PPagg<-data

#GISSgTKckLM_O18
data<-GISSgTKckLM_O18_PPsmooth_m
for (i in 1:length(otherIndex)) {
  new<-insertRow(data, otherIndex[i],GISSgTKckLM_O18_PPother[i,])
  data<-new
}
GISSgTKckLM_O18_PPagg<-data

#GISSgTcsLM_O18
data<-GISSgTcsLM_O18_PPsmooth_m
for (i in 1:length(otherIndex)) {
  new<-insertRow(data, otherIndex[i],GISSgTcsLM_O18_PPother[i,])
  data<-new
}
GISSgTcsLM_O18_PPagg<-data

#save aggreagate timeseries for binning
setwd("/home/spenn1/PSR_paleo/PSR_data/pseudoproxy/aggregate/")
write.csv(GISSgTckLM_O18_PPagg, file = "GISSgTckLM_O18_PPagg.csv", row.names = FALSE)
write.csv(GISSgTKckLM_O18_PPagg, file = "GISSgTKckLM_O18_PPagg.csv", row.names = FALSE)
write.csv(GISSgTcsLM_O18_PPagg, file = "GISSgTcsLM_O18_PPagg.csv", row.names = FALSE)
