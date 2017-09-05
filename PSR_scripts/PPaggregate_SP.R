#Stephanie Pennington | PSR summer research
#add bivalve and sclero sponge data back into O18 pseudoproxy timeseries after smoothing marine sediment data
#Created 7-17-17
#edited 8-7-17 to fix index issue and remove MgCa

library(miscTools)

#read in pseudoproxy timeseries with data
setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/noise")
HadpiC_O18_PP<-read.csv("HadpiC_O18_pp.csv")
GISSgCpiC_O18_PP<-read.csv("GISSgCpiC_O18_pp.csv")
GISSgy3piC_O18_PP<-read.csv("GISSgy3piC_O18_pp.csv")
GISSgTckLM_O18_PP<-read.csv("GISSgTckLM_O18_pp.csv")
GISSgTKckLM_O18_PP<-read.csv("GISSgTKckLM_O18_pp.csv")
GISSgTcsLM_O18_PP<-read.csv("GISSgTcsLM_O18_pp.csv")

#read in smoothed marine sediment pseudoproxy data
setwd("/home/spenn1/PSR_paleo/PSR_data/pseudoproxy/smooth")
HadpiC_O18_PPsmooth<-read.csv("HadpiC_O18_PPsmooth.csv")
GISSgCpiC_O18_PPsmooth<-read.csv("GISSgCpiC_O18_PPsmooth.csv")
GISSgy3piC_O18_PPsmooth<-read.csv("GISSgy3piC_O18_PPsmooth.csv")
GISSgTckLM_O18_PPsmooth<-read.csv("GISSgTckLM_O18_PPsmooth.csv")
GISSgTKckLM_O18_PPsmooth<-read.csv("GISSgTKckLM_O18_PPsmooth.csv")
GISSgTcsLM_O18_PPsmooth<-read.csv("GISSgTcsLM_O18_PPsmooth.csv")

#change from data frame to matrix 
HadpiC_O18_PPsmooth_m<-data.matrix(HadpiC_O18_PPsmooth)
GISSgCpiC_O18_PPsmooth_m<-data.matrix(GISSgCpiC_O18_PPsmooth)
GISSgy3piC_O18_PPsmooth_m<-data.matrix(GISSgy3piC_O18_PPsmooth)
GISSgTckLM_O18_PPsmooth_m<-data.matrix(GISSgTckLM_O18_PPsmooth)
GISSgTKckLM_O18_PPsmooth_m<-data.matrix(GISSgTKckLM_O18_PPsmooth)
GISSgTcsLM_O18_PPsmooth_m<-data.matrix(GISSgTcsLM_O18_PPsmooth)

otherIndex<- which(metadataO18$archiveType != "marine sediment" & metadataO18$archiveType != "marine sediments")

#pull out bivalve and sclerosponge archive type data from original PPs
HadpiC_O18_PPother<-HadpiC_O18_PP[otherIndex,]
GISSgCpiC_O18_PPother<-GISSgCpiC_O18_PP[otherIndex,]
GISSgy3piC_O18_PPother<-GISSgy3piC_O18_PP[otherIndex,]
GISSgTckLM_O18_PPother<-GISSgTckLM_O18_PP[otherIndex,]
GISSgTKckLM_O18_PPother<-GISSgTKckLM_O18_PP[otherIndex,]
GISSgTcsLM_O18_PPother<-GISSgTcsLM_O18_PP[otherIndex,]

#HadpiC_O18
data<-HadpiC_O18_PPsmooth_m

for (i in 1:length(otherIndex)) {
  new<-insertRow(data, otherIndex[i],HadpiC_O18_PPother[i,])
  data<-new
}

HadpiC_O18_PPagg<-data

#GISSgCpiC_O18
data<-GISSgCpiC_O18_PPsmooth_m

for (i in 1:length(otherIndex)) {
  new<-insertRow(data, otherIndex[i],GISSgCpiC_O18_PPother[i,])
  data<-new
}

GISSgCpiC_O18_PPagg<-data

#GISSgy3piC_O18
data<-GISSgy3piC_O18_PPsmooth_m

for (i in 1:length(otherIndex)) {
  new<-insertRow(data, otherIndex[i],GISSgy3piC_O18_PPother[i,])
  data<-new
}

GISSgy3piC_O18_PPagg<-data

#GISSgTckLM_O18
data<-GISSgTckLM_O18_PPsmooth_m

for (i in 1:length(otherIndex)) {
  new<-insertRow(data, otherIndex[i],GISSgTckLM_O18_PPother[i,])
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
setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/aggregate/")
write.csv(HadpiC_O18_PPagg, file = "HadpiC_O18_PPagg.csv", row.names = FALSE)
write.csv(GISSgCpiC_O18_PPagg, file = "GISSgCpiC_O18_PPagg.csv", row.names = FALSE)
write.csv(GISSgy3piC_O18_PPagg, file = "GISSgy3piC_O18_PPagg.csv", row.names = FALSE)
write.csv(GISSgTckLM_O18_PPagg, file = "GISSgTckLM_O18_PPagg.csv", row.names = FALSE)
write.csv(GISSgTKckLM_O18_PPagg, file = "GISSgTKckLM_O18_PPagg.csv", row.names = FALSE)
write.csv(GISSgTcsLM_O18_PPagg, file = "GISSgTcsLM_O18_PPagg.csv", row.names = FALSE)
