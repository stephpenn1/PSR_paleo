#Stephanie Pennington | PSR summer research
#Pseudoproxy timeseries smoothing
#apply weighed moving average with hanning filter
#Needs metadata from Proxybin_anomIND script
#Created 7-6-17

setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/pseudoproxy/model_data/")
library(signal)

#pull out only marine sediment archive type ------------------------------

sedIndex_O18<- which(metadataO18$archiveType == "marine sediment" | metadataO18$archiveType == "marine sediments")
res_O18<-metadataO18[,4]      #pull avg sample resolution from O18 metadata
resIndex_O18<-res_O18[sedIndex_O18]      #pull out only marine sediment sample resolutions
resIndex_MgCa<-MgCa[,5]     #pull out sample resolution from SST metadata


GISSgTckLM_O18<-read.csv(file = "GISSgTckLM_O18.csv")
GISSgTckLM_MgCa<-read.csv(file = "GISSgTckLM_MgCa.csv")
GISSgTKckLM_O18<-read.csv(file = "GISSgTKckLM_O18.csv")
GISSgTKckLM_MgCa<-read.csv(file = "GISSgTKckLM_MgCa.csv")
GISSgTcsLM_O18<-read.csv(file = "GISSgTcsLM_O18.csv")
GISSgTcsLM_MgCa<-read.csv(file = "GISSgTcsLM_MgCa.csv")

#pull out rows with marine sediment - O18
GISSgTckLM_O18_PPmarSed<-GISSgTckLM_O18[sedIndex_O18,]
GISSgTKckLM_O18_PPmarSed<-GISSgTKckLM_O18[sedIndex_O18,]
GISSgTcsLM_O18_PPmarSed<-GISSgTcsLM_O18[sedIndex_O18,]

#*note: O18 now has 40 rows*
setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/pseudoproxy/marine_sediments");

#save as CSV for smoothing
write.csv(GISSgTckLM_O18_PPmarSed, file = "GISSgTckLM_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTKckLM_O18_PPmarSed, file = "GISSgTKckLM_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTcsLM_O18_PPmarSed, file = "GISSgTcsLM_O18_PPmarSed.csv", row.names = FALSE)

# calculate weighted moving average ---------------------------------------
setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/pseudoproxy/marine_sediments");
GISSgTckLM_O18_PPmarSed<-read.csv("GISSgTckLM_O18_PPmarSed.csv")
GISSgTKckLM_O18_PPmarSed<-read.csv("GISSgTKckLM_O18_PPmarSed.csv")
GISSgTcsLM_O18_PPmarSed<-read.csv("GISSgTcsLM_O18_PPmarSed.csv")
setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/pseudoproxy/model_data/");
GISSgTckLM_MgCa<-read.csv("GISSgTckLM_MgCa.csv")
GISSgTKckLM_MgCa<-read.csv("GISSgTKckLM_MgCa.csv")
GISSgTcsLM_MgCa<-read.csv("GISSgTcsLM_MgCa.csv")


#GISSgTckLM_O18
data<-GISSgTckLM_O18_PPmarSed
s<-seq(1,ncol(data)) #set length of data

GISSgTckLM_O18_PPsmooth<-matrix(data = NA, nrow(data), ncol(data)) #create empty matrix for smoothed data

for(i in 1:nrow(data)) {
  sampleRes<-(resIndex_O18[i])/100
  winLength<- round(1/sampleRes)
  
  if (winLength %% 2 == 0) {  #make sure gaussian curve has odd number of weights
    winLength<-winLength+1
  }
  
  weight<-hamming(winLength)
  
  begin<-seq(1,ncol(data)-(length(weight)+1)) #begin location for running avg
  end<-seq(length(weight),ncol(data)) #end location for running avg
  shift<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)
  
  for (j in 1:length(begin)) {
    temp<-data[i,begin[j]:end[j]]
    avg<-sum(temp*weight)/sum(weight)
    GISSgTckLM_O18_PPsmooth[i,shift[j]]<-avg
  }
}

#GISSgTKckLM_O18
data<-GISSgTKckLM_O18_PPmarSed
s<-seq(1,ncol(data))

GISSgTKckLM_O18_PPsmooth<-matrix(data = NA, nrow(data), ncol(data))

for(i in 1:nrow(data)) {
  sampleRes<-(resIndex_O18[i])/100
  winLength<- round(1/sampleRes)
  
  if (winLength %% 2 == 0) {
    winLength<-winLength+1
  }
  
  weight<-hamming(winLength)
  
  begin<-seq(1,ncol(data)-(length(weight)+1))
  end<-seq(length(weight),ncol(data))
  shift<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)
  
  for (j in 1:length(begin)) {
    temp<-data[i,begin[j]:end[j]]
    avg<-sum(temp*weight)/sum(weight)
    GISSgTKckLM_O18_PPsmooth[i,shift[j]]<-avg
  }
}

#GISSgTcsLM_O18
data<-GISSgTcsLM_O18_PPmarSed
s<-seq(1,ncol(data))

GISSgTcsLM_O18_PPsmooth<-matrix(data = NA, nrow(data), ncol(data))

for(i in 1:nrow(data)) {
  sampleRes<-(resIndex_O18[i])/100
  winLength<- round(1/sampleRes)
  
  if (winLength %% 2 == 0) {
    winLength<-winLength+1
  }
  
  weight<-hamming(winLength)
  
  begin<-seq(1,ncol(data)-(length(weight)+1))
  end<-seq(length(weight),ncol(data))
  shift<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)
  
  for (j in 1:length(begin)) {
    temp<-data[i,begin[j]:end[j]]
    avg<-sum(temp*weight)/sum(weight)
    GISSgTcsLM_O18_PPsmooth[i,shift[j]]<-avg
  }
}

#GISSgTckLM_MgCa
data<-GISSgTckLM_MgCa
s<-seq(1,ncol(data))

GISSgTckLM_MgCa_PPsmooth<-matrix(data = NA, nrow(data), ncol(data))

for(i in 1:nrow(data)) {
  sampleRes<-(resIndex_MgCa[i])/100
  winLength<- round(1/sampleRes)
  
  if (winLength %% 2 == 0) {
    winLength<-winLength+1
  }
  
  weight<-hamming(winLength)
  
  begin<-seq(1,ncol(data)-(length(weight)+1))
  end<-seq(length(weight),ncol(data))
  shift<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)
  
  for (j in 1:length(begin)) {
    temp<-data[i,begin[j]:end[j]]
    avg<-sum(temp*weight)/sum(weight)
    GISSgTckLM_MgCa_PPsmooth[i,shift[j]]<-avg
  }
}

#GISSgTKckLM_MgCa
data<-GISSgTKckLM_MgCa
s<-seq(1,ncol(data))

GISSgTKckLM_MgCa_PPsmooth<-matrix(data = NA, nrow(data), ncol(data))

for(i in 1:nrow(data)) {
  sampleRes<-(resIndex_MgCa[i])/100
  winLength<- round(1/sampleRes)
  
  if (winLength %% 2 == 0) {
    winLength<-winLength+1
  }
  
  weight<-hamming(winLength)
  
  begin<-seq(1,ncol(data)-(length(weight)+1))
  end<-seq(length(weight),ncol(data))
  shift<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)
  
  for (j in 1:length(begin)) {
    temp<-data[i,begin[j]:end[j]]
    avg<-sum(temp*weight)/sum(weight)
    GISSgTKckLM_MgCa_PPsmooth[i,shift[j]]<-avg
  }
}

#GISSgTcsLM_MgCa
data<-GISSgTcsLM_MgCa
s<-seq(1,ncol(data))

GISSgTcsLM_MgCa_PPsmooth<-matrix(data = NA, nrow(data), ncol(data))

for(i in 1:nrow(data)) {
  sampleRes<-(resIndex_MgCa[i])/100
  winLength<- round(1/sampleRes)
  
  if (winLength %% 2 == 0) {
    winLength<-winLength+1
  }
  
  weight<-hamming(winLength)
  
  begin<-seq(1,ncol(data)-(length(weight)+1))
  end<-seq(length(weight),ncol(data))
  shift<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)
  
  for (j in 1:length(begin)) {
    temp<-data[i,begin[j]:end[j]]
    avg<-sum(temp*weight)/sum(weight)
    GISSgTcsLM_MgCa_PPsmooth[i,shift[j]]<-avg
  }
}

#save
setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/pseudoproxy/smooth/")
write.csv(GISSgTckLM_O18_PPsmooth, file = "GISSgTckLM_O18_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgTckLM_MgCa_PPsmooth, file = "GISSgTckLM_MgCa_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgTKckLM_O18_PPsmooth, file = "GISSgTKckLM_O18_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgTKckLM_MgCa_PPsmooth, file = "GISSgTKckLM_MgCa_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgTcsLM_O18_PPsmooth, file = "GISSgTcsLM_O18_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgTcsLM_MgCa_PPsmooth, file = "GISSgTcsLM_MgCa_PPsmooth.csv", row.names = FALSE)

#plot
x<-GISSgTckLM_O18_PPsmooth
plot(s,GISSgTckLM_O18_PPmarSed[40,], type = "l")
#title(main = "GISSgCpiC_O18 pseudoproxy - smooth")
lines(GISSgTckLM_O18_PPsmooth[40,],col="red")
