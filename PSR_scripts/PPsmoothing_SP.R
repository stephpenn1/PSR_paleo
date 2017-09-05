#Stephanie Pennington | PSR summer research
#Pseudoproxy timeseries smoothing
#apply weighed moving average with hanning filter
#Needs metadata from Proxybin_anomIND script
#Created 7-6-17

setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/noise")
library(signal)

#pull out only marine sediment archive type ------------------------------

sedIndex_O18<- which(metadataO18$archiveType == "marine sediment" | metadataO18$archiveType == "marine sediments")
res_O18<-metadataO18[,4]      #pull avg sample resolution from O18 metadata
resIndex_O18<-res_O18[sedIndex_O18]      #pull out only marine sediment sample resolutions

resIndex_MgCa<-MgCa[,5]     #pull out sample resolution from SST metadata

HadpiC_O18_PP<-read.csv(file = "HadpiC_O18_PP.csv")
HadpiC_MgCa_PP<-read.csv(file = "HadpiC_MgCa_PP.csv")
GISSgCpiC_O18_PP<-read.csv(file = "GISSgCpiC_O18_PP.csv")
GISSgCpiC_MgCa_PP<-read.csv(file = "GISSgCpiC_MgCa_PP.csv")
GISSgy3piC_O18_PP<-read.csv(file = "GISSgy3piC_O18_PP.csv")
GISSgy3piC_MgCa_PP<-read.csv(file = "GISSgy3piC_MgCa_PP.csv")
GISSgTckLM_O18_PP<-read.csv(file = "GISSgTckLM_O18_PP.csv")
GISSgTckLM_MgCa_PP<-read.csv(file = "GISSgTckLM_MgCa_PP.csv")
GISSgTKckLM_O18_PP<-read.csv(file = "GISSgTKckLM_O18_PP.csv")
GISSgTKckLM_MgCa_PP<-read.csv(file = "GISSgTKckLM_MgCa_PP.csv")
GISSgTcsLM_O18_PP<-read.csv(file = "GISSgTcsLM_O18_PP.csv")
GISSgTcsLM_MgCa_PP<-read.csv(file = "GISSgTcsLM_MgCa_PP.csv")

#pull out rows with marine sediment - O18
HadpiC_O18_PPmarSed<-HadpiC_O18_PP[sedIndex_O18,]
GISSgCpiC_O18_PPmarSed<-GISSgCpiC_O18_PP[sedIndex_O18,]
GISSgy3piC_O18_PPmarSed<-GISSgy3piC_O18_PP[sedIndex_O18,]
GISSgTckLM_O18_PPmarSed<-GISSgTckLM_O18_PP[sedIndex_O18,]
GISSgTKckLM_O18_PPmarSed<-GISSgTKckLM_O18_PP[sedIndex_O18,]
GISSgTcsLM_O18_PPmarSed<-GISSgTcsLM_O18_PP[sedIndex_O18,]

#*note: O18 now has 40 rows and MgCa has 31 rows*

setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/marine_sediments");

#save as CSV for smoothing
write.csv(HadpiC_O18_PPmarSed, file = "HadpiC_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgCpiC_O18_PPmarSed, file = "GISSgCpiC_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgy3piC_O18_PPmarSed, file = "GISSgy3piC_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTckLM_O18_PPmarSed, file = "GISSgTckLM_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTKckLM_O18_PPmarSed, file = "GISSgTKckLM_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTcsLM_O18_PPmarSed, file = "GISSgTcsLM_O18_PPmarSed.csv", row.names = FALSE)
write.csv(HadpiC_MgCa_PP, file = "HadpiC_MgCa_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgCpiC_MgCa_PP, file = "GISSgCpiC_MgCa_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgy3piC_MgCa_PP, file = "GISSgy3piC_MgCa_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTckLM_MgCa_PP, file = "GISSgTckLM_MgCa_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTKckLM_MgCa_PP, file = "GISSgTKckLM_MgCa_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTcsLM_MgCa_PP, file = "GISSgTcsLM_MgCa_PPmarSed.csv", row.names = FALSE)

# calculate weighted moving average ---------------------------------------

HadpiC_O18_PPmarSed<-read.csv("HadpiC_O18_PPmarSed.csv")
GISSgCpiC_O18_PPmarSed<-read.csv("GISSgCpiC_O18_PPmarSed.csv")
GISSgy3piC_O18_PPmarSed<-read.csv("GISSgy3piC_O18_PPmarSed.csv")
GISSgTckLM_O18_PPmarSed<-read.csv("GISSgTckLM_O18_PPmarSed.csv")
GISSgTKckLM_O18_PPmarSed<-read.csv("GISSgTKckLM_O18_PPmarSed.csv")
GISSgTcsLM_O18_PPmarSed<-read.csv("GISSgTcsLM_O18_PPmarSed.csv")
HadpiC_MgCa_PPmarSed<-read.csv("HadpiC_MgCa_PPmarSed.csv")
GISSgCpiC_MgCa_PPmarSed<-read.csv("GISSgCpiC_MgCa_PPmarSed.csv")
GISSgy3piC_MgCa_PPmarSed<-read.csv("GISSgy3piC_MgCa_PPmarSed.csv")
GISSgTckLM_MgCa_PPmarSed<-read.csv("GISSgTckLM_MgCa_PPmarSed.csv")
GISSgTKckLM_MgCa_PPmarSed<-read.csv("GISSgTKckLM_MgCa_PPmarSed.csv")
GISSgTcsLM_MgCa_PPmarSed<-read.csv("GISSgTcsLM_MgCa_PPmarSed.csv")

#HadpiC_O18

data<-HadpiC_O18_PPmarSed
s<-seq(1,ncol(data))            # set length of data

HadpiC_O18_PPsmooth<-matrix(nrow = nrow(data), ncol = ncol(data))   #create empty matrix for smoothed data

for(i in 1:nrow(data)) {
  sampleRes<-(resIndex_O18[i])/100
  winLength<- round(1/sampleRes)
  
  if (winLength %% 2 == 0) {
    winLength<-winLength+1
  }                              #to make sure gaussian curve has odd number of weights
  
  weight<-hamming(winLength)
  begin<-seq(1,ncol(data)-(length(weight)+1))       #begin location for running avg
  end<-seq(length(weight),ncol(data))          #end location for running avg
  shift<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)
  
  for (j in 1:length(begin)) {
    temp<-data[i,begin[j]:end[j]]
    avg<-sum(temp*weight)/sum(weight)
    HadpiC_O18_PPsmooth[i,shift[j]]<-avg
  }
}

#GISSgCpiC_O18

data<-GISSgCpiC_O18_PPmarSed
s<-seq(1,ncol(data))

GISSgCpiC_O18_PPsmooth<-matrix(data = NA, nrow(data), ncol(data))

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
    GISSgCpiC_O18_PPsmooth[i,shift[j]]<-avg
  }
}

#GISSgy3piC_O18

data<-GISSgy3piC_O18_PPmarSed
s<-seq(1,ncol(data))

GISSgy3piC_O18_PPsmooth<-matrix(data = NA, nrow(data), ncol(data))

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
    GISSgy3piC_O18_PPsmooth[i,shift[j]]<-avg
  }
}

#GISSgTckLM_O18

data<-GISSgTckLM_O18_PPmarSed
s<-seq(1,ncol(data))

GISSgTckLM_O18_PPsmooth<-matrix(data = NA, nrow(data), ncol(data))

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

#HadpiC_MgCa

data<-HadpiC_MgCa_PPmarSed
s<-seq(1,ncol(data))

HadpiC_MgCa_PPsmooth<-matrix(data = NA, nrow(data), ncol(data))

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
    HadpiC_MgCa_PPsmooth[i,shift[j]]<-avg
  }
}

#GISSgCpiC_MgCa

data<-GISSgCpiC_MgCa_PPmarSed
s<-seq(1,1100)

GISSgCpiC_MgCa_PPsmooth<-matrix(data = NA, nrow(data), ncol(data))

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
    GISSgCpiC_MgCa_PPsmooth[i,shift[j]]<-avg
  }
}

#GISSgy3piC_MgCa

data<-GISSgy3piC_MgCa_PPmarSed
s<-seq(1,160)

GISSgy3piC_MgCa_PPsmooth<-matrix(data = NA, nrow(data), ncol(data))

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
    GISSgy3piC_MgCa_PPsmooth[i,shift[j]]<-avg
  }
}

#GISSgTckLM_MgCa

data<-GISSgTckLM_MgCa_PPmarSed
s<-seq(1,1100)

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

data<-GISSgTKckLM_MgCa_PPmarSed
s<-seq(1,1100)

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

data<-GISSgTcsLM_MgCa_PPmarSed
s<-seq(1,999)

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
setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/smooth/")
write.csv(HadpiC_O18_PPsmooth, file = "HadpiC_O18_PPsmooth.csv", row.names = FALSE)
write.csv(HadpiC_MgCa_PPsmooth, file = "HadpiC_MgCa_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgCpiC_O18_PPsmooth, file = "GISSgCpiC_O18_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgCpiC_MgCa_PPsmooth, file = "GISSgCpiC_MgCa_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgy3piC_O18_PPsmooth, file = "GISSgy3piC_O18_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgy3piC_MgCa_PPsmooth, file = "GISSgy3piC_MgCa_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgTckLM_O18_PPsmooth, file = "GISSgTckLM_O18_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgTckLM_MgCa_PPsmooth, file = "GISSgTckLM_MgCa_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgTKckLM_O18_PPsmooth, file = "GISSgTKckLM_O18_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgTKckLM_MgCa_PPsmooth, file = "GISSgTKckLM_MgCa_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgTcsLM_O18_PPsmooth, file = "GISSgTcsLM_O18_PPsmooth.csv", row.names = FALSE)
write.csv(GISSgTcsLM_MgCa_PPsmooth, file = "GISSgTcsLM_MgCa_PPsmooth.csv", row.names = FALSE)

#plot
x<-HadpiC_O18_PPsmooth
x<-HadpiC_MgCa_PPsmooth
x<-GISSgCpiC_O18_PPsmooth
x<-GISSgCpiC_MgCa_PPsmooth
x<-GISSgy3piC_O18_PPsmooth
x<-GISSgy3piC_MgCa_PPsmooth
x<-GISSgTckLM_O18_PPsmooth
x<-GISSgTckLM_MgCa_PPsmooth
x<-GISSgTKckLM_O18_PPsmooth
x<-GISSgTKckLM_MgCa_PPsmooth
x<-GISSgTcsLM_O18_PPsmooth
x<-GISSgTcsLM_MgCa_PPsmooth

plot(s,GISSgTcsLM_MgCa_PPmarSed[1,], type = "l")
#title(main = "GISSgCpiC_O18 pseudoproxy - smooth")
lines(GISSgTcsLM_MgCa_PPsmooth[1,],col="red")