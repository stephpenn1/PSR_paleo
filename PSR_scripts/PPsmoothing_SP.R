#Stephanie Pennington | PSR summer research
#Pseudoproxy timeseries smoothing
#apply weighed moving average with hanning filter
#Created 7-6-17

setwd("/home/spenn1/PSR_paleo/PSR_data/pseudoproxy/");


# pull out only marine sediment archive type ------------------------------

sedIndex<- which(metadataO18$archiveType == "marine sediment" | metadataO18$archiveType == "marine sediments")

#O18
HadpiC_O18_PPmarSed<-HadpiC_O18_PP[sedIndex,]
GISSgCpiC_O18_PPmarSed<-GISSgCpiC_O18_PP[sedIndex,]
GISSgy3piC_O18_PPmarSed<-GISSgy3piC_O18_PP[sedIndex,]
GISSgTckLM_O18_PPmarSed<-GISSgTckLM_O18_PP[sedIndex,]
GISSgTKckLM_O18_PPmarSed<-GISSgTKckLM_O18_PP[sedIndex,]
GISSgTcsLM_O18_PPmarSed<-GISSgTcsLM_O18_PP[sedIndex,]

#MgCa
sedIndex_MgCa<-sedIndex[1:28]
HadpiC_MgCa_PPmarSed<-HadpiC_MgCa_PP[sedIndex_MgCa,]
GISSgCpiC_MgCa_PPmarSed<-GISSgCpiC_MgCa_PP[sedIndex_MgCa,]
GISSgy3piC_MgCa_PPmarSed<-GISSgy3piC_MgCa_PP[sedIndex_MgCa,]
GISSgTckLM_MgCa_PPmarSed<-GISSgTckLM_MgCa_PP[sedIndex_MgCa,]
GISSgTKckLM_MgCa_PPmarSed<-GISSgTKckLM_MgCa_PP[sedIndex_MgCa,]
GISSgTcsLM_MgCa_PPmarSed<-GISSgTcsLM_MgCa_PP[sedIndex_MgCa,]

#*note: O18 has 40 rows and MgCa has 28 rows*

setwd("/home/spenn1/PSR_paleo/PSR_data/pseudoproxy/marine sediments");

#save as CSV for smoothing
write.csv(HadpiC_O18_PPmarSed, file = "HadpiC_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgCpiC_O18_PPmarSed, file = "GISSgCpiC_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgy3piC_O18_PPmarSed, file = "GISSgy3piC_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTckLM_O18_PPmarSed, file = "GISSgTckLM_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTKckLM_O18_PPmarSed, file = "GISSgTKckLM_O18_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTcsLM_O18_PPmarSed, file = "GISSgTcsLM_O18_PPmarSed.csv", row.names = FALSE)
write.csv(HadpiC_MgCa_PPmarSed, file = "HadpiC_MgCa_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgCpiC_MgCa_PPmarSed, file = "GISSgCpiC_MgCa_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgy3piC_MgCa_PPmarSed, file = "GISSgy3piC_MgCa_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTckLM_MgCa_PPmarSed, file = "GISSgTckLM_MgCa_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTKckLM_MgCa_PPmarSed, file = "GISSgTKckLM_MgCa_PPmarSed.csv", row.names = FALSE)
write.csv(GISSgTcsLM_MgCa_PPmarSed, file = "GISSgTcsLM_MgCa_PPmarSed.csv", row.names = FALSE)

# calculate weighted moving average ---------------------------------------

setwd("/home/spenn1/PSR_paleo/PSR_data/pseudoproxy/marine sediments");
library(signal)

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

#HadpiC_O18_PPmarSed

data<-HadpiC_O18_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,39960)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
HadpiC_O18_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  HadpiC_O18_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}


#HadpiC_MgCa_PPmarSed

data<-HadpiC_MgCa_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,27972)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
HadpiC_MgCa_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  HadpiC_MgCa_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}


#GISSgCpiC_O18_PPmarSed

data<-GISSgCpiC_O18_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,44000)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
GISSgCpiC_O18_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  GISSgCpiC_O18_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}


#GISSgCpiC_MgCa_PPmarSed

data<-GISSgCpiC_MgCa_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,30800)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
GISSgCpiC_MgCa_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  GISSgCpiC_MgCa_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}


#GISSgy3piC_O18_PPmarSed

data<-GISSgy3piC_O18_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,6400)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
GISSgy3piC_O18_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  GISSgy3piC_O18_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}


#GISSgy3piC_MgCa_PPmarSed

data<-GISSgy3piC_MgCa_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,4448)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
GISSgy3piC_MgCa_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  GISSgy3piC_MgCa_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}


#GISSgTckLM_O18_PPmarSed

data<-GISSgTckLM_O18_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,44000)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
GISSgTckLM_O18_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  GISSgTckLM_O18_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}


#GISSgTckLM_MgCa_PPmarSed

data<-GISSgTckLM_MgCa_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,30800)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
GISSgTckLM_MgCa_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  GISSgTckLM_MgCa_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}


#GISSgTKckLM_O18_PPmarSed

data<-GISSgTKckLM_O18_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,44000)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
GISSgTKckLM_O18_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  GISSgTKckLM_O18_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}


#GISSgTKckLM_MgCa_PPmarSed

data<-GISSgTKckLM_MgCa_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,30800)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
GISSgTKckLM_MgCa_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  GISSgTKckLM_MgCa_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}


#GISSgTcsLM_O18_PPmarSed

data<-GISSgTcsLM_O18_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,39960)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
GISSgTcsLM_O18_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  GISSgTcsLM_O18_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}


#GISSgTcsLM_MgCa_PPmarSed

data<-GISSgTcsLM_MgCa_PPmarSed    #pseudoproxy timeseries to be smoothed

vectorData<-as.vector(t(data))     #to make matrix into vector

s<-seq(1,27972)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(vectorData)-(length(weight)+1))
end<-seq(length(weight),length(vectorData))
GISSgTcsLM_MgCa_PPsmooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-vectorData[begin[i]:end[i]]
  GISSgTcsLM_MgCa_PPsmooth[i]<-sum(temp*weight)/sum(weight)
}

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

plot(s,vectorData, type = "l")
#title(main = "GISSgCpiC_O18 pseudoproxy - smooth")
lines(s[((length(weight)+1)/2):(length(s)-(length(weight)-1)/2)],x,col="red")

#use model_XX_PPmarSed for vectorData