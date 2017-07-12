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

data<-HadpiC_O18_PPmarSed     #pseudoproxy timeseries to be smoothed

s<-seq(1,39960)
sampleRes<-10.85/100
winLength<- round(1/sampleRes)
weight<-hamming(winLength)

begin<-seq(1,length(testDat)-(length(w)+1))
end<-seq(length(weight),length(testDat))
smooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)

for (i in 1:length(begin)) {
  temp<-testDat[begin[i]:end[i]]
  smooth[i]<-sum(temp*weight)/sum(weight)
}

plot(s,testDat, type = "l")
lines(s[((length(weight)+1)/2):(length(s)-(length(weight)-1)/2)],smooth,col="red")

#to make matrix into vector
#vectorData<-as.vector(t(data))