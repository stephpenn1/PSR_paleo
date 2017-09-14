#Stephanie Pennington | PSR summer research
#Pseudoproxy binning according to sample resolution
#Created 7-18-17

sampleRes_O18<-metadataO18[,4]
sampleRes_MgCa<-MgCa[,5]

setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/aggregate/")
HadpiC_O18_PPagg<-read.csv("HadpiC_O18_PPagg.csv")
GISSgCpiC_O18_PPagg<-read.csv("GISSgCpiC_O18_PPagg.csv")
GISSgy3piC_O18_PPagg<-read.csv("GISSgy3piC_O18_PPagg.csv")
GISSgTckLM_O18_PPagg<-read.csv("GISSgTckLM_O18_PPagg.csv")
GISSgTKckLM_O18_PPagg<-read.csv("GISSgTKckLM_O18_PPagg.csv")
GISSgTcsLM_O18_PPagg<-read.csv("GISSgTcsLM_O18_PPagg.csv")

setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/smooth/")
HadpiC_MgCa_PPsmooth<-read.csv("HadpiC_MgCa_PPsmooth.csv")
GISSgCpiC_MgCa_PPsmooth<-read.csv("GISSgCpiC_MgCa_PPsmooth.csv")
GISSgy3piC_MgCa_PPsmooth<-read.csv("GISSgy3piC_MgCa_PPsmooth.csv")
GISSgTckLM_MgCa_PPsmooth<-read.csv("GISSgTckLM_MgCa_PPsmooth.csv")
GISSgTKckLM_MgCa_PPsmooth<-read.csv("GISSgTKckLM_MgCa_PPsmooth.csv")
GISSgTcsLM_MgCa_PPsmooth<-read.csv("GISSgTcsLM_MgCa_PPsmooth.csv")

data<-GISSgCpiC_MgCa_PPsmooth
time<-seq(1,ncol(data))
data.avg<-data.frame(matrix(data = NA, nrow(data),ncol(data))) #create empty matrix for time and dO18 data
data.time<-data.frame(matrix(data = NA, nrow(data),ncol(data)))
proxy.bounds<-matrix(data = NA,44,2) #create matrix with min and max of proxy data time
colnames(proxy.bounds)<-c("min", "max")
  
for (i in 1:max(SSTproxy$index)) {
  x<-which(SSTproxy$index == i)
  proxy.bounds[i,1]<-min(SSTproxy$year[x]) #find min and max for each site
  proxy.bounds[i,2]<-max(SSTproxy$year[x])
}
  
for (i in 1:nrow(data)) {
  if (sampleRes_MgCa[i] > 100) {
    sampleRes_MgCa[i] <- 100
  } else {
    s<-sampleRes_MgCa[i]/100    
  }
  
  binTotal<-ceiling(ncol(data) * s)
  
  binAsize<-floor(ncol(data)/binTotal)  
  binAlength<-binTotal - (ncol(data) %% binAsize)
  binBsize<-binAsize + 1
  binBlength<-ncol(data) %% binAsize
  
  binA<-rep(binAsize,binAlength)
  binB<-rep(binBsize,binBlength)
  binWidth<-sample(c(binA,binB))
  
  binEnd<-cumsum(binWidth)
  binStart<-(binEnd-binWidth)+1
    
  ##NEED to fiure out how to average sample size of 100##
  for (j in 1:length(binWidth)) {
    if (sampleRes_MgCa[i] == 100){
      data.avg[i,1]<-sum(data[i,], na.rm = TRUE)/ncol(data)
      data.time[i,1]<-sum(time[], na.rm = TRUE)/ncol(data)
    } else {
    data.avg[i,j]<-sum(data[i,binStart[j]:binEnd[j]], na.rm = TRUE)/binWidth[j]
    data.time[i,j]<-sum(time[binStart[j]:binEnd[j]], na.rm = TRUE)/binWidth[j]
    }
  }
}

index<-matrix(data = NA,nrow(data),ncol(data))  #create index corresponding to each bin size
for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    if (is.na(data.avg[i,j]) == FALSE) {
      index[i,j] = i
    } 
  }
}

#transpose matrix to vector
GISSgCpiC_MgCa_PPavg<-na.omit(as.vector(t(data.avg)))
GISSgCpiC_MgCa_PPtime<-na.omit(as.vector(t(data.time)))
GISSgCpiC_MgCa_PPindex<-na.omit(as.vector(t(index)))

index<-GISSgCpiC_MgCa_PPindex  #create data frame
year<-GISSgCpiC_MgCa_PPtime
dO18<-signif(GISSgCpiC_MgCa_PPavg, digits = 2)
test<-cbind(index,year,dO18)
colnames(test)<- c("index", "year", "dO18")
write.csv(test, file = "test.csv", row.names = FALSE)