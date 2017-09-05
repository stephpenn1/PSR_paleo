setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/aggregate/");

data<-read.csv("HadpiC_O18_PPagg.csv")
#data<-HadpiC_MgCa_PPsmooth
time<-seq(1,ncol(data))
data.avg<-data.frame(matrix(data = NA, nrow(data),ncol(data))) #create empty matrix for time and dO18 data
data.time<-data.frame(matrix(data = NA, nrow(data),ncol(data)))
sampleRes<-metadataO18[,4]
#sampleRes<-resIndex_MgCa
proxy.bounds<-matrix(data = NA,44,2) #create matrix with min and max of proxy data time
colnames(proxy.bounds)<-c("min", "max")

for (i in 1:max(O18proxy$index)) {
  x<-which(O18proxy$index == i)
  proxy.bounds[i,1]<-min(O18proxy$year[x]) #find min and max for each site
  proxy.bounds[i,2]<-max(O18proxy$year[x])
}

for (i in 1:nrow(data)) {
  if (sampleRes[i] > 100) {
    sampleRes[i] <- 100
  } else {
    s<-sampleRes[i]/100    
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
  
  for (j in 1:length(binWidth)) {
    data.avg[i,j]<-sum(data[i,binStart[j]:binEnd[j]], na.rm = TRUE)/binWidth[j]
    data.time[i,j]<-sum(time[binStart[j]:binEnd[j]], na.rm = TRUE)/binWidth[j]
    #if(data.time[i,j] < proxy.bounds[i,1] | data.time[i,j] > proxy.bounds[i,2]) {
    #  data.time[i,j]<-NA
    #  data.avg[i,j]<-NA
    #}
  }
}

#create index corresponding to each bin size
index<-matrix(data = NA,nrow(data),ncol(data))
for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    if (is.na(data.avg[i,j]) == FALSE) {
      index[i,j] = i
    } 
  }
}

#transpose matrix to vector
HadpiC_O18_PPavg<-na.omit(as.vector(t(data.avg)))
HadpiC_O18_PPtime<-na.omit(as.vector(t(data.time)))
HadpiC_O18_PPindex<-na.omit(as.vector(t(index)))

#create data frame
index<-HadpiC_O18_PPindex
year<-HadpiC_O18_PPtime
dO18<-signif(HadpiC_O18_PPavg, digits = 2)
test<-cbind(index,year,dO18)
colnames(test)<- c("index", "year", "dO18")
write.csv(test, file = "test.csv", row.names = FALSE)