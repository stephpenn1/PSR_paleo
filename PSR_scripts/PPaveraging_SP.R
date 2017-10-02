#Stephanie Pennington | PSR summer research
#Pseudoproxy averaging according to sample resolution
#Created 7-18-17

setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/aggregate/")
GISSgTckLM_O18_PPagg<-read.csv("GISSgTckLM_O18_PPagg.csv")
GISSgTKckLM_O18_PPagg<-read.csv("GISSgTKckLM_O18_PPagg.csv")
GISSgTcsLM_O18_PPagg<-read.csv("GISSgTcsLM_O18_PPagg.csv")

setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/smooth/")
GISSgTckLM_MgCa_PPsmooth<-read.csv("GISSgTckLM_MgCa_PPsmooth.csv")
GISSgTKckLM_MgCa_PPsmooth<-read.csv("GISSgTKckLM_MgCa_PPsmooth.csv")
GISSgTcsLM_MgCa_PPsmooth<-read.csv("GISSgTcsLM_MgCa_PPsmooth.csv")

sampleRes_O18<-metadataO18[,4]
sampleRes_MgCa<-MgCa[,5]

#functions to average MgCa and O18 data
average_MgCa<-function(data) {
  time<-seq(1,ncol(data))
  data.avg<-data.frame(matrix(data = NA, nrow(data),ncol(data))) #create empty matrix for time and dO18 data
  data.time<-data.frame(matrix(data = NA, nrow(data),ncol(data)))
  proxy.bounds<-matrix(data = NA,nrow(data),2) #create matrix with min and max of proxy data time
  colnames(proxy.bounds)<-c("min", "max")
    
  for (i in 1:nrow(data)) {
    x<-which(SSTproxy$index == i)
    proxy.bounds[i,1]<-min(SSTproxy$year[x]) #find min and max for each site
    proxy.bounds[i,2]<-max(SSTproxy$year[x])
  }
    
  for (i in 1:nrow(data)) {
    s<-sampleRes_MgCa[i]/100
    binTotal<-ceiling(ncol(data) * s)
  
    binAsize<-floor(ncol(data)/binTotal) #determine bin quantity and sizes 
    nbinA<-binTotal - (ncol(data) %% binAsize)
    binBsize<-binAsize + 1
    nbinB<-ncol(data) %% binAsize
    
    binA<-rep(binAsize,nbinA)
    binB<-rep(binBsize,nbinB)
    binWidth<-sample(c(binA,binB))
    
    binEnd<-cumsum(binWidth) #bin start location
    binStart<-(binEnd-binWidth)+1 #bin end location
      
    for (j in 1:length(binWidth)) {
      data.avg[i,j]<-sum(data[i,binStart[j]:binEnd[j]], na.rm = TRUE)/binWidth[j]
      data.time[i,j]<-sum(time[binStart[j]:binEnd[j]], na.rm = TRUE)/binWidth[j]
      if(data.time[i,j] < proxy.bounds[i,1] | data.time[i,j] > proxy.bounds[i,2]) { #truncate time
        data.time[i,j]<-NA
        data.avg[i,j]<-NA
      }
    }
  }
  
  index<-matrix(data = NA,nrow(data),ncol(data)) #create index corresponding to each bin size without NAs
  for (i in 1:nrow(data)) {
    for (j in 1:ncol(data)) {
      if (is.na(data.avg[i,j]) == FALSE) { #if number is not NA, add to final matrix
        index[i,j] = i
      } 
    }
  }
  
  #transpose matrix to vector to create data frame
  PPavg<-na.omit(as.vector(t(data.avg)))
  PPtime<-na.omit(as.vector(t(data.time)))
  PPindex<-na.omit(as.vector(t(index)))
  
  index<-PPindex  #create data frame
  year<-PPtime
  Tproxy.val<-signif(PPavg, digits = 2)
  avg<<-cbind(index,year,Tproxy.val)
  colnames(avg)<- c("index", "year", "Tproxy.val")
  return(avg)
}

average_O18<-function(data) {
  time<-seq(1,ncol(data))
  data.avg<-data.frame(matrix(data = NA, nrow(data),ncol(data))) #create empty matrix for time and dO18 data
  data.time<-data.frame(matrix(data = NA, nrow(data),ncol(data)))
  proxy.bounds<-matrix(data = NA,nrow(data),2) #create matrix with min and max of proxy data time
  colnames(proxy.bounds)<-c("min", "max")
  
  for (i in 1:nrow(data)) {
    x<-which(O18proxy$index == i)
    proxy.bounds[i,1]<-min(O18proxy$year[x]) #find min and max for each site
    proxy.bounds[i,2]<-max(O18proxy$year[x])
  }
  
  for (i in 1:nrow(data)) {     #
    s<-sampleRes_O18[i]/100
    binTotal<-ceiling(ncol(data) * s)
    
    binAsize<-floor(ncol(data)/binTotal)    
    nbinA<-binTotal - (ncol(data) %% binAsize)
    binBsize<-binAsize + 1
    nbinB<-ncol(data) %% binAsize
    
    binA<-rep(binAsize,nbinA)
    binB<-rep(binBsize,nbinB)
    binWidth<-sample(c(binA,binB))
    
    binEnd<-cumsum(binWidth)
    binStart<-(binEnd-binWidth)+1
    
    for (j in 1:length(binWidth)) {
      data.avg[i,j]<-sum(data[i,binStart[j]:binEnd[j]], na.rm = TRUE)/binWidth[j]
      data.time[i,j]<-sum(time[binStart[j]:binEnd[j]], na.rm = TRUE)/binWidth[j]
      if(data.time[i,j] < proxy.bounds[i,1] | data.time[i,j] > proxy.bounds[i,2]) {
        data.time[i,j]<-NA
        data.avg[i,j]<-NA
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
  PPavg<-na.omit(as.vector(t(data.avg)))
  PPtime<-na.omit(as.vector(t(data.time)))
  PPindex<-na.omit(as.vector(t(index)))
  
  index<-PPindex  #create data frame
  year<-PPtime
  dO18<-signif(PPavg, digits = 2)
  avg<<-cbind(index,year,dO18)
  colnames(avg)<- c("index", "year", "dO18")
  return(avg)
}

#run functions
#MgCa
average_MgCa(GISSgTckLM_MgCa_PPsmooth) 
GISSgTckLM_MgCa_PPavg<-avg #rename output file as model simulation
average_MgCa(GISSgTKckLM_MgCa_PPsmooth)
GISSgTKckLM_MgCa_PPavg<-avg
average_MgCa(GISSgTcsLM_MgCa_PPsmooth)
GISSgTcsLM_MgCa_PPavg<-avg

#O18
average_O18(GISSgTckLM_O18_PPagg)
GISSgTckLM_O18_PPavg<-avg
average_O18(GISSgTKckLM_O18_PPagg)
GISSgTKckLM_O18_PPavg<-avg
average_O18(GISSgTcsLM_O18_PPagg)
GISSgTcsLM_O18_PPavg<-avg

#save files
setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/average/")
write.csv(GISSgTckLM_O18_PPavg, file = "GISSgTckLM_O18_PPavg.csv", row.names = FALSE)
write.csv(GISSgTckLM_MgCa_PPavg, file = "GISSgTckLM_MgCa_PPavg.csv", row.names = FALSE)
write.csv(GISSgTKckLM_O18_PPavg, file = "GISSgTKckLM_O18_PPavg.csv", row.names = FALSE)
write.csv(GISSgTKckLM_MgCa_PPavg, file = "GISSgTKckLM_MgCa_PPavg.csv", row.names = FALSE)
write.csv(GISSgTcsLM_O18_PPavg, file = "GISSgTcsLM_O18_PPavg.csv", row.names = FALSE)
write.csv(GISSgTcsLM_MgCa_PPavg, file = "GISSgTcsLM_MgCa_PPavg.csv", row.names = FALSE)
