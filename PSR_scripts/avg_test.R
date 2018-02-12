SNR <- 0.5
#functions to average MgCa and O18 data
#average_MgCa<-function(data) {
setwd("~/Desktop/PSR_paleo/PSR_data/pseudoproxy/smooth2b/")
GISSgTckLM_MgCa_PPsmooth <- read.csv("GISSgTckLM_MgCa_PPsmooth.csv")
setwd("~/Desktop/PSR_paleo/PSR_data/pseudoproxy/model_data1/")
GISSgTckLM_MgCa <- read.csv("GISSgTckLM_MgCa.csv")
data <- GISSgTckLM_MgCa_PPsmooth
  time<-seq(850,ncol(data) + 850)
  data.avg<-data.frame(matrix(data = NA, nrow(data),ncol(data))) #create empty matrix for time and dO18 data
  data.time<-data.frame(matrix(data = NA, nrow(data),ncol(data)))
  proxy.bounds<-matrix(data = NA,nrow(data),2) #create matrix with min and max of proxy data time
  colnames(proxy.bounds)<-c("min", "max")
  
  for (i in 1:length(MgCa_index)) {
    x<-which(SSTproxy$index == MgCa_index[i])
    proxy.bounds[i,1]<-min(SSTproxy$year[x]) #find min and max for each site from SST proxy data
    proxy.bounds[i,2]<-max(SSTproxy$year[x])
  }
  for (i in 1:nrow(data)) {
    if (sampleRes_MgCa[i] > 70) { #set near annual resolution to 100 for simplicity
      sampleRes_MgCa[i] <- 100 #100 samples/century
    } else {
      s<-sampleRes_MgCa[i]/100 
    }
  } 
  
  for (i in 1:nrow(data)) {
    s<-sampleRes_MgCa[i]/100 #samples/century to samples/year
    binTotal<-ceiling(ncol(data) * s) #find total number of bins based on sample res
    
    binAsize<-floor(ncol(data)/binTotal) #determine number of data points in bin A  
    nbinA<-binTotal - (ncol(data) %% binAsize) #determine quantity of first set of bins
    binBsize<-binAsize + 1
    nbinB<-ncol(data) %% binAsize
    
    binA<-rep(binAsize,nbinA)
    binB<-rep(binBsize,nbinB)
    binWidth<-sample(c(binA,binB))
    
    binEnd<-cumsum(binWidth) #bin start location
    binStart<-(binEnd-binWidth)+1 #bin end location
    
    #take average
    for (j in 1:length(binWidth)) {
      data.avg[i,j]<-sum(data[i,binStart[j]:binEnd[j]], na.rm = TRUE)/binWidth[j]
      data.time[i,j]<-sum(time[binStart[j]:binEnd[j]], na.rm = TRUE)/binWidth[j]
      if(data.time[i,j] < proxy.bounds[i,1] | data.time[i,j] > proxy.bounds[i,2] | data.time[i,j] < 850 | data.time[i,j] > 1850) { #truncate time
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
  
  #add noise
  data.avg <- data.matrix(data.avg)
  for (i in 1:nrow(data)) {
    R<-nrow(data)
    C<-ncol(data)
    varNreq<-matrix(var(data.avg[i,], na.rm = TRUE)/(SNR^2),R,C)
    N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
    SNRcheck<-matrix(sqrt(var(data.avg[i,], na.rm = TRUE)/var(N[i,])),R,C)
  }
  data.avgNoise <-data.avg + N
  
  #transpose matrix to vector to create data frame
  PPavg<-na.omit(as.vector(t(data.avg)))
  PPtime<-na.omit(as.vector(t(data.time)))
  PPindex<-na.omit(as.vector(t(index)))
  
  index<-PPindex  #create data frame
  year<-round(PPtime)
  Tproxy.val<-signif(PPavg, digits = 2)
  avg<<-cbind(index,year,Tproxy.val)
  colnames(avg)<- c("index", "year", "Tproxy.val")
#  return(avg)
#}

average_MgCa(GISSgTckLM_MgCa_PPsmooth) 
GISSgTckLM_MgCa_PPavg<-avg #rename output file as model simulation
GISSgTckLM_MgCa <- as.matrix(GISSgTckLM_MgCa)

for (i in nrow(GISSgTckLM_MgCa)) {
SNRtest <- matrix(sqrt(var(GISSgTckLM_MgCa[i,])/var(data.avg[i,], na.rm = TRUE)), 31, 1100)
}

