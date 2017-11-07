
#set seed
set.seed(54)

#target SNR
SNR<-2

data <- GISSgTckLM_O18_PPagg
#average_O18<-function(data) {
  time<-seq(850,ncol(data) + 850)
  data.avg<-data.frame(matrix(data = NA, nrow(data),ncol(data))) #create empty matrix for time and dO18 data
  data.time<-data.frame(matrix(data = NA, nrow(data),ncol(data)))
  proxy.bounds<-matrix(data = NA,nrow(data),2) #create matrix with min and max of proxy data time
  colnames(proxy.bounds)<-c("min", "max")
  
  for (i in 1:nrow(data)) {
    x<-which(O18proxy$index == i)
    proxy.bounds[i,1]<-min(O18proxy$year[x]) #find min and max for each site
    proxy.bounds[i,2]<-max(O18proxy$year[x])
  }
  
  for (i in 1:nrow(data)) {
    if (sampleRes_O18[i] > 100) {
      sampleRes_O18[i] <- 100
    } else {
      s<-sampleRes_O18[i]/100 
    }
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
      if(data.time[i,j] <= proxy.bounds[i,1] | data.time[i,j] >= proxy.bounds[i,2] | data.time[i,j] <= 850 | data.time[i,j] >= 1850) {
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
  
  data.avg <- data.matrix(data.avg)
  for (i in 1:nrow(data)) {
    R<-nrow(data)
    C<-ncol(data)
    varNreq<-matrix(var(data.avg[i,], na.rm = TRUE)/(SNR^2),R,C)
    N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
    SNRcheck<-matrix(sqrt(var(data.avg[i,], na.rm = TRUE)/var(N[i,])),R,C)
  }
  GISSgTckLM_O18_SNR2<-data.avg + N
  
  setwd("/home/spenn1/PSR_paleo/PSR_data/pseudoproxy/SNR/")
  write.csv(GISSgTckLM_O18_SNR2, file = "GISSgTckLM_O18_SNR2.csv", row.names = FALSE)
  
  #transpose matrix to vector
#  PPavg<-na.omit(as.vector(t(data.avg)))
#  PPtime<-na.omit(as.vector(t(data.time)))
#  PPindex<-na.omit(as.vector(t(index)))
  
#  index<-PPindex  #create data frame
#  year<-round(PPtime)
#  dO18<-signif(PPavg, digits = 2)
#  avg<<-cbind(index,year,dO18)
#  colnames(avg)<- c("index", "year", "dO18")
#  return(avg)
#}