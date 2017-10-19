#Stephanie Pennington | PSR summer research
#Pseudoproxy test according to number of data points at each site/index
#Created 10-6-17

#pull out MgCa indices from SST proxy dataframe
MgCa_index<-MgCa$index

pseudo_MgCalength<-function(data) {
  psr_pseudo<-matrix(data = 0, 31, 1)
  for (i in 1:31) {
    z<-0 #set counter to zero
    x<-which(data[,1] == i) #create list of data with index i
    for (j in 1:length(x)) {
      if (data[,2][x[j]] >= 850 & data[,2][x[j]] <= 1850) { #if data is between 850-1850...
        z<-z+1 #...add to counter
      }
    }
    psr_pseudo[i] <- z
  }
  psr_pseudo<<-psr_pseudo #add variable to global environment
  return(psr_pseudo)
}

pseudo_O18length<-function(data) {
  psr_pseudo<-matrix(data = 0, 44, 1)
  for (i in 1:44) {
    z<-0
    x<-which(data[,1] == i)
    for (j in 1:length(x)) {
      if (data[,2][x[j]] >= 850 & data[,2][x[j]] <= 1850) {
        z<-z+1
      }
    }
    psr_pseudo[i] <- z
  }
  psr_pseudo<<-psr_pseudo
  return(psr_pseudo)
}

psr_O18proxy<-matrix(data = 0, 44, 1)
for (i in 1:44) {
  z<-0
  x<-which(O18proxy$index == i)
  for (j in 1:length(x)) {
    if (O18proxy$year[x[j]] >= 850 & O18proxy$year[x[j]] <= 1850) {
      z<-z+1
    }
  }
  psr_O18proxy[i] <- z
}

psr_SSTproxy<-matrix(data = 0, 31, 1)
for (i in 1:length(MgCa_index)) {
  z<-0
  x<-which(SSTproxy$index == MgCa_index[i])
  for (j in 1:length(x)) {
    if (SSTproxy$year[x[j]] >= 850 & SSTproxy$year[x[j]] <= 1850) {
      z<-z+1
    }
  }
  psr_SSTproxy[i] <- z
}


#run functions
#MgCa
pseudo_MgCalength(GISSgTckLM_MgCa_PPavg) 
GISSgTckLM_MgCa_PP<-psr_pseudo #rename output file as model simulation
pseudo_MgCalength(GISSgTKckLM_MgCa_PPavg)
GISSgTKckLM_MgCa_PPavg<-psr_pseudo
pseudo_MgCalength(GISSgTcsLM_MgCa_PPavg)
GISSgTcsLM_MgCa_PPavg<-psr_pseudo

#O18
pseudo_O18length(GISSgTckLM_O18_PPavg)
GISSgTckLM_O18_PPavg<-psr_pseudo
pseudo_O18length(GISSgTKckLM_O18_PPavg)
GISSgTKckLM_O18_PPavg<-psr_pseudo
pseudo_O18length(GISSgTcsLM_O18_PPavg)
GISSgTcsLM_O18_PPavg<-psr_pseudo

#plot
plot(psr_SSTproxy, psr_pseudo,
     xlab = "no. of data points - SST Proxy Data",
     ylab = "no. of data points - GISSgTckLM_MgCa")
title(main = "GISSgTckLM_MgCa")
abline(0,1, col = "red")

