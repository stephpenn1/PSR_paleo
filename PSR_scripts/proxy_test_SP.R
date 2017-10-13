#Stephanie Pennington | PSR summer research
#Pseudoproxy test
#Created 10-6-17

#pull out MgCa indices from SST proxy dataframe
MgCa_index<-MgCa$index

#pseudo_length<-function(data) {
  psr_pseudo<-matrix(data = 0, 44, 1)
  for (i in 1:44) {
    z<-0
    x<-which(GISSgTKckLM_O18_PPavg[,1] == i)
    for (j in 1:length(x)) {
      if (GISSgTKckLM_O18_PPavg[,2][x[j]] >= 850 & GISSgTKckLM_O18_PPavg[,2][x[j]] <= 1850) {
        z<-z+1
      }
    }
    psr_pseudo[i] <- z
  }
#  psr_pseudo<<-psr_pseudo
#  return(psr_pseudo)
#}

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
pseudo_length(GISSgTckLM_MgCa_PPavg) 
GISSgTckLM_MgCa_PP<-psr_pseudo #rename output file as model simulation
pseudo_length(GISSgTKckLM_MgCa_PPavg)
GISSgTKckLM_MgCa_PPavg<-psr_pseudo
pseudo_length(GISSgTcsLM_MgCa_PPavg)
GISSgTcsLM_MgCa_PPavg<-psr_pseudo

#O18
pseudo_length(GISSgTckLM_O18_PPavg)
GISSgTckLM_O18_PPavg<-psr_pseudo
pseudo_length(GISSgTKckLM_O18_PPavg)
GISSgTKckLM_O18_PPavg<-psr_pseudo
pseudo_length(GISSgTcsLM_O18_PPavg)
GISSgTcsLM_O18_PPavg<-psr_pseudo

plot(psr_SSTproxy[1:20], psr_pseudo[1:20],
     xlab = "SST Proxy Data",
     ylab = "GISSgTckLM_MgCa")
title(main = "GISSgTckLM_MgCa")
