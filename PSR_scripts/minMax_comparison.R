
proxy.boundsMgCa<-matrix(data = NA,31,2)
for (i in 1:length(MgCa_index)) {
  x<-which(SSTproxy$index == MgCa_index[i])
  proxy.boundsMgCa[i,1]<-min(SSTproxy$year[x]) #find min and max for each site from SST proxy data
  proxy.boundsMgCa[i,2]<-max(SSTproxy$year[x])
}

proxy.boundsO18<-matrix(data = NA,44,2)
for (i in 1:44) {
  x<-which(O18proxy$index == i)
  proxy.boundsO18[i,1]<-min(O18proxy$year[x], na.rm = TRUE) #find min and max for each site
  proxy.boundsO18[i,2]<-max(O18proxy$year[x], na.rm = TRUE)
}

GISSgTckLM_MgCa_PPavg <- read.csv("GISSgTckLM_MgCa_PPavg.csv")
pp.boundsMgCa <- matrix(data = NA, 31, 2)
data <- GISSgTckLM_MgCa_PPavg
for (i in 1:31) {
  x<-which(data$index == i)
  pp.boundsMgCa[i,1]<-min(data$year[x]) #find min and max for each site from SST proxy data
  pp.boundsMgCa[i,2]<-max(data$year[x])
}

GISSgTckLM_O18_PPavg <- read.csv("GISSgTckLM_O18_PPavg.csv")
pp.boundsO18 <- matrix(data = NA, 44, 2)
data <- GISSgTckLM_O18_PPavg
for (i in 1:44) {
  x<-which(data$index == i)
  pp.boundsO18[i,1]<-min(data$year[x]) #find min and max for each site from SST proxy data
  pp.boundsO18[i,2]<-max(data$year[x])
}
