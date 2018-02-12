
setwd("~/Desktop/PSR_paleo/PSR_data/pseudoproxy/model_data1/")
GISSgTckLM_MgCa <- as.matrix(read.csv("GISSgTckLM_MgCa.csv"))

setwd("~/Desktop/PSR_paleo/PSR_data/pseudoproxy/smooth2b/")
GISSgTckLM_MgCa_PPsmooth <- as.matrix(read.csv("GISSgTckLM_MgCa_PPsmooth.csv"))


setwd("~/Desktop/PSR_paleo/PSR_data/pseudoproxy/average3/")
GISSgTckLM_MgCa_PPavg <- read.csv("GISSgTckLM_MgCa_PPavg.csv")


time <- seq(850, 849+1100)
plot(time, GISSgTckLM_MgCa[2,], type = "l")
lines(time, GISSgTckLM_MgCa_PPsmooth[2,], col = "red")

GISSgTckLM_MgCa_beforeSNR <- data.avg
plot(data.time[2,], GISSgTckLM_MgCa_beforeSNR[2,], col = "blue")

m <- which(GISSgTckLM_MgCa_PPavg[,1] == 2)
lines(GISSgTckLM_MgCa_PPavg[,2][m], GISSgTckLM_MgCa_PPavg[,3][m], col = "purple")
