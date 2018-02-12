
x <- which(SSTproxy$index == 26)
plot(SSTproxy$year[x], SSTproxy$Tproxy.val[x], type = "l", ylim = c(0,8), main = "GISSgTckLM_MgCa - index 26")

w <- which(GISSgTckLM_MgCa_PPavg[,1] == 17)
lines(GISSgTckLM_MgCa_PPavg[,2][w], GISSgTckLM_MgCa_PPavg[,3][w], col = "purple")
points(GISSgTckLM_MgCa_PPavg[,2][w], GISSgTckLM_MgCa_PPavg[,3][w], col = "red", type = "l")
legend("topleft", 1, legend = c("proxy","signal", "SNR 2"), col = c("black", "red", "purple"), lty=1, cex=0.5)

y <- which(O18proxy$index == 1)
plot(O18proxy$year[y], O18proxy$d18O[y], type = "l", ylim = c(0,4))

z <- which(GISSgTckLM_O18_PPavg[,1] == 1)
points(GISSgTckLM_O18_PPavg[,2][z], GISSgTckLM_O18_PPavg[,3][z], type = "l", col = "purple")
lines(GISSgTckLM_O18_PPavg[,2][z], GISSgTckLM_O18_PPavg[,3][z], col = "red")

setwd("~/Desktop/PSR_paleo/PSR_data/pseudoproxy/model_data1/")
GISSgTckLM_MgCa <- read.csv("GISSgTckLm_MgCa.csv")
time <- seq(850, 849 + ncol(GISSgTckLM_MgCa), 1)
plot(time,GISSgTckLM_MgCa[1,], type = "l", xlim = c(0,2000),ylim = c(0,5), col = "red")
par(new = TRUE)
x <- which(SSTproxy$index == 1)
plot(SSTproxy$year[x], SSTproxy$Tproxy.val[x], type = "l",col = "black", xlim = c(0,2000),ylim = c(0,5))
