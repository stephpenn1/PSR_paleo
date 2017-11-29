
setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/proxy_data/")
O18proxy <- read.csv("O18ProxyDat.csv")

setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/pseudoproxy/SNR/")
SNR2 <- read.csv("GISSgTckLM_O18_SNR2.csv")
SNR25 <- read.csv("GISSgTckLM_O18_SNR25.csv")
SNR5 <- read.csv("GISSgTckLM_O18_SNR5.csv")

#loop through proxy data to find data at first site
proxy_data <- matrix(data = NA, 44, 1050)
for (i in 1:44) {
  x<-which(O18proxy$index == i)
  for (j in 1:length(x)) {
    if (O18proxy$year[x[j]] <= 1850 & O18proxy$year[x[j]] >= 850) {
      proxy_data[i,j] <- O18proxy$d18O[x[j]]
    }
  }
}

plot(seq(1,1050), proxy_data[1,], xlim = c(0,70), ylim = c(-1.5,2),type = "l",
     xlab = "year", ylab = "dO18")
#SNR of 2
points(seq(1,ncol(SNR2)), SNR2[1,], type = "l", xlim = c(0,70), col = "blue")
#SNR of .25
points(seq(1,ncol(SNR25)), SNR25[1,], type = "l", xlim = c(0,70), col = "red")
#SNR of .5
points(seq(1,ncol(SNR5)), SNR5[1,], type = "l", xlim = c(0,70), col = "green")
legend("bottomright", 1, legend=c("Proxy Data", "SNR 0.25", "SNR 0.5", "SNR 2"),
       col=c("black", "red", "green", "blue"), lty = 1, cex = 0.8)

setwd("/Users/sp/Desktop/PSR_paleo/PSR_data/pseudoproxy/SNR/plots")
for (i in 1:44) {
  mypath = file.path(paste("plot_", i, ".pdf", sep=""))
  pdf(mypath)
  if (i != 25 & i != 26 & i != 44) {
    max <- 200
  } else {
    max<-1000
  }
  
  #SNR of .25
  plot(seq(1,ncol(SNR25)), SNR25[i,], type = "l", col = "red", xlim = c(0,max),
       xlab = "year", ylab = "dO18")
  #proxy
  points(seq(1,1050), proxy_data[i,],type = "l", xlab = "year", ylab = "dO18")
  #SNR of .5
  points(seq(1,ncol(SNR5)), SNR5[i,], type = "l", col = "green")
  #SNR of 2
  points(seq(1,ncol(SNR2)), SNR2[i,], type = "l", col = "blue")
  
  legend("bottomright", 1, legend=c("Proxy Data", "SNR 0.25", "SNR 0.5", "SNR 2"),
         col=c("black", "red", "green", "blue"), lty = 1, cex = 0.8)
  dev.off()
  }

plot(seq(1,ncol(SNR25)), SNR25[1,], type = "l", xlim = c(-50,70), col = "red", xlab = "year", ylab = "dO18")
#proxy
points(seq(1,1050), proxy_data[1,],type = "l")
#SNR of 2
points(seq(1,ncol(SNR2)), SNR2[1,], type = "l", xlim = c(0,70), col = "blue")

#SNR of .5
points(seq(1,ncol(SNR5)), SNR5[1,], type = "l", xlim = c(0,70), col = "green")
legend("bottomright", 1, legend=c("Proxy Data", "SNR 0.25", "SNR 0.5", "SNR 2"),
       col=c("black", "red", "green", "blue"), lty = 1, cex = 0.5)
