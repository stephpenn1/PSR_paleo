#Stephanie Pennington | PSR summer research
#Pseudoproxy construction
#Created 6-16-17

setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/model_data/");

library(RColorBrewer)
library(hexbin)

#set seed
set.seed(54)

#target SNR
SNR<-0.5

#note: variance calculation works with matrix not data frame
for (i in 1:nrow(HadpiC_O18)) {
    R<-nrow(HadpiC_O18)
    C<-ncol(HadpiC_O18)
    varNreq<-matrix(var(HadpiC_O18[i,])/(SNR^2),R,C)  #find variance of noise needed to reach target SNR
    N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)  #create matrix of noise to match model dims
    SNRcheck<-matrix(sqrt(var(HadpiC_O18[i,])/var(N[i,])),R,C)  #verify SNR value
}
HadpiC_O18_PP<-HadpiC_O18+N  #add noise to model data

#high density scatter plot
#bin<-hexbin(HadpiC_O18_PP, HadpiC_O18, xbins=40)
#my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
#plot(bin, main="HadpiC_018" , colramp=my_colors , legend=F ) 

#scatter with regression/lowess
#plot(HadpiC_O18_PP,HadpiC_O18, xlab = "HadpiC_O18_PP", ylab = "HadpiC_O18") 
#title(main = "HadpiC_O18")
#abline(lm(HadpiC_O18~HadpiC_O18_PP, na.action = NULL), col="red")
#lines(lowess(HadpiC_O18_PP,HadpiC_O18), col="blue")

#plot
#plot(HadpiC_O18_PP[,1],HadpiC_O18[,1])
#plot(seq(1:length(X)),X)
#points(seq(1:length(X)),HadpiC_O18_PP,col="red")
#title(main = "HadpiC_018")

for (i in 1:nrow(HadpiC_MgCa)) {
  R<-nrow(HadpiC_MgCa)
  C<-ncol(HadpiC_MgCa)
  varNreq<-matrix(var(HadpiC_MgCa[i,])/(SNR^2),R,C)
  N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
  SNRcheck<-matrix(sqrt(var(HadpiC_MgCa[i,])/var(N[i,])),R,C)
}
HadpiC_MgCa_PP<-HadpiC_MgCa + N

for (i in 1:nrow(GISSgCpiC_O18)) {
  R<-nrow(GISSgCpiC_O18)
  C<-ncol(GISSgCpiC_O18)
  varNreq<-matrix(var(GISSgCpiC_O18[i,])/(SNR^2),R,C)
  N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
  SNRcheck<-matrix(sqrt(var(GISSgCpiC_O18[i,])/var(N[i,])),R,C)
}
GISSgCpiC_O18_PP<-GISSgCpiC_O18 + N

for (i in 1:nrow(GISSgCpiC_MgCa)) {
  R<-nrow(GISSgCpiC_MgCa)
  C<-ncol(GISSgCpiC_MgCa)
  varNreq<-matrix(var(GISSgCpiC_MgCa[i,])/(SNR^2),R,C)
  N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
  SNRcheck<-matrix(sqrt(var(GISSgCpiC_MgCa[i,])/var(N[i,])),R,C)
}
GISSgCpiC_MgCa_PP<-GISSgCpiC_MgCa + N

for (i in 1:nrow(GISSgy3piC_O18)) {
  R<-nrow(GISSgy3piC_O18)
  C<-ncol(GISSgy3piC_O18)
  varNreq<-matrix(var(GISSgy3piC_O18[i,])/(SNR^2),R,C)
  N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
  SNRcheck<-matrix(sqrt(var(GISSgy3piC_O18[i,])/var(N[i,])),R,C)
}
GISSgy3piC_O18_PP<-GISSgy3piC_O18 + N

for (i in 1:nrow(GISSgy3piC_MgCa)) {
  R<-nrow(GISSgy3piC_MgCa)
  C<-ncol(GISSgy3piC_MgCa)
  varNreq<-matrix(var(GISSgy3piC_MgCa[i,])/(SNR^2),R,C)
  N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
  SNRcheck<-matrix(sqrt(var(GISSgy3piC_MgCa[i,])/var(N[i,])),R,C)
}
GISSgy3piC_MgCa_PP<-GISSgy3piC_MgCa + N

for (i in 1:nrow(GISSgTckLM_O18)) {
  R<-nrow(GISSgTckLM_O18)
  C<-ncol(GISSgTckLM_O18)
  varNreq<-matrix(var(GISSgTckLM_O18[i,])/(SNR^2),R,C)
  N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
  SNRcheck<-matrix(sqrt(var(GISSgTckLM_O18[i,])/var(N[i,])),R,C)
}
GISSgTckLM_O18_PP<-GISSgTckLM_O18 + N

for (i in 1:nrow(GISSgTckLM_MgCa)) {
  R<-nrow(GISSgTckLM_MgCa)
  C<-ncol(GISSgTckLM_MgCa)
  varNreq<-matrix(var(GISSgTckLM_MgCa[i,])/(SNR^2),R,C)
  N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
  SNRcheck<-matrix(sqrt(var(GISSgTckLM_MgCa[i,])/var(N[i,])),R,C)
}
GISSgTckLM_MgCa_PP<-GISSgTckLM_MgCa + N

for (i in 1:nrow(GISSgTKckLM_O18)) {
  R<-nrow(GISSgTKckLM_O18)
  C<-ncol(GISSgTKckLM_O18)
  varNreq<-matrix(var(GISSgTKckLM_O18[i,])/(SNR^2),R,C)
  N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
  SNRcheck<-matrix(sqrt(var(GISSgTKckLM_O18[i,])/var(N[i,])),R,C)
}
GISSgTKckLM_O18_PP<-GISSgTKckLM_O18 + N

for (i in 1:nrow(GISSgTKckLM_MgCa)) {
  R<-nrow(GISSgTKckLM_MgCa)
  C<-ncol(GISSgTKckLM_MgCa)
  varNreq<-matrix(var(GISSgTKckLM_MgCa[i,])/(SNR^2),R,C)
  N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
  SNRcheck<-matrix(sqrt(var(GISSgTKckLM_MgCa[i,])/var(N[i,])),R,C)
}
GISSgTKckLM_MgCa_PP<-GISSgTKckLM_MgCa + N

for (i in 1:nrow(GISSgTcsLM_O18)) {
  R<-nrow(GISSgTcsLM_O18)
  C<-ncol(GISSgTcsLM_O18)
  varNreq<-matrix(var(GISSgTcsLM_O18[i,])/(SNR^2),R,C)
  N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
  SNRcheck<-matrix(sqrt(var(GISSgTcsLM_O18[i,])/var(N[i,])),R,C)
}
GISSgTcsLM_O18_PP<-GISSgTcsLM_O18 + N

for (i in 1:nrow(GISSgTcsLM_MgCa)) {
  R<-nrow(GISSgTcsLM_MgCa)
  C<-ncol(GISSgTcsLM_MgCa)
  varNreq<-matrix(var(GISSgTcsLM_MgCa[i,])/(SNR^2),R,C)
  N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
  SNRcheck<-matrix(sqrt(var(GISSgTcsLM_MgCa[i,])/var(N[i,])),R,C)
}
GISSgTcsLM_MgCa_PP<-GISSgTcsLM_MgCa + N

#save as CSV
setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/noise/")
write.csv(HadpiC_O18_PP, file = "HadpiC_O18_PP.csv", row.names = FALSE)
write.csv(HadpiC_MgCa_PP, file = "HadpiC_MgCa_PP.csv", row.names = FALSE)
write.csv(GISSgCpiC_O18_PP, file = "GISSgCpiC_O18_PP.csv", row.names = FALSE)
write.csv(GISSgCpiC_MgCa_PP, file = "GISSgCpiC_MgCa_PP.csv", row.names = FALSE)
write.csv(GISSgy3piC_O18_PP, file = "GISSgy3piC_O18_PP.csv", row.names = FALSE)
write.csv(GISSgy3piC_MgCa_PP, file = "GISSgy3piC_MgCa_PP.csv", row.names = FALSE)
write.csv(GISSgTckLM_O18_PP, file = "GISSgTckLM_O18_PP.csv", row.names = FALSE)
write.csv(GISSgTckLM_MgCa_PP, file = "GISSgTckLM_MgCa_PP.csv", row.names = FALSE)
write.csv(GISSgTKckLM_O18_PP, file = "GISSgTKckLM_O18_PP.csv", row.names = FALSE)
write.csv(GISSgTKckLM_MgCa_PP, file = "GISSgTKckLM_MgCa_PP.csv", row.names = FALSE)
write.csv(GISSgTcsLM_O18_PP, file = "GISSgTcsLM_O18_PP.csv", row.names = FALSE)
write.csv(GISSgTcsLM_MgCa_PP, file = "GISSgTcsLM_MgCa_PP.csv", row.names = FALSE)

