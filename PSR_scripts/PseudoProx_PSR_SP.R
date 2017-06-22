#Stephanie Pennington | PSR summer research
#Pseudoproxy construction
#Created 6-16-17

#target SNR
SNR<-0.5

##
X<-HadpiC_O18
#find variance of noise needed to reach target SNR
varNreq<-var(X)/(SNR^2)
#create matrix of noise to match model dims
R<-43
C<-999
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
#verify SNR value
SNRcheck<-sqrt(var(X)/var(N))
#add noise to model data
HadpiC_O18_PP<-X+N

plot(seq(1:length(X)),X)
points(seq(1:length(X)),HadpiC_O18_PP,col="red")
title(main = "HadpiC_018")

X<-HadpiC_MgCa
varNreq<-var(X)/(SNR^2)
R<-31
C<-999
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
SNRcheck<-sqrt(var(X)/var(N))
HadpiC_MgCa_PP<-X+N

X<-GISSgCpiC_O18
varNreq<-var(X)/(SNR^2)
R<-43
C<-1100
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
SNRcheck<-sqrt(var(X)/var(N))
GISSgCpiC_O18_PP<-X+N

X<-GISSgCpiC_MgCa
varNreq<-var(X)/(SNR^2)
R<-31
C<-1100
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
SNRcheck<-sqrt(var(X)/var(N))
GISSgCpiC_MgCa_PP<-X+N

X<-GISSgy3piC_O18
varNreq<-var(X)/(SNR^2)
R<-43
C<-160
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
SNRcheck<-sqrt(var(X)/var(N))
GISSgy3piC_O18_PP<-X+N

X<-GISSgy3piC_MgCa
varNreq<-var(X)/(SNR^2)
R<-31
C<-160
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
SNRcheck<-sqrt(var(X)/var(N))
GISSgy3piC_MgCa_PP<-X+N

X<-GISSgTckLM_O18
varNreq<-var(X)/(SNR^2)
R<-44
C<-1100
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
SNRcheck<-sqrt(var(X)/var(N))
GISSgTckLM_O18_PP<-X+N

X<-GISSgTckLM_MgCa
varNreq<-var(X)/(SNR^2)
R<-31
C<-1100
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
SNRcheck<-sqrt(var(X)/var(N))
GISSgTckLM_MgCa_PP<-X+N

X<-GISSgTKckLM_O18
varNreq<-var(X)/(SNR^2)
R<-44
C<-1100
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
SNRcheck<-sqrt(var(X)/var(N))
GISSgTKckLM_O18_PP<-X+N

X<-GISSgTKckLM_MgCa
varNreq<-var(X)/(SNR^2)
R<-31
C<-1100
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
SNRcheck<-sqrt(var(X)/var(N))
GISSgTKckLM_MgCa_PP<-X+N

X<-GISSgTcsLM_O18
varNreq<-var(X)/(SNR^2)
R<-44
C<-999
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
GISSgTcsLM_O18_PP<-X+N

X<-GISSgTcsLM_MgCa
varNreq<-var(X)/(SNR^2)
R<-31
C<-999
N<-matrix(rnorm(R*C,0,sqrt(varNreq)),R,C)
GISSgTcsLM_MgCa_PP<-X+N