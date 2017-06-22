SNR=sqrt(var(X)/var(N)), where X is proxy data and N is noise

SNR<-0.5						#target SNR
X<-rnorm(100,0,1) 				#pretend this is a record of proxy data. 
varNreq<-var(X)/sqrt(SNR)		#what does the variance of the noise need to be to give SNR = target
N<-rnorm(100,0,sqrt(varNreq))	#noise
SNRcheck<-sqrt(var(X)/var(N))	#check to see what SNR value
pp<-X+N							#pseudoproxy

plot(seq(1:length(X)),X)
points(seq(1:length(X)),pp,col="red")

