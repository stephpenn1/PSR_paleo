#The PSR itself. 
#All proxy data is included. Proxy data anomalies are calculated relative to their whole length.
#Model anomalies are calculated at the site level for the season relevant to each proxy data record
#piControl and LM simulations 
#PSR can accommodate bins of 50 or 100 years

setwd("/Users/Casey_JISAO/Google Drive/R_documents/PSR");
library(abind)
library(plyr)
library(maps);library(mapproj); library(mapdata);

#primary variables
n.it<-100																							#number of times to randomly sample calibration/validation
n<-3																								#number of model bins to keep
prct.cal<-0.66																						#proportion in calibration. Validation will be 1 minus this value
nn<-10 #n.it*0.005																						#number of best RE to keep
sc<-1																								#damping factor for proxy data
minres<-0.9																							#minimum number of samples per 100 years
bin.width<-100																						#width (in years) of bin. As of 3/21/17 only 100 or 50
proxy.type<-c("O18MgCa")																			#flag for proxy data type (this doesn't actually change a variable)		
						

#**********PROXY data, binned at same interval******
O18.proxy<-readRDS(paste0("proxy_data/O18_dt",bin.width/2,"wl",bin.width,"anom"))
O18.proxy[,7:ncol(O18.proxy)]<-O18.proxy[,7:ncol(O18.proxy)]*sc										#scale (if sc=1, no scaling)
O18.proxy<-O18.proxy[which(O18.proxy$resolution>=minres),]											#set min res (if minres=0.5, all are included)							
#O18.proxy<-O18.proxy[which(grepl("paired",O18.proxy$name)==TRUE),]									#MASK THIS if you want all data
#O18.proxy<-O18.proxy[which(O18.proxy$lat>=0 & O18.proxy$lon>=-110 & O18.proxy$lon<=30),]			#northern hemisphere, between 110W and 30E
O18.proxy.sd<-apply(O18.proxy[,7:ncol(O18.proxy)],1,sd,na.rm=TRUE)									#sd of O18 anomaly across all time bins for a given record
O18prxy.count<-colSums(!is.na(O18.proxy))[7:ncol(O18.proxy)]										#number of records through time
prxy.times<-as.numeric(colnames(O18.proxy)[7:ncol(O18.proxy)])

MgCa.proxy<-readRDS(paste0("proxy_data/MgCa_dt",bin.width/2,"wl",bin.width,"anom"))
MgCa.proxy[,7:ncol(MgCa.proxy)]<-MgCa.proxy[,7:ncol(MgCa.proxy)]*sc									#scale (if sc=1, no scaling)
MgCa.proxy<-MgCa.proxy[which(MgCa.proxy$resolution>=minres),]										#set min res (if minres=0.5, all are included)	
#MgCa.proxy<-MgCa.proxy[which(grepl("paired",MgCa.proxy$name)==TRUE),]								#MASK THIS if you want all data
MgCa.proxy.sd<-apply(MgCa.proxy[,7:ncol(MgCa.proxy)],1,sd,na.rm=TRUE)								#sd of MgCa across all time bins for a given record
MgCaprxy.count<-colSums(!is.na(MgCa.proxy))[7:ncol(MgCa.proxy)]										#number of records through time

proxy.all<-rbind(O18.proxy,MgCa.proxy)																#change if you want all or just one geochemical variable
proxy.sd<-c(O18.proxy.sd,MgCa.proxy.sd)
prxy.count<-O18prxy.count+MgCaprxy.count

#piControl simulation catalog and complementary AMOC/NAM 
O18PI<-readRDS(paste0("PSR/piC_O18.",bin.width))
O18PI<-O18PI[which(O18.proxy$resolution>=minres),]													#match min res of proxy data
O18PI.sd<-apply(O18PI,1,sd)																			#sd of simulation catalog at each site
MgCaPI<-readRDS(paste0("PSR/piC_MgCa.",bin.width))												#match min res of proxy data
MgCaPI<-MgCaPI[which(MgCa.proxy$resolution>=minres),]												#sd of simulation catalog at each site
MgCaPI.sd<-apply(MgCaPI,1,sd)
AMOCPI<-readRDS(paste0("PSR/piC_AMOC.",bin.width))													#AMOC
NAMPI<-readRDS(paste0("PSR/piC_NAM.",bin.width))													#NAM
NAOPI<-readRDS(paste0("PSR/piC_NAO.",bin.width))													#NAO

#LM simulation catalog and complementary AMOC/NAM 
O18LM<-readRDS(paste0("PSR/LM_O18.",bin.width))
O18LM<-O18LM[which(O18.proxy$resolution>=minres),]													#match min res of proxy data
O18LM.sd<-apply(O18LM,1,sd)																			#sd of simulation catalog at each site
MgCaLM<-readRDS(paste0("PSR/LM_MgCa.",bin.width))												
MgCaLM<-MgCaLM[which(MgCa.proxy$resolution>=minres),]												#match min res of proxy data
MgCaLM.sd<-apply(MgCaLM,1,sd)																		#sd of simulation catalog at each site
AMOCLM<-readRDS(paste0("PSR/LM_AMOC.",bin.width))													#AMOC
NAMLM<-readRDS(paste0("PSR/LM_NAM.",bin.width))														#NAM
NAOLM<-readRDS(paste0("PSR/LM_NAO.",bin.width))														#NAO

sim.pi<-rbind(O18PI,MgCaPI)
sim.lm<-rbind(O18LM,MgCaLM)
sim.all<-cbind(sim.pi,sim.lm)
sim.sd<-c(O18PI.sd,MgCaPI.sd,O18LM.sd,MgCaLM.sd)
AMOC.all<-c(AMOCPI,AMOCLM)	
NAO.all<-c(NAOPI,NAOLM)
NAM.all<-c(NAOPI,NAOLM)	

#********Initiate some empty matrices/arrays***********
d.site<-array(,dim=c(length(prxy.times),nrow(proxy.all),ncol(sim.all)))							#empty array with dimensions 1) each proxy timestep of past 2k 2) proxy-model distance at each site 3) each model bin
#d.site.med<-matrix(,nrow=nrow(proxy.all),ncol=length(prxy.times))								#empty matrix for median distances in each timestep REMOVE THIS?
d.site.wt<-array(,dim=dim(d.site))																#same as d.site but with weights applied

#**********Data distribution***********
col<-rainbow(101,start=0.45,end=0.9, alpha=0.7)
col1<-rainbow(101,start=0.45,end=0.9)															
dev.new(width=8, height=5)
par(mfrow=c(1,2))
plot(0,type="n",xlim=c(-2,2),ylim=c(-2,2),xlab="d18O proxy",ylab="d18O LM simulation",las=1,main=sc)
abline(0,1)
text(seq(-1,1,length.out=20)[c(1,4,7,10,13,16,19)],rep(-1.5,7),round_any(seq(0,101,length.out=20)[c(1,4,7,10,13,16,19)],10),pos=3,cex=0.8)
text(0,-1.5,"samples/100 yr",pos=1,cex=0.8)
points(seq(-1,1,length.out=101),rep(-1.5,101),pch=15,col=col1)
for (i in 1:nrow(O18.proxy)) {
	proxy.quant<-quantile(O18.proxy[i,7:45],probs=seq(0,1,by=0.1), na.rm=T)			
	sim.quant<-quantile(O18LM[i,],probs=seq(0,1,by=0.1), na.rm=T)
	points(proxy.quant,sim.quant,col=col[round(O18.proxy$resolution)[i]],cex=0.5)
	}
plot(0,type="n",xlim=c(-2,2),ylim=c(-2,2),xlab="MgCa proxy",ylab="MgCa LM simulation",las=1)
abline(0,1)
text(seq(-1,1,length.out=20)[c(1,4,7,10,13,16,19)],rep(-1.5,7),round_any(seq(0,101,length.out=20)[c(1,4,7,10,13,16,19)],10),pos=3,cex=0.8)
text(0,-1.5,"samples/100 yr",pos=1,cex=0.8)
points(seq(-1,1,length.out=101),rep(-1.5,101),pch=15,col=col1)
for (i in 1:nrow(MgCa.proxy)) {
	proxy.quant<-quantile(MgCa.proxy[i,7:45],probs=seq(0,1,by=0.1), na.rm=T)			
	sim.quant<-quantile(MgCaLM[i,],probs=seq(0,1,by=0.1), na.rm=T)
	points(proxy.quant,sim.quant,col=col[round(MgCa.proxy$resolution)[i]],cex=0.5)
	}
	
#*********DISTANCE calculations*******************
for (i in 1:length(prxy.times)) {																		#each i is a timestep
	d.site[i,,]<-(proxy.all[,i+6]-sim.all)^2															#(proxies in timestep i - all model analogs)^2 
	#d.site.sd<-apply(d.site[i,,],1,sd)																	#s.d. of distances in each timestep	WHY?
	d.site.wt[i,,]<-(1/proxy.sd)*d.site[i,,]															#WEIGHTED as 1/s.d. across all time bins for a given proxy record. 
																										#Gives same weight to each record regardless of timestep. 
																										#Helps account for different units (i.e. Mg/Ca vs. d18)O. Similar to Franke																				
	#d.site.wt[i,,]<-(1/sim.sd)*d.site[i,,]																#similar to above, but in terms of model sd. 
	}

#empty array for calibration results
stats.cal<-array(,dim=c(length(prxy.times),n,6,n.it))				#dimensions are timestep, n best surrogates, 7 metrics (index of best n iterations, (weighted) distance, r, p, intercept+SE, and intercept-SE)
cal.site<-array(,dim=c(length(prxy.times),max(prxy.count),n.it))	#dimensions are timestep, calibration proxy data, n iterations		

#empty array for validation results							
stats.val<-array(,dim=c(length(prxy.times),n,7,n.it))				#dimensions are timestep, n best surrogates, 7 metrics (distance, RE, CE, r, p, intercept+SE,intercept-SE), n iterations.										
best.val<-matrix(,nrow=length(prxy.times),ncol=nn)					#matrix to populate with which iterations that have the best RE
val.site<-array(,dim=c(length(prxy.times),max(prxy.count),n.it))	#dimensions are timestep, validation proxy data, n iterations	

#empty array for best PSR stats
PSRbest.stats<-matrix(,nrow=length(prxy.times),ncol=12)				#WHY ncol=12?

#*************FIND best analogs*********************
#i for loop is each timestep, j loop is each n.it iteration

for (i in 1:length(prxy.times)) {																			#each i is a timestep		
	nas<-which(is.na(proxy.all[,i+6]))																		#sites without proxy data in timestep i
	dat<-d.site.wt[i,-nas,]																					#******CHANGE THIS****** to d.site (d.site.wt) if you want unweighted (weighted) distance	
	proxy<-proxy.all[-nas,i+6]																				#skip sites without proxy data in a particular bin
	model<-sim.all[-nas,]																					#skip corresponding model sites
	n.cal<-round(nrow(dat)*prct.cal)																		#number of records in calibration period
	n.val<-nrow(dat)-n.cal																					#number of records in validation period
	
	for(j in 1:n.it) {
		rndm<-sample(seq(1:length(proxy)))																	#random sample from proxy data
		cal.dat<-dat[rndm[1:n.cal],]																		#select the first n.cal for calibration
		eud<-sqrt(colSums(cal.dat))/n.cal																	#euclidian distance per record in calibraiton period
		keep<-which(eud<=sort(eud)[n])[1:n]																	#keep the best n simulations with the shortest distance
		sur<-rowMeans(model[,keep])																			#mean O18 or MgCa anomaly of n kept					
		r<-cor(proxy[rndm[1:n.cal]],sur[rndm[1:n.cal]])														#pearson r	
		pval<-cor.test(proxy[rndm[1:n.cal]],sur[rndm[1:n.cal]])$p.value										#p value.
		int<-summary(lm(sur[rndm[1:n.cal]]~proxy[rndm[1:n.cal]]))$coef[1]									#intercept
		intSE<-summary(lm(sur[rndm[1:n.cal]]~proxy[rndm[1:n.cal]]))$coef[3]									#intercept standard error
		
		stats.cal[i,,1,j]<-keep																				#indicies of the best n
		stats.cal[i,,2,j]<-eud[keep]																		#distance of the best n
		stats.cal[i,,3,j]<-rep(r,n)																			#pearson r of best n
		stats.cal[i,,4,j]<-rep(pval,n)																		#p of best n
		stats.cal[i,,5,j]<-rep(int+intSE,n)																	#intercept+SE of best n
		stats.cal[i,,6,j]<-rep(int-intSE,n)																	#intercept-SE of best n
		cal.site[i,1:n.cal,j]<-rndm[1:n.cal]																#indicies of proxy data in calibration																	
		
		val.dat<-dat[rndm[(n.cal+1):nrow(dat)],keep]														#validation at the best n sites
		stats.val[i,,1,j]<-sqrt(colSums(val.dat))/n.val														#validation distance per record in validation
		num<-sum((proxy[rndm[(n.cal+1):nrow(dat)]]-sur[rndm[(n.cal+1):nrow(dat)]])^2)						#numerator of RE and CE calculation
		den.RE<-sum((proxy[rndm[(n.cal+1):nrow(dat)]]-mean(proxy[rndm[1:n.cal]]))^2)						#RE denom
		stats.val[i,,2,j]<-1-(num/den.RE)																	#RE
		den.CE<-sum((proxy[rndm[(n.cal+1):nrow(dat)]]-mean(proxy[rndm[(n.cal+1):nrow(dat)]]))^2)			#CE denom
		stats.val[i,,3,j]<-1-(num/den.CE)																	#CE
		stats.val[i,,4,j]<-cor(proxy[rndm[(n.cal+1):nrow(dat)]],sur[rndm[(n.cal+1):nrow(dat)]])				#pearson r
		stats.val[i,,5,j]<-cor.test(proxy[rndm[(n.cal+1):nrow(dat)]],sur[rndm[(n.cal+1):nrow(dat)]])$p.value	#p value
		int<-summary(lm(sur[rndm[(n.cal+1):nrow(dat)]]~proxy[rndm[(n.cal+1):nrow(dat)]]))$coef[1]	
		intSE<-summary(lm(sur[rndm[(n.cal+1):nrow(dat)]]~proxy[rndm[(n.cal+1):nrow(dat)]]))$coef[3]	
		stats.val[i,,6,j]<-rep(int+intSE)																	#intercept+SE
		stats.val[i,,7,j]<-rep(int-intSE)																	#intercept-SE
		val.site[i,(n.cal+1):nrow(dat),j]<-rndm[(n.cal+1):nrow(dat)]										#indicies of proxy data in validation	
	}
	best.val[i,]<-which((stats.val[i,1,2,]+stats.val[i,1,3,])>=sort((stats.val[i,1,2,]+stats.val[i,1,3,]),decreasing=T)[nn])[1:nn]					#iteration numbers with the best nn RE+CE values
}
	
#*************************Extract catalog info for best matches**************************
PSRbest.ind<-matrix(,nrow=length(prxy.times),ncol=n*nn)														#n model simulations that yield the best nn RE+CE, repeats are possible
PSRbest.proxy.val<-matrix(,nrow=length(prxy.times),ncol=ceiling(max(prxy.count)*(1-prct.cal))*nn)			#indicies of validation proxy data that yield the best RE

for (i in 1:length(prxy.times)) {
	PSRbest.ind[i,]<-stats.cal[i,,1,best.val[i,]]
	v<-val.site[i,,best.val[i,]][which(is.na(val.site[i,,best.val[i,]])==FALSE)]
	PSRbest.proxy.val[i,1:length(v)]<-v	
	}

#***********************************plots****************************
palb<-rainbow(n.it,start=0.5,end=0.55,alpha=0.1);
palg<-rainbow(n.it,start=0.08,end=0.13,alpha=0.1);

dev.new(width=6, height=6)
A<-length(which(colnames(O18PI)=="HadpiC")) 
B<-length(which(colnames(O18PI)=="GISSgCpiC")) 
C<-length(which(colnames(O18PI)=="GISSgy3piC"))
D<-length(which(colnames(O18LM)=="GISSgTckLM"))
E<-length(which(colnames(O18LM)=="GISSgTKckLM"))
F<-length(which(colnames(O18LM)=="GISSgTcsLM"))

#forcings
orig_forc<-readRDS("/Users/Casey_JISAO/Google Drive/R_documents/PSR/model_data/binned_forc")
HadpiC_forc<-readRDS("/Users/Casey_JISAO/Google Drive/R_documents/PSR/model_data/HadpiC_forc")
HadpiC_forc<-matrix(HadpiC_forc,nrow=A,ncol=length(HadpiC_forc),byrow=T)
GISSgy3piC_forc<-readRDS("/Users/Casey_JISAO/Google Drive/R_documents/PSR/model_data/GISSgy3piC_forc")
GISSgy3piC_forc<-matrix(GISSgy3piC_forc,nrow=B,ncol=length(GISSgy3piC_forc),byrow=T)
GISSgCpiC_forc<-readRDS("/Users/Casey_JISAO/Google Drive/R_documents/PSR/model_data/GISSgCpiC_forc")
GISSgCpiC_forc<-matrix(GISSgCpiC_forc,nrow=C,ncol=length(GISSgCpiC_forc),byrow=T)
GISSgTckLM_forc<-readRDS("/Users/Casey_JISAO/Google Drive/R_documents/PSR/model_data/GISSgTckLM_forc")[1:D,]
GISSgTKckLM_forc<-readRDS("/Users/Casey_JISAO/Google Drive/R_documents/PSR/model_data/GISSgTKckLM_forc")[1:E,]
GISSgTcsLM_forc<-readRDS("/Users/Casey_JISAO/Google Drive/R_documents/PSR/model_data/GISSgTcsLM_forc")[1:F,]
forc.all<-rbind(HadpiC_forc,GISSgy3piC_forc,GISSgCpiC_forc,GISSgTckLM_forc,GISSgTKckLM_forc,GISSgTcsLM_forc)

plot(0,type="n", xlab="proxy year (AD)",ylab="simulation catalog index",xlim=c(0,2000),ylim=c(0,ncol(sim.all)),las=1)
	for(i in 1:length(prxy.times)) {
		points(rep(prxy.times[i],ncol(PSRbest.ind)), PSRbest.ind[i,],pch=16,cex=0.5,col=palb[1])
		}
	abline(h=c(A,A+B,A+B+C,A+B+C+D,A+B+C+D+E,A+B+C+D+E+F),col="grey50")
	text(rep(0,6),c(A-10,A+B-10,A+B+C-4,A+B+C+D-10,A+B+C+D+E-10,A+B+C+D+E+F-10),c("HadpiC","GISSgCpiC","GISSgy3piC","GISSgTckLM","GISSgTKckLM","GISSgTcsLM"),pos=4,cex=0.7)
		
dev.new(width=6, height=5)
par(mfrow=c(2,2),mar=c(2,5,2,2))
plot(0,type="n",ylab="euclidean distance",xlim=c(0,2000),ylim=c(0,0.8),las=1,main="calibration")
	for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],nn),stats.cal[i,1,2,best.val[i,]],pch=16,cex=0.5,col=palb[1])
	PSRbest.stats[i,1]<-median(stats.cal[i,1,2,best.val[i,]])
	points(prxy.times[i], PSRbest.stats[i,1])
	}
	
plot(0,type="n",ylab="pearson r.",xlim=c(0,2000),ylim=c(-1,1),las=1)
	abline(h=0,col="grey50")
	for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],nn),stats.cal[i,1,3,best.val[i,]],pch=16,cex=0.5,col=palb[1])
	PSRbest.stats[i,2]<-median(stats.cal[i,1,3,best.val[i,]])
	points(prxy.times[i],PSRbest.stats[i,2])
	
	}

plot(0,type="n",ylab="log(p)",xlim=c(0,2000),ylim=c(-4,0),las=1)
	abline(h=-1,col="grey50")
	for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],nn),log10(stats.cal[i,1,4,best.val[i,]]),pch=16,cex=0.5,col=palb[1])
	PSRbest.stats[i,3]<-median(stats.cal[i,1,4,best.val[i,]])
	points(prxy.times[i],log10(as.numeric(PSRbest.stats[i,3])))
	}
	
plot(0,type="n",xlab="year",ylab="intercept",xlim=c(0,2000),ylim=c(-0.2,0.2),las=1)
	abline(h=0,col="grey50")
	for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],nn),stats.cal[i,1,5,best.val[i,]],pch=25,cex=0.5,col=palb[1])
	points(rep(prxy.times[i],nn),stats.cal[i,1,6,best.val[i,]],pch=24,cex=0.5,col=palb[1])
	PSRbest.stats[i,4]<-median(stats.cal[i,1,5,best.val[i,]])
	PSRbest.stats[i,5]<-median(stats.cal[i,1,6,best.val[i,]])
	points(prxy.times[i],PSRbest.stats[i,4],pch=6)
	points(prxy.times[i],PSRbest.stats[i,5],pch=2)
		}
	
dev.new(width=6, height=7.5)
par(mfrow=c(3,2),mar=c(2,5,2,2))
plot(0,type="n",ylab="RE",xlim=c(0,2000),ylim=c(-0.4,0.8),las=1, main="validation")
	abline(h=0,lty=3,col="grey50")
	for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],nn),stats.val[i,1,2,best.val[i,]],pch=16,cex=0.5,col=palg[1])
	PSRbest.stats[i,6]<-median(stats.val[i,1,2,best.val[i,]])
	points(prxy.times[i],PSRbest.stats[i,6])
	}

plot(0,type="n",ylab="CE",xlim=c(0,2000),ylim=c(-0.4,0.8),las=1)
	abline(h=0,col="grey50")
	for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],nn),stats.val[i,1,3,best.val[i,]],pch=16,cex=0.5,col=palg[1])
	PSRbest.stats[i,7]<-median(stats.val[i,1,3,best.val[i,]])
	points(prxy.times[i],PSRbest.stats[i,7])
	}
	
plot(0,type="n",ylab="euclidean distance",xlim=c(0,2000),ylim=c(0,0.8),las=1)
	for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],nn),stats.val[i,1,1,best.val[i,]],pch=16,cex=0.5,col=palg[1])
	PSRbest.stats[i,8]<-median(stats.val[i,1,1,best.val[i,]])
	points(prxy.times[i],PSRbest.stats[i,8])
	}
	
plot(0,type="n",ylab="pearson r",xlim=c(0,2000),ylim=c(-1,1),las=1)
	abline(h=0,col="grey50")
	for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],nn),stats.val[i,1,4,best.val[i,]],pch=16,cex=0.5,col=palg[1])
	PSRbest.stats[i,9]<-median(stats.val[i,1,4,best.val[i,]])
	points(prxy.times[i],PSRbest.stats[i,9])
	}
	
plot(0,type="n",ylab="log(p)",xlim=c(0,2000),ylim=c(-4,0),las=1)
	abline(h=-1,col="grey50")
	for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],nn),log10(stats.val[i,1,5,best.val[i,]]),pch=16,cex=0.5,col=palg[1])
	PSRbest.stats[i,10]<-median(stats.val[i,1,5,best.val[i,]])
	points(prxy.times[i],log10(as.numeric(PSRbest.stats[i,10])))
	}

plot(0,type="n",xlab="year",ylab="intercept",xlim=c(0,2000),ylim=c(-0.2,0.2),las=1)
	abline(h=0,col="grey50")
	for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],nn),stats.val[i,1,6,best.val[i,]],pch=25,cex=0.5,col=palg[1])
	points(rep(prxy.times[i],nn),stats.val[i,1,7,best.val[i,]],pch=24,cex=0.5,col=palg[1])
	PSRbest.stats[i,11]<-median(stats.val[i,1,6,best.val[i,]])
	PSRbest.stats[i,12]<-median(stats.val[i,1,7,best.val[i,]])
	points(prxy.times[i],PSRbest.stats[i,11],pch=6)
	points(prxy.times[i],PSRbest.stats[i,12],pch=2)
	}
	
PSRbest.stats<-cbind(rep(Sys.time(),nrow(PSRbest.stats)),rep(proxy.type, nrow(PSRbest.stats)), rep(n.it,nrow(PSRbest.stats)), rep(sc,nrow(PSRbest.stats)),rep(minres,nrow(PSRbest.stats)),rep(bin.width,nrow(PSRbest.stats)), prxy.times, PSRbest.stats)

dev.new(width=8, height=8)
par(mfrow=c(3,1))
palmask<-cm.colors(3,alpha=0.7);
plot(0,type="n",xlab="year",ylab="AMOC anomaly (Sv)",xlim=c(0,2000),ylim=c(-1,1),las=1)
tmp<-matrix(,ncol=3, nrow=length(prxy.times))
	abline(h=0,col="grey80")
	for(i in 1:length(prxy.times)) {
		points(rep(prxy.times[i],n*nn),AMOC.all[PSRbest.ind[i,]],col="grey50", pch=16,cex=0.5)
		tmp[i,1]<-mean(AMOC.all[PSRbest.ind[i,]])
		tmp[i,2]<-as.numeric(quantile((AMOC.all[PSRbest.ind[i,]]))[2])
		tmp[i,3]<-as.numeric(quantile((AMOC.all[PSRbest.ind[i,]]))[4])
		}
		lines(prxy.times,tmp[,1],lwd=2); lines(prxy.times,tmp[,2]); lines(prxy.times,tmp[,3])
		PSRbest.stats<-cbind(PSRbest.stats,tmp)
		for(i in 1:length(prxy.times)) {
			if (PSRbest.stats[i,14]<0 | PSRbest.stats[i,16]<0) {
			polygon(c(prxy.times[i]-50,prxy.times[i]-50,prxy.times[i]+50,prxy.times[i]+50),c(-1.2,1.2,1.2,-1.2),col=palmask[2],border=NA)
			}
		}
plot(0,type="n",xlab="year",ylab="NAM anomaly",xlim=c(0,2000),ylim=c(-0.4,0.4),las=1)
	abline(h=0,col="grey80")
	for(i in 1:length(prxy.times)) {
		points(rep(prxy.times[i],n*nn),NAM.all[PSRbest.ind[i,]],pch=16,cex=0.5,col="grey50")
		tmp[i,1]<-mean(NAM.all[PSRbest.ind[i,]])
		tmp[i,2]<-as.numeric(quantile((NAM.all[PSRbest.ind[i,]]))[2])
		tmp[i,3]<-as.numeric(quantile((NAM.all[PSRbest.ind[i,]]))[4])
		}
		lines(prxy.times,tmp[,1],lwd=2); lines(prxy.times,tmp[,2]); lines(prxy.times,tmp[,3])
		PSRbest.stats<-cbind(PSRbest.stats,tmp)
		for(i in 1:length(prxy.times)) {
			if (PSRbest.stats[i,14]<0 | PSRbest.stats[i,16]<0) {
			polygon(c(prxy.times[i]-50,prxy.times[i]-50,prxy.times[i]+50,prxy.times[i]+50),c(-0.5,0.5,0.5,-0.5),col=palmask[2],border=NA)
			}
		}
plot(0,type="n",xlab="year",ylab="NAO anomaly",xlim=c(0,2000),ylim=c(-0.4,0.4),las=1)
	abline(h=0,col="grey80")
	for(i in 1:length(prxy.times)) {
		points(rep(prxy.times[i],n*nn),NAO.all[PSRbest.ind[i,]],pch=16,cex=0.5,col="grey50")
		tmp[i,1]<-mean(NAO.all[PSRbest.ind[i,]])
		tmp[i,2]<-as.numeric(quantile(((NAO.all[PSRbest.ind[i,]])))[2])
		tmp[i,3]<-as.numeric(quantile(((NAO.all[PSRbest.ind[i,]])))[4])
		}
		lines(prxy.times,tmp[,1],lwd=2); lines(prxy.times,tmp[,2]); lines(prxy.times,tmp[,3])
		PSRbest.stats<-cbind(PSRbest.stats,tmp)
		for(i in 1:length(prxy.times)) {
			if (PSRbest.stats[i,14]<0 | PSRbest.stats[i,16]<0) {
			polygon(c(prxy.times[i]-50,prxy.times[i]-50,prxy.times[i]+50,prxy.times[i]+50),c(-0.5,0.5,0.5,-0.5),col=palmask[2],border=NA)
			}
		}

dev.new(width=8, height=10)
par(mfrow=c(4,1))
plot(0,type="n",xlab="year",ylab="forcing (W/m^2)",xlim=c(0,2000),ylim=c(-0.5,0.5),las=1,main="WMGHG")
lines(orig_forc[seq(1,101,by=5),1],orig_forc[seq(1,101,by=5),2])
for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],n*nn),forc.all[PSRbest.ind[i,],1],pch=16,cex=0.5,col="green3")
		tmp[i,1]<-mean(forc.all[PSRbest.ind[i,],1])
		tmp[i,2]<-as.numeric(quantile(((forc.all[PSRbest.ind[i,],1])))[2])
		tmp[i,3]<-as.numeric(quantile(((forc.all[PSRbest.ind[i,],1])))[4])
		}
		lines(prxy.times,tmp[,1],lwd=2,col="green3"); lines(prxy.times,tmp[,2],col="green3"); lines(prxy.times,tmp[,3],col="green3")
		PSRbest.stats<-cbind(PSRbest.stats,tmp)

plot(0,type="n",xlab="year",ylab="forcing (W/m^2)",xlim=c(0,2000),ylim=c(-1,0),las=1,main="volcanic")
lines(orig_forc[seq(1,101,by=5),1],orig_forc[seq(1,101,by=5),3],lty=2)
lines(orig_forc[seq(1,101,by=5),1],orig_forc[seq(1,101,by=5),4])
for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],n*nn),forc.all[PSRbest.ind[i,],2],pch=16,cex=0.5,col="red3")
		tmp[i,1]<-mean(forc.all[PSRbest.ind[i,],2])
		tmp[i,2]<-as.numeric(quantile(((forc.all[PSRbest.ind[i,],2])))[2])
		tmp[i,3]<-as.numeric(quantile(((forc.all[PSRbest.ind[i,],2])))[4])
		}
		lines(prxy.times,tmp[,1],lwd=2,col="red3"); lines(prxy.times,tmp[,2],col="red3"); lines(prxy.times,tmp[,3],col="red3")
		PSRbest.stats<-cbind(PSRbest.stats,tmp)

plot(0,type="n",xlab="year",ylab="forcing (W/m^2)",xlim=c(0,2000),ylim=c(-0.4,0.1),las=1,main="solar")
lines(orig_forc[seq(1,101,by=5),1],orig_forc[seq(1,101,by=5),5],lty=2)
lines(orig_forc[seq(1,101,by=5),1],orig_forc[seq(1,101,by=5),6])
for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],n*nn),forc.all[PSRbest.ind[i,],3],pch=16,cex=0.5,col="yellow2")
		tmp[i,1]<-mean(forc.all[PSRbest.ind[i,],3])
		tmp[i,2]<-as.numeric(quantile(((forc.all[PSRbest.ind[i,],3])))[2])
		tmp[i,3]<-as.numeric(quantile(((forc.all[PSRbest.ind[i,],3])))[4])
		}
		lines(prxy.times,tmp[,1],lwd=2,col="yellow2"); lines(prxy.times,tmp[,2],col="yellow2"); lines(prxy.times,tmp[,3],col="yellow2")
		PSRbest.stats<-cbind(PSRbest.stats,tmp)

plot(0,type="n",xlab="year",ylab="forcing (W/m^2)",xlim=c(0,2000),ylim=c(-0.4,0.1),las=1,main="land use")
lines(orig_forc[seq(1,101,by=5),1],orig_forc[seq(1,101,by=5),7],lty=2)
lines(orig_forc[seq(1,101,by=5),1],orig_forc[seq(1,101,by=5),8])
for(i in 1:length(prxy.times)) {
	points(rep(prxy.times[i],n*nn),forc.all[PSRbest.ind[i,],4],pch=16,cex=0.5,col="orange4")
		tmp[i,1]<-mean(forc.all[PSRbest.ind[i,],4])
		tmp[i,2]<-as.numeric(quantile(((forc.all[PSRbest.ind[i,],4])))[2])
		tmp[i,3]<-as.numeric(quantile(((forc.all[PSRbest.ind[i,],4])))[4])
		}
		lines(prxy.times,tmp[,1],lwd=2,col="orange4"); lines(prxy.times,tmp[,2],col="orange4"); lines(prxy.times,tmp[,3],col="orange4")
		PSRbest.stats<-cbind(PSRbest.stats,tmp)
		
colnames(PSRbest.stats)<-c("time","proxy.data.type","n.iterations","scaling","min.resolution","bin.width","bin.yr","med.distance.cal","med.cal.r","med.cal.p","med.cal.int.up","med.cal.int.dn","med.val.RE","med.val.CE","med.distance.val","med.val.r","med.val.p","med.val.int.up","med.val.int.dn","mean.AMOC","AMOC25","AMOC75","mean.NAM","NAM25","NAM75","mean.NAO","NAO25","NAO75","mean.WMGHG","WMGHG25","WMGHG75","mean.vol","vol25","vol75","mean.sol","sol25","sol75","mean.land","land25","land75")

tmp<-read.csv("PSR/best_indicies/PSRbest.stats"); tmp<-tmp[,-1]
PSRbest.stats<-rbind(tmp,PSRbest.stats)

write.csv(PSRbest.ind,paste0("PSR/best_indicies/indicies",proxy.type, n.it,"_", n, "_", prct.cal*100, (1-prct.cal)*100,"_",nn,"sc",sc,"minres",minres,bin.width))
write.csv(PSRbest.stats,"PSR/best_indicies/PSRbest.stats")





