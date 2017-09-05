#Script designed to do the following: 
#1. describe proxy data spatial/temporal distribution 
#2. calculate proxy anomalies RELATIVE TO THEIR WHOLE LENGTH 
#3. bin data for PSR
#records are interpolated to 1x1 lat/lon consistent with simulations, and an anomaly is calculated anomaly relative to point 2 above 
#Each row has an index, that is unique and the same as in the Excel .csv, lat, lon, and a series of columns
#The columns of each row following longitude are the mean proxy variable in years or each window length (wl) year bins, shifted every dt years

library(plyr)
library(maps);library(mapproj); library(mapdata);

#Main variables
wl<-50;																	#window length (in years) for binning
dt<-25;																		#timestep (in years) for binning

setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/proxydata/");
SSTproxy<-read.csv("/Users/SP/Desktop/PSR_paleo/PSR_data/proxy_data/SSTProxyDat.csv", header=T);
O18proxy<-read.csv("/Users/SP/Desktop/PSR_paleo/PSR_data/proxy_data/O18ProxyDat.csv", header=T);

metadataSST<-SSTproxy[which(is.na(SSTproxy$SampleRes)==FALSE),]				#make a matrix of just metadata
metadataSST<-droplevels(metadataSST)

metadataO18<-O18proxy[which(is.na(O18proxy$SampleRes)==FALSE),]				#make a matrix of just metadata
metadataO18<-droplevels(metadataO18)

#break down by proxy type
UK37<-which(metadataSST$Species=="UK'37 calibration to SST"); UK37<-metadataSST[UK37,]
diatomTF<-which(metadataSST$Species=="Diatom transfer function"); diatomTF<-metadataSST[diatomTF,]
TEX86<-which(metadataSST$Species=="TEX86 calibration to SST");TEX86<-metadataSST[TEX86,]
MgCa<-which(metadataSST$Species!="UK'37 calibration to SST" & metadataSST$Species!="Diatom transfer function" & metadataSST$Species!="TEX86 calibration to SST"); MgCa<-metadataSST[MgCa,]

#some fractions less than a given value
f100.UK37<-length(which(UK37$SampleRes>=1))/nrow(UK37)
f50.UK37<-length(which(UK37$SampleRes>=2))/nrow(UK37)
f100.MgCa<-length(which(MgCa$SampleRes>=1))/nrow(MgCa)
f50.MgCa<-length(which(MgCa$SampleRes>=2))/nrow(MgCa)
f100.TEX86<-length(which(TEX86$SampleRes>=1))/nrow(TEX86)
f50.TEX86<-length(which(TEX86$SampleRes>=2))/nrow(TEX86)
f100.O18<-length(which(metadataO18$SampleRes>=1))/nrow(metadataO18)
f50.O18<-length(which(metadataO18$SampleRes>=2))/nrow(metadataO18)

#some basic plots on temporal distribution
dev.new(width=7, height=7); par(mfrow=c(2,2)) 
hist(MgCa$SampleRes,breaks=seq(0,106,by=2),xlab="avg. sampling (samples/100yr)", main="Mg/Ca",col="grey50",las=1,ylim=c(0,14))		#sampling resolution
hist(UK37$SampleRes,breaks=seq(0,106,by=2),xlab="avg. sampling (samples/100yr)", main="Uk 37",col="grey50",las=1,ylim=c(0,14))		#sampling resolution
hist(TEX86$SampleRes,breaks=seq(0,106,by=2),xlab="avg. sampling (samples/100yr)", main="TEX 86",col="grey50",las=1,ylim=c(0,14))		#sampling resolution
hist(metadataO18$SampleRes,breaks=seq(0,106,by=2),xlab="avg. sampling (samples/100yr)", main="d18O",col="grey50",las=1,ylim=c(0,14))

#some maps of spatial distribution and proxy type
dev.new(width=12, height=6);
world<-map('worldHires',interior=FALSE,ylim=c(-65,88),xlim=c(-180,180))
#points(UK37$lon,UK37$lat,pch=1,col="blue",cex=1+log10(UK37$SampleRes),lwd=2)
points(MgCa$lon,MgCa$lat,pch=2,col="palegreen3",cex=1+log10(MgCa$SampleRes),lwd=2)
#points(diatomTF$lon,diatomTF$lat,pch=5,col="salmon3",cex=1+log10(diatomTF$SampleRes),lwd=2)
#points(TEX86$lon,TEX86$lat,pch=0,col="grey50",cex=1+log10(TEX86$SampleRes),lwd=2)
points(metadataO18$lon,metadataO18$lat,pch=6,col="magenta4",cex=1+log10(metadataO18$SampleRes),lwd=2)
legend(30,88,legend=c("Mg/Ca","d18O","1 sample/100 yrs","20 samples/100 yrs","100 samples/100 yrs"), pch=c(2,6,2,2,2), col=c("palegreen3","magenta4","grey50","grey50","grey50"), pt.cex=c(1,1,1,1+log10(20),3),bg="white")

#distill metadata 
infobinUK37<-UK37[,-c(2:4,6,8,9,12,14:17)]
names(infobinUK37)<-c("index","resolution","name","lat","lon","season")
infobinO18<-metadataO18[,-c(2,3,5,7,8,11,13)]
names(infobinO18)<-c("index","resolution","name","lat","lon","season")
infobinMgCa<-MgCa[,-c(2:4,6,8,9,12,14:17)]
names(infobinMgCa)<-c("index","resolution","name","lat","lon","season")
infobinTEX86<-TEX86[,-c(2:4,6,8,9,12,14:17)]
names(infobinTEX86)<-c("index","resolution","name","lat","lon","season")

#bin with anomalies relative to mean of each record
bin<-seq(2000-wl/2,wl/2,by=-dt)
MgCabin100.IND<-matrix(,nrow=nrow(MgCa),ncol=length(bin))					#empty matrix to populate with MgCa
O18bin100.IND<-matrix(,nrow=max(O18proxy$index),ncol=length(bin))			#empty matrix to populate with SST anomalies

dev.new(width=8, height=4);
par(mfrow=c(1,2))
plot(0,type="n",xlim=c(0,2000),ylim=c(-4,4),ylab="binned std. anomaly",xlab="year (AD)",main=paste0(wl, "yr bins"),las=1)

#O18 data
for (i in 1:max(O18proxy$index)) {
	dat<-which(O18proxy$index==i);
	for (j in 1:length(bin)) {												#the binning
		bb<-mean(O18proxy$d18O[dat][which(O18proxy$year[dat]<=bin[j]+wl/2 & O18proxy$year[dat]>bin[j]-wl/2)])
		O18bin100.IND[i,j]<-ifelse(is.nan(bb)==T,NA,bb)						#return NA if bin is empty
	}
}
O18bin100.IND<-O18bin100.IND-rowMeans(O18bin100.IND,na.rm=T) 						#anomalies	
normO18bin100.IND<-O18bin100.IND/(apply(O18bin100.IND,1,sd,na.rm=T))

#MgCa data
for (i in 1:nrow(MgCa)) {
	dat<-which(SSTproxy$index==MgCa$index[i]);
		for (j in 1:length(bin)) {												#the binning
		bb<-mean(SSTproxy$Tproxy.val[dat][which(SSTproxy$year[dat]<=bin[j]+wl/2 & SSTproxy$year[dat]>bin[j]-wl/2)])
		MgCabin100.IND[i,j]<-ifelse(is.nan(bb)==T,NA,bb)						#return NA if bin is empty
	}
}
MgCabin100.IND<-MgCabin100.IND-rowMeans(MgCabin100.IND,na.rm=T) 						#anomalies	
normMgCabin100.IND<-MgCabin100.IND/(apply(MgCabin100.IND,1,sd,na.rm=T))

for (i in 1:nrow(O18bin100.IND))	{	
	points(bin,normO18bin100.IND[i,],pch="-",col="magenta4") 
}
for (i in 1:nrow(MgCabin100.IND))	{	
	points(bin,normMgCabin100.IND[i,],pch="-",col="palegreen3") 
}

O18bin100.IND<-data.frame(O18bin100.IND); names(O18bin100.IND)<-bin
O18bin100.IND<-cbind(infobinO18,O18bin100.IND);
normO18bin100.IND<-data.frame(normO18bin100.IND); names(normO18bin100.IND)<-bin
normO18bin100.IND<-cbind(infobinO18,normO18bin100.IND);
MgCabin100.IND<-data.frame(MgCabin100.IND); names(MgCabin100.IND)<-bin
MgCabin100.IND<-cbind(infobinMgCa,MgCabin100.IND);
normMgCabin100.IND<-data.frame(normMgCabin100.IND); names(normMgCabin100.IND)<-bin
normMgCabin100.IND<-cbind(infobinMgCa,normMgCabin100.IND);

O18prxy.count<-colSums(!is.na(O18bin100.IND))[7:ncol(O18bin100.IND)]			#number of records through time
plot(bin,O18prxy.count,pch=20,las=1,xlab="year (AD)",ylab="number of records",ylim=c(10,45),col="magenta4");
segments(bin-wl/2,O18prxy.count,bin+wl/2,col="magenta4");
MgCa.count<-colSums(!is.na(MgCabin100.IND))[7:ncol(MgCabin100.IND)]			#number of MgCarecords through time
points(bin,MgCa.count,pch=20,col="palegreen3");
segments(bin-wl/2,MgCa.count,bin+wl/2,col="palegreen3")


saveRDS(MgCabin100.IND, file=paste0("proxy_data/MgCa_dt",dt,"wl",wl,"anom"))
saveRDS(normMgCabin100.IND, file=paste0("proxy_data/MgCa_dt",dt,"wl",wl,"std.anom"))
saveRDS(bin,file=paste0("proxy_data/MgCa_midYear_dt",dt,"wl",wl,"anom"))

saveRDS(O18bin100.IND, file=paste0("proxy_data/O18_dt",dt,"wl",wl,"anom"))
saveRDS(normO18bin100.IND, file=paste0("proxy_data/O18_dt",dt,"wl",wl,"std.anom"))
saveRDS(bin,file=paste0("proxy_data/O18_midYear_dt",dt,"wl",wl,"anom"))
