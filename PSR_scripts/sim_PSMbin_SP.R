# calculates forward models of Mg/Ca and d18O and then bins. Also bin AMOC, NAO and NAM

#----Stephanie Pennington version for generating pseudoproxy data----#
#Edited June 19, 2017

setwd("/Users/SP/Desktop/PSR_paleo/PSR_scripts/");

#constants
R18SMOW<-0.0020052																								#18R for VSMOW
MgCaTdep<-0.087																									#MgCa T dependence
width<-50																									#bin width
step<-10																										#bin step				

#read in annual simulation data
HadpiC_T_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/HadpiC_T_O18")
HadpiC_T_O18<-HadpiC_T_O18[,-1]																					#delete first column of NAs
HadpiC_H2O18_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/HadpiC_H2O18_O18")
HadpiC_H2O18_O18<-HadpiC_H2O18_O18[,-1]
HadpiC_H2O18_O18<-(HadpiC_H2O18_O18/R18SMOW-1)*1000																#convert to delta notation
HadpiC_T_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/HadpiC_T_MgCa")
HadpiC_T_MgCa<-HadpiC_T_MgCa[,-1]
HadpiC_salt_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/HadpiC_salt_MgCa")
HadpiC_salt_MgCa<-HadpiC_salt_MgCa[,-1]
HadpiC_AMOC<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/HadpiC_AMOC")
HadpiC_NAM<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/HadpiC_NAM")
HadpiC_NAO<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/HadpiC_NAO")

GISSgCpiC_T_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgCpiC_T_O18")
GISSgCpiC_T_O18[,848]<-rowMeans(GISSgCpiC_T_O18,na.rm=TRUE)														#replace an outlier in year 848 with mean
GISSgCpiC_T_O18<-GISSgCpiC_T_O18[,-1]																			
GISSgCpiC_H2O18_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgCpiC_H2O18_O18")
GISSgCpiC_H2O18_O18<-GISSgCpiC_H2O18_O18[,-1]
GISSgCpiC_T_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgCpiC_T_MgCa")
GISSgCpiC_T_MgCa[,848]<-rowMeans(GISSgCpiC_T_MgCa,na.rm=TRUE)													#replace an outlier in year 848 with mean
GISSgCpiC_T_MgCa<-GISSgCpiC_T_MgCa[,-1]
GISSgCpiC_salt_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgCpiC_salt_MgCa")
GISSgCpiC_salt_MgCa[,848]<-rowMeans(GISSgCpiC_salt_MgCa,na.rm=TRUE)												#replace an outlier in year 848 with mean
GISSgCpiC_salt_MgCa<-GISSgCpiC_salt_MgCa[,-1]
GISSgCpiC_AMOC<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgCpiC_AMOC")
GISSgCpiC_NAM<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgCpiC_NAM")
GISSgCpiC_NAO<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgCpiC_NAO")

GISSgy3piC_T_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgy3piC_T_O18")
GISSgy3piC_T_O18<-GISSgy3piC_T_O18[,-1]
GISSgy3piC_H2O18_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgy3piC_H2O18_O18")
GISSgy3piC_H2O18_O18<-GISSgy3piC_H2O18_O18[,-1]
GISSgy3piC_T_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgy3piC_T_MgCa")
GISSgy3piC_T_MgCa<-GISSgy3piC_T_MgCa[,-1]
GISSgy3piC_salt_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgy3piC_salt_MgCa")
GISSgy3piC_salt_MgCa<-GISSgy3piC_salt_MgCa[,-1]
GISSgy3piC_AMOC<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgy3piC_AMOC")
GISSgy3piC_NAM<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgy3piC_NAM")
GISSgy3piC_NAO<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgy3piC_NAO")

GISSgTckLM_T_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTckLM_T_O18")
GISSgTckLM_T_O18<-GISSgTckLM_T_O18[,-1]
GISSgTckLM_H2O18_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTckLM_H2O18_O18")
GISSgTckLM_H2O18_O18<-GISSgTckLM_H2O18_O18[,-1]
GISSgTckLM_T_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTckLM_T_MgCa")
GISSgTckLM_T_MgCa<-GISSgTckLM_T_MgCa[,-1]
GISSgTckLM_salt_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTckLM_salt_MgCa")
GISSgTckLM_salt_MgCa<-GISSgTckLM_salt_MgCa[,-1]
GISSgTckLM_AMOC<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTckLM_AMOC")
GISSgTckLM_NAM<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTckLM_NAM")
GISSgTckLM_NAO<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTckLM_NAO")

GISSgTKckLM_T_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTKckLM_T_O18")
GISSgTKckLM_T_O18<-GISSgTKckLM_T_O18[,-1]
GISSgTKckLM_H2O18_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTKckLM_H2O18_O18")
GISSgTKckLM_H2O18_O18<-GISSgTKckLM_H2O18_O18[,-1]
GISSgTKckLM_T_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTKckLM_T_MgCa")
GISSgTKckLM_T_MgCa<-GISSgTKckLM_T_MgCa[,-1]
GISSgTKckLM_salt_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTKckLM_salt_MgCa")
GISSgTKckLM_salt_MgCa<-GISSgTKckLM_salt_MgCa[,-1]
missing<-which(is.na(GISSgTKckLM_salt_MgCa)==TRUE)																	#some salt values seem to be missing. 
GISSgTKckLM_AMOC<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTKckLM_AMOC")
GISSgTKckLM_NAM<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTKckLM_NAM")
GISSgTKckLM_NAO<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTKckLM_NAO")

GISSgTcsLM_T_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTcsLM_T_O18")
GISSgTcsLM_T_O18<-GISSgTcsLM_T_O18[,-1]
GISSgTcsLM_H2O18_O18<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTcsLM_H2O18_O18")
GISSgTcsLM_H2O18_O18<-GISSgTcsLM_H2O18_O18[,-1]
GISSgTcsLM_T_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTcsLM_T_MgCa")
GISSgTcsLM_T_MgCa<-GISSgTcsLM_T_MgCa[,-1]
GISSgTcsLM_salt_MgCa<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTcsLM_salt_MgCa")
GISSgTcsLM_salt_MgCa<-GISSgTcsLM_salt_MgCa[,-1]
GISSgTcsLM_AMOC<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTcsLM_AMOC")
GISSgTcsLM_NAM<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTcsLM_NAM")
GISSgTcsLM_NAO<-readRDS("/Users/SP/Desktop/PSR_paleo/PSR_data/model_annualdata/GISSgTcsLM_NAO")

#for GISS detrend following Allegra's advice. Detrend all variables with the same parameters.
for (i in 1:nrow(GISSgCpiC_T_O18)) {
	tmp<-lowess(GISSgCpiC_T_O18[i,1:ncol(GISSgCpiC_T_O18)],f=0.3)
	GISSgCpiC_T_O18[i,]<-GISSgCpiC_T_O18[i,]-tmp$y+mean(tmp$y)
	tmp<-lowess(GISSgCpiC_H2O18_O18[i,1:ncol(GISSgCpiC_H2O18_O18)],f=0.3)
	GISSgCpiC_H2O18_O18[i,]<-GISSgCpiC_H2O18_O18[i,]-tmp$y+mean(tmp$y)
	}
for (i in 1:nrow(GISSgCpiC_T_MgCa)) {
	tmp<-lowess(GISSgCpiC_T_MgCa[i,1:ncol(GISSgCpiC_T_MgCa)],f=0.3)
	GISSgCpiC_T_MgCa[i,]<-GISSgCpiC_T_MgCa[i,]-tmp$y+mean(tmp$y)
	tmp<-lowess(GISSgCpiC_salt_MgCa[i,1:ncol(GISSgCpiC_salt_MgCa)],f=0.3)
	GISSgCpiC_salt_MgCa[i,]<-GISSgCpiC_salt_MgCa[i,]-tmp$y+mean(tmp$y)
	}
#detrend AMOC/NAM/NAO in the same way
tmp<-lowess(GISSgCpiC_AMOC,f=0.3); GISSgCpiC_AMOC<-GISSgCpiC_AMOC-tmp$y+mean(tmp$y)
tmp<-lowess(GISSgCpiC_NAM,f=0.3); GISSgCpiC_NAM<-GISSgCpiC_NAM-tmp$y+mean(tmp$y)
tmp<-lowess(GISSgCpiC_NAO,f=0.3); GISSgCpiC_NAO<-GISSgCpiC_NAO-tmp$y+mean(tmp$y)

for (i in 1:nrow(GISSgy3piC_T_O18)) {
	tmp<-lowess(GISSgy3piC_T_O18[i,1:ncol(GISSgy3piC_T_O18)],f=0.3)
	GISSgy3piC_T_O18[i,]<-GISSgy3piC_T_O18[i,]-tmp$y+mean(tmp$y)
	tmp<-lowess(GISSgy3piC_H2O18_O18[i,1:ncol(GISSgy3piC_H2O18_O18)],f=0.3)
	GISSgy3piC_H2O18_O18[i,]<-GISSgy3piC_H2O18_O18[i,]-tmp$y+mean(tmp$y)
	}
for (i in 1:nrow(GISSgy3piC_T_MgCa)) {
	tmp<-lowess(GISSgy3piC_T_MgCa[i,1:ncol(GISSgy3piC_T_MgCa)],f=0.3)
	GISSgy3piC_T_MgCa[i,]<-GISSgy3piC_T_MgCa[i,]-tmp$y+mean(tmp$y)
	tmp<-lowess(GISSgy3piC_salt_MgCa[i,1:ncol(GISSgy3piC_salt_MgCa)],f=0.3)
	GISSgy3piC_salt_MgCa[i,]<-GISSgy3piC_salt_MgCa[i,]-tmp$y+mean(tmp$y)
	}
tmp<-lowess(GISSgy3piC_AMOC,f=0.3); GISSgy3piC_AMOC<-GISSgy3piC_AMOC-tmp$y+mean(tmp$y)
tmp<-lowess(GISSgy3piC_NAM,f=0.3); GISSgy3piC_NAM<-GISSgy3piC_NAM-tmp$y+mean(tmp$y)
tmp<-lowess(GISSgy3piC_NAO,f=0.3); GISSgy3piC_NAO<-GISSgy3piC_NAO-tmp$y+mean(tmp$y)

for (i in 1:nrow(GISSgTckLM_T_O18)) {
	tmp<-lowess(GISSgTckLM_T_O18[i,1:ncol(GISSgTckLM_T_O18)],f=0.3)
	GISSgTckLM_T_O18[i,]<-GISSgTckLM_T_O18[i,]-tmp$y+mean(tmp$y)
	tmp<-lowess(GISSgTckLM_H2O18_O18[i,1:ncol(GISSgTckLM_H2O18_O18)],f=0.3)
	GISSgTckLM_H2O18_O18[i,]<-GISSgTckLM_H2O18_O18[i,]-tmp$y+mean(tmp$y)
	}
for (i in 1:nrow(GISSgTckLM_T_MgCa)) {
	tmp<-lowess(GISSgTckLM_T_MgCa[i,1:ncol(GISSgTckLM_T_MgCa)],f=0.3)
	GISSgTckLM_T_MgCa[i,]<-GISSgTckLM_T_MgCa[i,]-tmp$y+mean(tmp$y)
	tmp<-lowess(GISSgTckLM_salt_MgCa[i,1:ncol(GISSgTckLM_salt_MgCa)],f=0.3)
	GISSgTckLM_salt_MgCa[i,]<-GISSgTckLM_salt_MgCa[i,]-tmp$y+mean(tmp$y)
	}
tmp<-lowess(GISSgTckLM_AMOC,f=0.3); GISSgTckLM_AMOC<-GISSgTckLM_AMOC-tmp$y+mean(tmp$y)
tmp<-lowess(GISSgTckLM_NAM,f=0.3); GISSgTckLM_NAM<-GISSgTckLM_NAM-tmp$y+mean(tmp$y)
tmp<-lowess(GISSgTckLM_NAO,f=0.3); GISSgTckLM_NAO<-GISSgTckLM_NAO-tmp$y+mean(tmp$y)

for (i in 1:nrow(GISSgTKckLM_T_O18)) {
	tmp<-lowess(GISSgTKckLM_T_O18[i,1:ncol(GISSgTKckLM_T_O18)],f=0.3)
	GISSgTKckLM_T_O18[i,]<-GISSgTKckLM_T_O18[i,]-tmp$y+mean(tmp$y)
	tmp<-lowess(GISSgTKckLM_H2O18_O18[i,1:ncol(GISSgTKckLM_H2O18_O18)],f=0.3)
	GISSgTKckLM_H2O18_O18[i,]<-GISSgTKckLM_H2O18_O18[i,]-tmp$y+mean(tmp$y)
	}
for (i in 1:nrow(GISSgTKckLM_T_MgCa)) {
	tmp<-lowess(GISSgTKckLM_T_MgCa[i,1:ncol(GISSgTKckLM_T_MgCa)],f=0.3)
	GISSgTKckLM_T_MgCa[i,]<-GISSgTKckLM_T_MgCa[i,]-tmp$y+mean(tmp$y)
	missing<-which(is.na(GISSgTKckLM_salt_MgCa[i,1:ncol(GISSgTKckLM_salt_MgCa)])==TRUE | GISSgTKckLM_salt_MgCa[i,1:ncol(GISSgTKckLM_salt_MgCa)]==0) 	#some missing data
	if (length(missing)>0) {
		GISSgTKckLM_salt_MgCa[i,missing]<-rnorm(length(missing),mean=mean(GISSgTKckLM_salt_MgCa[i,-missing]),sd=sd(GISSgTKckLM_salt_MgCa[i,-missing]))	#replace with random at that site
	} 	
	tmp<-lowess(GISSgTKckLM_salt_MgCa[i,1:ncol(GISSgTKckLM_salt_MgCa)],f=0.3)
	GISSgTKckLM_salt_MgCa[i,]<-GISSgTKckLM_salt_MgCa[i,]-tmp$y+mean(tmp$y)
	}
	
tmp<-lowess(GISSgTKckLM_AMOC,f=0.3); GISSgTKckLM_AMOC<-GISSgTKckLM_AMOC-tmp$y+mean(tmp$y)
tmp<-lowess(GISSgTKckLM_NAM,f=0.3); GISSgTKckLM_NAM<-GISSgTKckLM_NAM-tmp$y+mean(tmp$y)
tmp<-lowess(GISSgTKckLM_NAO,f=0.3); GISSgTKckLM_NAO<-GISSgTKckLM_NAO-tmp$y+mean(tmp$y)

for (i in 1:nrow(GISSgTcsLM_T_O18)) {
	tmp<-lowess(GISSgTcsLM_T_O18[i,1:ncol(GISSgTcsLM_T_O18)],f=0.3)
	GISSgTcsLM_T_O18[i,]<-GISSgTcsLM_T_O18[i,]-tmp$y+mean(tmp$y)
	tmp<-lowess(GISSgTcsLM_H2O18_O18[i,1:ncol(GISSgTcsLM_H2O18_O18)],f=0.3)
	GISSgTcsLM_H2O18_O18[i,]<-GISSgTcsLM_H2O18_O18[i,]-tmp$y+mean(tmp$y)
	}
for (i in 1:nrow(GISSgTcsLM_T_MgCa)) {
	tmp<-lowess(GISSgTcsLM_T_MgCa[i,1:ncol(GISSgTcsLM_T_MgCa)],f=0.3)
	GISSgTcsLM_T_MgCa[i,]<-GISSgTcsLM_T_MgCa[i,]-tmp$y+mean(tmp$y)
	tmp<-lowess(GISSgTcsLM_salt_MgCa[i,1:ncol(GISSgTcsLM_salt_MgCa)],f=0.3)
	GISSgTcsLM_salt_MgCa[i,]<-GISSgTcsLM_salt_MgCa[i,]-tmp$y+mean(tmp$y)
	}
tmp<-lowess(GISSgTcsLM_AMOC,f=0.3); GISSgTcsLM_AMOC<-GISSgTcsLM_AMOC-tmp$y+mean(tmp$y)
tmp<-lowess(GISSgTcsLM_NAM,f=0.3); GISSgTcsLM_NAM<-GISSgTcsLM_NAM-tmp$y+mean(tmp$y)
tmp<-lowess(GISSgTcsLM_NAO,f=0.3); GISSgTcsLM_NAO<-GISSgTcsLM_NAO-tmp$y+mean(tmp$y)

#convert d18Ow to PDB
HadpiC_H2O18_O18<-(HadpiC_H2O18_O18-30.92)/1.03092
GISSgCpiC_H2O18_O18<-(GISSgCpiC_H2O18_O18-30.92)/1.03092
GISSgy3piC_H2O18_O18<-(GISSgy3piC_H2O18_O18-30.92)/1.03092
GISSgTckLM_H2O18_O18<-(GISSgTckLM_H2O18_O18-30.92)/1.03092
GISSgTKckLM_H2O18_O18<-(GISSgTKckLM_H2O18_O18-30.92)/1.03092
GISSgTcsLM_H2O18_O18<-(GISSgTcsLM_H2O18_O18-30.92)/1.03092

#forward model
HadpiC_alpha<-exp((17.99*(1000/(HadpiC_T_O18+273.15))-32.42)/1000)
HadpiC_O18<-HadpiC_alpha*(1000+HadpiC_H2O18_O18)-1000
HadpiC_MgCa<-exp(MgCaTdep*HadpiC_T_MgCa+0.039*HadpiC_salt_MgCa-2.25)

GISSgCpiC_alpha<-exp((17.99*(1000/(GISSgCpiC_T_O18+273.15))-32.42)/1000)
GISSgCpiC_O18<-GISSgCpiC_alpha*(1000+GISSgCpiC_H2O18_O18)-1000
GISSgCpiC_MgCa<-exp(MgCaTdep*GISSgCpiC_T_MgCa+0.039*GISSgCpiC_salt_MgCa-2.25)

GISSgy3piC_alpha<-exp((17.99*(1000/(GISSgy3piC_T_O18+273.15))-32.42)/1000)
GISSgy3piC_O18<-GISSgy3piC_alpha*(1000+GISSgy3piC_H2O18_O18)-1000
GISSgy3piC_MgCa<-exp(MgCaTdep*GISSgy3piC_T_MgCa+0.039*GISSgy3piC_salt_MgCa-2.25)

GISSgTckLM_alpha<-exp((17.99*(1000/(GISSgTckLM_T_O18+273.15))-32.42)/1000)
GISSgTckLM_O18<-GISSgTckLM_alpha*(1000+GISSgTckLM_H2O18_O18)-1000
GISSgTckLM_MgCa<-exp(MgCaTdep*GISSgTckLM_T_MgCa+0.039*GISSgTckLM_salt_MgCa-2.25)

GISSgTKckLM_alpha<-exp((17.99*(1000/(GISSgTKckLM_T_O18+273.15))-32.42)/1000)
GISSgTKckLM_O18<-GISSgTKckLM_alpha*(1000+GISSgTKckLM_H2O18_O18)-1000
GISSgTKckLM_MgCa<-exp(MgCaTdep*GISSgTKckLM_T_MgCa+0.039*GISSgTKckLM_salt_MgCa-2.25)

GISSgTcsLM_alpha<-exp((17.99*(1000/(GISSgTcsLM_T_O18+273.15))-32.42)/1000)
GISSgTcsLM_O18<-GISSgTcsLM_alpha*(1000+GISSgTcsLM_H2O18_O18)-1000
GISSgTcsLM_MgCa<-exp(MgCaTdep*GISSgTcsLM_T_MgCa+0.039*GISSgTcsLM_salt_MgCa-2.25)

#empty bin
wt.HadpiC<-floor((ncol(HadpiC_O18)-width)/step+1)
wt.GISSgCpiC<-floor((ncol(GISSgCpiC_O18)-width)/step+1)
wt.GISSgy3piC<-floor((ncol(GISSgy3piC_O18)-width)/step+1)
wt.GISSgTckLM<-floor((ncol(GISSgTckLM_O18)-width)/step+1)
wt.GISSgTKckLM<-floor((ncol(GISSgTKckLM_O18)-width)/step+1)
wt.GISSgTcsLM<-floor((ncol(GISSgTcsLM_O18)-width)/step+1)

HadpiC_O18.bin<-data.frame(matrix(,nrow=nrow(HadpiC_O18),ncol=wt.HadpiC))
HadpiC_MgCa.bin<-data.frame(matrix(,nrow=nrow(HadpiC_MgCa),ncol=wt.HadpiC))
HadpiC_AMOC.bin<-seq(,length.out=wt.HadpiC)
HadpiC_NAM.bin<-seq(,length.out=wt.HadpiC)
HadpiC_NAO.bin<-seq(,length.out=wt.HadpiC)

GISSgCpiC_O18.bin<-data.frame(matrix(,nrow=nrow(GISSgCpiC_O18),ncol=wt.GISSgCpiC))
GISSgCpiC_MgCa.bin<-data.frame(matrix(,nrow=nrow(GISSgCpiC_MgCa),ncol=wt.GISSgCpiC))
GISSgCpiC_AMOC.bin<-seq(,length.out=wt.GISSgCpiC)
GISSgCpiC_NAM.bin<-seq(,length.out=wt.GISSgCpiC)
GISSgCpiC_NAO.bin<-seq(,length.out=wt.GISSgCpiC)

GISSgy3piC_O18.bin<-data.frame(matrix(,nrow=nrow(GISSgy3piC_O18),ncol=wt.GISSgy3piC))
GISSgy3piC_MgCa.bin<-data.frame(matrix(,nrow=nrow(GISSgy3piC_MgCa),ncol=wt.GISSgy3piC))
GISSgy3piC_AMOC.bin<-seq(,length.out=wt.GISSgy3piC)
GISSgy3piC_NAO.bin<-seq(,length.out=wt.GISSgy3piC)
GISSgy3piC_NAM.bin<-seq(,length.out=wt.GISSgy3piC)

GISSgTckLM_O18.bin<-data.frame(matrix(,nrow=nrow(GISSgTckLM_O18),ncol=wt.GISSgTckLM))
GISSgTckLM_MgCa.bin<-data.frame(matrix(,nrow=nrow(GISSgTckLM_MgCa),ncol=wt.GISSgTckLM))
GISSgTckLM_AMOC.bin<-seq(,length.out=wt.GISSgTckLM)
GISSgTckLM_NAO.bin<-seq(,length.out=wt.GISSgTckLM)
GISSgTckLM_NAM.bin<-seq(,length.out=wt.GISSgTckLM)

GISSgTKckLM_O18.bin<-data.frame(matrix(,nrow=nrow(GISSgTKckLM_O18),ncol=wt.GISSgTKckLM))
GISSgTKckLM_MgCa.bin<-data.frame(matrix(,nrow=nrow(GISSgTKckLM_MgCa),ncol=wt.GISSgTKckLM))
GISSgTKckLM_AMOC.bin<-seq(,length.out=wt.GISSgTKckLM)
GISSgTKckLM_NAO.bin<-seq(,length.out=wt.GISSgTKckLM)
GISSgTKckLM_NAM.bin<-seq(,length.out=wt.GISSgTKckLM)

GISSgTcsLM_O18.bin<-data.frame(matrix(,nrow=nrow(GISSgTcsLM_O18),ncol=wt.GISSgTcsLM))
GISSgTcsLM_MgCa.bin<-data.frame(matrix(,nrow=nrow(GISSgTcsLM_MgCa),ncol=wt.GISSgTcsLM))
GISSgTcsLM_AMOC.bin<-seq(,length.out=wt.GISSgTcsLM)
GISSgTcsLM_NAO.bin<-seq(,length.out=wt.GISSgTcsLM)
GISSgTcsLM_NAM.bin<-seq(,length.out=wt.GISSgTcsLM)

#bin 
for (i in 1:wt.HadpiC) {
	HadpiC_O18.bin[,i]<-rowMeans(HadpiC_O18[,(1+(step*(i-1))):(width+(step*(i-1)))])
	HadpiC_MgCa.bin[,i]<-rowMeans(HadpiC_MgCa[,(1+(step*(i-1))):(width+(step*(i-1)))])
	HadpiC_AMOC.bin[i]<-mean(HadpiC_AMOC[(1+(step*(i-1))):(width+(step*(i-1)))])
	HadpiC_NAM.bin[i]<-mean(HadpiC_NAM[(1+(step*(i-1))):(width+(step*(i-1)))])
	HadpiC_NAO.bin[i]<-mean(HadpiC_NAO[(1+(step*(i-1))):(width+(step*(i-1)))])
}
for (i in 1:wt.GISSgCpiC) {
	GISSgCpiC_O18.bin[,i]<-rowMeans(GISSgCpiC_O18[,(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgCpiC_MgCa.bin[,i]<-rowMeans(GISSgCpiC_MgCa[,(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgCpiC_AMOC.bin[i]<-mean(GISSgCpiC_AMOC[(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgCpiC_NAM.bin[i]<-mean(GISSgCpiC_NAM[(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgCpiC_NAO.bin[i]<-mean(GISSgCpiC_NAO[(1+(step*(i-1))):(width+(step*(i-1)))])
}
for (i in 1:wt.GISSgy3piC) {
	GISSgy3piC_O18.bin[,i]<-rowMeans(GISSgy3piC_O18[,(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgy3piC_MgCa.bin[,i]<-rowMeans(GISSgy3piC_MgCa[,(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgy3piC_AMOC.bin[i]<-mean(GISSgy3piC_AMOC[(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgy3piC_NAM.bin[i]<-mean(GISSgy3piC_NAM[(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgy3piC_NAO.bin[i]<-mean(GISSgy3piC_NAO[(1+(step*(i-1))):(width+(step*(i-1)))])
}
for (i in 1:wt.GISSgTckLM) {
	GISSgTckLM_O18.bin[,i]<-rowMeans(GISSgTckLM_O18[,(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTckLM_MgCa.bin[,i]<-rowMeans(GISSgTckLM_MgCa[,(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTckLM_AMOC.bin[i]<-mean(GISSgTckLM_AMOC[(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTckLM_NAM.bin[i]<-mean(GISSgTckLM_NAM[(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTckLM_NAO.bin[i]<-mean(GISSgTckLM_NAO[(1+(step*(i-1))):(width+(step*(i-1)))])
}
for (i in 1:wt.GISSgTKckLM) {
	GISSgTKckLM_O18.bin[,i]<-rowMeans(GISSgTKckLM_O18[,(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTKckLM_MgCa.bin[,i]<-rowMeans(GISSgTKckLM_MgCa[,(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTKckLM_AMOC.bin[i]<-mean(GISSgTKckLM_AMOC[(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTKckLM_NAM.bin[i]<-mean(GISSgTKckLM_NAM[(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTKckLM_NAO.bin[i]<-mean(GISSgTKckLM_NAO[(1+(step*(i-1))):(width+(step*(i-1)))])
}
for (i in 1:wt.GISSgTcsLM) {
	GISSgTcsLM_O18.bin[,i]<-rowMeans(GISSgTcsLM_O18[,(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTcsLM_MgCa.bin[,i]<-rowMeans(GISSgTcsLM_MgCa[,(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTcsLM_AMOC.bin[i]<-mean(GISSgTcsLM_AMOC[(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTcsLM_NAM.bin[i]<-mean(GISSgTcsLM_NAM[(1+(step*(i-1))):(width+(step*(i-1)))])
	GISSgTcsLM_NAO.bin[i]<-mean(GISSgTcsLM_NAO[(1+(step*(i-1))):(width+(step*(i-1)))])
}

#calculate anomalies and std anomalies
HadpiC_O18.bin<-HadpiC_O18.bin-rowMeans(HadpiC_O18.bin)
HadpiC_O18.binstd<-HadpiC_O18.bin/apply(HadpiC_O18.bin,1,sd)
HadpiC_MgCa.bin<-HadpiC_MgCa.bin-rowMeans(HadpiC_MgCa.bin)
HadpiC_MgCa.binstd<-HadpiC_MgCa.bin/apply(HadpiC_MgCa.bin,1,sd)
HadpiC_AMOC.bin<-HadpiC_AMOC.bin-mean(HadpiC_AMOC.bin)
HadpiC_AMOC.binstd<-HadpiC_AMOC.bin/sd(HadpiC_AMOC.bin)
HadpiC_NAO.bin<-HadpiC_NAO.bin-mean(HadpiC_NAO.bin)
HadpiC_NAO.binstd<-HadpiC_NAO.bin/sd(HadpiC_NAO.bin)
HadpiC_NAM.bin<-HadpiC_NAM.bin-mean(HadpiC_NAM.bin)
HadpiC_NAM.binstd<-HadpiC_NAM.bin/sd(HadpiC_NAM.bin)

GISSgCpiC_O18.bin<-GISSgCpiC_O18.bin-rowMeans(GISSgCpiC_O18.bin)
GISSgCpiC_O18.binstd<-GISSgCpiC_O18.bin/apply(GISSgCpiC_O18.bin,1,sd)
GISSgCpiC_MgCa.bin<-GISSgCpiC_MgCa.bin-rowMeans(GISSgCpiC_MgCa.bin)
GISSgCpiC_MgCa.binstd<-GISSgCpiC_MgCa.bin/apply(GISSgCpiC_MgCa.bin,1,sd)
GISSgCpiC_AMOC.bin<-GISSgCpiC_AMOC.bin-mean(GISSgCpiC_AMOC.bin)
GISSgCpiC_AMOC.binstd<-GISSgCpiC_AMOC.bin/sd(GISSgCpiC_AMOC.bin)
GISSgCpiC_NAO.bin<-GISSgCpiC_NAO.bin-mean(GISSgCpiC_NAO.bin)
GISSgCpiC_NAO.binstd<-GISSgCpiC_NAO.bin/sd(GISSgCpiC_NAO.bin)
GISSgCpiC_NAM.bin<-GISSgCpiC_NAM.bin-mean(GISSgCpiC_NAM.bin)
GISSgCpiC_NAM.binstd<-GISSgCpiC_NAM.bin/sd(GISSgCpiC_NAM.bin)

GISSgy3piC_O18.bin<-GISSgy3piC_O18.bin-rowMeans(GISSgy3piC_O18.bin)
GISSgy3piC_O18.binstd<-GISSgy3piC_O18.bin/apply(GISSgy3piC_O18.bin,1,sd)
GISSgy3piC_MgCa.bin<-GISSgy3piC_MgCa.bin-rowMeans(GISSgy3piC_MgCa.bin)
GISSgy3piC_MgCa.binstd<-GISSgy3piC_MgCa.bin/apply(GISSgy3piC_MgCa.bin,1,sd)
GISSgy3piC_AMOC.bin<-GISSgy3piC_AMOC.bin-mean(GISSgy3piC_AMOC.bin)
GISSgy3piC_AMOC.binstd<-GISSgy3piC_AMOC.bin/sd(GISSgy3piC_AMOC.bin)
GISSgy3piC_NAO.bin<-GISSgy3piC_NAO.bin-mean(GISSgy3piC_NAO.bin)
GISSgy3piC_NAO.binstd<-GISSgy3piC_NAO.bin/sd(GISSgy3piC_NAO.bin)
GISSgy3piC_NAM.bin<-GISSgy3piC_NAM.bin-mean(GISSgy3piC_NAM.bin)
GISSgy3piC_NAM.binstd<-GISSgy3piC_NAM.bin/sd(GISSgy3piC_NAM.bin)

GISSgTckLM_O18.bin<-GISSgTckLM_O18.bin-rowMeans(GISSgTckLM_O18.bin)
GISSgTckLM_O18.binstd<-GISSgTckLM_O18.bin/apply(GISSgTckLM_O18.bin,1,sd)
GISSgTckLM_MgCa.bin<-GISSgTckLM_MgCa.bin-rowMeans(GISSgTckLM_MgCa.bin)
GISSgTckLM_MgCa.binstd<-GISSgTckLM_MgCa.bin/apply(GISSgTckLM_MgCa.bin,1,sd)
GISSgTckLM_AMOC.bin<-GISSgTckLM_AMOC.bin-mean(GISSgTckLM_AMOC.bin)
GISSgTckLM_AMOC.binstd<-GISSgTckLM_AMOC.bin/sd(GISSgTckLM_AMOC.bin)
GISSgTckLM_NAO.bin<-GISSgTckLM_NAO.bin-mean(GISSgTckLM_NAO.bin)
GISSgTckLM_NAO.binstd<-GISSgTckLM_NAO.bin/sd(GISSgTckLM_NAO.bin)
GISSgTckLM_NAM.bin<-GISSgTckLM_NAM.bin-mean(GISSgTckLM_NAM.bin)
GISSgTckLM_NAM.binstd<-GISSgTckLM_NAM.bin/sd(GISSgTckLM_NAM.bin)

GISSgTKckLM_O18.bin<-GISSgTKckLM_O18.bin-rowMeans(GISSgTKckLM_O18.bin)
GISSgTKckLM_O18.binstd<-GISSgTKckLM_O18.bin/apply(GISSgTKckLM_O18.bin,1,sd)
GISSgTKckLM_MgCa.bin<-GISSgTKckLM_MgCa.bin-rowMeans(GISSgTKckLM_MgCa.bin)
GISSgTKckLM_MgCa.binstd<-GISSgTKckLM_MgCa.bin/apply(GISSgTKckLM_MgCa.bin,1,sd)
GISSgTKckLM_AMOC.bin<-GISSgTKckLM_AMOC.bin-mean(GISSgTKckLM_AMOC.bin)
GISSgTKckLM_AMOC.binstd<-GISSgTKckLM_AMOC.bin/sd(GISSgTKckLM_AMOC.bin)
GISSgTKckLM_NAO.bin<-GISSgTKckLM_NAO.bin-mean(GISSgTKckLM_NAO.bin)
GISSgTKckLM_NAO.binstd<-GISSgTKckLM_NAO.bin/sd(GISSgTKckLM_NAO.bin)
GISSgTKckLM_NAM.bin<-GISSgTKckLM_NAM.bin-mean(GISSgTKckLM_NAM.bin)
GISSgTKckLM_NAM.binstd<-GISSgTKckLM_NAM.bin/sd(GISSgTKckLM_NAM.bin)

GISSgTcsLM_O18.bin<-GISSgTcsLM_O18.bin-rowMeans(GISSgTcsLM_O18.bin)
GISSgTcsLM_O18.binstd<-GISSgTcsLM_O18.bin/apply(GISSgTcsLM_O18.bin,1,sd)
GISSgTcsLM_MgCa.bin<-GISSgTcsLM_MgCa.bin-rowMeans(GISSgTcsLM_MgCa.bin)
GISSgTcsLM_MgCa.binstd<-GISSgTcsLM_MgCa.bin/apply(GISSgTcsLM_MgCa.bin,1,sd)
GISSgTcsLM_AMOC.bin<-GISSgTcsLM_AMOC.bin-mean(GISSgTcsLM_AMOC.bin)
GISSgTcsLM_AMOC.binstd<-GISSgTcsLM_AMOC.bin/sd(GISSgTcsLM_AMOC.bin)
GISSgTcsLM_NAO.bin<-GISSgTcsLM_NAO.bin-mean(GISSgTcsLM_NAO.bin)
GISSgTcsLM_NAO.binstd<-GISSgTcsLM_NAO.bin/sd(GISSgTcsLM_NAO.bin)
GISSgTcsLM_NAM.bin<-GISSgTcsLM_NAM.bin-mean(GISSgTcsLM_NAM.bin)
GISSgTcsLM_NAM.binstd<-GISSgTcsLM_NAM.bin/sd(GISSgTcsLM_NAM.bin)

#combine anomalies 
piC_O18<-cbind(HadpiC_O18.bin,GISSgCpiC_O18.bin,GISSgy3piC_O18.bin)
colnames(piC_O18)<-c(rep("HadpiC",wt.HadpiC),rep("GISSgCpiC",wt.GISSgCpiC),rep("GISSgy3piC",wt.GISSgy3piC))
piC_MgCa<-cbind(HadpiC_MgCa.bin,GISSgCpiC_MgCa.bin,GISSgy3piC_MgCa.bin)
colnames(piC_MgCa)<-c(rep("HadpiC",wt.HadpiC),rep("GISSgCpiC",wt.GISSgCpiC),rep("GISSgy3piC",wt.GISSgy3piC))
piC_AMOC<-c(HadpiC_AMOC.bin,GISSgCpiC_AMOC.bin,GISSgy3piC_AMOC.bin)
piC_NAM<-c(HadpiC_NAM.bin,GISSgCpiC_NAM.bin,GISSgy3piC_NAM.bin)
piC_NAO<-c(HadpiC_NAO.bin,GISSgCpiC_NAO.bin,GISSgy3piC_NAO.bin)

LM_O18<-cbind(GISSgTckLM_O18.bin,GISSgTKckLM_O18.bin,GISSgTcsLM_O18.bin)
colnames(LM_O18)<-c(rep("GISSgTckLM",wt.GISSgTckLM),rep("GISSgTKckLM",wt.GISSgTKckLM),rep("GISSgTcsLM",wt.GISSgTcsLM))
LM_MgCa<-cbind(GISSgTckLM_MgCa.bin,GISSgTKckLM_MgCa.bin,GISSgTcsLM_MgCa.bin)
colnames(LM_MgCa)<-c(rep("GISSgTckLM",wt.GISSgTckLM),rep("GISSgTKckLM",wt.GISSgTKckLM),rep("GISSgTcsLM",wt.GISSgTcsLM))
LM_AMOC<-c(GISSgTckLM_AMOC.bin, GISSgTKckLM_AMOC.bin, GISSgTcsLM_AMOC.bin)
LM_NAM<-c(GISSgTckLM_NAM.bin, GISSgTKckLM_NAM.bin, GISSgTcsLM_NAM.bin)
LM_NAO<-c(GISSgTckLM_NAO.bin, GISSgTKckLM_NAO.bin, GISSgTcsLM_NAO.bin)

#save anomalies
saveRDS(piC_O18,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/piC_O18.",width))
saveRDS(piC_MgCa,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/piC_MgCa.",width))
saveRDS(piC_AMOC,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/piC_AMOC.",width))
saveRDS(piC_NAM,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/piC_NAM.",width))
saveRDS(piC_NAO,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/piC_NAO.",width))

saveRDS(LM_O18,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/LM_O18.",width))
saveRDS(LM_MgCa,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/LM_MgCa.",width))
saveRDS(LM_AMOC,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/LM_AMOC.",width))
saveRDS(LM_NAM,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/LM_NAM.",width))
saveRDS(LM_NAO,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/LM_NAO.",width))

#combine std anomalies 
piC_O18std<-cbind(HadpiC_O18.binstd,GISSgCpiC_O18.binstd,GISSgy3piC_O18.binstd)
colnames(piC_O18std)<-c(rep("HadpiC",wt.HadpiC),rep("GISSgCpiC",wt.GISSgCpiC),rep("GISSgy3piC",wt.GISSgy3piC))
piC_MgCastd<-cbind(HadpiC_MgCa.binstd,GISSgCpiC_MgCa.binstd,GISSgy3piC_MgCa.binstd)
colnames(piC_MgCastd)<-c(rep("HadpiC",wt.HadpiC),rep("GISSgCpiC",wt.GISSgCpiC),rep("GISSgy3piC",wt.GISSgy3piC))
piC_AMOCstd<-c(HadpiC_AMOC.binstd,GISSgCpiC_AMOC.binstd,GISSgy3piC_AMOC.binstd)
piC_NAMstd<-c(HadpiC_NAM.binstd,GISSgCpiC_NAM.binstd,GISSgy3piC_NAM.binstd)
piC_NAOstd<-c(HadpiC_NAO.binstd,GISSgCpiC_NAO.binstd,GISSgy3piC_NAO.binstd)

LM_O18std<-cbind(GISSgTckLM_O18.binstd,GISSgTKckLM_O18.binstd,GISSgTcsLM_O18.binstd)  
colnames(LM_O18std)<-c(rep("GISSgTckLM",wt.GISSgTckLM),rep("GISSgTKckLM",wt.GISSgTKckLM),rep("GISSgTcsLM",wt.GISSgTcsLM))
LM_MgCastd<-cbind(GISSgTckLM_MgCa.binstd,GISSgTKckLM_MgCa.binstd,GISSgTcsLM_MgCa.binstd)
colnames(LM_MgCastd)<-c(rep("GISSgTckLM",wt.GISSgTckLM),rep("GISSgTKckLM",wt.GISSgTKckLM),rep("GISSgTcsLM",wt.GISSgTcsLM))
LM_AMOCstd<-c(GISSgTckLM_AMOC.binstd,GISSgTKckLM_AMOC.binstd,GISSgTcsLM_AMOC.binstd)
LM_NAMstd<-c(GISSgTckLM_NAM.binstd,GISSgTKckLM_NAM.binstd,GISSgTcsLM_NAM.binstd)
LM_NAOstd<-c(GISSgTckLM_NAO.binstd,GISSgTKckLM_NAO.binstd,GISSgTcsLM_NAO.binstd)

#save std anomalies
saveRDS(piC_O18std,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/piC_O18std.",width))
saveRDS(piC_MgCastd,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/piC_MgCastd.",width))
saveRDS(piC_AMOCstd,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/piC_AMOCstd.",width))
saveRDS(piC_NAMstd,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/piC_NAMstd.",width))
saveRDS(piC_NAOstd,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/piC_NAOstd.",width))

saveRDS(LM_O18std,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/LM_O18std.",width))
saveRDS(LM_MgCastd,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/LM_MgCastd.",width))
saveRDS(LM_AMOCstd,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/LM_AMOCstd.",width))
saveRDS(LM_NAMstd,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/LM_NAMstd.",width))
saveRDS(LM_NAOstd,file=paste0("/Users/Casey_JISAO/Google Drive/R_documents/PSR/PSR/LM_NAOstd.",width))



