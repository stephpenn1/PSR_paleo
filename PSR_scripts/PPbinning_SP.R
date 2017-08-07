#Stephanie Pennington | PSR summer research
#Pseudoproxy binning according to sample resolution
#Created 7-18-17

setwd("/Users/SP/Desktop/PSR_paleo/PSR_data/pseudoproxy/aggregate/")

sampleRes<-metadataO18[,4]

HadpiC_O18_PPagg<-read.csv("HadpiC_O18_PPagg.csv")
GISSgCpiC_O18_PPagg<-read.csv("GISSgCpiC_O18_PPagg.csv")
GISSgy3piC_O18_PPagg<-read.csv("GISSgy3piC_O18_PPagg.csv")
GISSgTckLM_O18_PPagg<-read.csv("GISSgTckLM_O18_PPagg.csv")
GISSgTKckLM_O18_PPagg<-read.csv("GISSgTKckLM_O18_PPagg.csv")
GISSgTcsLM_O18_PPagg<-read.csv("GISSgTcsLM_O18_PPagg.csv")
HadpiC_MgCa_PPagg<-read.csv("HadpiC_MgCa_PPagg.csv")
GISSgCpiC_MgCa_PPagg<-read.csv("GISSgCpiC_MgCa_PPagg.csv")
GISSgy3piC_MgCa_PPagg<-read.csv("GISSgy3piC_MgCa_PPagg.csv")
GISSgTckLM_MgCa_PPagg<-read.csv("GISSgTckLM_MgCa_PPagg.csv")
GISSgTKckLM_MgCa_PPagg<-read.csv("GISSgTKckLM_MgCa_PPagg.csv")
GISSgTcsLM_MgCa_PPagg<-read.csv("GISSgTcsLM_MgCa_PPagg.csv")

