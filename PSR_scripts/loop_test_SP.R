#test
library(signal)
testData<-HadpiC_O18_PPmarSed
s<-seq(1,999)

res<-metadataO18[,4]
resIndex<-res[sed]
smooth_matrix<-matrix(data = NA, 40, 999)

for(i in 1:nrow(testData)) {
  
  
  sampleRes<-(resIndex[i])/100
  winLength<- round(1/sampleRes)
  weight<-hamming(winLength)
  
  begin<-seq(1,ncol(testData)-(length(weight)+1))
  end<-seq(length(weight),ncol(testData))
  testData_smooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)
  
  for (j in 1:length(begin)) {
    temp<-testData[i,begin[j]:end[j]]
    avg<-sum(temp*weight)/sum(weight)
    smooth_matrix[i,testData_smooth[j]]<-avg #need to shift over
  }
  
}

plot(s,testData[1,], type = "l")
#title(main = "GISSgCpiC_O18 pseudoproxy - smooth")
lines(smooth_matrix[1,],col="red")
