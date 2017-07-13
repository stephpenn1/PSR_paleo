#test

testData<-HadpiC_O18_PPmarSed
s<-seq(1,999)

for(i in 1:nrow(testData)) {
  for (j in 1:ncol(testData)) {
    sampleRes<-(resIndex[i])/100
    winLength<- round(1/sampleRes)
    weight<-hamming(winLength)
    
    for (k in 1:length(begin)) {
      begin<-seq(1,length(testData)-(length(weight)+1))
      end<-seq(length(weight),length(testData))
      testData_smooth<-seq((length(weight)+1)/2,length(s) - (length(weight)-1)/2)
      
      temp<-testData[begin[k]:end[k]]
      testData_smooth[k]<-sum(temp*weight)/sum(weight)
    }
    
    
  }
}

plot(s,testData[1,], type = "l")
#title(main = "GISSgCpiC_O18 pseudoproxy - smooth")
lines(s[((length(weight)+1)/2):(length(s)-(length(weight)-1)/2)],testData_smooth,col="red")
