args <- commandArgs(trailingOnly = T) #pull in all arguments from the qsub file
k<-as.numeric(args[1])

load(file="/mnt/iusers01/jw01/mbmhsby2/scratch/TEST/splitteddata.RData")
load(file="/mnt/iusers01/jw01/mbmhsby2/scratch/TEST/stationsloc.RData")

library(data.table)
mydist<-function(row,df2){
  
  # check the dimension of df2 and if it is Null return a data frame with na binded to it 
  if(dim(df2)[1]==0){
    return(cbind(rbindlist(list(df2,as.list(c(rep(NA,11))))),distance=NA))
  }
  else {
    rad<-pi/180
    a1<-row[["Lat"]]*rad
    a2<-row[["Lng"]]*rad
    b1<-df2$LAT*rad
    b2<-df2$LON*rad
    dlat<-b1-a1
    dlon<-b2-a2
    a<-sin(dlat/2)^2+cos(a1)*cos(b1)*sin(dlon/2)^2
    c <- 2 * asin(pmin(1,sqrt(a)))
    R<-6378.145
    dists<-R*c  
    return(cbind(df2[which.min(dists),],distance=min(dists)))
  }
}

df1<-df11[[as.numeric(args[1])]]
z<-cbind(df1,do.call(rbind,lapply(1:nrow(df1), function(x) mydist(df1[x,], df2[df2$Date%in%df1[x,]$Date&df2$join_time%in%df1[x,]$join_time&df2$nonmissingcount_tot==4,]))))  
z
