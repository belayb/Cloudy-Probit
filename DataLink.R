#' This function computes the distance between each individuals location data and
#' weather station locations. The function then binds the individual location to the
#' nearest station. 
#' The input for the functions are 
#' @parm row - a data frame of the individual location data at each day/time. 
#'             the row data is expected to have the following column names: 
#'             id, Date, time, Lat, Lng
#' @parm df2 - a data frame of the weather station location data along with 
#'             the weather information. df2 is expected to have the following
#'             column names. USAF (sation id), Date, time,LAT,LON,
#'             and the column names for the weather infromation of interest. 

mydist<-function(row,df2){
  
  # check the dimension of df2 and if it is Null return a data frame with na binded to it. 
  # That is, there is no station to join at that time. 
  # Make sure to change the number of columns (11 here according to the dimension of your
  # weather data) 
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

z<-cbind(df1, 
           do.call(rbind,lapply(1:nrow(df1), function(x) 
           mydist(df1[x,], df2[df2$Date%in%df1[x,]$Date&df2$time%in%df1[x,]$time,]))))  

# Once you make the linkage, you can set the distance treshold to use. E.g you may say 
# stations more than 100km away are not trustworthy. 

