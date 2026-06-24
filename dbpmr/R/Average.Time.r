Average.Time<-function(species,time.lim){

  dt<-species@grid@tout
  tmin<-min(species@trange)
  tmax<-max(species@trange)

  #-------------#  
  # Error Steps #
  #-------------#
  
  if(missing(time.lim)){ time.lim=vector() ; time.lim[1]=tmin ; time.lim[2]=tmax}
  if(!is.numeric(time.lim)) stop("Time limits must be numeric")
  if(length(time.lim)>2) stop("Too many time limits entered")
  if(length(time.lim==1)) time.lim[2]=time.lim[1]

  time.start=time.lim[1]
  time.end=time.lim[2]
  if(time.start>tmax) stop("Start Time is too high")
  if(time.end>tmax) stop("End Time is too high")
  if(time.start<tmin) stop("Start Time is too low")
  if(time.end<tmin) stop("End Time is too low")
  if(time.start>time.end) stop("Start Time is greater than End Time")
  if(!sum(species@trange==time.start)) stop("Time.start value is not in discretisation")
  if(!sum(species@trange==time.end)) stop("Time.end value is not in discretisation")

  average<-new("timestep.data")

  average@spatial.dim<-as.integer(species@run[species@run=='spatial.dim',2])
  average@trange<-seq(time.start,time.end,dt)
  average@mrange<-species@mrange
  average@xrange<-species@xrange
  average@yrange<-species@yrange

  mnum<-length(species@mrange)
  jmin<-which(species@trange==time.start)
  jmax<-which(species@trange==time.end)
  hold<-matrix(0,length(species@xrange)*length(species@yrange),(mnum+2))
  for(i in jmin:jmax){
    temp<-as.matrix(species@data[species@data$t==species@trange[i],2:(mnum+3)])
    hold<-hold+temp
  }
  average@data<-as.data.frame(hold/length(average@trange))
  
  return(average)
}

  