Animate<-function(species,time.lim,time.step,mass.lim,x.lim,y.lim,...){

  dt<-species@grid@toutstep
  tmin<-min(species@trange)
  tmax<-max(species@trange)

  #-------------#
  # Error Steps #
  #-------------#

  if(missing(time.lim)) {time.lim=vector() ; time.lim[1]=tmin ; time.lim[2]=tmax}
  if(missing(time.step)) {time.step=dt}
  if(!is.numeric(time.lim)) stop("Time limits must be numeric")
  if(!is.numeric(time.step)) stop("Time limits must be numeric")
  if(length(time.lim)>2) stop("Too many time limits entered")
  if(length(time.step)>1) stop("Too many time steps entered")
  if(length(time.lim)==1) time.lim[2]=tmax

  time.start=time.lim[1]
  time.end=time.lim[2]
  if(time.start>tmax) stop("Start Time is too high")
  if(time.end>tmax) stop("End Time is too high")
  if(time.start<tmin) stop("Start Time is too low")
  if(time.end<tmin) stop("End Time is too low")
  if(time.start>time.end) stop("Start Time is greater than End Time")
  if(!sum(species@trange==time.start)) stop("Time.start value is not in discretisation")
  if(!sum(species@trange==time.end)) stop("Time.end value is not in discretisation")
  if((time.step/dt)-floor(time.step/dt)!=0) stop("Time.step value is not a multiple of grid tstep value")
  
  #Define time values to be animated
  tvals<-seq(time.start,time.end,time.step)
  
  #Plot graphs sequentially
  for(t in tvals){
    temp<-Extract.Time(species,t)
    Plot.Spectrum(temp,mass.lim,x.lim,y.lim,main=t,...)
  }
  
}
    
