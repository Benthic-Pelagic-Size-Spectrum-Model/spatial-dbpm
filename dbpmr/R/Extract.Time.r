Extract.Time<-function(species,time,section){

  tiny=1e-7
  
  if(missing(time)) stop("Time value is missing")
  if(!is.numeric(time)) stop("Time value must be numeric")
  if(!sum(abs(species@trange-time)<tiny)) stop("Time value is not in discretisation")
  
  extract<-new("timestep.data")
  mnum<-length(species@mrange)
  
  if(missing(section)){
    extract@data<-as.data.frame(lapply(species@uvals[abs(species@uvals$t-time)<tiny,2:(mnum+3)],as.numeric))
    names(extract@data)<-names(species@uvals[abs(species@uvals$t-time)<tiny,2:(mnum+3)])
    extract@spatial_dim<-species@run@spatial_dim
    extract@trange<-time
    extract@mrange<-species@mrange
    extract@xrange<-species@xrange
    extract@yrange<-species@yrange
  }

  return(extract)
}
