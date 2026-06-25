#' Extract the size spectrum at a single time
#'
#' Pulls the size spectrum for one output time from a species results object and
#' returns it as a [timestep.data] object ready for plotting.
#'
#' @param species A species results object (e.g. from [Read.In()]).
#' @param time The output time to extract; must lie on the output grid.
#' @param section Reserved for future use; currently ignored.
#'
#' @return A [timestep.data] object.
#' @seealso [Average.Time()], [Plot.Spectrum()]
#' @export
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
