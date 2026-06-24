.Last.lib <- function(libpath)
{
  
  removeMethod("show","run.params")
  removeMethod("show","grid.params")
  removeMethod("show","pelagic.params")
  removeMethod("show","benthic.params")
  removeMethod("show","plankton.params")
  removeMethod("show","detritus.params")  


  removeClass("run.params")
  removeClass("grid.params")
  removeClass("pelagic.params")
  removeClass("benthic.params")
  removeClass("plankton.params")
  removeClass("detritus.params")
  removeClass("pelagic.results")
  removeClass("benthic.results")
  removeClass("plankton.results")
  removeClass("detritus.results")
  removeClass("timestep.data")
  
  library.dynam.unload("dbpmr",libpath)
  cat("dbpmr 0.1-1 unloaded\n")
  
}