#' Set up the grid discretisation
#'
#' Builds a [grid.params] object describing the mass, time and spatial
#' discretisation of a run. Spatial dimensions not used by `run` are collapsed
#' to a single trivial cell.
#'
#' @param run A [run.params] object (used to decide the spatial dimension).
#' @param mmin,mmax,mstep,moutstep Minimum, maximum, step and output step of the
#'   log-mass grid.
#' @param t1,tmax,tstep Start, end and step of the time grid (in years).
#' @param toutmin,toutmax,toutstep Minimum, maximum and step of the time output
#'   grid (in years).
#' @param xmin,xmax,xstep,xoutstep Minimum, maximum, step and output step of the
#'   x-space grid.
#' @param ymin,ymax,ystep,youtstep Minimum, maximum, step and output step of the
#'   y-space grid.
#'
#' @return A [grid.params] object.
#' @seealso [Setup.Run()]
#' @export
Setup.Grid<-function(run, mmin=-28,mmax=14,mstep=0.2,moutstep=1,t1=0,tmax=1,tstep=(1/365),toutmin=0,toutmax=tmax,toutstep=(73/365),xmin=0,xmax=1000,xstep=50,xoutstep=50,ymin=0,ymax=1000,ystep=50,youtstep=50){

#------------------#  
# Create Grid List #
#------------------#
  grid.params<-new("grid.params")

  grid.params@mmin<-mmin
  grid.params@mmax<-mmax
  grid.params@mstep<-mstep
  grid.params@moutstep<-moutstep

  grid.params@t1<-t1
  grid.params@tmax<-tmax         #in years
  grid.params@tstep<-tstep       #in years
  grid.params@toutmin<-toutmin   #in years
  grid.params@toutmax<-toutmax   #in years
  grid.params@toutstep<-toutstep #in years
  
  grid.params@xmin<-xmin
  grid.params@xmax<-xmax
  grid.params@xstep<-xstep
  grid.params@xoutstep<-xoutstep

  grid.params@ymin<-ymin
  grid.params@ymax<-ymax
  grid.params@ystep<-ystep
  grid.params@youtstep<-youtstep
  
  if(run@spatial_dim==1){
    grid.params@ymin<-0
    grid.params@ymax<-0
    grid.params@ystep<-1
    grid.params@youtstep<-1
  }
  
  if(run@spatial_dim==0){
    grid.params@xmin<-0
    grid.params@xmax<-0
    grid.params@xstep<-1
    grid.params@xoutstep<-1

    grid.params@ymin<-0
    grid.params@ymax<-0
    grid.params@ystep<-1
    grid.params@youtstep<-1
  }  
  
  return(grid.params)
}