#' Set up the grid discretisation
#'
#' Builds a [grid.params] object describing the mass, time and spatial
#' discretisation of a run. Spatial dimensions not used by `run` are collapsed
#' to a single trivial cell.
#'
#' @param run A [run.params] object (used to decide the spatial dimension).
#' @param mmin,mmax,mstep,moutstep Minimum, maximum, step and output step of the
#'   log-mass grid. The defaults are the canonical FishMIP/DBPM log10 grid
#'   expressed in natural log: `10^-12` to `10^6` g in steps of `0.1` log10
#'   units (`0.1 * log(10)`), so the consumer minimum `10^-3` g lands exactly on
#'   a grid point.
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
Setup.Grid<-function(run, mmin=-12*log(10),mmax=6*log(10),mstep=0.1*log(10),moutstep=0.1*log(10),t1=0,tmax=1,tstep=(1/365),toutmin=0,toutmax=tmax,toutstep=(73/365),xmin=0,xmax=1000,xstep=50,xoutstep=50,ymin=0,ymax=1000,ystep=50,youtstep=50, land_mask=1){

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

#' Set up the grid discretisation
#'
#' Builds a [grid.params] object describing the mass, time and spatial
#' discretisation of a run. Spatial dimensions not used by `run` are collapsed
#' to a single trivial cell.
#'
#' @param run Terra SpatRaster object that contains 0 values (where land or invalid
#'   cells are present) and 1 values (where water cells are present for modelling).
#'
#' @return A [grid.params] object.
#' @seealso [Setup.Run()]
#' @export
Setup.Grid <- function (
    rast, 
    mmin=-12*log(10), 
    mmax=6*log(10), 
    mstep=0.1*log(10), 
    moutstep=0.1*log(10), 
    t1=0, 
    tmax=1, 
    tstep=1/365, 
    toutmin=0, 
    toutmax=tmax, 
    toutstep=73/365
){

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
    
    if (crs(rast) == "") {
        stop("rast terra SpatRaster must have a valid CRS attribute `crs(rast)`.")
    }

    if (!is.lonlat(rast)) {
        stop("rast terra SpatRaster must use a projected Coordinate Reference System")
    }

    # Check that all values in `rast` are 0s or 1s designating land and water cells.
    if (any(!values(rast) %in% c(0, 1))) {
        stop("rast terra SpatRaster must only contain 0 values (land cells) and 1 values (water cells).")
    }

    ny <- ncol(rast)
    grid_mask_matrix <- t(as.matrix(rast, wide=TRUE))[, ny:1]
    grid.params@spatial_grid_matrix <- grid_mask_matrix

    grid.params@x_res <- res(rast)[1]
    gird.params@y_res <- res(rast)[2]


    return(grid.params)
}