#' Write initial-condition or time-series input for a species
#'
#' Generates the `_ts.txt` input file used to seed a species' initial size
#' spectrum (when `ts_flag` is `FALSE`) or its full time series (when `ts_flag`
#' is `TRUE`). The values may be supplied as a function of mass/time/space, as
#' an R matrix, or read from a CSV data file. Exactly one of `func`, `mat` or
#' `dataname` must be given.
#'
#' A plankton time series can instead be supplied in memory as `intercept` and
#' `slope` arrays, with no input file written. In that case the function returns
#' the plankton object carrying the arrays, which must be assigned and passed on
#' to [SizeSpectrum()]: `pl <- Setup.ts(pl, run, grid, intercept = a, slope = b)`.
#'
#' Starting values can also be held in memory, with `initial`, for the plankton,
#' a pelagic or benthic species, or the detritus. Again no file is written and
#' the returned object must be assigned:
#' `pe <- Setup.ts(pe, run, grid, initial = u0)`.
#'
#' @param species A species parameter object (e.g. from [Setup.Pelagic()]).
#' @param run.in A [run.params] object.
#' @param grid.in A [grid.params] object.
#' @param func A function of `(m, t, x, y)` returning abundances.
#' @param mat A matrix of size-spectrum values (rows = spatial/time cells,
#'   columns = the species' mass bins).
#' @param dataname Path to a CSV file whose first line is the file type
#'   (`"intercept"`, `"interslope"` or `"spectrum"`).
#' @param intercept,slope Plankton only (`ts_flag = TRUE`): numeric arrays with
#'   dimensions `[x, y, layer]` giving, for each spatial cell and layer, the
#'   plankton spectrum `u(m) = intercept * exp(slope * m)` (`m` is natural-log
#'   body mass; `intercept` is the abundance density at `m = 0`, not its
#'   logarithm). For a 1-D run a matrix `[x, layer]` and for a 0-D run a vector
#'   with one value per layer are also accepted.
#' @param layer Used with `intercept`/`slope`: for each model time point
#'   (`seq(0, tmax, tstep)`) the layer it uses, as whole numbers from 1. May be
#'   omitted when there is one layer per time point, or a single layer.
#' @param time Used with `intercept`/`slope` instead of `layer`: the start time
#'   (in model years) of each layer; each time point uses the latest layer that
#'   has started.
#' @param initial Starting values held in memory (`initial_flag = TRUE`,
#'   `ts_flag = FALSE`). For the plankton or a pelagic or benthic species: the
#'   abundance density in each size class of that species
#'   (`seq(mmin, mmax, by = mstep)`), as a vector (the same in every cell) or an
#'   array `[x, y, size class]`. For the detritus: the biomass, as a single
#'   number or a matrix `[x, y]`.
#'
#' @details For an initial condition (`ts_flag = FALSE`) the supplied `mat` /
#'   data may have either a single row, applied uniformly to every spatial cell,
#'   or one row per spatial cell (`x * y` rows) ordered with x outermost and y
#'   innermost (`row = (k - 1) * ny + l`). For a full time series
#'   (`ts_flag = TRUE`) it must have `t * x * y` rows in the same x/y order
#'   nested within time.
#'
#' @return Invisibly `NULL` when writing an input file under the run's `Input`
#'   directory. When `intercept` and `slope`, or `initial`, are given, invisibly
#'   the object with the arrays attached (assign it).
#' @seealso [Setup.Rep()], [Setup.fishing()]
#' @export
`Setup.ts`<-
function(species, run.in, grid.in, func, mat, dataname, intercept, slope, layer, time, initial){
#func is a function of the variables m, x, y, and z

#Starting values held in memory: nothing is written
if(!missing(initial)){
  if(!missing(func) || !missing(mat) || !missing(dataname) || !missing(intercept) || !missing(slope) || !missing(layer) || !missing(time))
    stop("initial cannot be combined with func, mat, dataname, intercept, slope, layer or time")
  return(invisible(initial_array(species, grid.in, initial)))
}

#Plankton time series held in memory as intercept/slope arrays: nothing is written
if(!missing(intercept) || !missing(slope)){
  if(!missing(func) || !missing(mat) || !missing(dataname)) stop("intercept/slope cannot be combined with func, mat or dataname")
  return(invisible(plankton_ts_arrays(species, grid.in,
                                      if(missing(intercept)) NULL else intercept,
                                      if(missing(slope)) NULL else slope,
                                      if(missing(layer)) NULL else layer,
                                      if(missing(time)) NULL else time)))
}
if(!missing(layer) || !missing(time)) stop("layer and time are only used with intercept and slope")
#dataname is a filename containing either a vector of intercepts, two vectors of intercepts and slopes or a full matrix of ts values

#Input checking stuff
if(species@initial_flag==FALSE && species@ts_flag==FALSE) stop ("Starting values or time series are not to be specified since initial_flag==F or ts_flag=F")

if(missing(func) && missing(dataname) && missing(mat)) stop("A function, matrix or data for the initial step or time series must be given")
if(!missing(func) && !missing(dataname) && !missing(mat)) stop("Only one of a function or data for the initial step step or time series must be given")
if(!missing(func) && !missing(dataname) && missing(mat)) stop("Only one of a function or data for the initial step step or time series must be given")
if(!missing(func) && missing(dataname) && !missing(mat)) stop("Only one of a function or data for the initial step step or time series must be given")
if(missing(func) && !missing(dataname) && !missing(mat)) stop("Only one of a function or data for the initial step step or time series must be given")

#Grid stuff
mass<-grid_seq(grid.in@mmin,grid.in@mmax,grid.in@mstep)
xrange<-grid_seq(grid.in@xmin,grid.in@xmax,grid.in@xstep)
yrange<-grid_seq(grid.in@ymin,grid.in@ymax,grid.in@ystep)
trange<-seq(0,grid.in@tmax,grid.in@tstep)

#Species stuff
filename<-paste(species@filename,"_ts.txt",sep="")

if(species@speciestype!="detritus"){
  spmin<-species@mmin
  spmax<-species@mmax
}
if(species@speciestype=="detritus"){
  spmin<-0
  spmax<-0
  mass<-0
}

m<-length(mass)
x<-length(xrange)
y<-length(yrange)
t<-length(trange)

#Calculation of initial values /time series using function
if(!missing(func)){
  if(species@ts_flag==TRUE){
    temp<-matrix(0,nrow=(t*x*y),ncol=m)
    for(j in 1:t){
      for(k in 1:x){
        for(l in 1:y){
          for(i in mass_index(mass,spmin):mass_index(mass,spmax)){
            temp[((j-1)*x*y)+((k-1)*y)+l,i]=round(func(mass[i],trange[j],xrange[k],yrange[l]),16)
          }
        }
      }
    }
  }
  if(species@ts_flag==FALSE){
  temp<-matrix(0,nrow=(x*y),ncol=m)
    for(k in 1:x){
      for(l in 1:y){
        for(i in mass_index(mass,spmin):mass_index(mass,spmax)){
          temp[((k-1)*y)+l,i]=round(func(mass[i],trange[1],xrange[k],yrange[l]),16)
        }
      }
    }
  }
}

#Calcualtion of initial values / time series using data file
if(!missing(dataname)){
#read in csv files or intercepts and/or slopes
#Each csv file must contain a header row with 'intercept', 'interslope' or 'spectrum'
  if(!is.character(dataname)) stop("Please enter a valid filename string")
  if(!file.exists(dataname)) stop("Filename specified does not exist")

  filetype<-as.character(read.csv(dataname,nrows=1,header=FALSE,strip.white=TRUE,stringsAsFactors=FALSE))
  dat<-read.csv(dataname,skip=1,header=FALSE,strip.white=TRUE,stringsAsFactors=FALSE)
  dat<-as.matrix(dat)
  
  if(species@ts_flag==TRUE){
    #Check whether the number of lines in data file matches the number of time steps previously specified
    if(length(dat[,1])!=t*x*y){
      stop("Incorrect number of rows in data file")
    }
    
    #Setup input storage array
    temp<-matrix(0,nrow=(t*x*y),ncol=m)
    
    if(species@speciestype!="detritus"){
    
      if(filetype=="intercept"){
        for(i in mass_index(mass,spmin):mass_index(mass,spmax)){
          temp[,i]=round(dat[,1]*exp(species@lambda*mass[i]),16)
        }
      }
    
      if(filetype=="interslope"){
        if(length(dat[1,])<2) stop("Insufficient columns in data file")
        for(i in mass_index(mass,spmin):mass_index(mass,spmax)){
          temp[,i]=round(dat[,1]*exp(dat[,2]*mass[i]),16)
        }
      }
    
      if(filetype=="spectrum"){
        if(length(dat[1,])!=length(grid_seq(spmin,spmax,grid.in@mstep))) stop("Incorrect number of columns in data file")
        for(i in mass_index(mass,spmin):mass_index(mass,spmax)){
          temp[,i]=round(dat[,(i-mass_index(mass,spmin)+1)],16)
        }
      }
      if(!(filetype=="intercept" || filetype=="interslope" || filetype=="spectrum")) stop("Unknown filetype entered")
    
    }
    
    if(species@speciestype=="detritus"){
      temp[,1]=round(dat[,1],16)
    }
    
  }
  
  if(species@ts_flag==FALSE){
    #Initial conditions over the x*y spatial cells. The data file must hold
    #either one row (applied uniformly to every cell) or x*y rows, one per
    #spatial cell in (x outer, y inner) order: row = (k-1)*y + l.
    ncell<-x*y
    nr<-nrow(dat)
    if(nr!=1 && nr!=ncell) stop("Data file must have 1 row (uniform) or x*y rows (one per spatial cell)")
    src<-if(nr==1) rep(1,ncell) else seq_len(ncell)
    temp<-matrix(0,nrow=ncell,ncol=m)

    if(species@speciestype!="detritus"){

      if(filetype=="intercept"){
        for(cc in 1:ncell){
          for(i in mass_index(mass,spmin):mass_index(mass,spmax)){
            temp[cc,i]=round(dat[src[cc],1]*exp(species@lambda*mass[i]),16)
          }
        }
      }

      if(filetype=="interslope"){
        if(length(dat[1,])<2) stop("Insufficient columns in data file")
        for(cc in 1:ncell){
          for(i in mass_index(mass,spmin):mass_index(mass,spmax)){
            temp[cc,i]=round(dat[src[cc],1]*exp(dat[src[cc],2]*mass[i]),16)
          }
        }
      }

      if(filetype=="spectrum"){
        if(length(dat[1,])!=length(grid_seq(spmin,spmax,grid.in@mstep))) stop("Incorrect number of columns in data file")
        for(cc in 1:ncell){
          for(i in mass_index(mass,spmin):mass_index(mass,spmax)){
            temp[cc,i]=round(dat[src[cc],(i-mass_index(mass,spmin)+1)],16)
          }
        }
      }
      if(!(filetype=="intercept" || filetype=="interslope" || filetype=="spectrum")) stop("Unknown filetype entered")

    }

    if(species@speciestype=="detritus"){
      for(cc in 1:ncell){
        temp[cc,1]=round(dat[src[cc],1],16)
      }
    }
  }

}


#Calcualtion of initial values / time series using R matrix
if(!missing(mat)){
#read in matrix containing full size spectra
  if(!is.matrix(mat)) stop("mat must be of type matrix")
  
  if(species@ts_flag==TRUE){
    #Check whether the number of lines in the matrix matches the number of time steps previously specified
    if(length(mat[,1])!=t*x*y){
      stop("Incorrect number of rows in data file")
    }

    #Setup input storage array
    temp<-matrix(0,nrow=(t*x*y),ncol=m)
    
    if(species@speciestype!="detritus"){

      if(length(mat[1,])!=length(grid_seq(spmin,spmax,grid.in@mstep))) stop("Incorrect number of columns in matrix")
        for(i in mass_index(mass,spmin):mass_index(mass,spmax)){
          temp[,i]=round(mat[,(i-mass_index(mass,spmin)+1)],16)
        }
    }
    
    if(species@speciestype=="detritus"){
      temp[,1]=round(mat[,1],16)
    }
    
  }
  
  if(species@ts_flag==FALSE){
    #Initial conditions over the x*y spatial cells: mat must have either one row
    #(applied uniformly) or x*y rows, one per spatial cell in (x outer, y inner)
    #order: row = (k-1)*y + l.
    ncell<-x*y
    nr<-nrow(mat)
    if(nr!=1 && nr!=ncell) stop("mat must have 1 row (uniform) or x*y rows (one per spatial cell)")
    src<-if(nr==1) rep(1,ncell) else seq_len(ncell)

    #Setup input storage array
    temp<-matrix(0,nrow=ncell,ncol=m)

    if(species@speciestype!="detritus"){
      if(length(mat[1,])!=length(grid_seq(spmin,spmax,grid.in@mstep))) stop("Incorrect number of columns in data file")
      for(cc in 1:ncell){
        for(i in mass_index(mass,spmin):mass_index(mass,spmax)){
          temp[cc,i]=round(mat[src[cc],(i-mass_index(mass,spmin)+1)],16)
        }
      }
    }

    if(species@speciestype=="detritus"){
      for(cc in 1:ncell){
        temp[cc,1]=round(mat[src[cc],1],16)
      }
    }
  }

}

#Create Input directory
dir.create(paste(run.in@filename,"/Input",sep=""),showWarnings=FALSE)

#Write full table all at once
write.table(temp,file.path(run.in@filename,"Input",filename),append=FALSE,row.names=FALSE,col.names=FALSE,sep=",")

}
