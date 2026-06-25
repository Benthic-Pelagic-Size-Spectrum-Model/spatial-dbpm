#' Write fishing input for a species
#'
#' Generates the `_fishing_ts.txt` fishing-mortality time-series input file for a
#' pelagic or benthic species whose `fishing_flag` is `TRUE`. The values may be
#' supplied as a function of mass/time/space or read from a CSV data file.
#' Exactly one of `func` or `dataname` must be given.
#'
#' @param species A pelagic or benthic parameter object with `fishing_flag` set.
#' @param run.in A [run.params] object.
#' @param grid.in A [grid.params] object.
#' @param func A function of `(m, t, x, y)` returning fishing mortality.
#' @param dataname Path to a CSV file of fishing values (first line is the file
#'   type).
#'
#' @return Invisibly `NULL`; called for the side effect of writing the input
#'   file under the run's `Input` directory.
#' @seealso [Setup.ts()], [Setup.Rep()]
#' @export
`Setup.fishing`<-
function(species, run.in, grid.in, func, dataname){
#func is a function of the variables m, x, y, and z
#dataname is a filename containing either a vector of intercepts, two vectors of intercepts and slopes or a full matrix of ts values

#Input checking stuff
if(species@speciestype=='plankton' || species@speciestype=='detritus') stop ("Fishing can only be applied to pelagic or benthic systems")
if(species@fishing_flag==FALSE) stop ("Fishing values are not to be specified since fishing_flag==F")

if(missing(func) && missing(dataname)) stop("A function or data for the initial step must be given")
if(!missing(func) && !missing(dataname)) stop("Only one of a function or data for the initial step must be given")


#Grid stuff
mass<-seq(grid.in@mmin,grid.in@mmax,grid.in@mstep)
xrange<-seq(grid.in@xmin,grid.in@xmax,grid.in@xstep)
yrange<-seq(grid.in@ymin,grid.in@ymax,grid.in@ystep)
trange<-seq(0,grid.in@tmax,grid.in@tstep)

#Species stuff
filename<-paste(species@filename,"_fishing_ts.txt",sep="")
spmin<-species@mmin
spmax<-species@mmax

m<-length(mass)
x<-length(xrange)
y<-length(yrange)
t<-length(trange)

#Calculation of initial values /time series using function
if(missing(dataname)){
  temp<-matrix(0,nrow=(t*x*y),ncol=m)
  for(j in 1:t){
    for(k in 1:x){
      for(l in 1:y){
        for(i in which(mass==spmin):which(mass==spmax)){
          temp[((j-1)*x*y)+((k-1)*y)+l,i]=signif(func(mass[i],trange[j],xrange[k],yrange[l]))
        }
      }
    }
  }
}

#Calcualtion of fishing time series using data file
if(missing(func)){

  #Each csv file must contain a header row with 'fishing'
  if(!is.character(dataname)) stop("Please enter a valid filename string")
  if(!file.exists(dataname)) stop("Filename specified does not exist")

  filetype<-as.character(read.csv(dataname,nrows=1,header=FALSE,strip.white=TRUE,stringsAsFactors=FALSE))
  dat<-read.csv(dataname,skip=1,header=FALSE,strip.white=TRUE,stringsAsFactors=FALSE)
  dat<-as.matrix(dat)

  #Check whether the number of lines in data file matches the number of time steps previously specified
  if(length(dat[,1])!=t*x*y) stop("Incorrect number of rows in data file")
  #Checks whether the number of columns in data file matches the size discretisation for the species
  if(length(dat[1,])!=length(seq(spmin,spmax,grid.in@mstep))) stop("Incorrect number of columns in data file")

  temp<-matrix(0,nrow=(t*x*y),ncol=m)
  for(i in which(mass==spmin):which(mass==spmax)){
    temp[,i]=signif(dat[,(i-which(mass==spmin)+1)])
  }
}

#Create Input directory
dir.create(paste(run.in@filename,"/Input",sep=""),showWarnings=FALSE)

#Write full table all at once
write.table(temp,file.path(run.in@filename,"Input",filename),append=FALSE,row.names=FALSE,col.names=FALSE,sep=",")

}