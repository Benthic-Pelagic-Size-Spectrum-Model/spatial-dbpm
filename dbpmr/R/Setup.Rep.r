#' Write reproduction input for a species
#'
#' Generates the `_rep_ts.txt` reproduction time-series input file for a pelagic
#' or benthic species whose `rep_method` is `1`. The values may be supplied as a
#' function of time/space or read from a CSV data file. Exactly one of `func` or
#' `dataname` must be given.
#'
#' @param species A pelagic or benthic parameter object.
#' @param run.in A [run.params] object.
#' @param grid.in A [grid.params] object.
#' @param func A function of `(t, x, y)` returning reproduction values.
#' @param dataname Path to a CSV file of reproduction values (first line is the
#'   file type).
#'
#' @return Invisibly `NULL`; called for the side effect of writing the input
#'   file under the run's `Input` directory.
#' @seealso [Setup.ts()], [Setup.fishing()]
#' @export
`Setup.Rep`<-
function(species, run.in, grid.in, func, dataname){
#func is a function of the variables m, x, y, and z
#dataname is a filename containing either a vector of intercepts, two vectors of intercepts and slopes or a full matrix of ts values

#Input checking stuff
if(species@speciestype=='plankton' || species@speciestype=='detritus') stop ("Reproduction can only be applied to pelagic or benthic systems")
if(species@rep_method!=1) stop ("Reproduction values are not to be specified since rep_method!=1")

if(missing(func) && missing(dataname)) stop("A function or data for the reproduction values must be given")
if(!missing(func) && !missing(dataname)) stop("Only one of a function or data for the reproduction values must be given")


#Grid stuff
xrange<-seq(grid.in@xmin,grid.in@xmax,grid.in@xstep)
yrange<-seq(grid.in@ymin,grid.in@ymax,grid.in@ystep)
trange<-seq(0,grid.in@tmax,grid.in@tstep)

#Species stuff
filename<-paste(run.in@filename,"_",species@filename,"_rep_ts.txt",sep="")
spmin<-species@mmin

x<-length(xrange)
y<-length(yrange)
t<-length(trange)

#Calculation of reproduction values using function
if(missing(dataname)){
  temp<-matrix(0,nrow=(t*x*y),ncol=1)
  for(j in 1:t){
    for(k in 1:x){
      for(l in 1:y){
          temp[((j-1)*x*y)+((k-1)*y)+l,1]=signif(func(trange[j],xrange[k],yrange[l]))
      }
    }
  }
}

#Calcualtion of fishing time series using data file
if(missing(func)){

  if(!is.character(dataname)) stop("Please enter a valid filename string")
  if(!file.exists(dataname)) stop("Filename specified does not exist")

  #Each csv file must contain a header row with 'reproduction'
  filetype<-as.character(read.csv(dataname,nrows=1,header=FALSE,strip.white=TRUE,stringsAsFactors=FALSE))
  dat<-read.csv(dataname,skip=1,header=FALSE,strip.white=TRUE,stringsAsFactors=FALSE)
  dat<-as.matrix(dat)

  #Check whether the number of lines in data file matches the number of time steps previously specified
  if(length(dat[,1])!=t*x*y) stop("Incorrect number of rows in data file")
  #Checks whether the number of columns in data file matches the size discretisation for the species
  if(length(dat[1,])!=1) stop("Incorrect number of columns in data file")

  temp<-matrix(0,nrow=(t*x*y),ncol=1)
  temp=signif(dat)
}

#Create Input directory
dir.create(paste(run.in@filename,"/Input",sep=""),showWarnings=FALSE)

#Write full table all at once
write.table(temp,paste(run.in@filename,"/input/",filename,sep=""),append=FALSE,row.names=FALSE,col.names=FALSE,sep=",")

}
