#' Set up the detritus pool
#'
#' Builds a [detritus.params] object and creates the detritus output directory
#' under the run directory.
#'
#' @param run.in A [run.params] object.
#' @param w_0 Initial detritus biomass.
#' @param initial_flag,ts_flag Logical flags indicating whether
#'   initial-condition and time-series input files are supplied.
#' @param filename Character string naming the detritus pool (and its output
#'   directory). Required.
#'
#' @return A [detritus.params] object.
#' @seealso [Setup.Benthic()]
#' @export
Setup.Detritus<-function(run.in, w_0=0.6, initial_flag=FALSE, ts_flag=FALSE, filename){

  #Assign Default Values
  if(missing(filename)) stop("A Species filename must be given")
  if( !is.character(filename) ) stop("Please enter a valid value for filename: (character string)")
  if( !(initial_flag==0 || initial_flag==1) ) stop("Please enter a logical value for initial_flag: (T,F)")
  if( !(ts_flag==0 || ts_flag==1) ) stop("Please enter a logical value for ts_flag: (T,F)")

  #Create list
  species<-new("detritus.params")
  
  species@filename<-filename
  species@speciestype<-"detritus"
  
  species@w_0<-w_0

  species@initial_flag<-as.logical(initial_flag)
  species@ts_flag<-as.logical(ts_flag)

  dir.create(paste(run.in@filename,"/",filename,sep=""),showWarnings=FALSE)

  return(species)

}