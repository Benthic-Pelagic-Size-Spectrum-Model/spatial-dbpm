#' Set up the plankton resource
#'
#' Builds a [plankton.params] object and creates the species output directory
#' under the run directory.
#'
#' @param run.in A [run.params] object.
#' @param mmin,mmax Minimum and maximum log-mass of the plankton spectrum.
#' @param mu_0 Mortality constant.
#' @param beta Mortality mass-scaling exponent.
#' @param u_0 Abundance constant.
#' @param lambda Spectrum slope.
#' @param initial_flag,ts_flag Logical flags indicating whether initial-condition
#'   and time-series input files are supplied (see [Setup.ts()]).
#' @param filename Character string naming the species (and its output
#'   directory). Required.
#'
#' @return A [plankton.params] object.
#' @export
Setup.Plankton<-function(run.in, mmin=-28,mmax=-14,mu_0=0.2, beta=-0.25, u_0=0.01, lambda=-1, initial_flag=FALSE, ts_flag=FALSE, filename){

  #Assign Default Values
  if(missing(filename)) stop("A Species filename must be given")
  if( !is.character(filename) ) stop("Please enter a valid value for filename: (character string)")
  if( !(initial_flag==0 || initial_flag==1) ) stop("Please enter a logical value for initial_flag: (T,F)")
  if( !(ts_flag==0 || ts_flag==1) ) stop("Please enter a logical value for ts_flag: (T,F)")

  #Create list
  species<-new("plankton.params")

  species@filename<-filename
  species@speciestype<-"plankton"

  species@mmin<-mmin
  species@mmax<-mmax

  species@mu_0<-mu_0
  species@beta<-beta
  species@u_0<-u_0
  species@lambda<-lambda
  
  species@initial_flag<-as.logical(initial_flag)
  species@ts_flag<-as.logical(ts_flag)

  dir.create(paste(run.in@filename,"/",filename,sep=""),showWarnings=FALSE)

  return(species)

}