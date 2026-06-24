Setup.Plankton<-function(run.in, mmin=-28,mmax=-14,mu_0=0.2, beta=-0.25, u_0=0.01, lambda=-1, initial_flag=F, ts_flag=F, filename){

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

  dir.create(paste(run.in@filename,"/",filename,sep=""),showWarnings=F)

  return(species)

}