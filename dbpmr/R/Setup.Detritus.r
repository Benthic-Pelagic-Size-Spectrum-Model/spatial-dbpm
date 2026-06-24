Setup.Detritus<-function(run.in, w_0=0.6, initial_flag=F, ts_flag=F, filename){

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

  dir.create(paste(run.in@filename,"/",filename,sep=""),showWarnings=F)

  return(species)

}