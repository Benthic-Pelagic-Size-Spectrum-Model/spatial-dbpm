Setup.Benthic<-function(run.in, mmin=-17, mmat=-16, mmax=9, A=64, alpha=0.75, mu_0=0.2, beta=-0.25, mu_s=0.1, epsilon=0.1, u_0=0.01, lambda=-0.75, K_det=0.2, R_det=0.2, Ex_det=0.2, pref_det=1, rep_method=1, initial_flag=F, ts_flag=F, fishing_flag=F, filename){

  #Assign Default Values
  if( missing(filename) ) stop("A Species filename must be given")
  if( !is.character(filename) ) stop("Please enter a valid value for filename: (character string)")
  if( !(rep_method==0 || rep_method ==1 || rep_method==2 || rep_method ==3) ) stop("Please enter a valid value for rep_flag: (0,1)")
  if( !(initial_flag==0 || initial_flag==1) ) stop("Please enter a logical value for initial_flag: (T,F)")
  if( !(ts_flag==0 || ts_flag==1) ) stop("Please enter a logical value for ts_flag: (T,F)")
  if( !(fishing_flag==0 || fishing_flag==1) ) stop("Please enter a logical value for fishing_flag: (T,F)")

  #Create list
  species<-new("benthic.params")
  
  species@filename<-filename
  species@speciestype<-"benthic"

  species@mmin<-mmin
  species@mmat<-mmat
  species@mmax<-mmax
  
  species@A<-A
  species@alpha<-alpha
  species@mu_0<-mu_0
  species@beta<-beta
  species@mu_s<-mu_s
  species@epsilon<-epsilon
  species@u_0<-u_0
  species@lambda<-lambda
  
  species@K_det<-K_det
  species@R_det<-R_det
  species@Ex_det<-Ex_det
  
  species@pref_det<-pref_det
  
  species@rep_method<-as.integer(rep_method)
  species@initial_flag<-as.logical(initial_flag)
  species@ts_flag<-as.logical(ts_flag)
  species@fishing_flag<-as.logical(fishing_flag)

  dir.create(paste(run.in@filename,"/",filename,sep=""),showWarnings=F)
  
  return(species)
}