Setup.Pelagic<-function(run.in, mmin=-14, mmat=7, mmax=14, A=640, alpha=0.82, mu_0=0.2, beta=-0.25, mu_s=0.1, epsilon=0.1, u_0=0.01, lambda=-1, K_pla, R_pla, Ex_pla, K_pel=0.2, R_pel=0.2, Ex_pel=0.3, K_ben=0.1, R_ben=0.2, Ex_ben=0.4, pref_pla=1, pref_pel=1, pref_ben=1, q_0=log(100), sig=log(10), trunc=2, prey=0, pred=0, comp=0.1, gamma_prey=0.33, gamma_pred=0.33, gamma_comp=0.75, rep_method=1, initial_flag=F, ts_flag=F, fishing_flag=F, filename){

  #Assign Default Values
  if( missing(filename) ) stop("A Species filename must be given")
  if( !is.character(filename) ) stop("Please enter a valid value for filename: (character string)")
  if( !(rep_method==0 || rep_method ==1 || rep_method==2 || rep_method ==3) ) stop("Please enter a valid value for rep_method: (0,1)")
  if( !(initial_flag==0 || initial_flag==1) ) stop("Please enter a logical value for initial_flag: (T,F)")
  if( !(ts_flag==0 || ts_flag==1) ) stop("Please enter a logical value for ts_flag: (T,F)")
  if( !(fishing_flag==0 || fishing_flag==1) ) stop("Please enter a logical value for fishing_flag: (T,F)")

  if(missing(K_pla))  K_pla<-K_pel
  if(missing(R_pla))  R_pla<-R_pel
  if(missing(Ex_pla)) Ex_pla<-Ex_pel

  
  #Create list
  species<-new("pelagic.params")

  species@filename<-filename
  species@speciestype<-"pelagic"

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

  species@K_pla<-K_pla
  species@R_pla<-R_pla
  species@Ex_pla<-Ex_pla
  species@K_pel<-K_pel
  species@R_pel<-R_pel
  species@Ex_pel<-Ex_pel
  species@K_ben<-K_ben
  species@R_ben<-R_ben
  species@Ex_ben<-Ex_ben
  
  species@pref_pla<-pref_pla
  species@pref_pel<-pref_pel
  species@pref_ben<-pref_ben

  species@q_0<-q_0
  species@sig<-sig
  species@trunc<-trunc
  
  species@prey<-prey
  species@pred<-pred
  species@comp<-comp
  species@gamma_prey<-gamma_pred
  species@gamma_pred<-gamma_pred
  species@gamma_comp<-gamma_comp

  species@rep_method<-as.integer(rep_method)
  species@initial_flag<-as.logical(initial_flag)
  species@ts_flag<-as.logical(ts_flag)
  species@fishing_flag<-as.logical(fishing_flag)

  dir.create(paste(run.in@filename,"/",filename,sep=""),showWarnings=F)
  
  return(species)

}