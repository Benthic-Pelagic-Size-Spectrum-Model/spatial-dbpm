'Calculate.du'<-function(run.in,grid.in,plankton.in,pelagic.in,benthic.in,detritus.in,uvals){

  grid.in@tmax<-grid.in@tstep
  grid.in@toutstep<-grid.in@tstep
  grid.in@toutmax<-grid.in@tstep
  run.in@filename<-"Temp"

  
  #Ensure pelagic and benthic parameter inputs are lists
  pelagic.in<-c(pelagic.in)
  if(!missing(benthic.in)){
    benthic.in<-c(benthic.in)  
  }

  dir.create(run.in@filename,showWarnings=F)
    dir.create(paste(run.in@filename,"/",plankton.in@filename,sep=""),showWarnings=F)
   for(i in 1:length(pelagic.in)){
    dir.create(paste(run.in@filename,"/",pelagic.in[[i]]@filename,sep=""),showWarnings=F)
  }
  if(!missing(benthic.in)){
    for(i in 1:length(benthic.in)){
      dir.create(paste(run.in@filename,"/",benthic.in[[i]]@filename,sep=""),showWarnings=F)
    }
    dir.create(paste(run.in@filename,"/",detritus.in@filename,sep=""),showWarnings=F)
  }
  dir.create(paste(run.in@filename,"/Input",sep=""),showWarnings=F)

  plankton.in@initial_flag=F
  
  for(i in 1:length(pelagic.in)){
    pelagic.in[[i]]@initial_flag=T
  }
  
  if(!missing(benthic.in)){
    for(i in 1:length(benthic.in)){
      benthic.in[[i]]@initial_flag=T
    }
    detritus.in@initial_flag=T
  }

  plankton.in@ts_flag=F
  for(i in 1:length(pelagic.in)){
    pelagic.in[[i]]@ts_flag=F
  }
  if(!missing(benthic.in)){
    for(i in 1:length(benthic.in)){
      benthic.in[[i]]@ts_flag=F
    }
    detritus.in@ts_flag=F
  }


  #mat_pla<-t(as.matrix(uvals[1:lpla]))
  mat_pel<-t(as.matrix(c(pelagic.in[[1]]@u_0*exp(pelagic.in[[1]]@mmin)^-1,uvals[1:length(pelf)])))
  mat_spec<-t(as.matrix(c(pelagic.in[[2]]@u_0*exp(pelagic.in[[2]]@mmin)^-1,uvals[(length(pelf)+1):length(uvals)])))
  if(!missing(benthic.in)){
    mat_ben<-t(as.matrix(uvals[1:lben+lpla+lpel+lspec]))
    mat_det<-t(as.matrix(uvals[1+lben+lpel+lpla+lspec]))
  }


 # Setup.ts(plankton.in,run.in,grid.in,mat=mat_pla)
  Setup.ts(pelagic.in[[1]],run.in,grid.in,mat=mat_pel)
  Setup.ts(pelagic.in[[2]],run.in,grid.in,mat=mat_spec)
  if(!missing(benthic.in)){
    Setup.ts(benthic.in[[1]],run.in,grid.in,mat=mat_ben)
    Setup.ts(detritus.in,run.in,grid.in,mat=mat_det)
  }
  
  #Calling the Sizespectrum function which calls the C code #
  filenames<-SizeSpectrum(run.in,grid.in,plankton.in,pelagic.in,benthic.in,detritus.in)

  
#  if(run.in@no_benthic!=0){
#    fnames<-vector(length=(run.in@no_pelagic+run.in@no_benthic+2))
#  }
#  else{
#    fnames<-vector(length=(run.in@no_pelagic+1))
#  }
#  for(i in 1:length(fnames)){
#    fnames[i]<-Read.In(filenames[1],filenames[i+1])
#  }


  
 
 # pla<-Read.In(filenames[1],filenames[2])
  pel1<-Read.In(filenames[1],filenames[3])
  pel2<-Read.In(filenames[1],filenames[4])
 # ben<-Read.In(filenames[1],filenames[5])
 # det<-Read.In(filenames[1],filenames[6])
   
  #plaf<-as.vector(pla@finaluvals[1:lpla+which(mrange==plankton.in@mmin)+2],"numeric")
  pel1f<-as.vector(pel1@finaluvals[1:(lpel-1)+which(mrange==as.character(pelagic.in[[1]]@mmin))+3],"numeric")
  pel2f<-as.vector(pel2@finaluvals[1:(lspec-1)+which(mrange==as.character(pelagic.in[[2]]@mmin))+3],"numeric") 

 # benf<-as.vector(ben@finaluvals[1:lben+which(mrange==as.character(benthic.in@mmin))+2],"numeric") 
 # detf<-as.vector(det@biomass[2],"numeric")

  newuvals<-c(pel1f,pel2f)
  
  #newuvals<-c(plaf,pel1f,pel2f,benf,detf)

  du<-round((newuvals-uvals),signif=16)/grid.in@tstep

  return(du)
}
