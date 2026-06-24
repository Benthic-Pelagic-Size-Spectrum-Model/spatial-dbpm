Read.In<-function(run,species,filename){
#This is a generic function for reading in all files

  #Check directory step
  if(!is.character(run)) stop("Please enter a valid run string")
  if(!is.character(species)) stop("Please enter a valid speciesname string")
  if(!file.exists(run)) stop("Run directory specified does not exist")
  if(!file.exists(paste(run,"/",species,sep=""))) stop("Species directory specified does not exist")
  
  #Specify number of each parameter
  no_run<-6
  no_grid<-18
  no_plankton<-10
  no_pelagic<-38
  no_benthic<-21
  no_detritus<-5
  
  #Get length of data section (i.e. number of timesteps)
  if(!file.exists(paste(run,"/parameters.txt",sep=""))) stop("Parameters file does not exist")
  temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=(6+no_run+7),nrows=3,header=F,strip.white=T,stringsAsFactors=F)
  toutmin<-temp[temp=='toutmin',2]
  toutmax<-temp[temp=='toutmax',2]
  toutstep<-temp[temp=='toutstep',2]
  no_data<-length(seq(toutmin,toutmax,toutstep))


  #Read in everything about species
  if(missing(filename)){

    #Check what type of species it is
    if(!file.exists(paste(run,"/",species,"/summary.txt",sep=""))) stop("Species summary file does not exist")
    temp<-read.csv(paste(run,"/",species,"/summary.txt",sep=""),sep=':',skip=1,nrows=1,header=F,strip.white=T,stringsAsFactors=F)
    speciestype<-temp[temp=='speciestype',2]
    
    
    #--------------------------------------------------------------------------#
    
    if(speciestype=="plankton"){
    
      #Setup a new plankton results object
      ans<-new("plankton.results")
      
      
      #Run params
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=4,nrows=no_run,header=F,strip.white=T,stringsAsFactors=F)
      
      ans@run@filename<-temp[temp=='filename',2]
      ans@run@no_pelagic<-as.integer(temp[temp=='no_pelagic',2])
      ans@run@no_benthic<-as.integer(temp[temp=='no_benthic',2])
      ans@run@spatial_dim<-as.integer(temp[temp=='spatial_dim',2])
      ans@run@coupled_flag<-as.logical(as.numeric(temp[temp=='coupled_flag',2]))
      ans@run@diff_method<-as.integer(temp[temp=='diff_method',2])
      
      #Grid params
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=(6+no_run),nrows=no_grid,header=F,strip.white=T,stringsAsFactors=F)
      
      ans@grid@mmin<-temp[temp=='mmin',2]
      ans@grid@mmax<-temp[temp=='mmax',2]
      ans@grid@mstep<-temp[temp=='mstep',2]
      ans@grid@moutstep<-temp[temp=='moutstep',2]
      
      ans@grid@t1<-temp[temp=='t1',2]
      ans@grid@tmax<-temp[temp=='tmax',2]
      ans@grid@tstep<-temp[temp=='tstep',2]
      ans@grid@toutmin<-temp[temp=='toutmin',2]
      ans@grid@toutmax<-temp[temp=='toutmax',2]
      ans@grid@toutstep<-temp[temp=='toutstep',2]
      
      ans@grid@xmin<-temp[temp=='xmin',2]
      ans@grid@xmax<-temp[temp=='xmax',2]
      ans@grid@xstep<-temp[temp=='xstep',2]
      ans@grid@xoutstep<-temp[temp=='xoutstep',2]
      
      ans@grid@ymin<-temp[temp=='ymin',2]
      ans@grid@ymax<-temp[temp=='ymax',2]
      ans@grid@ystep<-temp[temp=='ystep',2]
      ans@grid@youtstep<-temp[temp=='youtstep',2]
      
      #Calculate ranges
      ans@mrange<-seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)
      ans@trange<-seq(ans@grid@toutmin,ans@grid@toutmax,ans@grid@toutstep)
      ans@xrange<-seq(ans@grid@xmin,ans@grid@xmax,ans@grid@xoutstep)
      ans@yrange<-seq(ans@grid@ymin,ans@grid@ymax,ans@grid@youtstep)
    
      #U values
      if(!file.exists(paste(run,"/",species,"/results.txt",sep=""))) stop("Species results file does not exist")
      ans@uvals<-read.csv(paste(run,"/",species,"/results.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
      names(ans@uvals)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)))
      
      #Final U values
      ans@finaluvals<-read.csv(paste(run,"/",species,"/results.txt",sep=""),skip=(7+no_data),nrows=1,header=T,stringsAsFactors=F)
      names(ans@finaluvals)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@mstep)))
      
      #Summary file values
      temp<-read.csv(paste(run,"/",species,"/summary.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
        #Biomass
        ans@biomass<-temp[,4]
    
      #Plankton params
      #Can do this because there can only ever be one plankton system
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=(8+no_run+no_grid),nrows=no_plankton,header=F,strip.white=T,stringsAsFactors=F)
    
      ans@species@filename<-temp[temp=='filename',2]
      ans@species@speciestype<-temp[temp=='speciestype',2]
      
      ans@species@mmin<-as.numeric(temp[temp=='mmin',2])
      ans@species@mmax<-as.numeric(temp[temp=='mmax',2])
      
      ans@species@mu_0<-as.numeric(temp[temp=='mu_0',2])
      ans@species@beta<-as.numeric(temp[temp=='beta',2])
      ans@species@u_0<-as.numeric(temp[temp=='u_0',2])
      ans@species@lambda<-as.numeric(temp[temp=='lambda',2])
      
      ans@species@initial_flag<-as.logical(as.numeric(temp[temp=='initial_flag',2]))
      ans@species@ts_flag<-as.logical(as.numeric(temp[temp=='ts_flag',2]))
    }
    
    
    #--------------------------------------------------------------------------#
    
    if(speciestype=="pelagic"){

      ans<-new("pelagic.results")
      
      
      #Run params
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=4,nrows=no_run,header=F,strip.white=T,stringsAsFactors=F)
      
      ans@run@filename<-temp[temp=='filename',2]
      ans@run@no_pelagic<-as.integer(temp[temp=='no_pelagic',2])
      ans@run@no_benthic<-as.integer(temp[temp=='no_benthic',2])
      ans@run@spatial_dim<-as.integer(temp[temp=='spatial_dim',2])
      ans@run@coupled_flag<-as.logical(as.numeric(temp[temp=='coupled_flag',2]))
      ans@run@diff_method<-as.integer(temp[temp=='diff_method',2])
      
      #Grid params
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=(6+no_run),nrows=no_grid,header=F,strip.white=T,stringsAsFactors=F)
      
      ans@grid@mmin<-temp[temp=='mmin',2]
      ans@grid@mmax<-temp[temp=='mmax',2]
      ans@grid@mstep<-temp[temp=='mstep',2]
      ans@grid@moutstep<-temp[temp=='moutstep',2]
      
      ans@grid@t1<-temp[temp=='t1',2]
      ans@grid@tmax<-temp[temp=='tmax',2]
      ans@grid@tstep<-temp[temp=='tstep',2]
      ans@grid@toutmin<-temp[temp=='toutmin',2]
      ans@grid@toutmax<-temp[temp=='toutmax',2]
      ans@grid@toutstep<-temp[temp=='toutstep',2]
      
      ans@grid@xmin<-temp[temp=='xmin',2]
      ans@grid@xmax<-temp[temp=='xmax',2]
      ans@grid@xstep<-temp[temp=='xstep',2]
      ans@grid@xoutstep<-temp[temp=='xoutstep',2]
      
      ans@grid@ymin<-temp[temp=='ymin',2]
      ans@grid@ymax<-temp[temp=='ymax',2]
      ans@grid@ystep<-temp[temp=='ystep',2]
      ans@grid@youtstep<-temp[temp=='youtstep',2]
      
      #Calculate ranges
      ans@mrange<-seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)
      ans@trange<-seq(ans@grid@toutmin,ans@grid@toutmax,ans@grid@toutstep)
      ans@xrange<-seq(ans@grid@xmin,ans@grid@xmax,ans@grid@xoutstep)
      ans@yrange<-seq(ans@grid@ymin,ans@grid@ymax,ans@grid@youtstep)
  
      #U values 
      if(!file.exists(paste(run,"/",species,"/results.txt",sep=""))) stop("Species results file does not exist")
      ans@uvals<-read.csv(paste(run,"/",species,"/results.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
      names(ans@uvals)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)))
                                                                                   
      #Final U values
      ans@finaluvals<-read.csv(paste(run,"/",species,"/results.txt",sep=""),skip=(7+no_data),nrows=1,header=T,stringsAsFactors=F)
      names(ans@finaluvals)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@mstep)))
      
      #Growth values
      if(!file.exists(paste(run,"/",species,"/growth.txt",sep=""))) stop("Species growth file does not exist")
      ans@growth<-read.csv(paste(run,"/",species,"/growth.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
      names(ans@growth)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)))
      
      #Mortality values
      if(!file.exists(paste(run,"/",species,"/mortality.txt",sep=""))) stop("Species mortality file does not exist")
      ans@mortality<-read.csv(paste(run,"/",species,"/mortality.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
      names(ans@mortality)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)))
      
      #Predation values
 #     if(!file.exists(paste(run,"/",species,"/predation.txt",sep=""))) stop("Species predation file does not exist")
 #     ans@predation<-read.csv(paste(run,"/",species,"/predation.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
 #     names(ans@predation)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)))
      
      #Fishing values
      if(!file.exists(paste(run,"/",species,"/fishing.txt",sep=""))) stop("Species fishing file does not exist")
      ans@fishing<-read.csv(paste(run,"/",species,"/fishing.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
      names(ans@fishing)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)))
      
      #Summary file values
      temp<-read.csv(paste(run,"/",species,"/summary.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
        #Biomass
        ans@biomass<-temp[,4]
        #Plankton in
        ans@plabio<-temp[,5]
        #Pelagic in
        ans@pelbio<-temp[,6]
        #Benthic in
        ans@benbio<-temp[,7]
        #Predation out
        ans@predbio<-temp[,8]
        #Fishing out
        ans@fishbio<-temp[,9]
        #Reproduction
        ans@reproduction<-temp[,10]
        
      #Pelagic params
      #Select parameter values for this particular species
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=(10+no_run+no_grid+no_plankton),header=F,strip.white=T,stringsAsFactors=F)
      i<-which(temp[,2]==species)
      temp<-temp[i:(i+no_pelagic-1),]

      ans@species@filename<-temp[temp=='filename',2]
      ans@species@speciestype<-temp[temp=='speciestype',2]
      
      ans@species@mmin<-as.numeric(temp[temp=='mmin',2])
      ans@species@mmat<-as.numeric(temp[temp=='mmat',2])
      ans@species@mmax<-as.numeric(temp[temp=='mmax',2])
      
      ans@species@A<-as.numeric(temp[temp=='A',2])
      ans@species@alpha<-as.numeric(temp[temp=='alpha',2])
      ans@species@mu_0<-as.numeric(temp[temp=='mu_0',2])
      ans@species@beta<-as.numeric(temp[temp=='beta',2])
      ans@species@mu_s<-as.numeric(temp[temp=='mu_s',2])
      ans@species@epsilon<-as.numeric(temp[temp=='epsilon',2])
      ans@species@u_0<-as.numeric(temp[temp=='u_0',2])
      ans@species@lambda<-as.numeric(temp[temp=='lambda',2])
      
      ans@species@K_pla<-as.numeric(temp[temp=='K_pla',2])
      ans@species@R_pla<-as.numeric(temp[temp=='R_pla',2])
      ans@species@Ex_pla<-as.numeric(temp[temp=='Ex_pla',2])
      ans@species@K_pel<-as.numeric(temp[temp=='K_pel',2])
      ans@species@R_pel<-as.numeric(temp[temp=='R_pel',2])
      ans@species@Ex_pel<-as.numeric(temp[temp=='Ex_pel',2])
      ans@species@K_ben<-as.numeric(temp[temp=='K_ben',2])
      ans@species@R_ben<-as.numeric(temp[temp=='R_ben',2])
      ans@species@Ex_ben<-as.numeric(temp[temp=='Ex_ben',2])
      
      ans@species@pref_pla<-as.numeric(temp[temp=='pref_pel',2])
      ans@species@pref_pel<-as.numeric(temp[temp=='pref_pla',2])
      ans@species@pref_ben<-as.numeric(temp[temp=='pref_ben',2])
      
      ans@species@q_0<-as.numeric(temp[temp=='q_0',2])
      ans@species@sig<-as.numeric(temp[temp=='sig',2])
      ans@species@trunc<-as.numeric(temp[temp=='trunc',2])
      
      ans@species@prey<-as.numeric(temp[temp=='prey',2])
      ans@species@pred<-as.numeric(temp[temp=='pred',2])
      ans@species@comp<-as.numeric(temp[temp=='comp',2])
      ans@species@gamma_prey<-as.numeric(temp[temp=='gamma_pred',2])
      ans@species@gamma_pred<-as.numeric(temp[temp=='gamma_pred',2])
      ans@species@gamma_comp<-as.numeric(temp[temp=='gamma_comp',2])
      
      ans@species@rep_method<-as.integer(as.numeric(temp[temp=='rep_method',2]))
      ans@species@initial_flag<-as.logical(as.numeric(temp[temp=='initial_flag',2]))
      ans@species@ts_flag<-as.logical(as.numeric(temp[temp=='ts_flag',2]))
      ans@species@fishing_flag<-as.logical(as.numeric(temp[temp=='fishing_flag',2]))
      
    }


    #--------------------------------------------------------------------------#
      
    if(speciestype=="benthic"){
    
      ans<-new("benthic.results")
      
      
      #Run params
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=4,nrows=no_run,header=F,strip.white=T,stringsAsFactors=F)
      
      ans@run@filename<-temp[temp=='filename',2]
      ans@run@no_pelagic<-as.integer(temp[temp=='no_pelagic',2])
      ans@run@no_benthic<-as.integer(temp[temp=='no_benthic',2])
      ans@run@spatial_dim<-as.integer(temp[temp=='spatial_dim',2])
      ans@run@coupled_flag<-as.logical(as.numeric(temp[temp=='coupled_flag',2]))
      ans@run@diff_method<-as.integer(temp[temp=='diff_method',2])
      
      #Grid params
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=(6+no_run),nrows=no_grid,header=F,strip.white=T,stringsAsFactors=F)
      
      ans@grid@mmin<-temp[temp=='mmin',2]
      ans@grid@mmax<-temp[temp=='mmax',2]
      ans@grid@mstep<-temp[temp=='mstep',2]
      ans@grid@moutstep<-temp[temp=='moutstep',2]
      
      ans@grid@t1<-temp[temp=='t1',2]
      ans@grid@tmax<-temp[temp=='tmax',2]
      ans@grid@tstep<-temp[temp=='tstep',2]
      ans@grid@toutmin<-temp[temp=='toutmin',2]
      ans@grid@toutmax<-temp[temp=='toutmax',2]
      ans@grid@toutstep<-temp[temp=='toutstep',2]
      
      ans@grid@xmin<-temp[temp=='xmin',2]
      ans@grid@xmax<-temp[temp=='xmax',2]
      ans@grid@xstep<-temp[temp=='xstep',2]
      ans@grid@xoutstep<-temp[temp=='xoutstep',2]
      
      ans@grid@ymin<-temp[temp=='ymin',2]
      ans@grid@ymax<-temp[temp=='ymax',2]
      ans@grid@ystep<-temp[temp=='ystep',2]
      ans@grid@youtstep<-temp[temp=='youtstep',2]
      
      #Calculate ranges
      ans@mrange<-seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)
      ans@trange<-seq(ans@grid@toutmin,ans@grid@toutmax,ans@grid@toutstep)
      ans@xrange<-seq(ans@grid@xmin,ans@grid@xmax,ans@grid@xoutstep)
      ans@yrange<-seq(ans@grid@ymin,ans@grid@ymax,ans@grid@youtstep)
      
      #U values 
      if(!file.exists(paste(run,"/",species,"/results.txt",sep=""))) stop("Species results file does not exist")
      ans@uvals<-read.csv(paste(run,"/",species,"/results.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
      names(ans@uvals)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)))
      
      #Final U values
      ans@finaluvals<-read.csv(paste(run,"/",species,"/results.txt",sep=""),skip=(7+no_data),nrows=1,header=T,stringsAsFactors=F)
      names(ans@finaluvals)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@mstep)))
      
      #Growth values
      if(!file.exists(paste(run,"/",species,"/growth.txt",sep=""))) stop("Species growth file does not exist")
      ans@growth<-read.csv(paste(run,"/",species,"/growth.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
      names(ans@growth)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)))
      
      #Mortality values
      if(!file.exists(paste(run,"/",species,"/mortality.txt",sep=""))) stop("Species mortality file does not exist")
      ans@mortality<-read.csv(paste(run,"/",species,"/mortality.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
      names(ans@mortality)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)))
      
      #Predation values
      if(!file.exists(paste(run,"/",species,"/predation.txt",sep=""))) stop("Species predation file does not exist")
      ans@predation<-read.csv(paste(run,"/",species,"/predation.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
      names(ans@predation)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)))
      
      #Fishing values
      if(!file.exists(paste(run,"/",species,"/fishing.txt",sep=""))) stop("Species fishing file does not exist")
      ans@fishing<-read.csv(paste(run,"/",species,"/fishing.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
      names(ans@fishing)<-c("t","x","y",as.character(seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)))
      
      #Summary file values
      temp<-read.csv(paste(run,"/",species,"/summary.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
        #Biomass
        ans@biomass<-temp[,4]
        #Detritus in
        ans@detbio<-temp[,5]
        #Predation out
        ans@predbio<-temp[,6]
        #Fishing out
        ans@fishbio<-temp[,7]
        #Reproduction
        ans@reproduction<-temp[,8]
        
      #Benthic params
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=(10+no_run+no_grid+no_plankton),header=F,strip.white=T,stringsAsFactors=F)
      i<-which(temp[,2]==species)
      temp<-temp[i:(i+no_benthic-1),]
      
      ans@species@filename<-temp[temp=='filename',2]
      ans@species@speciestype<-temp[temp=='speciestype',2]

      ans@species@mmin<-as.numeric(temp[temp=='mmin',2])
      ans@species@mmat<-as.numeric(temp[temp=='mmat',2])
      ans@species@mmax<-as.numeric(temp[temp=='mmax',2])
    
      ans@species@A<-as.numeric(temp[temp=='A',2])
      ans@species@alpha<-as.numeric(temp[temp=='alpha',2])
      ans@species@mu_0<-as.numeric(temp[temp=='mu_0',2])
      ans@species@beta<-as.numeric(temp[temp=='beta',2])
      ans@species@mu_s<-as.numeric(temp[temp=='mu_s',2])
      ans@species@epsilon<-as.numeric(temp[temp=='epsilon',2])
      ans@species@u_0<-as.numeric(temp[temp=='u_0',2])
      ans@species@lambda<-as.numeric(temp[temp=='lambda',2])

      ans@species@K_det<-as.numeric(temp[temp=='K_det',2])
      ans@species@R_det<-as.numeric(temp[temp=='R_det',2])
      ans@species@Ex_det<-as.numeric(temp[temp=='Ex_det',2])
    
      ans@species@pref_det<-as.numeric(temp[temp=='pref_det',2])
    
      ans@species@rep_method<-as.integer(as.numeric(temp[temp=='rep_method',2]))
      ans@species@initial_flag<-as.logical(as.numeric(temp[temp=='initial_flag',2]))
      ans@species@ts_flag<-as.logical(as.numeric(temp[temp=='ts_flag',2]))
      ans@species@fishing_flag<-as.logical(as.numeric(temp[temp=='fishing_flag',2]))
    }
    
    
    #--------------------------------------------------------------------------#
    
    if(speciestype=="detritus"){
    
      ans<-new("detritus.results")
      
      #Run params
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=4,nrows=no_run,header=F,strip.white=T,stringsAsFactors=F)
      
      ans@run@filename<-temp[temp=='filename',2]
      ans@run@no_pelagic<-as.integer(temp[temp=='no_pelagic',2])
      ans@run@no_benthic<-as.integer(temp[temp=='no_benthic',2])
      ans@run@spatial_dim<-as.integer(temp[temp=='spatial_dim',2])
      ans@run@coupled_flag<-as.logical(as.numeric(temp[temp=='coupled_flag',2]))
      ans@run@diff_method<-as.integer(temp[temp=='diff_method',2])
      
      #Grid params
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=(6+no_run),nrows=no_grid,header=F,strip.white=T,stringsAsFactors=F)
      
      ans@grid@mmin<-temp[temp=='mmin',2]
      ans@grid@mmax<-temp[temp=='mmax',2]
      ans@grid@mstep<-temp[temp=='mstep',2]
      ans@grid@moutstep<-temp[temp=='moutstep',2]
      
      ans@grid@t1<-temp[temp=='t1',2]
      ans@grid@tmax<-temp[temp=='tmax',2]
      ans@grid@tstep<-temp[temp=='tstep',2]
      ans@grid@toutmin<-temp[temp=='toutmin',2]
      ans@grid@toutmax<-temp[temp=='toutmax',2]
      ans@grid@toutstep<-temp[temp=='toutstep',2]
      
      ans@grid@xmin<-temp[temp=='xmin',2]
      ans@grid@xmax<-temp[temp=='xmax',2]
      ans@grid@xstep<-temp[temp=='xstep',2]
      ans@grid@xoutstep<-temp[temp=='xoutstep',2]
      
      ans@grid@ymin<-temp[temp=='ymin',2]
      ans@grid@ymax<-temp[temp=='ymax',2]
      ans@grid@ystep<-temp[temp=='ystep',2]
      ans@grid@youtstep<-temp[temp=='youtstep',2]
      
      #Calculate ranges
      ans@mrange<-seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)
      ans@trange<-seq(ans@grid@toutmin,ans@grid@toutmax,ans@grid@toutstep)
      ans@xrange<-seq(ans@grid@xmin,ans@grid@xmax,ans@grid@xoutstep)
      ans@yrange<-seq(ans@grid@ymin,ans@grid@ymax,ans@grid@youtstep)
    
      #Summary file values              
      temp<-read.csv(paste(run,"/",species,"/summary.txt",sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
        #Biomass
        ans@biomass<-temp[,4]
        #Detritus in
        ans@detin<-temp[,5]
        #Predation out
        ans@detout<-temp[,6]

      #Detritus params
      temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=(10+no_run+no_grid+no_plankton),header=F,strip.white=T,stringsAsFactors=F)
      i<-which(temp[,2]==species)
      temp<-temp[i:(i+no_detritus-1),]

      ans@species@filename<-temp[temp=='filename',2]
      ans@species@speciestype<-temp[temp=='speciestype',2]
      
      ans@species@w_0<-as.numeric(temp[temp=='w_0',2])
      
      ans@species@initial_flag<-as.logical(as.numeric(temp[temp=='initial_flag',2]))
      ans@species@ts_flag<-as.logical(as.numeric(temp[temp=='ts_flag',2]))
    }
    
    
    #--------------------------------------------------------------------------#
    
  }
  
  if(!missing(filename)){
  
    if(!file.exists(paste(run,"/",species,"/",filename,sep=""))) stop(paste("File: ",filename," does not exist"))
    ans=new("singlefile.results")
    
    #Run params
    temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=4,nrows=no_run,header=F,strip.white=T,stringsAsFactors=F)
      
    ans@run@filename<-temp[temp=='filename',2]
    ans@run@no_pelagic<-as.integer(temp[temp=='no_pelagic',2])
    ans@run@no_benthic<-as.integer(temp[temp=='no_benthic',2])
    ans@run@spatial_dim<-as.integer(temp[temp=='spatial_dim',2])
    ans@run@coupled_flag<-as.logical(as.numeric(temp[temp=='coupled_flag',2]))
    ans@run@diff_method<-as.integer(temp[temp=='diff_method',2])
      
    #Grid params
    temp<-read.csv(paste(run,"/parameters.txt",sep=""),sep=':',skip=(6+no_run),nrows=no_grid,header=F,strip.white=T,stringsAsFactors=F)
      
    ans@grid@mmin<-temp[temp=='mmin',2]
    ans@grid@mmax<-temp[temp=='mmax',2]
    ans@grid@mstep<-temp[temp=='mstep',2]
    ans@grid@moutstep<-temp[temp=='moutstep',2]
    
    ans@grid@t1<-temp[temp=='t1',2]
    ans@grid@tmax<-temp[temp=='tmax',2]
    ans@grid@tstep<-temp[temp=='tstep',2]
    ans@grid@toutmin<-temp[temp=='toutmin',2]
    ans@grid@toutmax<-temp[temp=='toutmax',2]
    ans@grid@toutstep<-temp[temp=='toutstep',2]
    
    ans@grid@xmin<-temp[temp=='xmin',2]
    ans@grid@xmax<-temp[temp=='xmax',2]
    ans@grid@xstep<-temp[temp=='xstep',2]
    ans@grid@xoutstep<-temp[temp=='xoutstep',2]
    
    ans@grid@ymin<-temp[temp=='ymin',2]
    ans@grid@ymax<-temp[temp=='ymax',2]
    ans@grid@ystep<-temp[temp=='ystep',2]
    ans@grid@youtstep<-temp[temp=='youtstep',2]
    
    #Calculate ranges
    ans@mrange<-seq(ans@grid@mmin,ans@grid@mmax,ans@grid@moutstep)
    ans@trange<-seq(ans@grid@toutmin,ans@grid@toutmax,ans@grid@toutstep)
    ans@xrange<-seq(ans@grid@xmin,ans@grid@xmax,ans@grid@xoutstep)
    ans@yrange<-seq(ans@grid@ymin,ans@grid@ymax,ans@grid@youtstep)
    
    #Filetype
    ans@filetype<-as.character(read.csv(paste(run,"/",species,"/",filename,sep=""),nrows=1,header=F,stringsAsFactors=F))
    
    #Speciestype
    ans@speciestype<-as.character(read.csv(paste(run,"/",species,"/",filename,sep=""),skip=1,nrows=1,header=F,stringsAsFactors=F))
    
    #Data
    ans@data<-read.csv(paste(run,"/",species,"/",filename,sep=""),skip=4,nrows=no_data,header=T,stringsAsFactors=F)
  }
   
  return(ans)
}