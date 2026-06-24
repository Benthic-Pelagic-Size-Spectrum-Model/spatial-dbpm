`SizeSpectrum` <-
function(run.in, grid.in, plankton.in, pelagic.in, benthic.in, detritus.in){

#run.in is an object of class of run.params (2 strings and 6 integers relating to the entire run process)

#grid.in is an object of class of grid.params (18 floats relating to the grid discretisation)

#plankton.in is an object of class plankton.params (2 strings and 6 floats)

#pelagic.in is a concatenated list of objects of class pelagic.params
#pelagic.params (2 strings, 31 floats and 1 integer)

#benthic.in is a concatenated list of objects of class benthic.params
#benthic.params (2 strings, 15 floats and 1 integer)

#detritus.in is an object of class detritus.params (2 strings and 1 float)

tiny<-1e-7


#------------------------------#
# Check systems to be modelled #
#------------------------------#
  
  #If benthic system is not to be modelled
  if(missing(benthic.in) && missing(detritus.in) && run.in@no_benthic==0) {benthic.in<-0 ; detritus.in<-0}
  else{
    if((missing(benthic.in) || missing(detritus.in)) && run.in@no_benthic!=0) stop("No benthic/detritus parameters entered")
    else{
      if(!(missing(benthic.in) && missing(detritus.in)) && run.in@no_benthic==0) {benthic.in<-0 ; detritus.in<-0 ; warning("Benthic/detritus parameters entered. No Benthic or Detrital system will be simulated")}
    }
  }
  
  
  #Ensure pelagic and benthic parameter inputs are lists
  pelagic.in<-c(pelagic.in)
  benthic.in<-c(benthic.in)  
  
#-------------------------------------------#
# Assign spatial null values if appropriate #
#-------------------------------------------#

  if(run.in@spatial_dim==0){
      grid.in@xmin<-0
      grid.in@xmax<-0
      grid.in@xstep<-1
      grid.in@xoutstep<-1
  
      grid.in@ymin<-0
      grid.in@ymax<-0
      grid.in@ystep<-1
      grid.in@youtstep<-1
  }
  if(run.in@spatial_dim==1){
      grid.in@ymin<-0
      grid.in@ymax<-0
      grid.in@ystep<-1
      grid.in@youtstep<-1
  }


#--------------------------------------------------------#
#  Construct vectors and matrices to be passed to C code #
#--------------------------------------------------------#

  #Run Parameters (Integer vector)

    #Setup run transfer vector
    run.params<-c(as.integer(run.in@no_pelagic), as.integer(run.in@no_benthic),
                  as.integer(run.in@spatial_dim),
                  as.integer(run.in@coupled_flag),
                  as.integer(run.in@diff_method))


  
  #Grid Parameters (Float vector)
  
    #Calculate maximum and minimum grid discretisation
    mmin<-min(plankton.in@mmin)
    mmax<-max(plankton.in@mmax)
    for(i in 1:length(pelagic.in)){
      mmin<-min(mmin,pelagic.in[[i]]@mmin)
      mmax<-max(mmax,pelagic.in[[i]]@mmax)
    }
    if(run.in@no_benthic!=0){
      for(i in 1:length(benthic.in)){
        mmin<-min(mmin,benthic.in[[i]]@mmin)
        mmax<-max(mmax,benthic.in[[i]]@mmax)
      }
    }
    grid.in@mmin<-mmin
    grid.in@mmax<-mmax
  
    #Setup grid transfer vector
    grid.params<-c(grid.in@mmin, grid.in@mmax, grid.in@mstep, grid.in@moutstep,
                  grid.in@t1, grid.in@tmax, grid.in@tstep, grid.in@toutmin, grid.in@toutmax, grid.in@toutstep,
                  grid.in@xmin, grid.in@xmax, grid.in@xstep, grid.in@xoutstep,
                  grid.in@ymin, grid.in@ymax, grid.in@ystep, grid.in@youtstep)

  
  #Plankton Parameters (Float vector)
  
    #Setup plankton transfer vector
    pla.params<-c(plankton.in@mmin, plankton.in@mmax,
                  plankton.in@mu_0, plankton.in@beta, plankton.in@u_0, plankton.in@lambda)
  
  #Pelagic Parameters (Float vector)
  
    #Setup pelagic transfer vector
    pel.params<-NULL
    for(i in 1:length(pelagic.in)){
      pel.params<-c(pel.params,
                    pelagic.in[[i]]@mmin, pelagic.in[[i]]@mmat, pelagic.in[[i]]@mmax,
                    pelagic.in[[i]]@A, pelagic.in[[i]]@alpha, pelagic.in[[i]]@mu_0, pelagic.in[[i]]@beta, pelagic.in[[i]]@mu_s, pelagic.in[[i]]@epsilon, pelagic.in[[i]]@u_0, pelagic.in[[i]]@lambda,
                    pelagic.in[[i]]@K_pla, pelagic.in[[i]]@R_pla, pelagic.in[[i]]@Ex_pla, pelagic.in[[i]]@K_pel, pelagic.in[[i]]@R_pel, pelagic.in[[i]]@Ex_pel, pelagic.in[[i]]@K_ben, pelagic.in[[i]]@R_ben, pelagic.in[[i]]@Ex_ben,
                    pelagic.in[[i]]@pref_pla, pelagic.in[[i]]@pref_pel, pelagic.in[[i]]@pref_ben,
                    pelagic.in[[i]]@q_0, pelagic.in[[i]]@sig, pelagic.in[[i]]@trunc,
                    pelagic.in[[i]]@prey, pelagic.in[[i]]@pred, pelagic.in[[i]]@comp, pelagic.in[[i]]@gamma_prey, pelagic.in[[i]]@gamma_pred, pelagic.in[[i]]@gamma_comp)
    }

  
  #Benthic & Detritus Parameters (Float vector)
  
    #Setup benthic transfer vector
    if(run.in@no_benthic!=0){
      ben.params<-NULL
      for(i in 1:length(benthic.in)){
        ben.params<-c(ben.params,
                      benthic.in[[i]]@mmin, benthic.in[[i]]@mmat, benthic.in[[i]]@mmax,
                      benthic.in[[i]]@A, benthic.in[[i]]@alpha, benthic.in[[i]]@mu_0, benthic.in[[i]]@beta, benthic.in[[i]]@mu_s, benthic.in[[i]]@epsilon, benthic.in[[i]]@u_0, benthic.in[[i]]@lambda,
                      benthic.in[[i]]@K_det, benthic.in[[i]]@R_det, benthic.in[[i]]@Ex_det,
                      benthic.in[[i]]@pref_det)
      }
  
      #Setup detritus transfer vector
      det.params<-c(detritus.in@w_0)  
    }
    if(run.in@no_benthic==0){
      ben.params<-0
      det.params<-0
    }
  
  #Filenames (Character vector)
  
    pla.names<-c(plankton.in@filename)
  
    pel.names<-NULL
    for(i in 1:length(pelagic.in)){
      pel.names<-c(pel.names,pelagic.in[[i]]@filename)
    }
  
    ben.names<-NULL
    det.names<-NULL
    if(run.in@no_benthic!=0){
      for(i in 1:length(benthic.in)){
        ben.names<-c(ben.names,benthic.in[[i]]@filename)
      }
      det.names<-c(detritus.in@filename)
    }
    
    #Setup filenames transfer vector
    names.params<-c(run.in@filename, pla.names, pel.names, ben.names, det.names)

  #Flag parameters (Integer vector)

    pla.flags<-c(as.integer(plankton.in@initial_flag), as.integer(plankton.in@ts_flag))
  
    pel.flags<-NULL
    for(i in 1:length(pelagic.in)){
      pel.flags<-c(pel.flags,as.integer(pelagic.in[[i]]@rep_method), as.integer(pelagic.in[[i]]@initial_flag), as.integer(pelagic.in[[i]]@ts_flag), as.integer(pelagic.in[[i]]@fishing_flag))
    }
  
    ben.flags<-NULL
    det.flags<-NULL
    if(run.in@no_benthic!=0){
      for(i in 1:length(benthic.in)){
        ben.flags<-c(ben.flags,as.integer(benthic.in[[i]]@rep_method), as.integer(benthic.in[[i]]@initial_flag), as.integer(benthic.in[[i]]@ts_flag), as.integer(benthic.in[[i]]@fishing_flag))
      }
      det.flags<-c(as.integer(detritus.in@initial_flag), as.integer(detritus.in@ts_flag))
    }
    
    #Setup flags transfer vector
    flags.params<-c(pla.flags, pel.flags, ben.flags, det.flags)  
    

#----------------------#
# Error Checking Steps #
#----------------------#

  # Check sufficient parameters have been inputed
  if(sum(is.infinite(grid.params))!=0) stop("Not enough Grid Values Inputed")
  if(sum(is.infinite(pla.params))!=0) stop("Not enough Plankton Parameter Values Inputed")
  if(sum(is.infinite(pel.params))!=0) stop("Not enough Pelagic Parameter Values Inputed")
  if(run.in@no_benthic!=0){
    if(sum(is.infinite(ben.params))!=0) stop("Not enough Benthic Parameter Values Inputed")
    if(sum(is.infinite(det.params))!=0) stop("Not enough Detritus Parameter Values Inputed")
  }

  if((length(pelagic.in))!=run.in@no_pelagic) stop("Incorrect number of Pelagic species Inputed")
  if(run.in@no_benthic!=0){
    if((length(benthic.in))!=run.in@no_benthic) stop("Incorrect number of Benthic species Inputed")
  }

  #Check uniqueness of filenames
  if(run.in@no_benthic!=0){
    if(length(unique(c(pla.names,pel.names,ben.names,det.names)))!=(run.in@no_pelagic+run.in@no_benthic+2)) stop("All Species must have different filenames")
  }
  if(run.in@no_benthic==0){
    if(length(unique(c(pla.names,pel.names)))!=(run.in@no_pelagic+1)) stop("All Species must have different filenames")
  }
  
  
  #Check initial input files exist
  if((plankton.in@initial_flag==1 || plankton.in@ts_flag==1) && !file.exists(paste(run.in@filename,"/Input/",plankton.in@filename,"_ts.txt",sep=""))) stop("A required Plankton initial/ts file does not exist")
  
  for(i in 1:length(pelagic.in)){
    if((pelagic.in[[i]]@initial_flag==1 || pelagic.in[[i]]@ts_flag==1) && !file.exists(paste(run.in@filename,"/Input/",pelagic.in[[i]]@filename,"_ts.txt",sep=""))) stop("A required Pelagic initial/ts file does not exist")
    if( pelagic.in[[i]]@fishing_flag==1 && !file.exists(paste(run.in@filename,"/Input/",pelagic.in[[i]]@filename,"_fishing_ts.txt",sep=""))) stop("A required Pelagic fishing file does not exist")
    if( pelagic.in[[i]]@rep_method==1 && !file.exists(paste(run.in@filename,"/Input/",pelagic.in[[i]]@filename,"_rep_ts.txt",sep=""))) stop("A required Pelagic reproduction file does not exist")
    
  }
  if(run.in@no_benthic!=0){
    for(i in 1:length(benthic.in)){
      if((benthic.in[[i]]@initial_flag==1 || benthic.in[[i]]@ts_flag==1) && !file.exists(paste(run.in@filename,"/Input/",benthic.in[[i]]@filename,"_ts.txt",sep=""))) stop("A required Benthic initial/ts file does not exist")
      if( benthic.in[[i]]@fishing_flag==1 && !file.exists(paste(run.in@filename,"/Input/",benthic.in[[i]]@filename,"_fishing_ts.txt",sep=""))) stop("A required Benthic fishing file does not exist")
      if( benthic.in[[i]]@rep_method==1 && !file.exists(paste(run.in@filename,"/Input/",benthic.in[[i]]@filename,"_rep_ts.txt",sep=""))) stop("A required Benthic reproduction file does not exist")
    }
    if((detritus.in@initial_flag==1 || detritus.in@ts_flag==1) && !file.exists(paste(run.in@filename,"/Input/",detritus.in@filename,"_ts.txt",sep=""))) stop("A required Detritus initial/ts file does not exist")
  }
  
  
  #Check consistency in input files
  for(i in 1:length(pelagic.in)){
    if(pelagic.in[[i]]@ts_flag==1 && pelagic.in[[i]]@rep_method==1) warning("Cannot specify both time series and reproduction values. Rep values will be ignored")
  }
  if(run.in@no_benthic!=0){
    for(i in 1:length(benthic.in)){
      if(benthic.in[[i]]@ts_flag==1 && benthic.in[[i]]@rep_method==1) warning("Cannot specify both time series and reproduction values. Rep values will be ignored")
    }
  }
  
  
  #Check Grid Discretisation
  if(grid.in@toutstep < grid.in@tstep) stop("Cannot output results. Time Print increments are too small")
  if(grid.in@tmax < grid.in@tstep || grid.in@tmax < grid.in@toutstep || (grid.in@toutmax-grid.in@toutmin) < grid.in@toutstep) stop("Cannot output results. Time Range is too small")
  if(abs((grid.in@toutstep/grid.in@tstep)-round(grid.in@toutstep/grid.in@tstep)) > tiny) stop("Cannot output results. Time Print increment is not a multiple of Time Step")

  if(grid.in@moutstep < grid.in@mstep) stop("Cannot output results. Mass Print increments are too small")
  if((grid.in@mmax-grid.in@mmin) < grid.in@mstep || (grid.in@mmax-grid.in@mmin) < grid.in@moutstep) stop("Cannot output results. Mass Range is too small")
  if(abs((grid.in@moutstep/grid.in@mstep)-round(grid.in@moutstep/grid.in@mstep)) > tiny) stop("Cannot output results. Mass Print increment is not a multiple of Mass Step")

  if(run.in@spatial_dim==1 || run.in@spatial_dim==2){
    if(grid.in@xoutstep < grid.in@xstep) stop("Cannot output results. X-Space Print increments are too small")
    if((grid.in@xmax-grid.in@xmin) < grid.in@xstep || (grid.in@xmax-grid.in@xmin) < grid.in@xoutstep) stop("Cannot output results. X-Space Range is too small")
    if(abs((grid.in@xoutstep/grid.in@xstep)-round(grid.in@xoutstep/grid.in@xstep)) > tiny) stop("Cannot output results. X-Space Print increment is not a multiple of X-Space Step")

    if(run.in@spatial_dim==2){
      if(grid.in@youtstep < grid.in@ystep) stop("Cannot output results. Y-Space Print increments are too small")
      if((grid.in@ymax-grid.in@ymin) < grid.in@ystep || (grid.in@ymax-grid.in@ymin) < grid.in@youtstep) stop("Cannot output results. Y-Space Range is too small")
      if(abs((grid.in@youtstep/grid.in@ystep)-round(grid.in@youtstep/grid.in@ystep)) >tiny ) stop("Cannot output results. Y-Space Print increment is not a multiple of X-Space Step")
    }
  }

#-------------#
# Call C code #
#-------------#

  x<-.C('SizeSpectrum', run.params, grid.params, pla.params, pel.params, ben.params, det.params, names.params, flags.params, PACKAGE='dbpmr')

#-----------------------------------------------------------------#
# Return filenames of files that have been produced by the C code #
#-----------------------------------------------------------------#
  
  if(run.in@no_benthic!=0){
    filenames<-vector(mode='character',length=(run.in@no_pelagic+run.in@no_benthic+3))
  }
  else{
    filenames<-vector(mode='character',length=(run.in@no_pelagic+2))
  }

  filenames[1]<-paste(getwd(),"/",run.in@filename,sep="")  
  filenames[2]<-plankton.in@filename
  for(i in 1:run.in@no_pelagic){
    filenames[i+2]<-pel.names[i]
  }
  if(run.in@no_benthic!=0){
    for(i in 1:run.in@no_benthic){
      filenames[2+run.in@no_pelagic+i]<-ben.names[i]
    }
    filenames[run.in@no_pelagic+run.in@no_benthic+3]<-detritus.in@filename
  }
  
  return(filenames)
}

