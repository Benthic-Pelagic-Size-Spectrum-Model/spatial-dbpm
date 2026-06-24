`Setup.ts`<-
function(species, run.in, grid.in, func, mat, dataname){
#func is a function of the variables m, x, y, and z
#dataname is a filename containing either a vector of intercepts, two vectors of intercepts and slopes or a full matrix of ts values

#Input checking stuff
if(species@initial_flag==F && species@ts_flag==F) stop ("Starting values or time series are not to be specified since initial_flag==F or ts_flag=F")

if(missing(func) && missing(dataname) && missing(mat)) stop("A function, matrix or data for the initial step or time series must be given")
if(!missing(func) && !missing(dataname) && !missing(mat)) stop("Only one of a function or data for the initial step step or time series must be given")
if(!missing(func) && !missing(dataname) && missing(mat)) stop("Only one of a function or data for the initial step step or time series must be given")
if(!missing(func) && missing(dataname) && !missing(mat)) stop("Only one of a function or data for the initial step step or time series must be given")
if(missing(func) && !missing(dataname) && !missing(mat)) stop("Only one of a function or data for the initial step step or time series must be given")

#Grid stuff
mass<-seq(grid.in@mmin,grid.in@mmax,grid.in@mstep)
xrange<-seq(grid.in@xmin,grid.in@xmax,grid.in@xstep)
yrange<-seq(grid.in@ymin,grid.in@ymax,grid.in@ystep)
trange<-seq(0,grid.in@tmax,grid.in@tstep)

#Species stuff
filename<-paste(species@filename,"_ts.txt",sep="")

if(species@speciestype!="detritus"){
  spmin<-species@mmin
  spmax<-species@mmax
}
if(species@speciestype=="detritus"){
  spmin<-0
  spmax<-0
  mass<-0
}

m<-length(mass)
x<-length(xrange)
y<-length(yrange)
t<-length(trange)

#Calculation of initial values /time series using function
if(!missing(func)){
  if(species@ts_flag==T){
    temp<-matrix(0,nrow=(t*x*y),ncol=m)
    for(j in 1:t){
      for(k in 1:x){
        for(l in 1:y){
          for(i in which(mass==spmin):which(mass==spmax)){
            temp[((j-1)*x*y)+((k-1)*y)+l,i]=round(func(mass[i],trange[j],xrange[k],yrange[l]),16)
          }
        }
      }
    }
  }
  if(species@ts_flag==F){
  temp<-matrix(0,nrow=(x*y),ncol=m)
    for(k in 1:x){
      for(l in 1:y){
        for(i in which(mass==spmin):which(mass==spmax)){
          temp[((k-1)*y)+l,i]=round(func(mass[i],trange[1],xrange[k],yrange[l]),16)
        }
      }
    }
  }
}

#Calcualtion of initial values / time series using data file
if(!missing(dataname)){
#read in csv files or intercepts and/or slopes
#Each csv file must contain a header row with 'intercept', 'interslope' or 'spectrum'
  if(!is.character(dataname)) stop("Please enter a valid filename string")
  if(!file.exists(dataname)) stop("Filename specified does not exist")

  filetype<-as.character(read.csv(dataname,nrows=1,header=F,strip.white=T,stringsAsFactors=F))
  dat<-read.csv(dataname,skip=1,header=F,strip.white=T,stringsAsFactors=F)
  dat<-as.matrix(dat)
  
  if(species@ts_flag==T){
    #Check whether the number of lines in data file matches the number of time steps previously specified
    if(length(dat[,1])!=t*x*y){
      stop("Incorrect number of rows in data file")
    }
    
    #Setup input storage array
    temp<-matrix(0,nrow=(t*x*y),ncol=m)
    
    if(species@speciestype!="detritus"){
    
      if(filetype=="intercept"){
        for(i in which(mass==spmin):which(mass==spmax)){
          temp[,i]=round(dat[,1]*exp(species@lambda*mass[i]),16)
        }
      }
    
      if(filetype=="interslope"){
        if(length(dat[1,])<2) stop("Insufficient columns in data file")
        for(i in which(mass==spmin):which(mass==spmax)){
          temp[,i]=round(dat[,1]*exp(dat[,2]*mass[i]),16)
        }
      }
    
      if(filetype=="spectrum"){
        if(length(dat[1,])!=length(seq(spmin,spmax,grid.in@mstep))) stop("Incorrect number of columns in data file")
        for(i in which(mass==spmin):which(mass==spmax)){
          temp[,i]=round(dat[,(i-which(mass==spmin)+1)],16)
        }
      }
      if(!(filetype=="intercept" || filetype=="interslope" || filetype=="spectrum")) stop("Unknown filetype entered")
    
    }
    
    if(species@speciestype=="detritus"){
      temp[,1]=round(dat[,1],16)
    }
    
  }
  
  if(species@ts_flag==F){
    if(length(dat[,1])<1) stop("Insufficient number of rows in data file")
    temp<-matrix(0,nrow=(x*y),ncol=m)
    
    if(species@speciestype!="detritus"){    
    
      if(filetype=="intercept"){
        for(i in which(mass==spmin):which(mass==spmax)){
          temp[1,i]=round(dat[1,1]*exp(species@lambda*mass[i]),16)
        }
      }
    
      if(filetype=="interslope"){
        if(length(dat[1,])<2) stop("Insufficient columns in data file")
        for(i in which(mass==spmin):which(mass==spmax)){
          temp[1,i]=round(dat[1,1]*exp(dat[1,2]*mass[i]),16)
        }
      }
    
      if(filetype=="spectrum"){
        if(length(dat[1,])!=length(seq(spmin,spmax,grid.in@mstep))) stop("Incorrect number of columns in data file")
        for(i in which(mass==spmin):which(mass==spmax)){
          temp[1,i]=round(dat[1,(i-which(mass==spmin)+1)],16)
        }
      }
      if(!(filetype=="intercept" || filetype=="interslope" || filetype=="spectrum")) stop("Unknown filetype entered")
      
    }

    if(species@speciestype=="detritus"){
      temp[1,1]=round(dat[1,1],16)
    }
  }

}


#Calcualtion of initial values / time series using R matrix
if(!missing(mat)){
#read in matrix containing full size spectra
  if(!is.matrix(mat)) stop("mat must be of type matrix")
  
  if(species@ts_flag==T){
    #Check whether the number of lines in the matrix matches the number of time steps previously specified
    if(length(mat[,1])!=t*x*y){
      stop("Incorrect number of rows in data file")
    }

    #Setup input storage array
    temp<-matrix(0,nrow=(t*x*y),ncol=m)
    
    if(species@speciestype!="detritus"){

      if(length(mat[1,])!=length(seq(spmin,spmax,grid.in@mstep))) stop("Incorrect number of columns in matrix")
        for(i in which(mass==spmin):which(mass==spmax)){
          temp[,i]=round(mat[,(i-which(mass==spmin)+1)],16)
        }
    }
    
    if(species@speciestype=="detritus"){
      temp[,1]=round(mat[,1],16)
    }
    
  }
  
  if(species@ts_flag==F){
    if(length(mat[,1])<1) stop("Insufficient number of rows in data file")
    
    #Setup input storage array
    temp<-matrix(0,nrow=(x*y),ncol=m)
    
    if(species@speciestype!="detritus"){    
      if(length(mat[1,])!=length(seq(spmin,spmax,grid.in@mstep))) stop("Incorrect number of columns in data file")
      for(i in which(mass==spmin):which(mass==spmax)){
        temp[1,i]=round(mat[1,(i-which(mass==spmin)+1)],16)
      }
    }

    if(species@speciestype=="detritus"){
      temp[1,1]=round(mat[1,1],16)
    }
  }

}

#Create Input directory
dir.create(paste(run.in@filename,"/Input",sep=""),showWarnings=F)

#Write full table all at once
write.table(temp,paste(run.in@filename,"/input/",filename,sep=""),append=F,row.names=F,col.names=F,sep=",")

}
