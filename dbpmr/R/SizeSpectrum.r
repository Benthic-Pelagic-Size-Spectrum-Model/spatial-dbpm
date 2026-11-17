#' Run the dynamic size-spectrum simulation
#'
#' Validates the supplied parameter objects, assembles the transfer vectors and
#' calls the underlying C engine to integrate the dynamic benthic-pelagic size
#' spectrum. Output text files are written under the run directory.
#'
#' @param run.in A [run.params] object.
#' @param grid.in A [grid.params] object.
#' @param plankton.in A [plankton.params] object.
#' @param pelagic.in A [pelagic.params] object, or a list of them (one per
#'   pelagic species).
#' @param benthic.in A [benthic.params] object, or a list of them. Optional when
#'   `run.in@no_benthic == 0`.
#' @param detritus.in A [detritus.params] object. Optional when
#'   `run.in@no_benthic == 0`.
#' @param forcing Optional environmental forcing held in memory, from
#'   [Setup.Forcing()]. When given, it replaces the `forcing_ts.txt` input file
#'   (it is an error to supply both).
#' @param gravity Optional allocation of fishing effort among cells during the
#'   run, from [Setup.Gravity()]. The rates given to [Setup.fishing()] are then
#'   the fishing mortality at evenly spread effort.
#'
#' @param output Where the model output goes. `"files"` (the default) writes the
#'   text result files under the run directory. `"arrays"` returns the output to
#'   R as arrays and writes no files at all; the record of the run's parameters
#'   that `parameters.txt` would hold comes back in the output object. `"both"`
#'   does both.
#' @param variables Used with `output = "arrays"` or `"both"`: which variables
#'   to return. Abundance is always returned; `"growth"`, `"mortality"`,
#'   `"predation"`, `"fishing"` and `"catch"` are returned by default and can be
#'   left out to save memory on large grids.
#'
#' @return With `output = "files"`, a character vector of the paths and species
#'   filenames produced by the run, suitable for passing to [Read.In()].
#'
#'   With `output = "arrays"` or `"both"`, an object of class `dbpmr_output`, a
#'   list with:
#'
#'   * `t` (years), `x`, `y` and `m` (natural-log body mass): the output
#'     coordinates;
#'   * `abundance`: a named list with one array per group (the plankton and each
#'     pelagic and benthic species, named by their file names), each with
#'     dimensions `[x, y, m, t]`;
#'   * `growth`, `mortality`, `fishing`: named lists with one `[x, y, m, t]`
#'     array per pelagic and benthic species, holding the growth rate, the total
#'     mortality rate and the fishing mortality rate;
#'   * `predation`: predation mortality rate. For a pelagic species, a list with
#'     one array per pelagic predator species; for a benthic species, one array
#'     (predation by all pelagic species);
#'   * `catch`: a named list with one `[x, y, m, t]` array per pelagic and
#'     benthic species, holding the catch rate density: fishing mortality times
#'     abundance times body mass, per unit natural-log body mass. Summed over
#'     size classes and multiplied by the size step (`diff(m)[1]`) it is the
#'     total catch rate. Like every other variable it is the rate at the output
#'     time, not a total over the interval between outputs. (In a run without
#'     movement the total equals the `fishing` column of the species'
#'     `summary.txt`; with movement that column is computed before the
#'     movement step and differs slightly.);
#'   * `detritus`: detritus biomass as `[x, y, t]` (when the run has benthos);
#'   * `parameters`: the parameters of the run, as a list with `run`, `grid`,
#'     `groups` (one entry per group, named by its file name, holding its
#'     parameter values and, in `in_memory`, the names of the inputs it was given
#'     in memory), `forcing` (the forcing channels supplied in memory and
#'     `temp_on_mort`), `gravity`, and `environment` (the settings that were
#'     taken from environment variables). [Write.NetCDF()] saves it with the
#'     file.
#'
#'   The values are those the text files would hold (`results.txt`,
#'   `growth.txt`, `mortality.txt`, `predation*.txt`, `fishing.txt`): the same
#'   output times, cells and size classes (set by `toutstep`, `xoutstep`,
#'   `youtstep` and `moutstep` in [Setup.Grid()]), with zeros outside a group's
#'   size range. With `"both"` the character vector described above is attached
#'   as the attribute `"files"`. See [Write.NetCDF()] to save the object.
#' @seealso [Read.In()], [Write.NetCDF()], [Setup.Run()]
#' @export
`SizeSpectrum` <-
function(run.in, grid.in, plankton.in, pelagic.in, benthic.in, detritus.in, forcing = NULL, output = c("files", "arrays", "both"),
         variables = c("abundance", "growth", "mortality", "predation", "fishing", "catch"), gravity = NULL){
  output <- match.arg(output)
  variables <- match.arg(variables, several.ok = TRUE)

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

  # Check that every value handed to the C code is finite. .C() used to reject NA/NaN/Inf
  # itself; .Call() does not, so a NaN parameter would otherwise run to completion and
  # fill the output files with NaN.
  check_finite <- function(v, what, nms = NULL){
    bad <- which(!is.finite(v))
    if(length(bad)!=0){
      where <- if(is.null(nms)) paste("position", paste(bad, collapse=", ")) else paste(nms[bad], collapse=", ")
      stop("Non-finite value (NA, NaN or Inf) in the ", what, " parameters: ", where)
    }
  }
  check_finite(grid.params, "grid")
  check_finite(pla.params, "plankton", c("mmin", "mmax", "mu_0", "beta", "u_0", "lambda"))
  check_finite(pel.params, "pelagic")
  if(run.in@no_benthic!=0){
    check_finite(ben.params, "benthic")
    check_finite(det.params, "detritus")
  }
  if(anyNA(run.params)) stop("Missing value (NA) in the run parameters")
  if(anyNA(flags.params)) stop("Missing value (NA) in the species flags")

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
  #(a plankton time series held in memory as intercept/slope arrays needs no file)
  if(is.null(plankton.in@arrays$ts) && is.null(plankton.in@arrays$initial) && is.null(plankton.in@arrays$series) && (plankton.in@initial_flag==1 || plankton.in@ts_flag==1) && !file.exists(paste(run.in@filename,"/Input/",plankton.in@filename,"_ts.txt",sep=""))){
    stop("A required Plankton initial/ts file does not exist. Write one with Setup.ts(), or, for values held in memory, ",
         "assign the result: pl <- Setup.ts(pl, run, grid, intercept = , slope = ) or pl <- Setup.ts(pl, run, grid, initial = )")
  }
  
  for(i in 1:length(pelagic.in)){
    if(is.null(pelagic.in[[i]]@arrays$initial) && is.null(pelagic.in[[i]]@arrays$series) && (pelagic.in[[i]]@initial_flag==1 || pelagic.in[[i]]@ts_flag==1) && !file.exists(paste(run.in@filename,"/Input/",pelagic.in[[i]]@filename,"_ts.txt",sep=""))) stop("A required Pelagic initial/ts file does not exist. Write one with Setup.ts(), or, for starting values held in memory, assign the result: pe <- Setup.ts(pe, run, grid, initial = )")
    if( pelagic.in[[i]]@fishing_flag==1 && is.null(pelagic.in[[i]]@arrays$fishing) && !file.exists(paste(run.in@filename,"/Input/",pelagic.in[[i]]@filename,"_fishing_ts.txt",sep=""))) stop("A required Pelagic fishing file does not exist. Write one with Setup.fishing(), or, for a rate held in memory, assign the result: pe <- Setup.fishing(pe, run, grid, rate = )")
    if( pelagic.in[[i]]@rep_method==1 && is.null(pelagic.in[[i]]@arrays$reproduction) && !file.exists(paste(run.in@filename,"/Input/",pelagic.in[[i]]@filename,"_rep_ts.txt",sep=""))) stop("A required Pelagic reproduction file does not exist. Write one with Setup.Rep(), or, for values held in memory, assign the result: pe <- Setup.Rep(pe, run, grid, eggs = )")
    
  }
  if(run.in@no_benthic!=0){
    for(i in 1:length(benthic.in)){
      if(is.null(benthic.in[[i]]@arrays$initial) && is.null(benthic.in[[i]]@arrays$series) && (benthic.in[[i]]@initial_flag==1 || benthic.in[[i]]@ts_flag==1) && !file.exists(paste(run.in@filename,"/Input/",benthic.in[[i]]@filename,"_ts.txt",sep=""))) stop("A required Benthic initial/ts file does not exist. Write one with Setup.ts(), or, for starting values held in memory, assign the result: be <- Setup.ts(be, run, grid, initial = )")
      if( benthic.in[[i]]@fishing_flag==1 && is.null(benthic.in[[i]]@arrays$fishing) && !file.exists(paste(run.in@filename,"/Input/",benthic.in[[i]]@filename,"_fishing_ts.txt",sep=""))) stop("A required Benthic fishing file does not exist. Write one with Setup.fishing(), or, for a rate held in memory, assign the result: be <- Setup.fishing(be, run, grid, rate = )")
      if( benthic.in[[i]]@rep_method==1 && is.null(benthic.in[[i]]@arrays$reproduction) && !file.exists(paste(run.in@filename,"/Input/",benthic.in[[i]]@filename,"_rep_ts.txt",sep=""))) stop("A required Benthic reproduction file does not exist. Write one with Setup.Rep(), or, for values held in memory, assign the result: be <- Setup.Rep(be, run, grid, eggs = )")
    }
    if(is.null(detritus.in@arrays$initial) && is.null(detritus.in@arrays$series) && (detritus.in@initial_flag==1 || detritus.in@ts_flag==1) && !file.exists(paste(run.in@filename,"/Input/",detritus.in@filename,"_ts.txt",sep=""))) stop("A required Detritus initial/ts file does not exist. Write one with Setup.ts(), or, for a starting value held in memory, assign the result: de <- Setup.ts(de, run, grid, initial = )")
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

  # Inputs passed to C as arrays instead of text files (NULL = none)
  inputs <- NULL
  if(!is.null(plankton.in@arrays$ts)){
    ta <- plankton.in@arrays$ts
    inputs <- list(plankton_intercept=ta$intercept, plankton_slope=ta$slope,
                   plankton_layer=as.integer(ta$layer-1L), plankton_nlayer=as.integer(dim(ta$intercept)[3]))
  }
  # Fishing mortality held in memory: one entry per species (NULL where it comes from a file or the species is not fished)
  fish_in <- function(sp.list, what){
    lapply(sp.list, function(sp){
      fa <- sp@arrays$fishing
      if(is.null(fa)) return(NULL)
      if(!isTRUE(sp@fishing_flag)) stop("fishing rate supplied for ", sp@filename, " but its fishing_flag is FALSE")
      if(file.exists(file.path(run.in@filename, "Input", paste0(sp@filename, "_fishing_ts.txt"))))
        stop("Both a fishing rate held in memory and the file ", sp@filename, "_fishing_ts.txt were supplied for ", what, " species ", sp@filename, "; use one or the other")
      fgrid <- c(length(grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep)), length(grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep)))
      if(!identical(as.integer(dim(fa$rate)[1:2]), as.integer(fgrid)) || length(fa$layer)!=length(seq(0, grid.in@tmax, grid.in@tstep)))
        stop("the fishing rate of ", sp@filename, " was built for a different grid; rebuild it with Setup.fishing() using this grid")
      list(rate=fa$rate, layer=as.integer(fa$layer-1L), nlayer=as.integer(dim(fa$rate)[3]), range=as.double(fa$range))
    })
  }
  fpel <- fish_in(pelagic.in, "pelagic")
  fben <- if(run.in@no_benthic!=0) fish_in(benthic.in, "benthic") else list()
  if(any(!vapply(fpel, is.null, TRUE)) || any(!vapply(fben, is.null, TRUE))){
    if(is.null(inputs)) inputs <- list()
    inputs$fishing_pelagic <- fpel
    if(run.in@no_benthic!=0) inputs$fishing_benthic <- fben
  }
  # Other inputs held in memory, one entry per group (NULL where a group has none): starting values, a series
  # forced through time, reproduction, refuge and per-cell parameters. Each replaces its text file; supplying
  # both is an error.
  ncell <- c(length(grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep)), length(grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep)))
  ntime <- length(seq(0, grid.in@tmax, grid.in@tstep))
  in_file <- function(nm, suffix) file.exists(file.path(run.in@filename, "Input", paste0(nm, suffix)))
  not_both <- function(nm, suffix, what) if(in_file(nm, suffix))
    stop("Both ", what, " held in memory and the file ", nm, suffix, " were supplied for ", nm, "; use one or the other")
  same_grid <- function(a, nm, what, setup) if(!identical(as.integer(dim(a)[1:2]), as.integer(ncell)))
    stop("the ", what, " of ", nm, " was built for a different grid; rebuild it with ", setup, " using this grid")
  initial_in <- function(sp){
    a <- sp@arrays$initial
    if(is.null(a)) return(NULL)
    if(!isTRUE(sp@initial_flag) || isTRUE(sp@ts_flag)) stop("starting values were supplied for ", sp@filename, " but its initial_flag is FALSE (or its ts_flag is TRUE)")
    not_both(sp@filename, "_ts.txt", "starting values")
    same_grid(a, sp@filename, "starting values", "Setup.ts()")
    a
  }
  sp_in <- function(sp.list){
    lapply(sp.list, function(sp){
      a <- sp@arrays; out <- list()
      out$initial <- initial_in(sp)
      if(!is.null(a$series)){
        if(!isTRUE(sp@ts_flag)) stop("a forced series was supplied for ", sp@filename, " but its ts_flag is FALSE")
        if(!is.null(a$ts)) stop("Both a series and intercept/slope arrays were supplied for ", sp@filename, "; use one or the other")
        not_both(sp@filename, "_ts.txt", "a forced series")
        same_grid(a$series$values, sp@filename, "forced series", "Setup.ts()")
        if(length(a$series$layer)!=ntime) stop("the forced series of ", sp@filename, " was built for a different grid; rebuild it with Setup.ts() using this grid")
        nd <- dim(a$series$values)
        out$series <- a$series$values; out$series_layer <- as.integer(a$series$layer-1L); out$series_nlayer <- as.integer(nd[length(nd)])
      }
      if(!is.null(a$reproduction)){
        if(sp@rep_method!=1) stop("reproduction values were supplied for ", sp@filename, " but its rep_method is not 1")
        not_both(sp@filename, "_rep_ts.txt", "reproduction values")
        same_grid(a$reproduction$values, sp@filename, "reproduction values", "Setup.Rep()")
        if(length(a$reproduction$layer)!=ntime) stop("the reproduction values of ", sp@filename, " were built for a different grid; rebuild them with Setup.Rep() using this grid")
        out$rep <- a$reproduction$values; out$rep_layer <- as.integer(a$reproduction$layer-1L); out$rep_nlayer <- as.integer(dim(a$reproduction$values)[3])
      }
      if(!is.null(a$refuge)){
        not_both(sp@filename, "_refuge.txt", "a refuge")
        same_grid(a$refuge$values, sp@filename, "refuge", "Setup.Refuge()")
        if(dim(a$refuge$values)[3]!=length(grid_seq(grid.in@mmin, grid.in@mmax, grid.in@mstep)))
          stop("the refuge of ", sp@filename, " has ", dim(a$refuge$values)[3], " mass classes but the run has ", length(grid_seq(grid.in@mmin, grid.in@mmax, grid.in@mstep)),
               "; build it with a grid whose mmin and mmax span all the species")
        out$refuge <- a$refuge$values; out$refuge_min_a <- as.double(a$refuge$min_A)
      }
      if(methods::is(sp, "pelagic.params") && !is.na(sp@imm_frac)) out$imm_frac <- as.double(sp@imm_frac)   # recruitment floor given as a parameter
      if(!is.null(a$pref_ben)){
        not_both(sp@filename, "_cellpar.txt", "per-cell parameters")
        same_grid(a$pref_ben, sp@filename, "per-cell parameters", "Setup.CellParams()")
        out$pref_ben <- a$pref_ben
      }
      if(length(out)) out else NULL
    })
  }
  spel <- sp_in(pelagic.in)
  sben <- if(run.in@no_benthic!=0) sp_in(benthic.in) else list()
  spla <- sp_in(list(plankton.in))
  sdet <- if(run.in@no_benthic!=0) sp_in(list(detritus.in)) else list()
  given <- function(l) any(!vapply(l, is.null, TRUE))
  if(given(spel) || given(sben) || given(spla) || given(sdet)){
    if(is.null(inputs)) inputs <- list()
    if(given(spel)) inputs$species_pelagic <- spel
    if(given(sben)) inputs$species_benthic <- sben
    if(given(spla)) inputs$species_plankton <- spla
    if(given(sdet)) inputs$species_detritus <- sdet
  }
  if(!is.null(forcing)){
    if(!inherits(forcing, "dbpmr_forcing")) stop("forcing must be created with Setup.Forcing()")
    if(file.exists(file.path(run.in@filename, "Input", "forcing_ts.txt")))
      stop("Both forcing arrays and the file ", file.path(run.in@filename, "Input", "forcing_ts.txt"),
           " were supplied; use one or the other")
    fgrid <- c(length(grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep)), length(grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep)))
    if(!identical(as.integer(dim(forcing$channels[[1]])[1:2]), as.integer(fgrid)) ||
       length(forcing$layer)!=length(seq(0, grid.in@tmax, grid.in@tstep)))
      stop("forcing was built for a different grid; rebuild it with Setup.Forcing() using this grid")
    fnames <- c("pel_tempeff", "ben_tempeff", "sinking_rate", "depth", "residence_time",
                "export_attn", "export_gamma", "export_magg", "gravity")
    fin <- forcing$channels
    names(fin) <- paste0("forcing_", names(fin))
    fnl <- vapply(fnames, function(nm) if(is.null(forcing$channels[[nm]])) 0L else dim(forcing$channels[[nm]])[3], integer(1))
    if(is.null(inputs)) inputs <- list()
    inputs <- c(inputs, fin, list(forcing_layer=as.integer(forcing$layer-1L), forcing_nlayer=as.integer(fnl)))
    if(!is.null(forcing$temp_on_mort)) inputs$forcing_temp_on_mort <- as.integer(forcing$temp_on_mort)
  }

  if(!is.null(gravity)){
    if(!inherits(gravity, "dbpmr_gravity")) stop("gravity must be created with Setup.Gravity()")
    fgrid <- c(length(grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep)), length(grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep)))
    if(!identical(as.integer(dim(gravity$access)), as.integer(fgrid)))
      stop("gravity was built for a different grid; rebuild it with Setup.Gravity() using this grid")
    fished <- c(vapply(pelagic.in, function(sp) isTRUE(sp@fishing_flag), TRUE),
                if(run.in@no_benthic!=0) vapply(benthic.in, function(sp) isTRUE(sp@fishing_flag), TRUE))
    if(!any(fished)) stop("gravity was supplied but no species is fished (fishing_flag = TRUE and a rate from Setup.fishing())")
    if(is.null(inputs)) inputs <- list()
    inputs$gravity <- list(access=as.double(gravity$access), thick_pel=as.double(gravity$pelagic_thickness),
                           thick_ben=as.double(gravity$benthic_thickness), nupdate=as.integer(gravity$nupdate))
  }

  # The text output goes in the run folder and one folder per group
  group.names <- c(pla.names, pel.names, ben.names, det.names)
  if(output != "arrays"){
    for(d in file.path(run.in@filename, group.names)) dir.create(d, showWarnings=FALSE, recursive=TRUE)
  }

  x <- .Call(
    C_SizeSpectrum, 
    as.integer(run.params), 
    as.double(grid.params), 
    as.double(pla.params), 
    as.double(pel.params), 
    as.double(ben.params), 
    as.double(det.params), 
    as.character(names.params), 
    as.integer(flags.params),
    inputs,
    c(c(files = 0L, both = 1L, arrays = 2L)[[output]],
      as.integer(c("growth", "mortality", "predation", "fishing", "catch") %in% variables))
  )

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
  
  if(output == "files") return(filenames)

  # Output returned from C as arrays: give them their dimensions and names
  dims <- c(length(x$x), length(x$y), length(x$m), length(x$t))
  shape <- function(a){ dim(a) <- dims; a }
  has_ben <- run.in@no_benthic!=0
  spnames <- c(pel.names, if(has_ben) ben.names)
  # one array per pelagic then per benthic species, named by species
  by_species <- function(pel, ben){ l <- c(lapply(pel, shape), lapply(ben, shape)); names(l) <- spnames; l }
  out <- list(t = x$t, x = x$x, y = x$y, m = x$m,
              abundance = c(structure(list(shape(x$plankton)), names = plankton.in@filename), by_species(x$pelagic, x$benthic)))
  if("growth" %in% variables)    out$growth    <- by_species(x$pelagic_growth, x$benthic_growth)
  if("mortality" %in% variables) out$mortality <- by_species(x$pelagic_mortality, x$benthic_mortality)
  if("predation" %in% variables){
    pp <- lapply(x$pelagic_predation, function(l){ l <- lapply(l, shape); names(l) <- pel.names; l })   # by predator species
    out$predation <- c(pp, lapply(x$benthic_predation, shape)); names(out$predation) <- spnames
  }
  if("fishing" %in% variables)   out$fishing   <- by_species(x$pelagic_fishing, x$benthic_fishing)
  if("catch" %in% variables)     out$catch     <- by_species(x$pelagic_catch, x$benthic_catch)
  if(has_ben){
    det <- x$detritus; dim(det) <- dims[c(1, 2, 4)]
    out$detritus <- det
  }
  # The parameters of the run, kept with the output (what parameters.txt records on the text-file route)
  record <- function(obj){
    sl <- setdiff(methods::slotNames(obj), "arrays")
    rec <- stats::setNames(lapply(sl, function(n) methods::slot(obj, n)), sl)
    if("arrays" %in% methods::slotNames(obj)) rec$in_memory <- names(obj@arrays)
    rec
  }
  groups <- c(list(plankton.in), pelagic.in, if(has_ben) c(benthic.in, list(detritus.in)))
  out$parameters <- list(run = record(run.in), grid = record(grid.in),
                         groups = stats::setNames(lapply(groups, record), vapply(groups, function(g) g@filename, "")),
                         forcing = if(!is.null(forcing)) list(channels = names(forcing$channels), temp_on_mort = forcing$temp_on_mort),
                         gravity = if(!is.null(gravity)) list(update = gravity$update))
  # settings still taken from the environment (where the matching argument was left at its default)
  env <- Sys.getenv(c("PEL_IMM_FRAC", "TEMP_ON_MORT", "REFUGE_MIN_A", "FORCING_PERCELL"), unset = NA)
  out$parameters$environment <- env[!is.na(env)]
  class(out) <- "dbpmr_output"
  if(output == "both") attr(out, "files") <- filenames
  # Arrays-only output writes nothing: remove the folders made for this run if they are still empty
  if(output == "arrays"){
    is_empty <- function(d) dir.exists(d) && length(list.files(d, all.files=TRUE, no..=TRUE))==0
    for(d in c(file.path(run.in@filename, c(group.names, "Input")), run.in@filename)) if(is_empty(d)) unlink(d, recursive=TRUE)
  }
  return(out)
}

