#' Set up a benthic species
#'
#' Builds a [benthic.params] object and creates the species output directory
#' under the run directory.
#'
#' @param run.in A [run.params] object.
#' @param mmin,mmat,mmax Minimum, maturation and maximum log-mass of the species.
#'   The defaults are the canonical FishMIP detritivore range `10^-3` to `10^4` g
#'   (natural log, `-3*log(10)` to `4*log(10)`), so the consumer minimum is
#'   `10^-3` g. Require `mmin < mmat < mmax`: a maturation size at or above `mmax`
#'   leaves no mature size classes, which zeroes endogenous reproduction
#'   (`rep_method` 2/3) and silently collapses the species — this is warned.
#' @param A Search/encounter rate constant.
#' @param alpha Search-rate mass-scaling exponent.
#' @param mu_0 Background mortality constant.
#' @param beta Mortality mass-scaling exponent.
#' @param mu_s Senescence mortality constant.
#' @param epsilon Senescence size offset: added to the maximum log10 size in the
#'   senescence denominator
#'   `mu_s * (log10 w - log10 w_min) / ((log10 w_max + epsilon) - log10 w)`, so
#'   senescence mortality stays finite at `w_max`. This is **not** a growth
#'   efficiency (the growth efficiency is `K_det`).
#' @param u_0 Abundance constant.
#' @param lambda Spectrum slope.
#' @param K_det,R_det,Ex_det Growth (assimilation-to-growth), reproduction
#'   (retention) and excretion fractions of intake when feeding on detritus.
#'   Together with defecation they partition each unit of intake, so
#'   `K + R + Ex = 1 - defecation`.
#' @param pref_det Feeding preference for detritus.
#' @param rep_method Integer reproduction method. `0` = fixed recruitment held
#'   at the initial density of the smallest size class; `1` = a prescribed
#'   reproduction time series read from an input file (see [Setup.Rep()]);
#'   `2` (the default) = reproduction allocated from assimilated energy of
#'   mature individuals; `3` = reproduction proportional to mature biomass.
#'   Methods `0`, `2` and `3` require no input file.
#' @param initial_flag,ts_flag,fishing_flag Logical flags indicating whether
#'   initial-condition, time-series and fishing input files are supplied.
#' @param filename Character string naming the species (and its output
#'   directory). Required.
#'
#' @return A [benthic.params] object.
#' @note `rep_method = 1` requires a reproduction time-series input file (see
#'   [Setup.Rep()]). These input data must be supplied by the user and are
#'   **not** provided with the package; without the file [SizeSpectrum()] stops
#'   with an error. The other methods (`0`, `2`, `3`) need no input file.
#' @seealso [Setup.Pelagic()], [Setup.Detritus()]
#' @export
Setup.Benthic<-function(run.in, mmin=-3*log(10), mmat=0, mmax=4*log(10), A=64, alpha=0.75, mu_0=0.2, beta=-0.25, mu_s=0.1, epsilon=0.1, u_0=0.01, lambda=-0.75, K_det=0.2, R_det=0.2, Ex_det=0.2, pref_det=1, rep_method=2, initial_flag=FALSE, ts_flag=FALSE, fishing_flag=FALSE, filename){

  #Assign Default Values
  if( missing(filename) ) stop("A Species filename must be given")
  if( !is.character(filename) ) stop("Please enter a valid value for filename: (character string)")
  if( !(rep_method==0 || rep_method ==1 || rep_method==2 || rep_method ==3) ) stop("Please enter a valid value for rep_flag: (0,1)")
  if( !(initial_flag==0 || initial_flag==1) ) stop("Please enter a logical value for initial_flag: (T,F)")
  if( !(ts_flag==0 || ts_flag==1) ) stop("Please enter a logical value for ts_flag: (T,F)")
  if( !(fishing_flag==0 || fishing_flag==1) ) stop("Please enter a logical value for fishing_flag: (T,F)")

  #Maturation-size sanity: mmin < mmat < mmax. If mmat >= mmax there are no mature
  #size classes, so endogenous reproduction (rep_method 2/3) produces no eggs and
  #the species collapses silently with no error - warn rather than fail.
  if( mmin >= mmax ) stop("mmin must be less than mmax")
  if( mmat >= mmax ) warning("mmat (", signif(mmat, 4), ") >= mmax (", signif(mmax, 4),
    "): no mature size classes. With rep_method 2 or 3 this gives zero reproduction and ",
    "the species will collapse silently - set mmat below mmax.", call.=FALSE)
  if( mmat <= mmin ) warning("mmat (", signif(mmat, 4), ") <= mmin (", signif(mmin, 4),
    "): every size class is treated as mature.", call.=FALSE)

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

  dir.create(paste(run.in@filename,"/",filename,sep=""),showWarnings=FALSE)
  
  return(species)
}