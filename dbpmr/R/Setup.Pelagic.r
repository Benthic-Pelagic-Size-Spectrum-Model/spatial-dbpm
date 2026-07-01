#' Set up a pelagic species
#'
#' Builds a [pelagic.params] object and creates the species output directory
#' under the run directory.
#'
#' @param run.in A [run.params] object.
#' @param mmin,mmat,mmax Minimum, maturation and maximum log-mass of the species.
#'   The defaults are the canonical FishMIP consumer range `10^-3` to `10^6` g
#'   (natural log, `-3*log(10)` to `6*log(10)`), so the consumer minimum is
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
#'   efficiency (the growth efficiencies are `K_pla`/`K_pel`/`K_ben`).
#' @param u_0 Abundance constant.
#' @param lambda Spectrum slope.
#' @param K_pla,R_pla,Ex_pla Growth (assimilation-to-growth), reproduction
#'   (retention) and excretion fractions of intake when feeding on plankton.
#'   Together with defecation they partition each unit of intake, so
#'   `K + R + Ex = 1 - defecation`. Default to the corresponding pelagic values.
#' @param K_pel,R_pel,Ex_pel Growth, reproduction and excretion fractions for
#'   feeding on other pelagics.
#' @param K_ben,R_ben,Ex_ben Growth, reproduction and excretion fractions for
#'   feeding on benthos.
#' @param pref_pla,pref_pel,pref_ben Feeding preferences for plankton, pelagic
#'   and benthic prey.
#' @param q_0 Preferred predator-prey log-mass ratio.
#' @param sig Width of the feeding kernel.
#' @param trunc Truncation of the feeding kernel.
#' @param prey,pred,comp Density-dependence coefficients for prey, predator and
#'   competition effects.
#' @param gamma_prey,gamma_pred,gamma_comp Mass-scaling exponents for the prey,
#'   predator and competition effects.
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
#' @return A [pelagic.params] object.
#' @note `rep_method = 1` requires a reproduction time-series input file (see
#'   [Setup.Rep()]). These input data must be supplied by the user and are
#'   **not** provided with the package; without the file [SizeSpectrum()] stops
#'   with an error. The other methods (`0`, `2`, `3`) need no input file.
#' @seealso [Setup.Benthic()], [Setup.Rep()], [Setup.fishing()]
#' @export
Setup.Pelagic<-function(run.in, mmin=-3*log(10), mmat=2*log(10), mmax=6*log(10), A=640, alpha=0.82, mu_0=0.2, beta=-0.25, mu_s=0.1, epsilon=0.1, u_0=0.01, lambda=-1, K_pla, R_pla, Ex_pla, K_pel=0.2, R_pel=0.2, Ex_pel=0.3, K_ben=0.1, R_ben=0.2, Ex_ben=0.4, pref_pla=1, pref_pel=1, pref_ben=1, q_0=log(100), sig=log(10), trunc=2, prey=0, pred=0, comp=0.1, gamma_prey=0.33, gamma_pred=0.33, gamma_comp=0.75, rep_method=2, initial_flag=FALSE, ts_flag=FALSE, fishing_flag=FALSE, filename){

  #Assign Default Values
  if( missing(filename) ) stop("A Species filename must be given")
  if( !is.character(filename) ) stop("Please enter a valid value for filename: (character string)")
  if( !(rep_method==0 || rep_method ==1 || rep_method==2 || rep_method ==3) ) stop("Please enter a valid value for rep_method: (0,1)")
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
  species@gamma_prey<-gamma_prey
  species@gamma_pred<-gamma_pred
  species@gamma_comp<-gamma_comp

  species@rep_method<-as.integer(rep_method)
  species@initial_flag<-as.logical(initial_flag)
  species@ts_flag<-as.logical(ts_flag)
  species@fishing_flag<-as.logical(fishing_flag)

  dir.create(paste(run.in@filename,"/",filename,sep=""),showWarnings=FALSE)
  
  return(species)

}