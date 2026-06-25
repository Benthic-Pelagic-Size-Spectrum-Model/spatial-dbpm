#' Set up the run-level parameters
#'
#' Validates the high-level settings for a DBPM run, creates the output
#' directory named `filename`, and returns a [run.params] object.
#'
#' @param filename Character string naming the run; also used as the output
#'   directory that is created.
#' @param no_pelagic Number of pelagic species (positive integer).
#' @param no_benthic Number of benthic species (non-negative integer).
#' @param spatial_dim Spatial dimension of the run: `0`, `1` or `2`.
#' @param coupled_flag Logical; whether the benthic and pelagic systems are
#'   coupled.
#' @param diff_method Integer flagging the differencing/integration method
#'   (`0`, `1` or `2`).
#'
#' @return A [run.params] object.
#' @seealso [Setup.Grid()], [SizeSpectrum()]
#' @export
Setup.Run<-function(filename, no_pelagic, no_benthic, spatial_dim, coupled_flag, diff_method){

#---------------------#
# Error checking step #
#---------------------#

  if( !is.character(filename) ) stop("Please enter a valid value for filename: (character string)")
  if( no_pelagic<=0 || (no_pelagic%%1)!=0 ) stop("Please enter a valid value for no_pelagic: (positive integer)")
  if( no_benthic<0 || (no_benthic%%1)!=0 ) stop("Please enter a valid value for no_benthic: (positive integer)")
  if( no_pelagic==0 && no_benthic==0) stop("Run must contain at least one species")
  if( !(spatial_dim==0 || spatial_dim==1 || spatial_dim==2) ) stop("Please enter a valid value for spatial_dim: (0,1,2)")
  if( !(coupled_flag==0 || coupled_flag==1) ) stop("Please enter a logical value for coupled_flag: (T,F)")
  if( !(diff_method==0 || diff_method==1 || diff_method==2) ) stop("Please enter an integer value for the differencing method: (0,1,2)")


#-------------------------#
# Setup run.params object #
#-------------------------#

run.params<-new("run.params")

run.params@filename<-filename
run.params@no_pelagic<-as.integer(no_pelagic)
run.params@no_benthic<-as.integer(no_benthic)
run.params@spatial_dim<-as.integer(spatial_dim)
run.params@coupled_flag<-as.logical(coupled_flag)
run.params@diff_method<-as.integer(diff_method)

dir.create(filename)

return(run.params)

}
