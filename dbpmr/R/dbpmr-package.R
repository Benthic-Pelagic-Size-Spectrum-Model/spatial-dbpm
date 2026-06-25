#' dbpmr: Dynamic Benthic and Pelagic Size Spectra for Marine Ecosystems
#'
#' An R implementation of the Dynamic Benthic Pelagic Model (DBPM) with spatial
#' size-spectrum capabilities. The package provides helpers to configure a model
#' run ([Setup.Run()], [Setup.Grid()], [Setup.Plankton()], [Setup.Pelagic()],
#' [Setup.Benthic()], [Setup.Detritus()]), an interface to the underlying C
#' simulation engine ([SizeSpectrum()]), and tools to read in, summarise and
#' plot the resulting size spectra ([Read.In()], [Extract.Time()],
#' [Plot.Spectrum()]).
#'
#' The model is based on the dynamic, coupled benthic-pelagic size-spectrum
#' model described in the references below. Run `citation("dbpmr")` for these in
#' BibTeX form.
#'
#' @references
#' Blanchard, J.L., Jennings, S., Law, R., Castle, M.D., McCloghrie, P.,
#' Rochet, M.-J. & Benoit, E. (2009). How does abundance scale with body size in
#' coupled size-structured food webs? \emph{Journal of Animal Ecology},
#' \strong{78}(1), 270-280. \doi{10.1111/j.1365-2656.2008.01466.x}
#'
#' Blanchard, J.L., Law, R., Castle, M.D. & Jennings, S. (2011). Coupled energy
#' pathways and the resilience of size-structured food webs. \emph{Theoretical
#' Ecology}, \strong{4}(3), 289-300. \doi{10.1007/s12080-010-0078-9}
#'
#' Castle, M.D., Blanchard, J.L. & Jennings, S. (2011). Predicted effects of
#' behavioural movement and passive transport on individual growth and community
#' size structure in marine ecosystems. \emph{Advances in Ecological Research},
#' \strong{45}, 41-66. \doi{10.1016/B978-0-12-386475-8.00002-2}
#'
#' @keywords internal
#' @useDynLib dbpmr, .registration = TRUE, .fixes = "C_"
#' @import methods
#' @importFrom utils read.csv write.table
#' @importFrom graphics persp plot points
"_PACKAGE"
