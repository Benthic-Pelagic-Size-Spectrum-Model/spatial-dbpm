#' Supply a size-structured prey refuge (reef structural complexity)
#'
#' Writes the optional per-size *refuge* input that activates the reef
#' prey-refuge mechanism of Rogers, Blanchard & Mumby (2014) for a pelagic
#' species. The engine converts the refuge into a per-size *availability*
#' multiplier
#' \deqn{a(w) = \max(\mathrm{min\_A},\; 1 - \mathrm{refuge}(w)/\mathrm{competitor}(w)),}
#' where the competitor pool is the species' own abundance density. `a(w)`
#' scales both the predation mortality on, and predator intake of, prey of body
#' size `w`, so sheltered prey are both eaten less and, as predators, hunt less
#' --- reproducing the paper's structural-complexity effect on the size spectrum.
#'
#' With no refuge supplied (this function not called, or an all-zero refuge) the
#' mechanism is inactive and the engine is bit-identical to standard `dbpmr`.
#'
#' @param species A pelagic species parameter object (from [Setup.Pelagic()]).
#' @param run.in A [run.params] object.
#' @param grid.in A [grid.params] object.
#' @param refuge Either a numeric vector of refuge-space density (count m^-2),
#'   one value per mass class on the model's size grid, or a function of log10...
#'   actually of natural-log body mass `m` returning the refuge density at `m`.
#' @param min_A The availability floor \eqn{\mathrm{min\_A}} (the minimum
#'   proportion of a size class that remains vulnerable to predation; `0`–`1`).
#'   Default `0`. Sets the `REFUGE_MIN_A` environment variable read by the
#'   engine (process-global; affects every species with a refuge in this run).
#'
#' @return Invisibly `NULL`; called for the side effect of writing
#'   `<run>/Input/<species>_refuge.txt` and setting `REFUGE_MIN_A`.
#' @references Rogers A, Blanchard JL, Mumby PJ (2014) Vulnerability of coral
#'   reef fisheries to a loss of structural complexity. *Current Biology*
#'   24(9): 1000--1005. \doi{10.1016/j.cub.2014.03.026}
#' @seealso [Setup.Pelagic()], [Setup.ts()]
#' @export
`Setup.Refuge` <-
function(species, run.in, grid.in, refuge, min_A = 0){
  if(species@speciestype != "pelagic")
    stop("Setup.Refuge currently applies to a pelagic species only")
  if(min_A < 0 || min_A > 1) stop("min_A must be in [0, 1]")
  mass <- grid_seq(grid.in@mmin, grid.in@mmax, grid.in@mstep)
  if(is.function(refuge)){
    rv <- vapply(mass, refuge, numeric(1))
  } else {
    rv <- as.numeric(refuge)
    if(length(rv) != length(mass))
      stop(sprintf("refuge must have one value per mass class (%d), got %d", length(mass), length(rv)))
  }
  rv[!is.finite(rv) | rv < 0] <- 0
  dir.create(file.path(run.in@filename, "Input"), showWarnings = FALSE, recursive = TRUE)
  writeLines(formatC(rv, format = "g", digits = 10),
             file.path(run.in@filename, "Input", paste0(species@filename, "_refuge.txt")))
  Sys.setenv(REFUGE_MIN_A = as.character(min_A))
  invisible(NULL)
}
