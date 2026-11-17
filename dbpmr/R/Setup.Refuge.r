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
#' With `file = FALSE` the refuge is held in memory: no file is written, the
#' floor `min_A` travels with it (the environment variable is not used), and the
#' function returns the species object, which must be assigned and passed on to
#' [SizeSpectrum()]: `pe <- Setup.Refuge(pe, run, grid, refuge, file = FALSE)`.
#' Held in memory, the refuge may differ between cells.
#'
#' @param refuge A numeric vector of refuge-space density (count m^-2), one value
#'   per mass class on the model's size grid, or a function of natural-log body
#'   mass `m` returning the refuge density at `m`. With `file = FALSE` it may also
#'   be an array `[x, y, mass class]`, giving each cell its own refuge.
#' @param min_A The availability floor \eqn{\mathrm{min\_A}} (the minimum
#'   proportion of a size class that remains vulnerable to predation; `0`–`1`).
#'   Default `0`. Sets the `REFUGE_MIN_A` environment variable read by the
#'   engine (process-global; affects every species with a refuge in this run).
#'   With `file = FALSE` it applies to this species only.
#' @param file `TRUE` (the default) writes the refuge to its input file.
#'   `FALSE` holds it in memory and returns the species object.
#'
#' @return With `file = TRUE`, invisibly `NULL`; called for the side effect of
#'   writing `<run>/Input/<species>_refuge.txt` and setting `REFUGE_MIN_A`. With
#'   `file = FALSE`, invisibly the species object with the refuge attached
#'   (assign it).
#' @references Rogers A, Blanchard JL, Mumby PJ (2014) Vulnerability of coral
#'   reef fisheries to a loss of structural complexity. *Current Biology*
#'   24(9): 1000--1005. \doi{10.1016/j.cub.2014.03.026}
#' @seealso [Setup.Pelagic()], [Setup.ts()]
#' @export
`Setup.Refuge` <-
function(species, run.in, grid.in, refuge, min_A = 0, file = TRUE){
  if(species@speciestype != "pelagic")
    stop("Setup.Refuge currently applies to a pelagic species only")
  if(min_A < 0 || min_A > 1) stop("min_A must be in [0, 1]")
  mass <- grid_seq(grid.in@mmin, grid.in@mmax, grid.in@mstep)
  if(length(dim(refuge)) == 3L){
    # a refuge for each cell, [x, y, mass class]: held in memory only
    if(isTRUE(file)) stop("a refuge that differs between cells is held in memory: use file = FALSE and assign the result")
    nx <- length(grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep))
    ny <- length(grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep))
    if(!is.numeric(refuge) || !identical(as.integer(dim(refuge)), c(nx, ny, length(mass))))
      stop(sprintf("a refuge array must be [x = %d, y = %d, mass classes = %d]; got [%s]", nx, ny, length(mass), paste(dim(refuge), collapse = ", ")))
    rv <- refuge
    rv[!is.finite(rv) | rv < 0] <- 0
    storage.mode(rv) <- "double"
    species@arrays$refuge <- list(values = rv, min_A = as.numeric(min_A))
    return(invisible(species))
  }
  if(is.function(refuge)){
    rv <- vapply(mass, refuge, numeric(1))
  } else {
    rv <- as.numeric(refuge)
    if(length(rv) != length(mass))
      stop(sprintf("refuge must have one value per mass class (%d), got %d", length(mass), length(rv)))
  }
  rv[!is.finite(rv) | rv < 0] <- 0
  if(!isTRUE(file)){
    # held in memory, the same in every cell
    nx <- length(grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep))
    ny <- length(grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep))
    species@arrays$refuge <- list(values = array(rep(rv, each = nx * ny), c(nx, ny, length(mass))), min_A = as.numeric(min_A))
    return(invisible(species))
  }
  dir.create(file.path(run.in@filename, "Input"), showWarnings = FALSE, recursive = TRUE)
  writeLines(formatC(rv, format = "g", digits = 10),
             file.path(run.in@filename, "Input", paste0(species@filename, "_refuge.txt")))
  Sys.setenv(REFUGE_MIN_A = as.character(min_A))
  invisible(NULL)
}
