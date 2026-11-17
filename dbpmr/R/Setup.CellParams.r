#' Supply per-cell model parameters for a spatial run
#'
#' By default every model parameter of a species is a single number that applies
#' to the whole spatial grid. `Setup.CellParams()` lets a parameter vary from cell
#' to cell, so that a gridded configuration that is normally run as many
#' independent 0-D cells (each with its own parameter value) can be run as one
#' spatial grid.
#'
#' Currently one parameter can be supplied per cell: `pref_ben`, the preference
#' of a pelagic species for benthic prey, which gridded runs typically make a
#' function of water depth. Where it is supplied, the per-cell value replaces the
#' scalar `pref_ben` given to [Setup.Pelagic()], both in the pelagic species'
#' intake of benthic prey and in the predation mortality it imposes on the
#' benthos.
#'
#' The values are static (read once at the start of the run). If this function is
#' not called the scalar parameter is used in every cell and results are
#' unchanged.
#'
#' With `file = FALSE` the values are held in memory: no file is written and the
#' function returns the species object, which must be assigned and passed on to
#' [SizeSpectrum()]: `pe <- Setup.CellParams(pe, run, grid, pref_ben = p, file = FALSE)`.
#'
#' @param species A pelagic species parameter object (from [Setup.Pelagic()]).
#' @param run.in A [run.params] object.
#' @param grid.in A [grid.params] object.
#' @param pref_ben Either a numeric matrix with one row per x cell and one column
#'   per y cell (`pref_ben[i, j]` is the value at the `i`-th x and `j`-th y grid
#'   position; for a 1-D or 0-D run it has a single column), or a function of
#'   `(x, y)` returning the value at that grid position. Values must be finite and
#'   non-negative.
#' @param file `TRUE` (the default) writes the values to their input file.
#'   `FALSE` holds them in memory and returns the species object.
#'
#' @return With `file = FALSE`, invisibly the species object with the values
#'   attached (assign it). With `file = TRUE`,
#'   invisibly the `x` by `y` matrix that was written; called for the side
#'   effect of writing `<run>/Input/<species>_cellpar.txt` (a header line naming
#'   the columns, then one row per cell with x as the outer and y as the inner
#'   index, the same cell order as [Setup.ts()] and per-cell forcing).
#' @seealso [Setup.Pelagic()], [Setup.ts()]
#' @examples
#' \dontrun{
#' run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
#' grid <- Setup.Grid(run, tmax = 5, tstep = 1/48, toutstep = 1,
#'                    xmin = 0, xmax = 2, xstep = 1, xoutstep = 1,
#'                    ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
#' pe   <- Setup.Pelagic(run, filename = "fish")
#' depth <- matrix(c(200, 800, 3000, 150, 1200, 4000), nrow = 3)   # [x, y]
#' Setup.CellParams(pe, run, grid, pref_ben = 0.8 * exp(-depth / 1500))
#' }
#' @export
`Setup.CellParams` <-
function(species, run.in, grid.in, pref_ben, file = TRUE){
  if(species@speciestype != "pelagic")
    stop("Setup.CellParams currently applies to a pelagic species only")
  if(missing(pref_ben)) stop("pref_ben must be given")
  xrange <- grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep)
  yrange <- grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep)
  nx <- length(xrange); ny <- length(yrange)
  if(is.function(pref_ben)){
    pb <- outer(xrange, yrange, Vectorize(pref_ben))
  } else {
    pb <- pref_ben
    if(!is.matrix(pb)) pb <- matrix(as.numeric(pb), ncol = 1)
    if(nrow(pb) != nx || ncol(pb) != ny)
      stop(sprintf("pref_ben must be a %d x %d matrix (x cells by y cells), got %d x %d",
                   nx, ny, nrow(pb), ncol(pb)))
  }
  if(any(!is.finite(pb)) || any(pb < 0)) stop("pref_ben values must be finite and non-negative")
  if(!isTRUE(file)){
    storage.mode(pb) <- "double"
    species@arrays$pref_ben <- pb
    return(invisible(species))
  }
  dir.create(file.path(run.in@filename, "Input"), showWarnings = FALSE, recursive = TRUE)
  # one row per cell, x outer / y inner
  writeLines(c("pref_ben", formatC(as.vector(t(pb)), format = "g", digits = 17)),
             file.path(run.in@filename, "Input", paste0(species@filename, "_cellpar.txt")))
  invisible(pb)
}
