#' Supply environmental forcing as arrays
#'
#' Builds the environmental forcing of a run (temperature effects, export to the
#' sea floor, depth, ...) from R arrays held in memory. The result is passed to
#' [SizeSpectrum()] through its `forcing` argument; no `forcing_ts.txt` file is
#' written. In a spatial run every channel can differ between cells.
#'
#' @details
#' Channels are given as named arguments. The recognised names are those of the
#' forcing file:
#'
#' * `pel_tempeff`, `ben_tempeff`: temperature multipliers on pelagic and benthic
#'   feeding (and on background mortality unless the `TEMP_ON_MORT` environment
#'   variable is `0`);
#' * `sinking_rate`: fraction of surface-origin detritus reaching the sea floor;
#' * `depth`: water-column depth (m);
#' * `export_attn`, `export_gamma`, `export_magg`: size-dependent export;
#' * `residence_time`, `gravity`: must be the same in every cell.
#'
#' Each channel may be
#'
#' * a single number (constant in space and time);
#' * a vector with one value per layer (the same in every cell, varying in time);
#' * a matrix `[x, y]` (varying between cells, constant in time);
#' * an array `[x, y, layer]` (varying in both).
#'
#' All channels that vary in time must have the same number of layers. `layer`
#' or `time` says which layer each model time point uses, exactly as for the
#' plankton arrays in [Setup.ts()], so the same index can be used for both. The
#' forcing applied while the model advances to a time point is that time point's
#' layer.
#'
#' @param run.in A [run.params] object.
#' @param grid.in A [grid.params] object.
#' @param ... Forcing channels, as named arguments (see Details).
#' @param layer For each model time point (`seq(0, tmax, tstep)`) the layer it
#'   uses, as whole numbers from 1. May be omitted when there is one layer per
#'   time point, or when no channel varies in time.
#' @param time Instead of `layer`: the start time (in model years) of each
#'   layer; each time point uses the latest layer that has started.
#'
#' @return An object of class `dbpmr_forcing` to pass to [SizeSpectrum()] as
#'   `forcing`.
#' @seealso [SizeSpectrum()], [Setup.ts()], [Setup.CellParams()]
#' @examples
#' \dontrun{
#' run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
#' grid <- Setup.Grid(run, tmax = 5, tstep = 1/48, toutstep = 1,
#'                    xmin = 0, xmax = 2, xstep = 1, xoutstep = 1,
#'                    ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
#' month <- (floor(seq(0, 5, 1/48) * 12) %% 12) + 1            # layer used at each time point
#' temp  <- array(runif(3 * 2 * 12, 0.8, 1.2), c(3, 2, 12))    # [x, y, month]
#' fo <- Setup.Forcing(run, grid, pel_tempeff = temp, ben_tempeff = 1,
#'                     depth = matrix(200, 3, 2), layer = month)
#' # SizeSpectrum(run, grid, pl, pe, be, de, forcing = fo)
#' }
#' @export
`Setup.Forcing` <-
function(run.in, grid.in, ..., layer = NULL, time = NULL){
  chan <- list(...)
  known <- c("pel_tempeff", "ben_tempeff", "sinking_rate", "depth", "residence_time",
             "export_attn", "export_gamma", "export_magg", "gravity")
  if(length(chan) == 0) stop("at least one forcing channel must be given")
  if(is.null(names(chan)) || any(names(chan) == "") || anyDuplicated(names(chan)))
    stop("forcing channels must be given as uniquely named arguments")
  bad <- setdiff(names(chan), known)
  if(length(bad)) stop("unknown forcing channel: ", paste(bad, collapse = ", "),
                       " (recognised: ", paste(known, collapse = ", "), ")")
  nx <- length(grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep))
  ny <- length(grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep))

  # layers of each channel: 3-D arrays and plain vectors vary in time; scalars and [x, y] matrices do not
  nlay <- vapply(names(chan), function(nm){
    a <- chan[[nm]]
    if(!is.numeric(a)) stop(nm, " must be numeric")
    if(any(!is.finite(a))) stop(nm, " must not contain NA, NaN or Inf")
    d <- dim(a)
    if(is.null(d)) return(length(a))
    if(length(d) == 2L){
      if(d[1] != nx || d[2] != ny) stop(sprintf("%s: a matrix must be [x = %d, y = %d]; got [%s]", nm, nx, ny, paste(d, collapse = ", ")))
      return(1L)
    }
    if(length(d) != 3L || d[1] != nx || d[2] != ny)
      stop(sprintf("%s: an array must be [x = %d, y = %d, layers]; got [%s]", nm, nx, ny, paste(d, collapse = ", ")))
    d[3]
  }, numeric(1))
  nl <- unique(nlay[nlay > 1])
  if(length(nl) > 1) stop("channels that vary in time must all have the same number of layers; got ",
                          paste(sprintf("%s: %d", names(nlay)[nlay > 1], nlay[nlay > 1]), collapse = ", "))
  if(length(nl) == 0) nl <- 1L
  lay <- layer_index(grid.in, nl, layer, time)

  # every channel as [x, y, 1] or [x, y, layers]
  arrs <- lapply(names(chan), function(nm){
    a <- chan[[nm]]; d <- dim(a)
    out <- if(is.null(d)) array(rep(as.numeric(a), each = nx * ny), c(nx, ny, length(a)))
           else if(length(d) == 2L) array(as.numeric(a), c(nx, ny, 1L))
           else a
    storage.mode(out) <- "double"
    out
  })
  names(arrs) <- names(chan)
  for(nm in intersect(c("residence_time", "gravity"), names(arrs))){
    a <- arrs[[nm]]
    if(any(apply(a, 3, function(z) diff(range(z))) != 0)) stop(nm, " must be the same in every cell")
  }
  structure(list(channels = arrs, layer = lay), class = "dbpmr_forcing")
}

#' @export
print.dbpmr_forcing <- function(x, ...){
  d <- dim(x$channels[[1]])
  cat(sprintf("dbpmr forcing arrays on a %d x %d grid\n", d[1], d[2]))
  for(nm in names(x$channels)){
    a <- x$channels[[nm]]
    cat(sprintf("  %-15s %s, range %.4g to %.4g\n", nm,
                if(dim(a)[3] == 1) "constant in time" else sprintf("%d layers", dim(a)[3]), min(a), max(a)))
  }
  invisible(x)
}
