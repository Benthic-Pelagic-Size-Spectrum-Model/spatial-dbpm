#' Build a grid axis robustly (internal)
#'
#' Like `seq(from, to, by)` but uses an explicit, rounded point count so that
#' floating-point error in `by` cannot drop the final point. This matters when
#' the step is irrational (e.g. a log10 grid expressed in natural log, where
#' `by = 0.1 * log(10)`): plain `seq()` can accumulate enough error that the last
#' value just exceeds `to` and is dropped, yielding one fewer bin than the C
#' engine, which sizes its grid as `round((to - from) / by) + 1`.
#'
#' Built on R's own [seq()] so the floating-point values match what the rest of
#' the package expects (several functions locate sizes with exact `==`); the only
#' change is to append the final point(s) when `seq()` has dropped them.
#'
#' @param from,to,by Axis start, end and step.
#' @return Numeric vector of `round((to - from) / by) + 1` evenly spaced points.
#' @keywords internal
#' @noRd
grid_seq <- function(from, to, by) {
  s <- seq(from, to, by)
  n <- round((to - from) / by) + 1L
  if (length(s) < n) s <- c(s, s[length(s)] + by * seq_len(n - length(s)))
  s
}

#' Locate a size on the mass grid robustly (internal)
#'
#' Returns the index of the grid point nearest `target`. Replaces the fragile
#' `which(mass == target)`, which relies on exact floating-point equality: on an
#' irrational grid (a log10 grid expressed in natural log) a size computed
#' directly (e.g. `-3 * log(10)`) and the same size accumulated by [seq()] can
#' differ in the last bit, so `==` finds nothing. The match is the nearest grid
#' point, validated to lie within a small fraction of the target (well inside one
#' bin), so a genuinely off-grid size still errors rather than snapping silently.
#'
#' @param mass Numeric mass-grid vector.
#' @param target Size to locate.
#' @return Integer index into `mass`.
#' @keywords internal
#' @noRd
mass_index <- function(mass, target) {
  i <- which.min(abs(mass - target))
  if (length(i) == 0L || abs(mass[i] - target) > 1e-6 * max(1, abs(target)))
    stop("size ", target, " does not align with the mass grid", call. = FALSE)
  i
}

#' Attach an in-memory plankton time series to a plankton object (internal)
#'
#' Validates intercept/slope arrays and the time-point-to-layer mapping and
#' stores them in the `ts_arrays` slot. Used by [Setup.ts()].
#'
#' @param species A plankton parameter object with `ts_flag = TRUE`.
#' @param grid.in A [grid.params] object.
#' @param intercept,slope Arrays `[x, y, layer]`.
#' @param layer,time See [Setup.ts()].
#' @return The plankton object with `ts_arrays` filled.
#' @keywords internal
#' @noRd
plankton_ts_arrays <- function(species, grid.in, intercept, slope, layer = NULL, time = NULL) {
  if (species@speciestype != "plankton")
    stop("intercept and slope arrays are currently supported for plankton only")
  if (!isTRUE(species@ts_flag))
    stop("intercept and slope arrays supply a time series: create the plankton with ts_flag = TRUE")
  if (is.null(intercept) || is.null(slope)) stop("both intercept and slope must be given")
  nx <- length(grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep))
  ny <- length(grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep))
  as3 <- function(a, nm) {
    if (!is.numeric(a)) stop(nm, " must be numeric")
    if (is.null(dim(a))) a <- array(a, c(1L, 1L, length(a)))                       # 0-D: one value per layer
    else if (length(dim(a)) == 2L) a <- array(a, c(dim(a)[1], 1L, dim(a)[2]))      # 1-D: [x, layer]
    if (length(dim(a)) != 3L || dim(a)[1] != nx || dim(a)[2] != ny)
      stop(sprintf("%s must have dimensions [x = %d, y = %d, layers]; got [%s]", nm, nx, ny, paste(dim(a), collapse = ", ")))
    if (any(!is.finite(a))) stop(nm, " must not contain NA, NaN or Inf")
    storage.mode(a) <- "double"
    a
  }
  intercept <- as3(intercept, "intercept"); slope <- as3(slope, "slope")
  if (!identical(dim(intercept), dim(slope))) stop("intercept and slope must have the same dimensions")
  if (any(intercept < 0))
    stop("intercept must be non-negative: it is the abundance density at m = 0, not its logarithm")
  nl <- dim(intercept)[3]
  layer <- layer_index(grid.in, nl, layer, time)
  species@ts_arrays <- list(intercept = intercept, slope = slope, layer = as.integer(layer))
  species
}

#' Layer used at each model time point (internal)
#'
#' Shared by the array inputs: resolves `layer` (one entry per time point) or
#' `time` (start time of each layer) into a whole-number layer index, 1-based,
#' for each time point `seq(0, tmax, tstep)`.
#'
#' @param grid.in A [grid.params] object.
#' @param nl Number of layers in the arrays.
#' @param layer,time See [Setup.ts()].
#' @return Integer vector, one layer per model time point.
#' @keywords internal
#' @noRd
layer_index <- function(grid.in, nl, layer = NULL, time = NULL) {
  trange <- seq(0, grid.in@tmax, grid.in@tstep)
  nb <- length(trange)
  if (!is.null(layer) && !is.null(time)) stop("give only one of layer or time")
  if (!is.null(time)) {
    if (length(time) != nl || is.unsorted(time, strictly = TRUE)) stop("time must give one strictly increasing start time per layer")
    if (time[1] > 0) stop("time[1] must be at or before the start of the run (0)")
    # a time point that coincides with a layer start (to within rounding) takes that layer
    layer <- findInterval(trange + grid.in@tstep * 1e-6, time)
  }
  if (is.null(layer)) {
    if (nl == nb) layer <- seq_len(nb)
    else if (nl == 1L) layer <- rep(1L, nb)
    else stop(sprintf(paste0("the arrays have %d layers but the run has %d time points (0 to tmax by tstep): ",
                             "give `layer` (the layer each time point uses) or `time` (the start time of each layer)"), nl, nb))
  }
  if (length(layer) != nb) stop(sprintf("layer must have one entry per time point (%d), got %d", nb, length(layer)))
  if (any(!is.finite(layer)) || any(layer != round(layer)) || any(layer < 1) || any(layer > nl))
    stop(sprintf("layer must be whole numbers between 1 and the number of layers (%d)", nl))
  as.integer(layer)
}
