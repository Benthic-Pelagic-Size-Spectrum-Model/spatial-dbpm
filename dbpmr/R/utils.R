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
#' stores them as `ts` in the `arrays` slot. Used by [Setup.ts()].
#'
#' @param species A plankton parameter object with `ts_flag = TRUE`.
#' @param grid.in A [grid.params] object.
#' @param intercept,slope Arrays `[x, y, layer]`.
#' @param layer,time See [Setup.ts()].
#' @return The plankton object with `arrays$ts` filled.
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
  species@arrays$ts <- list(intercept = intercept, slope = slope, layer = as.integer(layer))
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

#' Attach in-memory fishing mortality to a species object (internal)
#'
#' Validates a fishing-mortality rate given as a scalar, per-layer vector,
#' `[x, y]` matrix or `[x, y, layer]` array, the fished size range and the
#' time-point-to-layer mapping, and stores them as `fishing` in the `arrays` slot.
#' Used by [Setup.fishing()].
#'
#' @param species A pelagic or benthic parameter object with `fishing_flag = TRUE`.
#' @param grid.in A [grid.params] object.
#' @param rate,min_size,max_size,layer,time See [Setup.fishing()].
#' @return The species object with `arrays$fishing` filled.
#' @keywords internal
#' @noRd
fishing_arrays <- function(species, grid.in, rate, min_size = -Inf, max_size = Inf, layer = NULL, time = NULL) {
  a <- xyl_array(rate, "rate", grid.in)
  if (!is.numeric(min_size) || !is.numeric(max_size) || length(min_size) != 1 || length(max_size) != 1 || is.na(min_size) || is.na(max_size) || min_size > max_size)
    stop("min_size and max_size must be single numbers with min_size <= max_size")
  species@arrays$fishing <- list(rate = a, layer = layer_index(grid.in, dim(a)[3], layer, time),
                                 range = c(as.numeric(min_size), as.numeric(max_size)))
  species
}

#' A non-negative input as an `[x, y, layer]` array (internal)
#'
#' Shared by the array inputs that vary in space and time: a single number, one
#' value per layer, an `[x, y]` matrix or an `[x, y, layer]` array becomes an
#' `[x, y, layer]` array of doubles.
#'
#' @param a The values.
#' @param nm Name used in error messages.
#' @param grid.in A [grid.params] object.
#' @return A numeric array `[x, y, layer]`.
#' @keywords internal
#' @noRd
xyl_array <- function(a, nm, grid.in) {
  nx <- length(grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep))
  ny <- length(grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep))
  if (!is.numeric(a)) stop(nm, " must be numeric")
  if (any(!is.finite(a))) stop(nm, " must not contain NA, NaN or Inf")
  if (any(a < 0)) stop(nm, " must be non-negative")
  d <- dim(a)
  out <- if (is.null(d)) array(rep(as.numeric(a), each = nx * ny), c(nx, ny, length(a)))    # scalar or one value per layer
         else if (length(d) == 2L) {
           if (d[1] != nx || d[2] != ny) stop(sprintf("%s: a matrix must be [x = %d, y = %d]; got [%s]", nm, nx, ny, paste(d, collapse = ", ")))
           array(as.numeric(a), c(nx, ny, 1L))
         } else {
           if (length(d) != 3L || d[1] != nx || d[2] != ny)
             stop(sprintf("%s: an array must be [x = %d, y = %d, layers]; got [%s]", nm, nx, ny, paste(d, collapse = ", ")))
           a
         }
  storage.mode(out) <- "double"
  out
}

#' Attach an in-memory starting state to a species object (internal)
#'
#' Validates a starting abundance (plankton, pelagic or benthic species: one
#' value per size class of the species, the same in every cell, or an
#' `[x, y, size class]` array) or starting biomass (detritus: one number or an
#' `[x, y]` matrix) and stores it as `initial` in the `arrays` slot. Used by
#' [Setup.ts()].
#'
#' @param species A parameter object with `initial_flag = TRUE` and `ts_flag = FALSE`.
#' @param grid.in A [grid.params] object.
#' @param initial See [Setup.ts()].
#' @return The object with `arrays$initial` filled.
#' @keywords internal
#' @noRd
initial_array <- function(species, grid.in, initial) {
  if (!isTRUE(species@initial_flag)) stop("initial supplies starting values: create ", species@filename, " with initial_flag = TRUE")
  if (isTRUE(species@ts_flag)) stop("initial supplies starting values only; ", species@filename, " has ts_flag = TRUE (a full time series)")
  nx <- length(grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep))
  ny <- length(grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep))
  if (!is.numeric(initial)) stop("initial must be numeric")
  if (any(!is.finite(initial))) stop("initial must not contain NA, NaN or Inf")
  if (any(initial < 0)) stop("initial must be non-negative")
  d <- dim(initial)
  if (species@speciestype == "detritus") {
    a <- if (length(initial) == 1L) matrix(as.numeric(initial), nx, ny)
         else if (is.null(d) && ny == 1L && length(initial) == nx) matrix(as.numeric(initial), nx, 1L)
         else if (length(d) == 2L && d[1] == nx && d[2] == ny) matrix(as.numeric(initial), nx, ny)
         else stop(sprintf("initial for detritus must be a single number or a matrix [x = %d, y = %d]", nx, ny))
  } else {
    ns <- length(grid_seq(species@mmin, species@mmax, grid.in@mstep))
    if (length(d) == 2L && d[1] == 1L) { initial <- as.vector(initial); d <- NULL }             # one row: the same in every cell
    a <- if (is.null(d)) {
           if (length(initial) != ns) stop(sprintf("initial must have one value per size class of %s (%d); got %d", species@filename, ns, length(initial)))
           array(rep(as.numeric(initial), each = nx * ny), c(nx, ny, ns))
         } else if (length(d) == 2L && ny == 1L && d[1] == nx && d[2] == ns) array(as.numeric(initial), c(nx, 1L, ns))   # 1-D run: [x, size class]
         else if (length(d) == 3L && d[1] == nx && d[2] == ny && d[3] == ns) initial
         else stop(sprintf("initial must be a vector with one value per size class of %s (%d) or an array [x = %d, y = %d, %d]; got [%s]",
                           species@filename, ns, nx, ny, ns, paste(d, collapse = ", ")))
  }
  storage.mode(a) <- "double"
  species@arrays$initial <- a
  species
}
