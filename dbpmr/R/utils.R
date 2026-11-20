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
