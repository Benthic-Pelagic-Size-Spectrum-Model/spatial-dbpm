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
