#' Allocate fishing effort among cells by gravity
#'
#' In a spatial run, moves fishing effort between cells during the run, towards
#' the cells where the accessible catch is highest. The result is passed to
#' [SizeSpectrum()] through its `gravity` argument.
#'
#' @details
#' Fishing has two parts: how much there is at each time, and where it falls.
#' The first is always set with [Setup.fishing()]. Without `gravity`, the second
#' is set there too, by giving a rate that differs between cells. With `gravity`,
#' the rate given to [Setup.fishing()] is taken as the fishing mortality when
#' effort is spread evenly over the grid (typically the same in every cell:
#' catchability times the regional effort at that time), and the model decides
#' where it falls.
#'
#' The share of effort of a cell is proportional to
#'
#' \deqn{access \times \sum_{species} thickness \times \int F\,u\,w\;dm,}
#'
#' the accessible catch rate per unit area that evenly spread effort would take
#' there, summed over all fished species (`F` the current rate from
#' [Setup.fishing()], `u` the abundance at the latest update, `w` body mass). The
#' fishing mortality of every fished species in the cell is the rate from
#' [Setup.fishing()] times the number of cells times that share. The shares sum
#' to one, so the multiplier averages one over the grid: effort is moved between
#' cells, not added.
#'
#' Effort follows changes in abundance at each update (every `update` years). A
#' change in the rates themselves takes effect straight away: a cell where the
#' rate is zero (a closure, such as a marine protected area) or `access` is zero
#' gets no effort, and the effort it would have had goes to the other cells from
#' the moment the closure starts. Closures that vary in space and time are
#' therefore given as zeros in the `[x, y, layer]` rate of [Setup.fishing()].
#' While there is no catch to be had anywhere (for example before fishing starts)
#' effort stays evenly spread.
#'
#' All cells count equally, so the grid cells are assumed to have similar areas.
#'
#' This is not the `gravity` channel of [Setup.Forcing()], which divides fishing
#' between the pelagic and benthic groups within a cell.
#'
#' @param run.in A [run.params] object.
#' @param grid.in A [grid.params] object.
#' @param access Accessibility of each cell: a single number or a matrix `[x, y]`
#'   of non-negative weights (for example falling with depth and with distance
#'   from port). Only its relative size between cells matters.
#' @param pelagic_thickness,benthic_thickness Thickness (m) of the layer each
#'   group occupies, a single number or a matrix `[x, y]`. They turn the model's
#'   densities per volume into densities per area so that the catch of pelagic
#'   and benthic species can be added.
#' @param update Time (years) between updates of the abundances the shares are
#'   based on, counted from the start of the run. `0` updates at every time
#'   step.
#'
#' @return An object of class `dbpmr_gravity` to pass to [SizeSpectrum()] as
#'   `gravity`.
#' @seealso [Setup.fishing()], [SizeSpectrum()]
#' @examples
#' \dontrun{
#' run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
#' grid <- Setup.Grid(run, tmax = 5, tstep = 1/48, toutstep = 1,
#'                    xmin = 0, xmax = 2, xstep = 1, xoutstep = 1,
#'                    ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
#' pe <- Setup.Pelagic(run, filename = "fish", fishing_flag = TRUE)
#' pe <- Setup.fishing(pe, run, grid, rate = 0.3, min_size = log(100))   # the same rate in every cell
#' gv <- Setup.Gravity(run, grid, access = matrix(c(1, 1, 0.5, 0.5, 0.1, 0.1), 3, 2), update = 1)
#' # SizeSpectrum(run, grid, pl, pe, be, de, gravity = gv)
#' }
#' @export
`Setup.Gravity` <-
function(run.in, grid.in, access = 1, pelagic_thickness = 1, benthic_thickness = 1, update = 1){
  if(run.in@spatial_dim == 0) stop("gravity moves effort between cells and needs a spatial run (spatial_dim of 1 or 2)")
  nx <- length(grid_seq(grid.in@xmin, grid.in@xmax, grid.in@xstep))
  ny <- if(run.in@spatial_dim == 2) length(grid_seq(grid.in@ymin, grid.in@ymax, grid.in@ystep)) else 1L
  field <- function(a, nm, positive){
    if(!is.numeric(a)) stop(nm, " must be numeric")
    if(any(!is.finite(a))) stop(nm, " must not contain NA, NaN or Inf")
    if(positive && any(a <= 0)) stop(nm, " must be positive")
    if(!positive && any(a < 0)) stop(nm, " must be non-negative")
    d <- dim(a)
    if(is.null(d) && length(a) == 1L) return(matrix(as.numeric(a), nx, ny))
    if(is.null(d) && ny == 1L && length(a) == nx) return(matrix(as.numeric(a), nx, 1L))
    if(length(d) != 2L || d[1] != nx || d[2] != ny)
      stop(sprintf("%s must be a single number or a matrix [x = %d, y = %d]; got %s", nm, nx, ny,
                   if(is.null(d)) sprintf("a vector of length %d", length(a)) else sprintf("[%s]", paste(d, collapse = ", "))))
    matrix(as.numeric(a), nx, ny)
  }
  acc <- field(access, "access", FALSE)
  if(all(acc == 0)) stop("access must be greater than zero in at least one cell")
  if(!is.numeric(update) || length(update) != 1L || !is.finite(update) || update < 0) stop("update must be a single non-negative number (years)")
  nup <- max(1L, as.integer(round(update / grid.in@tstep)))
  if(update > 0 && abs(update / grid.in@tstep - nup) > 1e-7) stop("update must be a whole number of time steps (tstep = ", format(grid.in@tstep), ")")
  structure(list(access = acc, pelagic_thickness = field(pelagic_thickness, "pelagic_thickness", TRUE),
                 benthic_thickness = field(benthic_thickness, "benthic_thickness", TRUE),
                 nupdate = nup, update = nup * grid.in@tstep),
            class = "dbpmr_gravity")
}

#' @export
print.dbpmr_gravity <- function(x, ...){
  d <- dim(x$access)
  cat(sprintf("dbpmr gravity allocation of fishing effort on a %d x %d grid\n", d[1], d[2]))
  cat(sprintf("  effort follows abundance every %s\n", if(x$nupdate == 1L) "time step" else sprintf("%.4g years", x$update)))
  cat(sprintf("  access: range %.4g to %.4g (%d cells with no access)\n", min(x$access), max(x$access), sum(x$access == 0)))
  invisible(x)
}
