#' Write model output arrays to a NetCDF file
#'
#' Saves the output returned by [SizeSpectrum()] with `output = "arrays"` (or
#' `"both"`) as a compressed NetCDF-4 file with one variable per group.
#'
#' @details
#' The file has the dimensions `x`, `y`, `size` (natural-log body mass in grams)
#' and `time` (model years). Every variable in the output object is written with
#' dimensions `[x, y, size, time]` and named `<group>_<variable>`, where the
#' group is the plankton or species file name: `<group>_abundance` (abundance
#' density per unit natural-log body mass), `<group>_growth`,
#' `<group>_mortality`, `<group>_fishing`, `<group>_catch` and, for predation mortality,
#' `<prey>_predation_by_<predator>` for a pelagic prey species and
#' `<group>_predation` for a benthic one. Detritus biomass, when present, is
#' written as `detritus` with dimensions `[x, y, time]`.
#'
#' By default `x` and `y` are the model grid coordinates. Supply `x` and `y`
#' (and `xname`, `yname`, `xunits`, `yunits`) to label the cells differently,
#' for example with longitude and latitude when the model grid is in metres.
#'
#' Requires the \pkg{ncdf4} package.
#'
#' @param output A `dbpmr_output` object from [SizeSpectrum()].
#' @param file Path of the NetCDF file to create (overwritten if it exists).
#' @param x,y Coordinate values for the output cells; default to the model grid
#'   coordinates in `output`.
#' @param xname,yname Names of the two horizontal dimensions.
#' @param xunits,yunits Units of the two horizontal dimensions.
#' @param compression Deflate level, `0` (none) to `9`.
#'
#' @return Invisibly, `file`.
#' @seealso [SizeSpectrum()]
#' @examples
#' \dontrun{
#' out <- SizeSpectrum(run, grid, pl, pe, be, de, output = "arrays")
#' Write.NetCDF(out, "run.nc")
#' }
#' @export
`Write.NetCDF` <-
function(output, file, x = output$x, y = output$y, xname = "x", yname = "y",
         xunits = "", yunits = "", compression = 4){
  if(!inherits(output, "dbpmr_output")) stop("output must come from SizeSpectrum(output = \"arrays\") or \"both\"")
  if(!requireNamespace("ncdf4", quietly = TRUE)) stop("Write.NetCDF() needs the ncdf4 package: install.packages(\"ncdf4\")")
  if(length(x) != length(output$x) || length(y) != length(output$y))
    stop(sprintf("x and y must have %d and %d values, one per output cell", length(output$x), length(output$y)))
  dx <- ncdf4::ncdim_def(xname, xunits, as.numeric(x))
  dy <- ncdf4::ncdim_def(yname, yunits, as.numeric(y))
  dm <- ncdf4::ncdim_def("size", "ln(g)", output$m, longname = "natural-log body mass")
  dt <- ncdf4::ncdim_def("time", "years", output$t, longname = "model time")
  # flatten to one named [x, y, size, time] array per NetCDF variable, with its units and description
  arrs <- list(); meta <- list()
  add <- function(name, a, units, longname){ arrs[[name]] <<- a; meta[[name]] <<- c(units, longname) }
  for(g in names(output$abundance)) add(paste0(g, "_abundance"), output$abundance[[g]], "", paste(g, "abundance density per unit natural-log body mass"))
  for(g in names(output$growth))    add(paste0(g, "_growth"), output$growth[[g]], "yr-1", paste(g, "growth rate (rate of change of natural-log body mass)"))
  for(g in names(output$mortality)) add(paste0(g, "_mortality"), output$mortality[[g]], "yr-1", paste(g, "total mortality rate"))
  for(g in names(output$predation)){
    p <- output$predation[[g]]
    if(is.list(p)) for(pr in names(p)) add(paste0(g, "_predation_by_", pr), p[[pr]], "yr-1", paste(g, "predation mortality rate from", pr))
    else add(paste0(g, "_predation"), p, "yr-1", paste(g, "predation mortality rate"))
  }
  for(g in names(output$fishing))   add(paste0(g, "_fishing"), output$fishing[[g]], "yr-1", paste(g, "fishing mortality rate"))
  for(g in names(output$catch))     add(paste0(g, "_catch"), output$catch[[g]], "", paste(g, "catch rate density per year (fishing mortality x abundance x body mass) per unit natural-log body mass"))
  if(anyDuplicated(names(arrs))) stop("variable names must be unique: ", paste(names(arrs)[duplicated(names(arrs))], collapse = ", "))
  vars <- lapply(names(arrs), function(nm)
    ncdf4::ncvar_def(nm, meta[[nm]][1], list(dx, dy, dm, dt), missval = NA_real_, prec = "double", compression = compression,
                     longname = meta[[nm]][2]))
  names(vars) <- names(arrs)
  if(!is.null(output$detritus))
    vars$detritus <- ncdf4::ncvar_def("detritus", "", list(dx, dy, dt), missval = NA_real_, prec = "double",
                                      compression = compression, longname = "detritus biomass")
  nc <- ncdf4::nc_create(file, vars, force_v4 = TRUE)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  for(nm in names(arrs)) ncdf4::ncvar_put(nc, vars[[nm]], arrs[[nm]])
  if(!is.null(output$detritus)) ncdf4::ncvar_put(nc, vars$detritus, output$detritus)
  ncdf4::ncatt_put(nc, 0, "title", "dbpmr model output")
  ncdf4::ncatt_put(nc, 0, "source", paste("dbpmr", as.character(utils::packageVersion("dbpmr"))))
  invisible(file)
}

#' @export
print.dbpmr_output <- function(x, ...){
  cat(sprintf("dbpmr output arrays: %d x %d cells, %d size classes, %d output times (t = %g to %g)\n",
              length(x$x), length(x$y), length(x$m), length(x$t), min(x$t), max(x$t)))
  cat("  groups:", paste(names(x$abundance), collapse = ", "), "\n")
  cat("  variables [x, y, m, t]:", paste(intersect(c("abundance", "growth", "mortality", "predation", "fishing", "catch"), names(x)), collapse = ", "), "\n")
  if(!is.null(x$detritus)) cat("  detritus [x, y, t]\n")
  invisible(x)
}
