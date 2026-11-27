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
    function(
        output,
        file,
        x = output@x,
        y = output@y,
        xname = "x",
        yname = "y",
        xunits = "",
        yunits = "",
        compression = 4
    ) {
        if (!inherits(output, "memory.outputs")) {
            stop("output must come from SizeSpectrum(memory.return=TRUE)")
        }

        dm <- ncdf4::ncdim_def(
            "size",
            "ln(g)",
            output@m,
            longname = "natural-log body mass"
        )
        dx <- ncdf4::ncdim_def(xname, xunits, as.numeric(x))
        dy <- ncdf4::ncdim_def(yname, yunits, as.numeric(y))
        dt <- ncdf4::ncdim_def(
            "time",
            "years",
            output@time,
            longname = "model time"
        )

        d_pel <- ncdf4::ncdim_def("pelagic_species", "", output@pelagic_species)
        d_ben <- ncdf4::ncdim_def("benthic_species", "", output@benthic_species)
        # pelagic_mu_pred contains prey and predator dimensions that are just the pelagic species:
        d_pelprey <- ncdf4::ncdim_def(
            "pelagic_prey_species",
            "",
            output@pelagic_species
        )
        d_pelpred <- ncdf4::ncdim_def(
            "pelagic_pred_species",
            "",
            output@pelagic_species
        )

        # flatten to one named [x, y, size, time] array per NetCDF variable, with its units and description
        add_variable <- function(name, dims, longname, units = "") {
            ncdf4::ncvar_def(
                name,
                units,
                dims,
                missval = NA_real_,
                longname = longname,
                prec = "double",
                compression = compression
            )
        }

        vars <- list(
            plankton_u = add_variable(
                "plankton_u",
                list(dm, dx, dy, dt),
                longname = "abundance density per unit natural-log body mass"
            ),
            pelagic_u = add_variable(
                "pelagic_u",
                list(dm, dx, dy, d_pel, dt),
                longname = "abundance density per unit natural-log body mass"
            ),
            pelagic_mu_pred = add_variable(
                "pelagic_mu_pred",
                list(dm, dx, dy, d_pelprey, d_pelpred, dt),
                longname = "predation mortality rate",
                units = "yr-1"
            ),
            pelagic_mu_fish = add_variable(
                "pelagic_mu_fish",
                list(dm, dx, dy, d_pel, dt),
                longname = "fishing mortality rate",
                units = "yr-1"
            ),
            pelagic_g = add_variable(
                "pelagic_g",
                list(dm, dx, dy, d_pel, dt),
                longname = "growth rate (rate of change of natural-log body mass)",
                units = "yr-1"
            ),
            benthic_u = add_variable(
                "benthic_u",
                list(dm, dx, dy, d_ben, dt),
                longname = "abundance density per unit natural-log body mass"
            ),
            benthic_mu_pred = add_variable(
                "benthic_mu_pred",
                list(dm, dx, dy, d_ben, dt),
                longname = "predation mortality rate",
                units = "yr-1"
            ),
            benthic_mu_fish = add_variable(
                "benthic_mu_fish",
                list(dm, dx, dy, d_ben, dt),
                longname = "fishing mortality rate",
                units = "yr-1"
            ),
            benthic_g = add_variable(
                "benthic_g",
                list(dm, dx, dy, d_ben, dt),
                longname = "growth rate (rate of change of natural-log body mass)",
                units = "yr-1"
            ),
            detritus_w = add_variable(
                "detritus_w",
                list(dx, dy, dt),
                longname = "detritus_biomass"
            ),
            detritus_g = add_variable(
                "detritus_g",
                list(dx, dy, dt),
                longname = "detritus increase rate",
                units = "yr-1"
            ),
            detritus_mu = add_variable(
                "detritus_mu",
                list(dx, dy, dt),
                longname = "detritus_decrease_rate",
                units = "yr-1"
            ),
            pelagic_catch = add_variable(
                "pelagic_catch",
                list(dm, dx, dy, d_pel, dt),
                longname = "catch rate density per year (fishing mortality x abundance x body mass) per unit natural-log body mass"
            ),
            benthic_catch = add_variable(
                "benthic_catch",
                list(dm, dx, dy, d_ben, dt),
                longname = "catch rate density per year (fishing mortality x abundance x body mass) per unit natural-log body mass"
            )
        )

        nc <- ncdf4::nc_create(file, vars, force_v4 = TRUE)
        for (nm in names(vars)) {
            ncdf4::ncvar_put(nc, vars[[nm]], slot(out, nm))
        }

        ncdf4::ncatt_put(nc, 0, "title", "dbpmr model output")
        ncdf4::ncatt_put(
            nc,
            0,
            "source",
            paste("dbpmr", as.character(utils::packageVersion("dbpmr")))
        )

        ncdf4::nc_close(nc)

        return(invisible(file))
    }
