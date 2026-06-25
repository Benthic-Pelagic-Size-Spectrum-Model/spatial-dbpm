# spatial-dbpm

`dbpmr` is an R package implementing the **Dynamic Benthic Pelagic Model
(DBPM)** with spatial size-spectrum capabilities for marine ecosystems. It
couples a pelagic predator size spectrum to a benthic detritivore spectrum and a
detritus pool, and integrates them forward in time using a C simulation engine.

The model implemented here is based on the coupled size-spectrum models
described in Blanchard et al. (2009, 2011) and the spatial extension of
Castle et al. (2011). Note that the spatial movement (behavioural and passive
transport) of Castle et al. applies to the **pelagic** size spectrum only —
benthic organisms are treated as sedentary. See [References](#references).

The R package lives in the [`dbpmr/`](dbpmr) directory.

## Installation

Installing from source compiles the bundled C engine, so a C toolchain is
required (see [System requirements](#system-requirements) below).

From a local clone, install the package sub-directory:

```r
# install.packages("remotes")
remotes::install_local("dbpmr")
```

or install straight from GitHub:

```r
remotes::install_github(
  "Benthic-Pelagic-Size-Spectrum-Model/spatial-dbpm",
  subdir = "dbpmr"
)
```

or, with the working directory at the repository root, during development:

```r
devtools::install("dbpmr")
```

### System requirements

`dbpmr` contains C code (the simulation engine in `src/`) that is compiled when
the package is installed from source, so you need a **C toolchain**. This is
only needed to *compile* the package — not to *use* an already-installed build.

| OS | Install | How |
|---|---|---|
| Windows | Rtools (matching your R version) | [CRAN Rtools](https://cran.r-project.org/bin/windows/Rtools/) |
| macOS | Xcode command-line tools | `xcode-select --install` |
| Linux (Debian/Ubuntu) | `r-base-dev` | `sudo apt-get install r-base-dev` |

If installation fails with an error such as `compilation failed`,
`make: command not found`, or a missing compiler, the toolchain is not yet
installed. Most users who have previously installed any R package containing
C/C++/Fortran code already have it.

## A simple example

This runs a minimal, aspatial (0-D) coupled benthic–pelagic model with one
pelagic and one benthic species, then reads the results back and plots the final
pelagic size spectrum. Output files are written under a new `MyRun/` directory
in the current working directory.

```r
library(dbpmr)

# Work in a clean directory (the model writes output files here)
setwd(tempdir())

# 1. Configure the run and the mass/time grid
run  <- Setup.Run("MyRun", no_pelagic = 1, no_benthic = 1,
                  spatial_dim = 0, coupled_flag = TRUE, diff_method = 1)
grid <- Setup.Grid(run, tmax = 10)          # integrate for 10 years

# 2. Configure the plankton resource, the two species and the detritus pool.
#    By default reproduction is computed from assimilated energy
#    (rep_method = 2), so no extra input files are needed. See ?Setup.Pelagic
#    for the other reproduction methods.
plankton <- Setup.Plankton(run, filename = "plankton")
pelagic  <- Setup.Pelagic(run,  filename = "fish")
benthic  <- Setup.Benthic(run,  filename = "benthos")
detritus <- Setup.Detritus(run, filename = "detritus")

# 3. Run the simulation (calls the C engine; writes output under MyRun/)
files <- SizeSpectrum(run, grid, plankton, pelagic, benthic, detritus)

# 4. Read the pelagic results back and plot the final size spectrum
fish <- Read.In("MyRun", "fish")
snap <- Extract.Time(fish, time = max(fish@trange))
Plot.Spectrum(snap, type = "l",
              xlab = "log body mass", ylab = "log abundance",
              main = "Pelagic size spectrum")
```

To explore the dynamics over time, `Animate(fish)` steps through the output
times, and `Average.Time(fish)` returns the time-averaged spectrum.

See the function help (`?Setup.Run`, `?Setup.Pelagic`, `?SizeSpectrum`,
`?Read.In`, `?Plot.Spectrum`) for the full set of parameters.

## Customising the grid

`Setup.Grid()` has sensible defaults; for a simple run you usually only set how
long to run and how often to save output:

```r
grid <- Setup.Grid(run,
  mstep    = 0.2,     # mass: computation step (log body mass)
  moutstep = 1,       # mass: output step (must be a multiple of mstep)
  tmax     = 20,      # time: run length (years)
  tstep    = 1/365,   # time: integration step (daily)
  toutstep = 73/365   # time: output step (must be a multiple of tstep)
)
```

The mass range (`mmin`, `mmax`) is widened automatically to cover the species
you define, so you rarely set it. For a `spatial_dim = 0` run the spatial
arguments are ignored (x and y collapse to a single cell).

**The rule that trips people up:** the *output* steps must be exact whole-number
multiples of the *computation* steps, or `SizeSpectrum()` errors. With the
default daily `tstep = 1/365`:

| Desired output | `toutstep` | whole steps? |
|---|---|---|
| every ~73 days (default) | `73/365` | 73 ✓ |
| weekly | `7/365` | 7 ✓ |
| ~monthly | `30/365` | 30 ✓ |
| annual | `1` | 365 ✓ |
| calendar-monthly | `1/12` | 30.4 ✗ (errors) |

The same applies to `moutstep`/`mstep`, and to `xoutstep`/`xstep` and
`youtstep`/`ystep` for spatial (`spatial_dim = 1` or `2`) runs.

## Spatial runs and external data

Spatial runs (`spatial_dim = 1` or `2`) use a **regular rectangular x/y grid**
and plain-text inputs. The package does **not** read shapefiles, rasters or
netCDF directly, and there is no land mask or environmental forcing of vital
rates. You can still drive the spatial inputs (`Setup.ts()` for initial
conditions / time series, `Setup.fishing()` for fishing) from such data by
sampling it onto the model grid with R's geospatial packages (`terra`, `sf`,
`ncdf4`). For example, a raster of fishing mortality:

```r
library(terra)
Fmap <- rast("fishing_mortality.tif")               # in the model's x/y units
Ffun <- function(m, t, x, y) terra::extract(Fmap, cbind(x, y))[1, 1]
Setup.fishing(pelagic, run, grid, func = Ffun)
```

The model grid cells are `seq(xmin, xmax, xstep)` × `seq(ymin, ymax, ystep)`,
and input rows are ordered time-outermost, then x, then y. The vignette
(`vignette("dbpmr")`, section *Driving the model with spatial data*) has the
full worked example, including netCDF time series and shapefiles.

## Development

```r
devtools::document("dbpmr")   # regenerate man/ and NAMESPACE from roxygen2
devtools::test("dbpmr")       # run the testthat suite
devtools::check("dbpmr")      # R CMD check
```

## References

`dbpmr` implements and extends the dynamic, coupled benthic–pelagic
size-spectrum model described in:

- Blanchard, J.L., Jennings, S., Law, R., Castle, M.D., McCloghrie, P.,
  Rochet, M.-J. & Benoît, E. (2009). How does abundance scale with body size in
  coupled size-structured food webs? *Journal of Animal Ecology*, **78**(1),
  270–280. <https://doi.org/10.1111/j.1365-2656.2008.01466.x>
- Blanchard, J.L., Law, R., Castle, M.D. & Jennings, S. (2011). Coupled energy
  pathways and the resilience of size-structured food webs. *Theoretical
  Ecology*, **4**(3), 289–300. <https://doi.org/10.1007/s12080-010-0078-9>
- Castle, M.D., Blanchard, J.L. & Jennings, S. (2011). Predicted effects of
  behavioural movement and passive transport on individual growth and community
  size structure in marine ecosystems. *Advances in Ecological Research*,
  **45**, 41–66. <https://doi.org/10.1016/B978-0-12-386475-8.00002-2>
  (spatial movement of the pelagic size spectrum)

Run `citation("dbpmr")` for these references in BibTeX form.

## Authors

Julia Blanchard and Matthew Castle.

## License

GPL-3.
