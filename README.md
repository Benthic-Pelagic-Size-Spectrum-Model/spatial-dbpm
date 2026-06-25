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

A C toolchain is required to compile the simulation engine (Rtools on Windows,
Xcode command-line tools on macOS, or `r-base-dev` on Linux).

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
#    rep_method = 0 uses a fixed reproduction amount (no input file needed);
#    the default rep_method = 1 reads a reproduction time series (see ?Setup.Rep).
plankton <- Setup.Plankton(run, filename = "plankton")
pelagic  <- Setup.Pelagic(run,  filename = "fish",    rep_method = 0)
benthic  <- Setup.Benthic(run,  filename = "benthos", rep_method = 0)
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
