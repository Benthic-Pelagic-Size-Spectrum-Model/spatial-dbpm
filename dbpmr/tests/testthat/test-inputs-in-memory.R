# The remaining text inputs held in memory instead (issue #58, step 2a): starting values of every group,
# reproduction, the refuge and the per-cell parameters. Each must give the same run as its text file, and a
# run that uses them needs no files under <run>/Input/.

# One run of nx cells in a row (a single cell when nx = 1), movement off. `memory` selects the input route for
# the same inputs. Starting abundance of the fish and benthos, and the detritus, differ between cells by `scale`.
#   pref:   per-cell pref_ben (vector) or NULL
#   eggs:   function(t, x) giving the egg density of the fish (rep_method 1) or NULL
#   refuge: per-size vector, or [x, 1, size] array (memory only), or NULL; min_A its floor
run_inputs <- function(memory, nx = 3, scale = 1 + 0.5 * (seq_len(nx) - 1), pref = NULL, eggs = NULL, refuge = NULL, min_A = 0,
                       tmax = 2, assign = TRUE, also_files = FALSE) {
  on.exit(Sys.unsetenv("REFUGE_MIN_A"), add = TRUE)
  wd <- tempfile("mem"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = if (nx > 1) 1 else 0, coupled_flag = TRUE, diff_method = 1)
  grid <- if (nx > 1) Setup.Grid(run, tmax = tmax, tstep = 1/48, toutstep = 0.5, xmin = 0, xmax = nx - 1, xstep = 1, xoutstep = 1)
          else Setup.Grid(run, tmax = tmax, tstep = 1/48, toutstep = 0.5)
  pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE)
  pe <- Setup.Pelagic(run, filename = "fish", comp = 0, prey = 0, pred = 0, initial_flag = TRUE, rep_method = if (is.null(eggs)) 2 else 1)
  be <- Setup.Benthic(run, filename = "benthos", initial_flag = TRUE)
  de <- Setup.Detritus(run, filename = "detritus", initial_flag = TRUE)
  sizes <- function(sp) dbpmr:::grid_seq(sp@mmin, sp@mmax, grid@mstep)
  spec  <- function(sp) sp@u_0 * exp(sp@lambda * sizes(sp))                    # default spectrum over the species' sizes
  cells <- function(sp, f) outer(f, spec(sp))                                  # [cell, size]
  files <- function() {
    Setup.ts(pl, run, grid, mat = cells(pl, rep(1, nx))); Setup.ts(pe, run, grid, mat = cells(pe, scale))
    Setup.ts(be, run, grid, mat = cells(be, scale)); Setup.ts(de, run, grid, func = function(m, t, x, y) 0.6 * scale[round(x) + 1])
    if (!is.null(pref)) Setup.CellParams(pe, run, grid, pref_ben = pref)
    if (!is.null(eggs)) Setup.Rep(pe, run, grid, func = function(t, x, y) eggs(t, x))
    if (!is.null(refuge)) Setup.Refuge(pe, run, grid, refuge, min_A = min_A)
  }
  if (!memory || also_files) files()
  if (memory) {
    p2 <- Setup.ts(pl, run, grid, initial = spec(pl))                          # one spectrum: the same in every cell
    f2 <- Setup.ts(pe, run, grid, initial = if (nx > 1) cells(pe, scale) else spec(pe) * scale)
    b2 <- Setup.ts(be, run, grid, initial = if (nx > 1) cells(be, scale) else spec(be) * scale)
    d2 <- Setup.ts(de, run, grid, initial = 0.6 * scale)
    if (!is.null(pref)) f2 <- Setup.CellParams(f2, run, grid, pref_ben = pref, file = FALSE)
    if (!is.null(eggs)) f2 <- Setup.Rep(f2, run, grid, eggs = array(outer(seq_len(nx) - 1, 0:(tmax - 1), function(x, t) eggs(t, x)), c(nx, 1, tmax)),
                                        time = 0:(tmax - 1))
    if (!is.null(refuge)) f2 <- Setup.Refuge(f2, run, grid, refuge, min_A = min_A, file = FALSE)
    if (assign) { pl <- p2; pe <- f2; be <- b2; de <- d2 }
  }
  out <- NULL
  invisible(capture.output(out <- SizeSpectrum(run, grid, pl, pe, be, de, output = "arrays")))
  attr(out, "input_files") <- list.files(file.path("R", "Input"))
  out
}
same_run <- function(a, b, tol = 1e-10) {
  for (v in c("abundance", "growth", "mortality")) expect_equal(a[[v]], b[[v]], tolerance = tol, info = v)
  expect_equal(a$detritus, b$detritus, tolerance = tol)
}
band <- function(lo, hi, strength = 1e20) {                                    # refuge on a band of sizes of the model grid
  m <- dbpmr:::grid_seq(-12 * log(10), 6 * log(10), 0.1 * log(10)) / log(10)
  ifelse(m >= lo & m <= hi, strength, 0)
}

test_that("starting values held in memory give the same run as the _ts files, and write no files", {
  for (nx in c(1, 3)) {
    f <- run_inputs(FALSE, nx); m <- run_inputs(TRUE, nx)
    expect_gt(length(attr(f, "input_files")), 0); expect_length(attr(m, "input_files"), 0)
    same_run(m, f)
    u0 <- matrix(m$abundance$fish[, 1, , 1], nx); has <- u0[1, ] > 0; sc <- 1 + 0.5 * (seq_len(nx) - 1)
    expect_equal(u0[, has, drop = FALSE], outer(sc, u0[1, has]), tolerance = 1e-12)   # the starting state is what was supplied
    expect_equal(as.vector(m$detritus[, 1, 1]), 0.6 * sc)
  }
})

test_that("per-cell pref_ben held in memory gives the same run as the _cellpar file", {
  pref <- c(0.2, 0.4, 0.8)
  f <- run_inputs(FALSE, pref = pref); m <- run_inputs(TRUE, pref = pref)
  same_run(m, f)
  expect_false(isTRUE(all.equal(m$abundance$benthos, run_inputs(TRUE)$abundance$benthos)))   # and it does something
  expect_length(attr(m, "input_files"), 0)
})

test_that("reproduction held in memory gives the same run as the _rep_ts file", {
  eggs <- function(t, x) 10 * (1 + 0.5 * x) * (1 + 0.25 * floor(t + 1e-9))     # changes between cells and between years
  f <- run_inputs(FALSE, eggs = eggs); m <- run_inputs(TRUE, eggs = eggs)
  same_run(m, f)
  first <- which(m$abundance$fish[1, 1, , 2] > 0)[1]                           # the smallest fish are the eggs supplied
  expect_equal(m$abundance$fish[, 1, first, which(abs(m$t - 1.5) < 1e-9)], eggs(1.5, 0:2), tolerance = 1e-12)
  expect_length(attr(m, "input_files"), 0)
})

test_that("a refuge held in memory gives the same run as the _refuge file", {
  rf <- band(0, 1.5)
  f <- run_inputs(FALSE, refuge = rf, min_A = 0.2); m <- run_inputs(TRUE, refuge = rf, min_A = 0.2)
  same_run(m, f)
  expect_false(isTRUE(all.equal(m$abundance$fish, run_inputs(TRUE)$abundance$fish)))
  expect_length(attr(m, "input_files"), 0)
})

test_that("a refuge that differs between cells gives each cell the result of its own single-cell run", {
  rfs <- list(band(0, 1.5), band(-1, 0.5), band(0, 0))                         # a different refuge in each cell (none in the third)
  nm <- length(rfs[[1]])
  grid_run <- run_inputs(TRUE, scale = rep(1, 3), refuge = aperm(array(unlist(rfs), c(nm, 3, 1)), c(2, 3, 1)), min_A = 0.2)
  for (k in 1:3) {
    one <- run_inputs(TRUE, nx = 1, refuge = rfs[[k]], min_A = 0.2)
    for (v in c("fish", "benthos")) expect_equal(grid_run$abundance[[v]][k, 1, , ], one$abundance[[v]][1, 1, , ], tolerance = 1e-10, info = paste(v, k))
  }
  expect_false(isTRUE(all.equal(grid_run$abundance$fish[1, 1, , ], grid_run$abundance$fish[2, 1, , ])))
})

test_that("forgetting to assign, or supplying both memory and file, is an error", {
  expect_error(run_inputs(TRUE, assign = FALSE), "assign the result")
  expect_error(run_inputs(TRUE, also_files = TRUE), "use one or the other")
})

test_that("the in-memory inputs are validated", {
  wd <- tempfile("v"); dir.create(wd); old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  run  <- Setup.Run("R", 1, 1, spatial_dim = 2, coupled_flag = TRUE, diff_method = 1)
  grid <- Setup.Grid(run, tmax = 1, tstep = 1/48, toutstep = 1, xmin = 0, xmax = 2, xstep = 1, xoutstep = 1, ymin = 0, ymax = 1, ystep = 1, youtstep = 1)
  pe <- Setup.Pelagic(run, filename = "fish", initial_flag = TRUE, rep_method = 1); ns <- length(dbpmr:::grid_seq(pe@mmin, pe@mmax, grid@mstep))
  de <- Setup.Detritus(run, filename = "detritus", initial_flag = TRUE)
  out <- Setup.ts(pe, run, grid, initial = rep(1, ns))
  expect_s4_class(out, "pelagic.params"); expect_equal(dim(out@arrays$initial), c(3L, 2L, ns))
  expect_equal(dim(Setup.ts(pe, run, grid, initial = array(1, c(3, 2, ns)))@arrays$initial), c(3L, 2L, ns))
  expect_equal(dim(Setup.ts(de, run, grid, initial = 0.5)@arrays$initial), c(3L, 2L))
  expect_error(Setup.ts(pe, run, grid, initial = rep(1, 5)), "one value per size class")
  expect_error(Setup.ts(pe, run, grid, initial = array(1, c(2, 3, ns))), "x = 3, y = 2")
  expect_error(Setup.ts(pe, run, grid, initial = rep(-1, ns)), "non-negative")
  expect_error(Setup.ts(pe, run, grid, initial = rep(NA_real_, ns)), "NA, NaN or Inf")
  expect_error(Setup.ts(pe, run, grid, initial = rep(1, ns), mat = matrix(1, 1, ns)), "cannot be combined")
  expect_error(Setup.ts(Setup.Pelagic(run, filename = "fish"), run, grid, initial = rep(1, ns)), "initial_flag = TRUE")
  expect_error(Setup.ts(de, run, grid, initial = matrix(1, 2, 3)), "x = 3, y = 2")

  expect_equal(dim(Setup.Rep(pe, run, grid, eggs = c(1, 2), time = c(0, 0.5))@arrays$reproduction$values), c(3L, 2L, 2L))
  expect_error(Setup.Rep(pe, run, grid, eggs = -1), "non-negative")
  expect_error(Setup.Rep(pe, run, grid, eggs = 1, func = function(t, x, y) 1), "cannot be combined")
  expect_error(Setup.Rep(pe, run, grid, func = function(t, x, y) 1, time = 0), "only used with eggs")

  nm <- length(dbpmr:::grid_seq(grid@mmin, grid@mmax, grid@mstep))
  expect_equal(dim(Setup.Refuge(pe, run, grid, rep(0, nm), file = FALSE)@arrays$refuge$values), c(3L, 2L, nm))
  expect_equal(Setup.Refuge(pe, run, grid, array(0, c(3, 2, nm)), min_A = 0.3, file = FALSE)@arrays$refuge$min_A, 0.3)
  expect_error(Setup.Refuge(pe, run, grid, array(0, c(3, 2, nm))), "use file = FALSE")
  expect_error(Setup.Refuge(pe, run, grid, array(0, c(2, 3, nm)), file = FALSE), "x = 3, y = 2")

  expect_equal(Setup.CellParams(pe, run, grid, pref_ben = matrix(0.5, 3, 2), file = FALSE)@arrays$pref_ben, matrix(0.5, 3, 2))
  expect_false(file.exists(file.path("R", "Input", "fish_cellpar.txt")))
})
