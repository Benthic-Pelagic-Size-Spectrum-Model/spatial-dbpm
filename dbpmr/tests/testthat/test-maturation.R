# Maturation-size validation (Setup.Pelagic / Setup.Benthic): mmat >= mmax leaves
# no mature size classes, so endogenous reproduction (rep_method 2/3) silently
# produces zero eggs and the species collapses with no error. These setups should
# warn (and an inverted mmin/mmax should stop).

# TRUE if `expr` emits a maturation warning (ignoring unrelated warnings such as
# dir.create's "already exists").
mat_warned <- function(expr) {
  w <- character()
  withCallingHandlers(expr,
    warning = function(x) { w[[length(w) + 1]] <<- conditionMessage(x)
                            invokeRestart("muffleWarning") })
  any(grepl("mature size classes", w))
}

test_that("Setup.Pelagic warns only when mmat >= mmax", {
  LN10 <- log(10)
  old <- setwd(tempdir()); on.exit({ unlink("mrun", recursive = TRUE); setwd(old) })
  run <- Setup.Run("mrun", 1, 0, 0, TRUE, 1)
  # trap: default mmat (100 g) equals a 100 g mmax
  expect_true(mat_warned(Setup.Pelagic(run, filename = "f1", mmin = -3 * LN10, mmax = 2 * LN10)))
  # fixed: maturation below max
  expect_false(mat_warned(Setup.Pelagic(run, filename = "f2", mmin = -3 * LN10,
                                        mmat = 1 * LN10, mmax = 2 * LN10)))
  # inverted range is an error
  expect_error(Setup.Pelagic(run, filename = "f3", mmin = 2 * LN10, mmax = -3 * LN10),
               "mmin must be less than mmax")
})

test_that("Setup.Benthic warns only when mmat >= mmax", {
  LN10 <- log(10)
  old <- setwd(tempdir()); on.exit({ unlink("mrun", recursive = TRUE); setwd(old) })
  run <- Setup.Run("mrun", 1, 1, 0, TRUE, 1)
  expect_true(mat_warned(Setup.Benthic(run, filename = "b1", mmat = 4 * LN10, mmax = 4 * LN10)))
  expect_false(mat_warned(Setup.Benthic(run, filename = "b2", mmin = -3 * LN10,
                                        mmat = 0, mmax = 4 * LN10)))
})
