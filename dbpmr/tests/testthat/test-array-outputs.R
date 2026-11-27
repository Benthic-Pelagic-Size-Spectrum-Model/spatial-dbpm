# Model output returned to R as arrays (issue #58, step 1) instead of, or as well as, the text result
# files, and written to NetCDF. The arrays must hold exactly what the text files hold.

# a small 2-D run with movement and fishing; returns whatever SizeSpectrum() returns and leaves the run directory in place
run_out <- function(
    output,
    xoutstep = 1,
    moutstep = 0.1 * log(10),
    toutstep = 0.5,
    move = TRUE,
    ...
) {
    run <- Setup.Run(
        "R",
        1,
        1,
        spatial_dim = 2,
        coupled_flag = TRUE,
        diff_method = 1
    )
    grid <- Setup.Grid(
        run,
        tmax = 2,
        tstep = 1 / 48,
        toutstep = toutstep,
        moutstep = moutstep,
        xmin = 0,
        xmax = 4,
        xstep = 1,
        xoutstep = xoutstep,
        ymin = 0,
        ymax = 1,
        ystep = 1,
        youtstep = 1
    )
    uni <- function(sp) {
        m <- dbpmr:::grid_seq(sp@mmin, sp@mmax, grid@mstep)
        matrix(sp@u_0 * exp(sp@lambda * m), nrow = 1)
    }
    pl <- Setup.Plankton(run, filename = "plankton", initial_flag = TRUE)
    Setup.ts(pl, run, grid, func = function(m, t, x, y) {
        0.01 * exp(-m) * (1 + 0.2 * x + 0.1 * y)
    })
    pe <- Setup.Pelagic(
        run,
        filename = "pelagic",
        prey = if (move) 0.3 else 0,
        comp = if (move) 0.2 else 0,
        initial_flag = TRUE,
        fishing_flag = TRUE
    )
    Setup.ts(pe, run, grid, mat = uni(pe))
    Setup.fishing(pe, run, grid, func = function(m, t, x, y) {
        if (m > log(10)) 0.2 * (1 + 0.1 * x) else 0
    })
    be <- Setup.Benthic(run, filename = "benthic", initial_flag = TRUE)
    Setup.ts(be, run, grid, mat = uni(be))
    de <- Setup.Detritus(run, filename = "detritus")
    res <- NULL
    if (output == "files") {
        invisible(capture.output(
            res <- SizeSpectrum(
                run,
                grid,
                pl,
                pe,
                be,
                de,
                write.files = TRUE,
                memory.return = FALSE,
                ...
            )
        ))
    }
    if (output == "arrays") {
        invisible(capture.output(
            res <- SizeSpectrum(
                run,
                grid,
                pl,
                pe,
                be,
                de,
                write.files = FALSE,
                memory.return = TRUE,
                ...
            )
        ))
    }
    if (output == "both") {
        invisible(capture.output(
            res <- SizeSpectrum(
                run,
                grid,
                pl,
                pe,
                be,
                de,
                write.files = TRUE,
                memory.return = TRUE,
                ...
            )
        ))
    }
    if (output == "none") {
        invisible(capture.output(
            res <- SizeSpectrum(
                run,
                grid,
                pl,
                pe,
                be,
                de,
                write.files = FALSE,
                memory.return = FALSE,
                ...
            )
        ))
    }

    res
}
num <- function(z) {
    suppressWarnings(as.numeric(strsplit(z, ",", fixed = TRUE)[[1]]))
}
# a per-timestep text result file as an [x, y, m, t] array (rows are x-outer / y-inner within each output
# time; the final time is written a second time at the end of results.txt, so only the first block is used)
# read_txt <- function(group, file, nx, ny) {
#   L <- readLines(file.path("R", group, paste0(file, ".txt")))
#   m <- num(L[5]); m <- m[!is.na(m)]
#   rows <- lapply(L[-(1:5)], num); rows <- rows[vapply(rows, function(v) length(v) == length(m) + 3 && !anyNA(v[1:3]), TRUE)]
#   tt <- vapply(rows, `[`, 0, 1); ut <- unique(tt)
#   a <- array(NA_real_, c(nx, ny, length(m), length(ut)))
#   for (j in seq_along(ut)) { blk <- rows[tt == ut[j]][seq_len(nx * ny)]
#     for (r in seq_along(blk)) a[((r - 1) %/% ny) + 1, ((r - 1) %% ny) + 1, , j] <- blk[[r]][-(1:3)] }
#   list(a = a, m = m, t = ut)
# }
read_txt <- function(group, file, nx, ny, nt) {
    df <- read.csv(
        file.path("R", group, paste0(file, ".txt")),
        skip = 4,
        nrows = nx * ny * nt,
        check.names = FALSE
    )
    mat <- as.matrix(df[, -(1:3)]) # Size class values without dimension values

    ti <- match(df$t, sort(unique(df$t))) # Indices for dimensions to allocate to array
    xi <- match(df$x, sort(unique(df$x)))
    yi <- match(df$y, sort(unique(df$y)))

    res <- array(NA, dim = c(ncol(df[, -(1:3)]), nx, ny, nt))
    for (r in seq_len(nrow(df))) {
        res[, xi[r], yi[r], ti[r]] <- mat[r, ]
    }

    dimnames(res) <- list(
        m = dimnames(mat)[[2]],
        x = sort(unique(df$x)),
        y = sort(unique(df$y)),
        time = sort(unique(df$t))
    )
    return(res)
}

same <- function(a, b) {
    identical(dim(a), dim(b)) && all(abs(a - b) <= 1e-15 + 1e-12 * abs(b))
} # the files hold 16 decimal places

test_that("error if both output flags are unselected", {
    wd <- tempfile("out")
    dir.create(wd)
    old <- setwd(wd)
    on.exit(setwd(old), add = TRUE)
    expect_error(run_out("none"))
})

test_that("output arrays hold exactly what the text result files hold", {
    wd <- tempfile("out")
    dir.create(wd)
    old <- setwd(wd)
    on.exit(setwd(old), add = TRUE)
    out <- run_out("both")
    # expect_s3_class(out, "dbpmr_output")
    # expect_named(out, c("t", "x", "y", "m", "abundance", "growth", "mortality", "predation", "fishing", "catch", "detritus"))
    # expect_named(out$abundance, c("plankton", "fish", "benthos")); expect_named(out$growth, c("fish", "benthos"))
    # expect_equal(dim(out$abundance$plankton), c(5L, 2L, 181L, 5L))
    # expect_equal(out$t, seq(0, 2, 0.5)); expect_equal(out$x, 0:4); expect_equal(out$y, 0:1)
    # expect_type(attr(out, "files"), "character")
    # Read in pelagic results (u) values from file
    f <- read_txt("pelagic", "results", 5, 2, 5)

    # Check dimensions are consistent with memory outputs
    expect_equal(out@m, as.numeric(dimnames(f)$m), tolerance = 1e-6)
    expect_equal(out@time, as.numeric(dimnames(f)$t))
    expect_equal(out@x, as.numeric(dimnames(f)$x))
    expect_equal(out@y, as.numeric(dimnames(f)$y))

    # Check the data values that are underlying each output (m values significant figures do not match, so just removing dimensions)
    for (group in c("plankton", "pelagic", "benthic")) {
        memory_output <- slot(out, paste0(group, "_u"))
        if (group %in% c("pelagic", "benthic")) {
            memory_output <- memory_output[,,, "1", ] # Remove the species dimension from the pelagic array
        }
        file_output <- read_txt(group, "results", 5, 2, 5)
        expect_equal(as.vector(memory_output), as.vector(file_output))
    }
    for (group in c("pelagic", "benthic")) {
        for (vari in c("mu_pred", "mu_fish", "g")) {
            memory_output <- slot(out, paste0(group, "_", vari))
            if (vari == "mu_pred") {
                if (group == "pelagic") {
                    fn <- "predation0"
                    memory_output <- memory_output[,,, "1", , ] # mu_pred has an extra dimension for prey and predator species (!Only for the pelagic groups)
                } else {
                    fn <- "predation"
                    memory_output <- memory_output[,,, "1", ]
                }
            }
            if (vari == "mu_fish") {
                fn <- "fishing"
            }
            if (vari == "g") {
                fn <- "growth"
            }
            if (vari != "mu_pred") {
                memory_output <- memory_output[,,, "1", ]
            }

            file_output <- read_txt(group, fn, 5, 2, 5)
            expect_equal(as.vector(memory_output), as.vector(file_output))
        }
    }

    expect_gt(max(out@pelagic_catch), 0)
    expect_equal(max(out@benthic_catch), 0)
})

test_that("Test the NetCDF writing capability", {
    skip_if_not_installed("ncdf4")
    wd <- tempfile("out")
    dir.create(wd)
    old <- setwd(wd)
    on.exit(setwd(old), add = TRUE)
    out <- run_out("arrays")
    f <- Write.NetCDF(
        out,
        "run.nc",
        x = 70 + out@x,
        y = -55 + out@y,
        xname = "lon",
        yname = "lat",
        xunits = "degrees_east",
        yunits = "degrees_north"
    )
    nc <- ncdf4::nc_open(f)
    on.exit(ncdf4::nc_close(nc), add = TRUE)
    expect_setequal(
        names(nc$var),
        names(getSlots(class(out)))[
            !(names(getSlots(class(out))) %in%
                c("m", "x", "y", "time", "pelagic_species", "benthic_species"))
        ]
    )
    expect_equal(
        vapply(nc$var$pelagic_u$dim, `[[`, "", "name"),
        c("size", "lon", "lat", "pelagic_species", "time")
    )
    expect_equal(as.vector(nc$dim$lon$vals), 70 + out@x)
    expect_equal(as.vector(nc$dim$time$vals), out@time)
    expect_equal(as.vector(nc$dim$size$vals), out@m)
    get <- function(v) ncdf4::ncvar_get(nc, v, collapse_degen = FALSE)
    expect_equal(get("pelagic_u"), out@pelagic_u, ignore_attr = TRUE)
    expect_equal(get("plankton_u"), out@plankton_u, ignore_attr = TRUE)
    expect_equal(get("benthic_g"), out@benthic_g, ignore_attr = TRUE)
    expect_equal(
        get("pelagic_mu_pred"),
        out@pelagic_mu_pred,
        ignore_attr = TRUE
    )
    expect_equal(
        get("benthic_mu_pred"),
        out@benthic_mu_pred,
        ignore_attr = TRUE
    )
    expect_equal(
        get("pelagic_mu_fish"),
        out@pelagic_mu_fish,
        ignore_attr = TRUE
    )
    expect_equal(
        get("pelagic_catch"),
        out@pelagic_catch,
        ignore_attr = TRUE
    )
    expect_equal(get("detritus_w"), out@detritus_w, ignore_attr = TRUE)
})
