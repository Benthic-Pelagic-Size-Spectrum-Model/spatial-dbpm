#' Set up species from a table of parameters
#'
#' Builds the pelagic and benthic species of a run from a data frame with one row
#' per species, in the style of the `species_params` table used by the
#' \pkg{mizer} package. It is an alternative to calling [Setup.Pelagic()] and
#' [Setup.Benthic()] once per species: each row is passed to the matching
#' function, so the same defaults, checks and warnings apply.
#'
#' @details
#' The table needs two columns:
#'
#' * `species`: the species name, used as its file name (the `filename`
#'   argument of [Setup.Pelagic()] / [Setup.Benthic()]);
#' * `type`: `"pelagic"` or `"benthic"`.
#'
#' Any other column must be named after an argument of [Setup.Pelagic()] or
#' [Setup.Benthic()] (for example `mmin`, `mmat`, `mmax`, `A`, `mu_0`,
#' `pref_ben`, `comp`, `fishing_flag`). A column that is absent, or a cell that
#' is `NA`, leaves that parameter at the function default. A value in a column
#' that does not apply to the row's type (for example `pref_ben` for a benthic
#' species) is an error. Body masses are natural-log grams, as in the functions.
#'
#' The number of pelagic and benthic rows must match the run
#' ([Setup.Run()]'s `no_pelagic` and `no_benthic`).
#'
#' @param run.in A [run.params] object.
#' @param species_params A data frame with one row per species, or the path of
#'   a CSV file holding one.
#'
#' @return A list with elements `pelagic` and `benthic`, each a named list of
#'   species parameter objects, ready to pass to [SizeSpectrum()] as
#'   `pelagic.in` and `benthic.in`.
#' @seealso [Setup.Pelagic()], [Setup.Benthic()], [SizeSpectrum()]
#' @examples
#' \dontrun{
#' run <- Setup.Run("R", 1, 1, 0, TRUE, 1)
#' sp  <- data.frame(species = c("fish", "benthos"),
#'                   type    = c("pelagic", "benthic"),
#'                   mmat    = c(2, 1) * log(10),
#'                   mmax    = c(6, 4) * log(10),
#'                   mu_0    = c(0.1, 0.1),
#'                   pref_ben = c(0.5, NA))
#' groups <- Setup.Species(run, sp)
#' # SizeSpectrum(run, grid, pl, groups$pelagic, groups$benthic, de)
#' }
#' @export
`Setup.Species` <-
function(run.in, species_params){
  if(is.character(species_params) && length(species_params) == 1){
    if(!file.exists(species_params)) stop("species_params file not found: ", species_params)
    species_params <- utils::read.csv(species_params, stringsAsFactors = FALSE)
  }
  if(!is.data.frame(species_params)) stop("species_params must be a data frame or the path of a CSV file")
  sp <- species_params
  if(!all(c("species", "type") %in% names(sp))) stop("species_params needs the columns 'species' and 'type'")
  sp$species <- as.character(sp$species); sp$type <- as.character(sp$type)
  if(nrow(sp) == 0) stop("species_params has no rows")
  if(anyNA(sp$species) || any(sp$species == "") || anyDuplicated(sp$species)) stop("species names must be given and unique")
  if(!all(sp$type %in% c("pelagic", "benthic"))) stop("type must be \"pelagic\" or \"benthic\"")

  fun  <- list(pelagic = Setup.Pelagic, benthic = Setup.Benthic)
  args <- lapply(fun, function(f) setdiff(names(formals(f)), c("run.in", "filename")))
  cols <- setdiff(names(sp), c("species", "type"))
  unknown <- setdiff(cols, unlist(args))
  if(length(unknown)) stop("unknown column in species_params: ", paste(unknown, collapse = ", "),
                           " (columns must be arguments of Setup.Pelagic() or Setup.Benthic())")
  flags <- c("initial_flag", "ts_flag", "fishing_flag")

  n_pel <- sum(sp$type == "pelagic"); n_ben <- sum(sp$type == "benthic")
  if(n_pel != run.in@no_pelagic || n_ben != run.in@no_benthic)
    stop(sprintf("species_params has %d pelagic and %d benthic species but the run was set up for %d and %d",
                 n_pel, n_ben, run.in@no_pelagic, run.in@no_benthic))

  out <- list(pelagic = list(), benthic = list())
  for(r in seq_len(nrow(sp))){
    type <- sp$type[r]
    call_args <- list(run.in = run.in, filename = sp$species[r])
    for(col in cols){
      v <- sp[[col]][r]
      if(is.na(v)) next
      if(!col %in% args[[type]])
        stop(sprintf("column '%s' does not apply to a %s species (row %d, %s)", col, type, r, sp$species[r]))
      call_args[[col]] <- if(col %in% flags) as.logical(v) else v
    }
    out[[type]][[sp$species[r]]] <- do.call(fun[[type]], call_args)
  }
  out
}
