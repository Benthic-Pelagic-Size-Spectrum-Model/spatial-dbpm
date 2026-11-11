#include <SizeSpectra.h>

// This is the new R interface. SizeSpectrum() is now a wrapper function that
// takes the objects from R (.Call()). SizeSpectrum() checks the SEXP objects from
// R are the correct SXPtypes, sets C pointers to the R data, and then runs the
// internal function SizeSpectrum_internals() which performs the calculations and
// writes outputs. SizeSpectrum() returns NULL to R.
// TODO: return to R will be changed to be optional so that users can toggle
// file output, memory output or both.
SEXP SizeSpectrum(SEXP run_params_r, SEXP grid_params_r, SEXP plankton_params_r, SEXP pelagic_params_r,
                  SEXP benthic_params_r, SEXP detritus_params_r, SEXP names_params_r, SEXP flags_params_r,
                  SEXP output_params_r, SEXP land_mask_r)
{
    // Ensure that all of the SEXP objects are of the correct type after
    // passing in from R.
    if (!isInteger(run_params_r) || !isInteger(flags_params_r) || !isInteger(output_params_r) ||
        !isInteger(land_mask_r)
    )
    {
        error("run, flags, output *_params and land_mask R objects must be integer vectors");
    }
    if (!isReal(grid_params_r) || !isReal(plankton_params_r) || !isReal(pelagic_params_r) ||
        !isReal(benthic_params_r) || !isReal(detritus_params_r))
    {
        error("grid, plankton, pelagic, benthic and detritus params must be numeric vectors");
    }
    if (!isString(names_params_r))
    {
        error("names_params_r R object must be character");
    }

    // Checking the length of the parameter vectors that were passed from R.
    // The values here are taken from hardcoded indices throughout the code.
    if (LENGTH(run_params_r) != 5)
    {
        error("run_params_r must have 5 elements");
    }
    if (LENGTH(grid_params_r) != 18)
    {
        error("grid_params_r must have 18 elements");
    }
    
    int *land_mask = INTEGER(land_mask_r);
    R_xlen_t n_mask = XLENGTH(land_mask_r);
    for (R_xlen_t i=0; i<n_mask; i++) {
        if (land_mask[i] != 0 && land_mask[i] != 1) {
            error("land_mask must only contain values of 0 or 1.");
        }
    }

    if (LENGTH(plankton_params_r) != 6)
    {
        error("plankton_params_r must have 6 elements");
    }
    if (LENGTH(output_params_r) != 2)
    {
        error("output_params_r must have 2 elements");
    }

    int npel = INTEGER(run_params_r)[0];
    if (LENGTH(pelagic_params_r) != 32 * npel)
    {
        error("pelagic_params_r must have 32 * n_pel_species elements");
    }

    // It is possible to include no benthic species, in which case benthic_params_r
    // and detritus_params_r do not need to be checked for length.
    int nben = INTEGER(run_params_r)[1];
    if (nben > 0)
    {
        if (LENGTH(benthic_params_r) != 15 * nben)
        {
            error("benthic_params_r must have 15 * n_ben_species elements");
        }
        if (LENGTH(detritus_params_r) != 1)
        {
            error("detritus_params_r must have 1 element");
        }
        if (LENGTH(names_params_r) != (npel + nben + 1) + 2)
        {
            error("names_params_r not correct length");
        }
    }
    else if (LENGTH(names_params_r) != npel + 2)
    {
        error("names_params_r not correct length");
    }

    int nflags = 2 + 4 * (npel + nben) + (nben > 0 ? 2 : 0);
    // nflags is 2 for plankton, 4 for each pelagic / benthic species and 2
    // for detritus if there is at least 1 benthic species
    if (LENGTH(flags_params_r) != nflags)
    {
        error("flags_params_r not correct length");
    }

    // Assign pointers for the SEXP data to correct names for function use
    int *run_params = INTEGER(run_params_r);
    int *flags_params = INTEGER(flags_params_r);
    int *output_params = INTEGER(output_params_r);
    double *grid_params = REAL(grid_params_r);
    double *pla_params = REAL(plankton_params_r);
    double *pel_params = REAL(pelagic_params_r);
    double *ben_params = REAL(benthic_params_r);
    double *det_params = REAL(detritus_params_r);

    char **names_params = (char **)R_alloc(LENGTH(names_params_r), sizeof(char *));
    int name_idx;
    for (name_idx = 0; name_idx < LENGTH(names_params_r); name_idx++)
    {
        names_params[name_idx] = (char *)CHAR(STRING_ELT(names_params_r, name_idx));
    }

    SEXP result = SizeSpectrum_internals(run_params, grid_params, pla_params, pel_params, ben_params, det_params,
                                         names_params, flags_params, output_params, land_mask, n_mask);
    return (result);
}

//-----------------------------------------------------------------------------//
// This is the internal function that does the work of the program called by R //
//-----------------------------------------------------------------------------//

static SEXP SizeSpectrum_internals(int *run_params, double *grid_params, double *pla_params, double *pel_params,
                                   double *ben_params, double *det_params, char **names_params, int *flags_params,
                                   int *output_params, int *land_mask, R_xlen_t mask_length)
// run_params is an array of run parameters
// grid_params is an array of grid discretisation values
// pla_params is plankton params
// pel_params is pelagic params
// ben_params is benthic parameters
// det_params is detritus parameters
// names_params is a list of filenames
// flags_params is a list of selection 'flag' parameters
{
    int pel_species_idx, ben_species_idx;
    RUN run;
    GRID grid;
    COMMUNITY community;
    OUTPUT out = {0};
    MATRIX *pelmatrix;
    MATRIX *benmatrix = NULL;
    MATRIX xmatrix, ymatrix;
    double **temp_pel_params;
    double **temp_ben_params;

    /*Setup run*/
    setup_run(&run, run_params, names_params[0], output_params);

    /*Setup grid*/
    setup_grid(&grid, grid_params, land_mask, mask_length);

    SEXP result = R_NilValue;
    if (run.return_memory)
        result = PROTECT(allocVector(VECSXP, 18));

    /*Setup community*/
    community.plankton = (PLANKTON *)safe_malloc(sizeof(PLANKTON), __LINE__, __FILE__);
    community.pelagic = (PELAGIC *)safe_malloc(run.no_pelagic * sizeof(PELAGIC), __LINE__, __FILE__);
    if (run.no_benthic != 0)
    {
        community.benthic = (BENTHIC *)safe_malloc(run.no_benthic * sizeof(BENTHIC), __LINE__, __FILE__);
        community.detritus = (DETRITUS *)safe_malloc(sizeof(DETRITUS), __LINE__, __FILE__);
    }

    /*Setup plankton*/
    setup_plankton(&run, &grid, community.plankton, pla_params, names_params[1], flags_params[0], flags_params[1]);

    /*Setup pelagic*/
    temp_pel_params = setup_pelagic_params(run.no_pelagic, pel_params);
    for (pel_species_idx = 0; pel_species_idx < run.no_pelagic; pel_species_idx++)
    {
        setup_pelagic(&run, &grid, &(community.pelagic[pel_species_idx]), temp_pel_params[pel_species_idx],
                      names_params[pel_species_idx + 2], flags_params[4 * pel_species_idx + 2],
                      flags_params[4 * pel_species_idx + 3], flags_params[4 * pel_species_idx + 4],
                      flags_params[4 * pel_species_idx + 5]);
    }

    if (run.no_benthic != 0)
    {
        /*Setup benthic*/
        temp_ben_params = setup_benthic_params(run.no_benthic, ben_params);
        for (ben_species_idx = 0; ben_species_idx < run.no_benthic; ben_species_idx++)
        {
            setup_benthic(&run, &grid, &(community.benthic[ben_species_idx]), temp_ben_params[ben_species_idx],
                          names_params[ben_species_idx + run.no_pelagic + 2],
                          flags_params[4 * ben_species_idx + (2 + (4 * run.no_pelagic))],
                          flags_params[4 * ben_species_idx + (3 + (4 * run.no_pelagic))],
                          flags_params[4 * ben_species_idx + (4 + (4 * run.no_pelagic))],
                          flags_params[4 * ben_species_idx + (5 + (4 * run.no_pelagic))]);
        }

        /*Setup detritus*/
        setup_detritus(&run, &grid, community.detritus, det_params, names_params[run.no_pelagic + run.no_benthic + 2],
                       flags_params[2 + 4 * (run.no_pelagic + run.no_benthic)],
                       flags_params[3 + 4 * (run.no_pelagic + run.no_benthic)]);
    }

    /*Setup differencing matrices*/
    pelmatrix = (MATRIX *)safe_malloc(run.no_pelagic * sizeof(MATRIX), __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < run.no_pelagic; pel_species_idx++)
    {
        setup_matrix((community.pelagic[pel_species_idx].ipelmax - community.pelagic[pel_species_idx].ipelmin + 1),
                     &(pelmatrix[pel_species_idx]));
    }

    if (run.no_benthic != 0)
    {
        benmatrix = (MATRIX *)safe_malloc(run.no_benthic * sizeof(MATRIX), __LINE__, __FILE__);
        for (ben_species_idx = 0; ben_species_idx < run.no_benthic; ben_species_idx++)
        {
            setup_matrix((community.benthic[ben_species_idx].ibenmax - community.benthic[ben_species_idx].ibenmin + 1),
                         &(benmatrix[ben_species_idx]));
        }
    }

    setup_matrix(grid.xnum, &xmatrix);
    setup_matrix(grid.ynum, &ymatrix);

    setup_output(&run, &grid, &out, result);

    // Set land cells to zero incase they have incoming values from input files
    set_land_to_zero(&run, &grid, &community);

    /*Perform all finite differencing*/
    calculate_results(&run, &grid, &community, &out, pelmatrix, benmatrix, &xmatrix, &ymatrix);

    /*Free all allocated memory*/
    free_mem(&run, &grid, &community, pelmatrix, benmatrix, &xmatrix, &ymatrix, temp_pel_params, temp_ben_params);

    if (run.return_memory)
        UNPROTECT(1);
    return (result);
}

//-----------------------------------//
// Output Function. Not to be edited //
//-----------------------------------//

void calculate_results(RUN *run, GRID *grid, COMMUNITY *community, OUTPUT *out, MATRIX *pelmatrix, MATRIX *benmatrix,
                       MATRIX *xmatrix, MATRIX *ymatrix)
/* Sets up output txt files and implements time stepping process*/
{
    // Indices for each variable with readable names
    int sc_idx, ts_idx, x_idx, y_idx, pel_species_idx, ben_species_idx;
    // Counts for each variable with readable names
    int n_m, n_t, n_x, n_y, n_pel_species, n_ben_species;
    double time1 = 0, time2 = 0, timed = 0; // variables for countdown timer

    n_m = grid->mnum;
    n_t = grid->tnum;
    n_x = grid->xnum;
    n_y = grid->ynum;
    n_pel_species = run->no_pelagic;
    n_ben_species = run->no_benthic;

    /*Optional environmental-forcing plugin: auto-detect <run>/Input/forcing_ts.txt, whose
      first line names the channels present (comma-separated) and each subsequent line gives
      one timestep's values. We map the recognised channel names to column indices here; the
      per-timestep read + apply happens in mass_solver(). Absent file/channels -> neutral.*/
    run->forcing_flag = 0;
    run->forcing_ncol = 0;
    run->col_pel_tempeff = -1;
    run->col_ben_tempeff = -1;
    run->col_sinking_rate = -1;
    run->col_depth = -1;
    run->col_residence_time = -1;
    run->col_export_attn = -1;
    run->col_export_gamma = -1;
    run->col_export_magg = -1;
    run->col_gravity = -1;
    run->pel_tempeff = 1.0;
    run->ben_tempeff = 1.0;
    run->sinking_rate = 1.0;
    run->depth = 1.0;
    run->residence_time = 0.0;
    run->export_attn = 0.0;
    run->export_gamma = 0.0;
    run->export_magg = -1e30;
    run->gravity = 0.0;
    run->fptr_forcing = NULL;
    run->forcing_percell = 0;
    run->pel_tempeff_c = NULL;
    run->ben_tempeff_c = NULL;
    run->sinking_rate_c = NULL;
    run->depth_c = NULL;
    run->export_attn_c = NULL;
    run->export_gamma_c = NULL;
    run->export_magg_c = NULL;
    {
        char *e = getenv("TEMP_ON_MORT");
        run->temp_on_mort = e ? atoi(e) : 1;
    }
    {
        char fname_forcing[256];
        snprintf(fname_forcing, sizeof(fname_forcing), "%s/Input/forcing_ts.txt", run->filename);
        run->fptr_forcing = fopen(fname_forcing, "r");
        if (run->fptr_forcing != NULL)
        {
            run->forcing_flag = 1;
            char hdr[512];
            if (fgets(hdr, sizeof(hdr), run->fptr_forcing) != NULL)
            {
                /*strip trailing newline*/
                size_t hl = strlen(hdr);
                while (hl > 0 && (hdr[hl - 1] == '\n' || hdr[hl - 1] == '\r' || hdr[hl - 1] == ' '))
                {
                    hdr[hl - 1] = '\0';
                    hl--;
                }
                int col = 0;
                char *name = strtok(hdr, ",");
                while (name != NULL)
                {
                    /*trim leading spaces*/
                    while (*name == ' ')
                        name++;
                    if (strcmp(name, "pel_tempeff") == 0)
                        run->col_pel_tempeff = col;
                    else if (strcmp(name, "ben_tempeff") == 0)
                        run->col_ben_tempeff = col;
                    else if (strcmp(name, "sinking_rate") == 0)
                        run->col_sinking_rate = col;
                    else if (strcmp(name, "depth") == 0)
                        run->col_depth = col;
                    else if (strcmp(name, "residence_time") == 0)
                        run->col_residence_time = col;
                    else if (strcmp(name, "export_attn") == 0)
                        run->col_export_attn = col;
                    else if (strcmp(name, "export_gamma") == 0)
                        run->col_export_gamma = col;
                    else if (strcmp(name, "export_magg") == 0)
                        run->col_export_magg = col;
                    else if (strcmp(name, "gravity") == 0)
                        run->col_gravity = col;
                    col++;
                    name = strtok(NULL, ",");
                }
                run->forcing_ncol = col;
            }
        }
    }
    /* Per-cell forcing (opt-in): FORCING_PERCELL=1 + spatial_dim>=1 => forcing_ts.txt has one row per
       CELL per timestep (x*y rows, ordered x-outer y-inner, like the plankton _ts file), so temperature
       and export can vary across space. Allocate the per-cell multiplier arrays (init neutral=1). When
       off (default), the single global row is used and these arrays are unused -> behaviour unchanged. */
    if (run->forcing_flag == 1)
    {
        char *e = getenv("FORCING_PERCELL");
        run->forcing_percell = (e && atoi(e) != 0 && (run->spatial_dim == 1 || run->spatial_dim == 2)) ? 1 : 0;

        run->pel_tempeff_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->ben_tempeff_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->sinking_rate_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->depth_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->export_attn_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->export_gamma_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->export_magg_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            run->pel_tempeff_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->ben_tempeff_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->sinking_rate_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->depth_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->export_attn_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->export_gamma_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->export_magg_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                run->pel_tempeff_c[x_idx][y_idx] = 1.0;
                run->ben_tempeff_c[x_idx][y_idx] = 1.0;
                run->sinking_rate_c[x_idx][y_idx] = 1.0;
                run->depth_c[x_idx][y_idx] = 1.0;
                run->export_attn_c[x_idx][y_idx] = 0.0;
                run->export_gamma_c[x_idx][y_idx] = 0.0;
                run->export_magg_c[x_idx][y_idx] = -1e30;
            }
        }
    }

    /*Create temporary arrays for timestepping process*/
    double ***temp_plankton;
    double ****temp_pelagic;
    double ****temp_benthic;
    double **temp_detritus;

    temp_plankton = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        temp_plankton[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            temp_plankton[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                temp_plankton[sc_idx][x_idx][y_idx] = 0;
            }
        }
    }
    temp_pelagic = (double ****)safe_malloc(n_pel_species * sizeof(double ***), __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        temp_pelagic[pel_species_idx] = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
        for (sc_idx = 0; sc_idx < n_m; sc_idx++)
        {
            temp_pelagic[pel_species_idx][sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                temp_pelagic[pel_species_idx][sc_idx][x_idx] =
                    (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    temp_pelagic[pel_species_idx][sc_idx][x_idx][y_idx] = 0;
                }
            }
        }
    }

    if (run->no_benthic != 0)
    {
        temp_benthic = (double ****)safe_malloc(n_ben_species * sizeof(double ***), __LINE__, __FILE__);
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            temp_benthic[ben_species_idx] = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
            for (sc_idx = 0; sc_idx < n_m; sc_idx++)
            {
                temp_benthic[ben_species_idx][sc_idx] =
                    (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    temp_benthic[ben_species_idx][sc_idx][x_idx] =
                        (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        temp_benthic[ben_species_idx][sc_idx][x_idx][y_idx] = 0;
                    }
                }
            }
        }
        temp_detritus = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            temp_detritus[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                temp_detritus[x_idx][y_idx] = 0;
            }
        }
    }

    /*Calculate initial fishing pressures*/
    calculate_fishing(run, grid, community);

    /*Calculate starting g and mu values*/
    calculate_g_and_mu(run, grid, community);

    /*Calculate starting reproduction allocation*/
    calculate_reproduction(run, grid, community);

    /*Print run info and starting output values*/
    /*plankton*/
    if (grid->toutmin == 0 && run->return_memory)
    {
        record_output(run, grid, community, out, 0);
    }

    if (run->write_files)
    {
        if (grid->toutmin == 0)
        {
            print_timestep_plankton(grid, community->plankton, 0);
            print_plankton_summary(grid, community->plankton, 0);

            /*Pelagic*/
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                print_timestep_pelagic(run, grid, &(community->pelagic[pel_species_idx]), 0);
                print_pelagic_summary(grid, &(community->pelagic[pel_species_idx]), 0);
            }

            if (run->no_benthic != 0)
            {
                /*Benthic*/
                for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
                {
                    print_timestep_benthic(grid, &(community->benthic[ben_species_idx]), 0);
                    print_benthic_summary(grid, &(community->benthic[ben_species_idx]), 0);
                }
                /*Detritus*/
                print_detritus_summary(grid, community->detritus, 0);
            }
        }
    }

    /*Parameter File*/
    run->fptr_summ = safe_fopen(run->fname_summ, "a", __LINE__, __FILE__);

    print_run(run, run->fptr_summ);
    print_grid(grid, run->fptr_summ);
    print_plankton(community->plankton, run->fptr_summ);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        print_pelagic(&(community->pelagic[pel_species_idx]), run->fptr_summ);
    }
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            print_benthic(&(community->benthic[ben_species_idx]), run->fptr_summ);
        }
        print_detritus(community->detritus, run->fptr_summ);
    }
    fclose(run->fptr_summ);

    /*Time Stepping Loop*/
    for (ts_idx = 1; ts_idx < n_t; ts_idx++)
    {

        /*Countdown timer for the program*/

        // Checks if user has terminated function in R (e.g. by pressing escape)
        R_CheckUserInterrupt();

        /*Countdown timer for the program*/
        if (ts_idx == 1)
        {
            time2 = (double)clock() / CLOCKS_PER_SEC;
        }
        else
        {
            if (is_output_timestep(grid, ts_idx))
            {
                time1 = time2;
                time2 = (double)clock() / CLOCKS_PER_SEC;
                timed = (time2 - time1) * (grid->tstep / grid->toutstep);
            }
            Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            R_FlushConsole();
            Rprintf("%6.1f   %d", timed * (grid->tnum - ts_idx), grid->tnum - ts_idx);
            R_FlushConsole();
        }
        if (ts_idx == (n_t - 1))
        {
            Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            R_FlushConsole();
            Rprintf("\n");
            R_FlushConsole();
        }

        /*Store previous time step in temporary array for each species*/
        for (sc_idx = 0; sc_idx < n_m; sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    temp_plankton[sc_idx][x_idx][y_idx] = community->plankton->u_values[sc_idx][x_idx][y_idx];
                }
            }
        }
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            for (sc_idx = 0; sc_idx < n_m; sc_idx++)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        temp_pelagic[pel_species_idx][sc_idx][x_idx][y_idx] =
                            community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx];
                    }
                }
            }
        }

        if (run->no_benthic != 0)
        {
            for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
            {
                for (sc_idx = 0; sc_idx < n_m; sc_idx++)
                {
                    for (x_idx = 0; x_idx < n_x; x_idx++)
                    {
                        for (y_idx = 0; y_idx < n_y; y_idx++)
                        {
                            temp_benthic[ben_species_idx][sc_idx][x_idx][y_idx] =
                                community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx];
                        }
                    }
                }
            }
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    temp_detritus[x_idx][y_idx] = community->detritus->w_values[x_idx][y_idx];
                }
            }
        }

        /*Calculate u(m,t) using previous timestep values*/
        mass_solver(run, grid, community, pelmatrix, benmatrix, temp_plankton, temp_pelagic, temp_benthic,
                    temp_detritus);

        /*Calculate fishing pressures*/
        calculate_fishing(run, grid, community);

        /*Calculate g and mu using current timestep values for use in movement and next time step*/
        calculate_g_and_mu(run, grid, community);

        /*Calculate starting reproduction allocation*/
        calculate_reproduction(run, grid, community);

        /*If Spatial Model is being used (1D = x movement only; 2D = x then y).
          NOTE: the outer guard must admit spatial_dim==2 as well — the inner spatial_dim==2
          block below only runs when this outer test passes, so a `==1`-only guard left ALL
          movement dead for 2D runs (neither x nor y solver ran). */
        if (run->spatial_dim == 1 || run->spatial_dim == 2)
        {

            /*Store itermediate time step as temp_pelagic*/

            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                for (sc_idx = 0; sc_idx < n_m; sc_idx++)
                {
                    for (x_idx = 0; x_idx < n_x; x_idx++)
                    {
                        for (y_idx = 0; y_idx < n_y; y_idx++)
                        {
                            temp_pelagic[pel_species_idx][sc_idx][x_idx][y_idx] =
                                community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx];
                        }
                    }
                }
            }

            /*Calculate spatial movements using intermediate values*/
            xmove_solver(run, grid, community, xmatrix, temp_pelagic, ts_idx);

            if (run->spatial_dim == 2)
            {

                /*Store itermediate time step as temp_pelagic*/
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
                    {
                        for (x_idx = 0; x_idx < n_x; x_idx++)
                        {
                            for (y_idx = 0; y_idx < n_y; y_idx++)
                            {
                                temp_pelagic[pel_species_idx][sc_idx][x_idx][y_idx] =
                                    community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx];
                            }
                        }
                    }
                }
                /*Calculate spatial movements using intermediate values*/
                ymove_solver(run, grid, community, ymatrix, temp_pelagic, ts_idx);
            }
        }

        /*Print current timestep values to files*/
        if (is_output_timestep(grid, ts_idx))
        {
            if (run->write_files)
            {
                print_timestep_plankton(grid, community->plankton, ts_idx);
                print_plankton_summary(grid, community->plankton, ts_idx);
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    print_timestep_pelagic(run, grid, &(community->pelagic[pel_species_idx]), ts_idx);
                    print_pelagic_summary(grid, &(community->pelagic[pel_species_idx]), ts_idx);
                }
                if (run->no_benthic != 0)
                {
                    for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
                    {
                        print_timestep_benthic(grid, &(community->benthic[ben_species_idx]), ts_idx);
                        print_benthic_summary(grid, &(community->benthic[ben_species_idx]), ts_idx);
                    }
                    print_detritus_summary(grid, community->detritus, ts_idx);
                }
            }

            if (run->return_memory)
            {
                record_output(run, grid, community, out, ts_idx);
            }
        }
    }
    if (run->write_files)
    {
        /*Print final detailed timestep*/
        print_detailed_plankton(grid, community->plankton, n_t - 1);
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            print_detailed_pelagic(run, grid, &(community->pelagic[pel_species_idx]), n_t - 1);
        }
        if (run->no_benthic != 0)
        {
            for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
            {
                print_detailed_benthic(grid, &(community->benthic[ben_species_idx]), n_t - 1);
            }
        }
    }

    /*Free pointers*/
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            free(temp_plankton[sc_idx][x_idx]);
        }
        free(temp_plankton[sc_idx]);
    }
    free(temp_plankton);

    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        for (sc_idx = 0; sc_idx < n_m; sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                free(temp_pelagic[pel_species_idx][sc_idx][x_idx]);
            }
            free(temp_pelagic[pel_species_idx][sc_idx]);
        }
        free(temp_pelagic[pel_species_idx]);
    }
    free(temp_pelagic);

    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            for (sc_idx = 0; sc_idx < n_m; sc_idx++)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    free(temp_benthic[ben_species_idx][sc_idx][x_idx]);
                }
                free(temp_benthic[ben_species_idx][sc_idx]);
            }
            free(temp_benthic[ben_species_idx]);
        }
        free(temp_benthic);

        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            free(temp_detritus[x_idx]);
        }
        free(temp_detritus);
    }

    /*Close ts files*/
    if (community->plankton->ts_flag == 1)
    {
        fclose(community->plankton->fptr_ts);
    }
    if (run->forcing_flag == 1)
    {
        fclose(run->fptr_forcing);
        if (run->pel_tempeff_c != NULL)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                free(run->pel_tempeff_c[x_idx]);
                free(run->ben_tempeff_c[x_idx]);
                free(run->sinking_rate_c[x_idx]);
                free(run->depth_c[x_idx]);
                free(run->export_attn_c[x_idx]);
                free(run->export_gamma_c[x_idx]);
                free(run->export_magg_c[x_idx]);
            }
            free(run->pel_tempeff_c);
            free(run->ben_tempeff_c);
            free(run->sinking_rate_c);
            free(run->depth_c);
            free(run->export_attn_c);
            free(run->export_gamma_c);
            free(run->export_magg_c);
        }
    }
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (community->pelagic[pel_species_idx].ts_flag == 1)
        {
            fclose(community->pelagic[pel_species_idx].fptr_ts);
        }
    }
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            if (community->benthic[ben_species_idx].ts_flag == 1)
            {
                fclose(community->benthic[ben_species_idx].fptr_ts);
            }
        }
        if (community->detritus->ts_flag == 1)
        {
            fclose(community->detritus->fptr_ts);
        }
    }

    /*Close fishing files*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (community->pelagic[pel_species_idx].fishing_flag == 1)
        {
            fclose(community->pelagic[pel_species_idx].fptr_fish_ts);
        }
    }
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            if (community->benthic[ben_species_idx].fishing_flag == 1)
            {
                fclose(community->benthic[ben_species_idx].fptr_fish_ts);
            }
        }
    }

    /*Close reproduction files*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (community->pelagic[pel_species_idx].rep_method == 1)
        {
            if (community->pelagic[pel_species_idx].ts_flag == 0)
            {
                fclose(community->pelagic[pel_species_idx].fptr_rep_ts);
            }
        }
    }
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            if (community->benthic[ben_species_idx].rep_method == 1)
            {
                if (community->benthic[ben_species_idx].ts_flag == 0)
                {
                    fclose(community->benthic[ben_species_idx].fptr_rep_ts);
                }
            }
        }
    }
}

//-------------------------------------------------//
// Differencing Scheme Functions. Not to be edited //
//-------------------------------------------------//

void mass_solver(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *pelmatrix, MATRIX *benmatrix,
                 double ***temp_plankton, double ****temp_pelagic, double ****temp_benthic, double **temp_detritus)
/*Impliments the differencing scheme to obtain the values at the next time step*/
{
    int sc_idx, x_idx, y_idx, pel_species_idx, ben_species_idx;
    int n_x, n_y, n_pel_species, n_ben_species;
    double dt, dm;

    dm = grid->mstep;
    dt = grid->tstep;

    n_x = grid->xnum;
    n_y = grid->ynum;
    n_pel_species = run->no_pelagic;
    n_ben_species = run->no_benthic;

    /*Advance the environmental-forcing plugin one timestep (if active), in lockstep with the
      plankton re-read below. Parse this timestep's row into the recognised channels, then:
        - rescale every species' feeding (A) and background mortality (mu_0) from their stored
          base values by pel/ben_tempeff - the growth, predation-mortality and reproduction
          integrals computed in calculate_g_and_mu() (after this solve) then carry the
          temperature effect, matching sizemodel(): feeding & other-mortality are temperature-
          scaled; senescence (mu_s) and fishing are not;
        - store sinking_rate on run for use in g_det() (detritus input).*/
    if (run->forcing_flag == 1 && run->forcing_percell)
    {
        /* PER-CELL forcing: read x*y rows this timestep (one per cell, x-outer y-inner) into the
           per-cell temperature/sinking arrays. A and mu_0 are held at BASE here; the per-cell
           temperature is applied in the feeding/mortality integrals (via run->*_tempeff_c[x_idx][l]).
           Only pel_tempeff/ben_tempeff/sinking_rate are made per-cell in this mode. */
        char fline[1024];
        double vals[64];
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                if (fgets(fline, sizeof(fline), run->fptr_forcing) == NULL)
                    break;
                int nv = 0;
                char *tok = strtok(fline, ",");
                while (tok != NULL && nv < 64)
                {
                    vals[nv++] = atof(tok);
                    tok = strtok(NULL, ",");
                }
                run->pel_tempeff_c[x_idx][y_idx] =
                    (run->col_pel_tempeff >= 0 && run->col_pel_tempeff < nv) ? vals[run->col_pel_tempeff] : 1.0;
                run->ben_tempeff_c[x_idx][y_idx] =
                    (run->col_ben_tempeff >= 0 && run->col_ben_tempeff < nv) ? vals[run->col_ben_tempeff] : 1.0;
                run->sinking_rate_c[x_idx][y_idx] =
                    (run->col_sinking_rate >= 0 && run->col_sinking_rate < nv) ? vals[run->col_sinking_rate] : 1.0;
                run->depth_c[x_idx][y_idx] = (run->col_depth >= 0 && run->col_depth < nv) ? vals[run->col_depth] : 1.0;
                run->export_attn_c[x_idx][y_idx] =
                    (run->col_export_attn >= 0 && run->col_export_attn < nv) ? vals[run->col_export_attn] : 0.0;
                run->export_gamma_c[x_idx][y_idx] =
                    (run->col_export_gamma >= 0 && run->col_export_gamma < nv) ? vals[run->col_export_gamma] : 0.0;
                run->export_magg_c[x_idx][y_idx] =
                    (run->col_export_magg >= 0 && run->col_export_magg < nv) ? vals[run->col_export_magg] : -1e30;
            }
        }
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            community->pelagic[pel_species_idx].A = community->pelagic[pel_species_idx].A_base;
            community->pelagic[pel_species_idx].mu_0 = community->pelagic[pel_species_idx].mu_0_base;
        }
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            community->benthic[ben_species_idx].A = community->benthic[ben_species_idx].A_base;
            community->benthic[ben_species_idx].mu_0 = community->benthic[ben_species_idx].mu_0_base;
        }
    }
    else if (run->forcing_flag == 1)
    {
        char fline[1024];
        if (fgets(fline, sizeof(fline), run->fptr_forcing) != NULL)
        {
            double vals[64];
            int nv = 0;
            char *tok = strtok(fline, ",");
            while (tok != NULL && nv < 64)
            {
                vals[nv++] = atof(tok);
                tok = strtok(NULL, ",");
            }
            if (run->col_pel_tempeff >= 0 && run->col_pel_tempeff < nv)
                run->pel_tempeff = vals[run->col_pel_tempeff];
            if (run->col_ben_tempeff >= 0 && run->col_ben_tempeff < nv)
                run->ben_tempeff = vals[run->col_ben_tempeff];
            if (run->col_sinking_rate >= 0 && run->col_sinking_rate < nv)
                run->sinking_rate = vals[run->col_sinking_rate];
            if (run->col_depth >= 0 && run->col_depth < nv)
                run->depth = vals[run->col_depth];
            if (run->col_residence_time >= 0 && run->col_residence_time < nv)
                run->residence_time = vals[run->col_residence_time];
            if (run->col_export_attn >= 0 && run->col_export_attn < nv)
                run->export_attn = vals[run->col_export_attn];
            if (run->col_export_gamma >= 0 && run->col_export_gamma < nv)
                run->export_gamma = vals[run->col_export_gamma];
            if (run->col_export_magg >= 0 && run->col_export_magg < nv)
                run->export_magg = vals[run->col_export_magg];
            if (run->col_gravity >= 0 && run->col_gravity < nv)
                run->gravity = vals[run->col_gravity];
        }
        /* TEMP_ON_MORT (env, default 1): whether temperature scales background mortality.
           =0 -> temperature drives FEEDING only; mu_0 held at base. Rationale: the model does
           not resolve metabolic maintenance costs (respiration) or their temperature response,
           so mu_0*tempeff is a lone, unbalanced "cost of warmth"; feeding has a clean
           encounter-rate justification, mortality does not. Read once. */
        static int tom = -1;
        if (tom < 0)
        {
            char *e = getenv("TEMP_ON_MORT");
            tom = e ? atoi(e) : 1;
        }
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            community->pelagic[pel_species_idx].A = community->pelagic[pel_species_idx].A_base * run->pel_tempeff;
            community->pelagic[pel_species_idx].mu_0 =
                tom ? community->pelagic[pel_species_idx].mu_0_base * run->pel_tempeff
                    : community->pelagic[pel_species_idx].mu_0_base;
        }
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            community->benthic[ben_species_idx].A = community->benthic[ben_species_idx].A_base * run->ben_tempeff;
            community->benthic[ben_species_idx].mu_0 =
                tom ? community->benthic[ben_species_idx].mu_0_base * run->ben_tempeff
                    : community->benthic[ben_species_idx].mu_0_base;
        }
    }

    /*Used for calculating detritus updates*/
    double **temp_det_g;
    double **temp_det_mu;

    temp_det_g = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    temp_det_mu = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        temp_det_g[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        temp_det_mu[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            temp_det_g[x_idx][y_idx] = 0;
            temp_det_mu[x_idx][y_idx] = 0;
        }
    }
    /*******************/
    /* Plankton System */
    /*******************/
    if (community->plankton->ts_flag == 0)
    {
        // Calculate values for Plankton Size Spectra using internal dynamics //
        for (sc_idx = community->plankton->iplamin; sc_idx < (community->plankton->iplamax + 1); sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    community->plankton->u_values[sc_idx][x_idx][y_idx] =
                        temp_plankton[sc_idx][x_idx][y_idx]; // All plankton systems are invariant
                }
            }
        }
    }
    if (community->plankton->ts_flag == 1)
    {
        int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
        char *input;                   // input string to store line as string
        input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
        char sep[] = ",";    // identifying the ',' char as seperating values
        char *result = NULL; // declaring a pointer for storing each individual value as a char

        /*Read in a single timestep*/
        /* each line one at a time (corresponding to spectra at a spatial point)*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                fgets(input, MAX_LEN, community->plankton->fptr_ts);
                input[strlen(input) - 1] =
                    '\0'; // replacing the newline \n char of the line of text with a null \0 char

                result = (char *)strtok(input, sep); // this is the 1st value
                sc_idx = 0;
                while (result != NULL)
                {
                    community->plankton->u_values[sc_idx][x_idx][y_idx] = atof(result);
                    result = (char *)strtok(NULL, sep);
                    sc_idx++;
                }
            }
        }
        free(input);
    }

    /******************/
    /* Pelagic System */
    /******************/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (community->pelagic[pel_species_idx].ts_flag == 0)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    SKIP_IF_LAND(grid, x_idx, y_idx);

                    for (sc_idx = 0; sc_idx < pelmatrix[pel_species_idx].size; sc_idx++)
                    {
                        pelmatrix[pel_species_idx].a[sc_idx] =
                            -(dt / dm) *
                            community->pelagic[pel_species_idx]
                                .g_values[sc_idx - 1 + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx];
                        pelmatrix[pel_species_idx].b[sc_idx] =
                            (1 +
                             (dt * community->pelagic[pel_species_idx]
                                       .mu_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx]) +
                             (dt / dm) *
                                 community->pelagic[pel_species_idx]
                                     .g_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx]);
                        pelmatrix[pel_species_idx].c[sc_idx] = 0;
                        pelmatrix[pel_species_idx].r[sc_idx] =
                            temp_pelagic[pel_species_idx][sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx]
                                        [y_idx];
                    }
                    pelmatrix[pel_species_idx].a[0] = 0;
                    tridag(&(pelmatrix[pel_species_idx]));
                    for (sc_idx = 0; sc_idx < pelmatrix[pel_species_idx].size; sc_idx++)
                    {
                        community->pelagic[pel_species_idx]
                            .u_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx] =
                            pelmatrix[pel_species_idx].u[sc_idx];
                    }
                }
            }
        }
        if (community->pelagic[pel_species_idx].ts_flag == 1)
        {
            int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
            char *input;                   // input string to store line as string
            input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
            char sep[] = ",";    // identifying the ',' char as seperating values
            char *result = NULL; // declaring a pointer for storing each individual value as a char

            /*Read in a single timestep*/
            /* each line one at a time (corresponding to spectra at a spatial point)*/
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    fgets(input, MAX_LEN, community->pelagic[pel_species_idx].fptr_ts);
                    input[strlen(input) - 1] =
                        '\0'; // replacing the newline \n char of the line of text with a null \0 char

                    result = (char *)strtok(input, sep); // this is the 1st value
                    sc_idx = 0;
                    while (result != NULL)
                    {
                        community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] = atof(result);
                        result = (char *)strtok(NULL, sep);
                        sc_idx++;
                    }
                }
            }
            free(input);
        }
    }

    /******************/
    /* Benthic System */
    /******************/
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            if (community->benthic[ben_species_idx].ts_flag == 0)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        SKIP_IF_LAND(grid, x_idx, y_idx); 

                        for (sc_idx = 0; sc_idx < benmatrix[ben_species_idx].size; sc_idx++)
                        {
                            benmatrix[ben_species_idx].a[sc_idx] =
                                -(dt / dm) *
                                community->benthic[ben_species_idx]
                                    .g_values[sc_idx - 1 + community->benthic[ben_species_idx].ibenmin][x_idx][y_idx];
                            benmatrix[ben_species_idx].b[sc_idx] =
                                (1 +
                                 (dt *
                                  community->benthic[ben_species_idx]
                                      .mu_values[sc_idx + community->benthic[ben_species_idx].ibenmin][x_idx][y_idx]) +
                                 (dt / dm) *
                                     community->benthic[ben_species_idx]
                                         .g_values[sc_idx + community->benthic[ben_species_idx].ibenmin][x_idx][y_idx]);
                            benmatrix[ben_species_idx].c[sc_idx] = 0;
                            benmatrix[ben_species_idx].r[sc_idx] =
                                temp_benthic[ben_species_idx][sc_idx + community->benthic[ben_species_idx].ibenmin]
                                            [x_idx][y_idx];
                        }
                        benmatrix[ben_species_idx].a[0] = 0;
                        tridag(&(benmatrix[ben_species_idx]));
                        for (sc_idx = 0; sc_idx < benmatrix[ben_species_idx].size; sc_idx++)
                        {
                            community->benthic[ben_species_idx]
                                .u_values[sc_idx + community->benthic[ben_species_idx].ibenmin][x_idx][y_idx] =
                                benmatrix[ben_species_idx].u[sc_idx];
                        }
                    }
                }
            }
            if (community->benthic[ben_species_idx].ts_flag == 1)
            {
                int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
                char *input;                   // input string to store line as string
                input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
                char sep[] = ",";    // identifying the ',' char as seperating values
                char *result = NULL; // declaring a pointer for storing each individual value as a char

                /*Read in a single timestep*/
                /* each line one at a time (corresponding to spectra at a spatial point)*/
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        fgets(input, MAX_LEN, community->benthic[ben_species_idx].fptr_ts);
                        input[strlen(input) - 1] =
                            '\0'; // replacing the newline \n char of the line of text with a null \0 char

                        result = (char *)strtok(input, sep); // this is the 1st value
                        sc_idx = 0;
                        while (result != NULL)
                        {
                            community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] = atof(result);
                            result = (char *)strtok(NULL, sep);
                            sc_idx++;
                        }
                    }
                }
                free(input);
            }
        }

        /*******************/
        /* Detritus System */
        /*******************/
        if (community->detritus->ts_flag == 0)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    SKIP_IF_LAND(grid, x_idx, y_idx);

                    /* Analytical (exponential) integrator for the detritus pool.
                       Over a step the fish/plankton state is fixed, so detritus obeys a
                       LINEAR ODE  dW/dt = I - kappa*W , with input I = g_values (constant)
                       and loss mu_values = kappa*W (detritivory is proportional to W).
                       The exact solution  W(t+dt) = W* + (W0-W*)*exp(-kappa*dt),
                       W* = I/kappa,  is UNCONDITIONALLY STABLE and provably non-negative
                       (a convex combination of W0>=0 and W*>=0) - unlike the explicit
                       Euler step W += dt*(I-kappa*W), which overshoots to negative when
                       kappa*dt > 1 and drives the NaN/collapse cascade. Removes the
                       timestep stability limit on the detritus pool. */
                    {
                        double W0 = temp_detritus[x_idx][y_idx];
                        double In = community->detritus->g_values[x_idx][y_idx]; /* gross input rate I */
                        /* Dunne et al. (2007) eq.3 burial: a fraction of the incoming detrital
                           flux is buried (permanently lost) before joining the pool, EXACTLY as
                           sizemodel does: burial = input_w*(0.013+0.53*iw^2/(7+iw)^2) applied to
                           the model's native input_w. dbpmr's densities (like sizemodel's) are
                           per-VOLUME - anchored to the depth-averaged plankton concentration - so
                           input_w is a per-volume flux and the Dunne fraction is applied to it
                           DIRECTLY, with no areal (x depth) rescale, matching sizemodel. Burial is
                           OPT-IN: only when a `depth` forcing channel is supplied (col_depth>=0),
                           so existing runs without the forcing plugin keep their original dynamics. */
                        if (run->col_depth >= 0)
                        {
                            double bfrac = 0.013 + 0.53 * In * In / ((7.0 + In) * (7.0 + In));
                            In = In * (1.0 - bfrac); /* net of burial   */
                        }
                        double loss = community->detritus->mu_values[x_idx][y_idx]; /* = kappa*W0   */
                        double kappa = (W0 > 1e-30) ? loss / W0 : 0.0;
                        /* Alternative/opt-in residence-time closure: a first-order loss W/tau on the
                           detritus POOL (remineralisation/export from the benthic layer), added to the
                           detritivory rate. Dimensionally clean in dbpmr's per-volume units - no areal
                           ambiguity - and the exact integrator below absorbs it (K = kappa + 1/tau, so
                           the equilibrium is W* = In/(kappa+1/tau); with negligible detritivory that is
                           the classic W* = In*tau). Supplied via the residence_time channel (yr); an
                           alternative to Dunne burial (Julia's commented CMIP5 variant made this loss
                           proportional to W). Off when the channel is absent or tau<=0. */
                        if (run->col_residence_time >= 0 && run->residence_time > 0.0)
                        {
                            kappa += 1.0 / run->residence_time;
                        }
                        if (kappa > 1e-30)
                        {
                            double Wstar = In / kappa;
                            community->detritus->w_values[x_idx][y_idx] = Wstar + (W0 - Wstar) * exp(-kappa * dt);
                        }
                        else
                        {
                            community->detritus->w_values[x_idx][y_idx] = W0 + In * dt; /* kappa -> 0 limit */
                        }
                        /* defensive backstop only against a non-finite upstream rate */
                        if (!isfinite(community->detritus->w_values[x_idx][y_idx]) ||
                            community->detritus->w_values[x_idx][y_idx] < 0)
                        {
                            community->detritus->w_values[x_idx][y_idx] = 0;
                        }
                    }
                    temp_det_g[x_idx][y_idx] = community->detritus->g_values[x_idx][y_idx];
                    temp_det_mu[x_idx][y_idx] = community->detritus->mu_values[x_idx][y_idx];
                }
            }
        }
        if (community->detritus->ts_flag == 1)
        {
            int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
            char *input;                   // input string to store line as string
            input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
            char sep[] = ",";    // identifying the ',' char as seperating values
            char *result = NULL; // declaring a pointer for storing each individual value as a char

            /*Read in a single timestep*/
            /* each line one at a time (corresponding to spectra at a spatial point)*/
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    fgets(input, MAX_LEN, community->detritus->fptr_ts);
                    input[strlen(input) - 1] =
                        '\0'; // replacing the newline \n char of the line of text with a null \0 char

                    result = (char *)strtok(input, sep); // this is the 1st value
                    while (result != NULL)
                    {
                        community->detritus->w_values[x_idx][y_idx] = atof(result);
                        result = (char *)strtok(NULL, sep);
                    }
                }
            }
            free(input);
        }
    }

    // Calculate fully implicit method if required
    if (run->diff_method == 1)
    {
        // Use intermediate u values to calculate new g and mu for full implicit methods
        calculate_g_and_mu(run, grid, community);

        // Use new g and mu to claculate final u values
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            if (community->pelagic[pel_species_idx].ts_flag == 0)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        SKIP_IF_LAND(grid, x_idx, y_idx);

                        for (sc_idx = 0; sc_idx < pelmatrix[pel_species_idx].size; sc_idx++)
                        {
                            pelmatrix[pel_species_idx].a[sc_idx] =
                                -(dt / dm) *
                                community->pelagic[pel_species_idx]
                                    .g_values[sc_idx - 1 + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx];
                            pelmatrix[pel_species_idx].b[sc_idx] =
                                (1 +
                                 (dt *
                                  community->pelagic[pel_species_idx]
                                      .mu_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx]) +
                                 (dt / dm) *
                                     community->pelagic[pel_species_idx]
                                         .g_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx]);
                            pelmatrix[pel_species_idx].c[sc_idx] = 0;
                            pelmatrix[pel_species_idx].r[sc_idx] =
                                temp_pelagic[pel_species_idx][sc_idx + community->pelagic[pel_species_idx].ipelmin]
                                            [x_idx][y_idx];
                        }
                        pelmatrix[pel_species_idx].a[0] = 0;
                        tridag(&(pelmatrix[pel_species_idx]));
                        for (sc_idx = 0; sc_idx < pelmatrix[pel_species_idx].size; sc_idx++)
                        {
                            community->pelagic[pel_species_idx]
                                .u_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx] =
                                pelmatrix[pel_species_idx].u[sc_idx];
                        }
                    }
                }
            }
        }

        if (run->no_benthic != 0)
        {
            for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
            {
                if (community->benthic[ben_species_idx].ts_flag == 0)
                {
                    for (x_idx = 0; x_idx < n_x; x_idx++)
                    {
                        for (y_idx = 0; y_idx < n_y; y_idx++)
                        {
                            SKIP_IF_LAND(grid, x_idx, y_idx);

                            for (sc_idx = 0; sc_idx < benmatrix[ben_species_idx].size; sc_idx++)
                            {
                                benmatrix[ben_species_idx].a[sc_idx] =
                                    -(dt / dm) * community->benthic[ben_species_idx]
                                                     .g_values[sc_idx - 1 + community->benthic[ben_species_idx].ibenmin]
                                                              [x_idx][y_idx];
                                benmatrix[ben_species_idx].b[sc_idx] =
                                    (1 +
                                     (dt * community->benthic[ben_species_idx]
                                               .mu_values[sc_idx + community->benthic[ben_species_idx].ibenmin][x_idx]
                                                         [y_idx]) +
                                     (dt / dm) * community->benthic[ben_species_idx]
                                                     .g_values[sc_idx + community->benthic[ben_species_idx].ibenmin]
                                                              [x_idx][y_idx]);
                                benmatrix[ben_species_idx].c[sc_idx] = 0;
                                benmatrix[ben_species_idx].r[sc_idx] =
                                    temp_benthic[ben_species_idx][sc_idx + community->benthic[ben_species_idx].ibenmin]
                                                [x_idx][y_idx];
                            }
                            benmatrix[ben_species_idx].a[0] = 0;
                            tridag(&(benmatrix[ben_species_idx]));
                            for (sc_idx = 0; sc_idx < benmatrix[ben_species_idx].size; sc_idx++)
                            {
                                community->benthic[ben_species_idx]
                                    .u_values[sc_idx + community->benthic[ben_species_idx].ibenmin][x_idx][y_idx] =
                                    benmatrix[ben_species_idx].u[sc_idx];
                            }
                        }
                    }
                }
            }

            if (community->detritus->ts_flag == 0)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        SKIP_IF_LAND(grid, x_idx, y_idx);

                        community->detritus->w_values[x_idx][y_idx] =
                            temp_detritus[x_idx][y_idx] +
                            (dt / 2) * ((community->detritus->g_values[x_idx][y_idx] + temp_det_g[x_idx][y_idx]) -
                                        (community->detritus->mu_values[x_idx][y_idx] + temp_det_mu[x_idx][y_idx]));
                    }
                }
            }
        }
    }

    /***************************/
    /* Add reproduction values */
    /***************************/
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            SKIP_IF_LAND(grid, x_idx, y_idx);

            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                /*If values aren't given by timeseries*/
                if (community->pelagic[pel_species_idx].ts_flag == 0)
                {
                    /*Specifying Total Eggs*/
                    if (community->pelagic[pel_species_idx].rep_method == 0 ||
                        community->pelagic[pel_species_idx].rep_method == 1)
                    {
                        community->pelagic[pel_species_idx]
                            .u_values[community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx] =
                            community->pelagic[pel_species_idx].reproduction[x_idx][y_idx];
                    }
                    /*Adding Eggs*/
                    if (community->pelagic[pel_species_idx].rep_method == 2 ||
                        community->pelagic[pel_species_idx].rep_method == 3)
                    {
                        community->pelagic[pel_species_idx]
                            .u_values[community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx] +=
                            community->pelagic[pel_species_idx].reproduction[x_idx][y_idx];
                    }
                }
            }

            if (run->no_benthic != 0)
            {
                for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
                {
                    /*If values aren't given by timeseries*/
                    if (community->benthic[ben_species_idx].ts_flag == 0)
                    {
                        /*Constant Egg amount*/
                        if (community->benthic[ben_species_idx].rep_method == 0)
                        {
                            community->benthic[ben_species_idx]
                                .u_values[community->benthic[ben_species_idx].ibenmin][x_idx][y_idx] =
                                community->benthic[ben_species_idx].u_0 *
                                exp(community->benthic[ben_species_idx].lambda *
                                    community->benthic[ben_species_idx].benmin);
                        }
                        /*Eggs produced according to biomass intake*/
                        if (community->benthic[ben_species_idx].rep_method == 1 ||
                            community->benthic[ben_species_idx].rep_method == 2)
                        {
                            community->benthic[ben_species_idx]
                                .u_values[community->benthic[ben_species_idx].ibenmin][x_idx][y_idx] +=
                                community->benthic[ben_species_idx].reproduction[x_idx][y_idx];
                        }
                    }
                }
            }
        }
    }

    /*Free pointers*/
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        free(temp_det_g[x_idx]);
        free(temp_det_mu[x_idx]);
    }
    free(temp_det_g);
    free(temp_det_mu);
}