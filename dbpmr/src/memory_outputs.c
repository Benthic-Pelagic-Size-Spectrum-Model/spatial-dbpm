#include <memory_outputs.h>

void fill_na(SEXP x) // Helper function to fill R outputs with NA values (this covers land cells)
{
    double *p = REAL(x);
    R_xlen_t n = XLENGTH(x);
    for (R_xlen_t i = 0; i < n; i++) {p[i] = NA_REAL;}
}

//--------------------------------//
// Functions for handling outputs //
//--------------------------------//

void setup_output(RUN *run, GRID *grid, OUTPUT *out, SEXP result)
{
    int ts_idx, sc_idx, x_idx, y_idx, i, pel_species_idx, ben_species_idx;
    out->next_time = 0;
    out->nout_m = out->nout_t = out->nout_x = out->nout_y = 0;
    if (!run->return_memory)
    {
        out->m_sel = out->x_sel = out->y_sel = NULL;
        out->t_out = out->pla_u = NULL;

        return;
    }

    // Output intervals in terms of step sizes (e.g. How many msteps are between each m output step).
    int m_int = output_interval(grid->moutstep, grid->mstep);
    int x_int = output_interval(grid->xoutstep, grid->xstep);
    int y_int = output_interval(grid->youtstep, grid->ystep);

    // Count the number of output steps and add the values to *_sel values
    // (This has to be part of the loop as this is how output steps are already
    // set up in the calculations). Loops iterate, increasing the index by the interval
    // so only interval-indices are added
    for (ts_idx = 0; ts_idx < grid->tnum; ts_idx++)
    {
        if (is_output_timestep(grid, ts_idx))
            out->nout_t++;
    }
    for (sc_idx = 0; sc_idx < grid->mnum - 1; sc_idx += m_int)
    {
        out->nout_m++;
    };
    out->nout_m++;
    for (x_idx = 0; x_idx < grid->xnum; x_idx += x_int)
    {
        out->nout_x++;
    }
    for (y_idx = 0; y_idx < grid->ynum; y_idx += y_int)
    {
        out->nout_y++;
    }

    // Get all of the indices for output steps (these are internal to C - only the dimension values will be returned to
    // R)
    out->m_sel = (int *)R_alloc(out->nout_m, sizeof(int));
    out->x_sel = (int *)R_alloc(out->nout_x, sizeof(int));
    out->y_sel = (int *)R_alloc(out->nout_y, sizeof(int));

    for (i = 0, sc_idx = 0; sc_idx < grid->mnum - 1; sc_idx += m_int)
    {
        out->m_sel[i++] = sc_idx;
    }
    out->m_sel[i] = grid->mnum - 1;
    for (i = 0, x_idx = 0; x_idx < grid->xnum; x_idx += x_int)
    {
        out->x_sel[i++] = x_idx;
    }
    for (i = 0, y_idx = 0; y_idx < grid->ynum; y_idx += y_int)
    {
        out->y_sel[i++] = y_idx;
    }

    // Allocate all of the dimension outputs in R
    SEXP t = PROTECT(allocVector(REALSXP, out->nout_t));
    SET_VECTOR_ELT(result, 0, t);
    out->t_out = REAL(t);

    SEXP m = PROTECT(allocVector(REALSXP, out->nout_m));
    SEXP x = PROTECT(allocVector(REALSXP, out->nout_x));
    SEXP y = PROTECT(allocVector(REALSXP, out->nout_y));

    for (sc_idx = 0; sc_idx < out->nout_m; sc_idx++)
        REAL(m)[sc_idx] = grid->m_values[out->m_sel[sc_idx]];
    for (x_idx = 0; x_idx < out->nout_x; x_idx++)
        REAL(x)[x_idx] = grid->x_values[out->x_sel[x_idx]];
    for (y_idx = 0; y_idx < out->nout_y; y_idx++)
        REAL(y)[y_idx] = grid->y_values[out->y_sel[y_idx]];
    SET_VECTOR_ELT(result, 1, m);
    SET_VECTOR_ELT(result, 2, x);
    SET_VECTOR_ELT(result, 3, y);

    // Allocate holders for pelagic and benthic species dimensions
    SEXP pel_sp = PROTECT(allocVector(REALSXP, run->no_pelagic));
    for (pel_species_idx = 0; pel_species_idx < run->no_pelagic; pel_species_idx++)
        REAL(pel_sp)[pel_species_idx] = pel_species_idx + 1;
    SET_VECTOR_ELT(result, 4, pel_sp);
    SEXP ben_sp = PROTECT(allocVector(REALSXP, run->no_benthic));
    for (ben_species_idx = 0; ben_species_idx < run->no_benthic; ben_species_idx++)
        REAL(ben_sp)[ben_species_idx] = ben_species_idx + 1;
    SET_VECTOR_ELT(result, 5, ben_sp);

    // Allocate the plankton mass density array in R (pla_u)
    SEXP pla = PROTECT(allocVector(REALSXP, (R_xlen_t)out->nout_m * out->nout_x * out->nout_y * out->nout_t));
    SEXP pla_dim = PROTECT(allocVector(INTSXP, 4));
    INTEGER(pla_dim)[0] = out->nout_m;
    INTEGER(pla_dim)[1] = out->nout_x;
    INTEGER(pla_dim)[2] = out->nout_y;
    INTEGER(pla_dim)[3] = out->nout_t;
    setAttrib(pla, R_DimSymbol, pla_dim);
    SET_VECTOR_ELT(result, 6, pla);
    out->pla_u = REAL(pla);

    // Pelagic species //
    // Allocate pelagic species array dimensions (once for the group)
    SEXP pel_dim = PROTECT(allocVector(INTSXP, 5));
    INTEGER(pel_dim)[0] = out->nout_m;
    INTEGER(pel_dim)[1] = out->nout_x;
    INTEGER(pel_dim)[2] = out->nout_y;
    INTEGER(pel_dim)[3] = run->no_pelagic;
    INTEGER(pel_dim)[4] = out->nout_t;

    // Allocate the pelagic mass density array in R (pel_u)
    SEXP pel_u = PROTECT(
        allocVector(REALSXP, (R_xlen_t)out->nout_m * out->nout_x * out->nout_y * run->no_pelagic * out->nout_t));
    setAttrib(pel_u, R_DimSymbol, pel_dim);
    SET_VECTOR_ELT(result, 7, pel_u);
    out->pel_u = REAL(pel_u);

    // Allocate the pelagic mortality (predation) array in R (pel_mu_pred)
    // Predation mortality requires specific dimensions for prey and predator species combinations
    SEXP pel_pred_dim = PROTECT(allocVector(INTSXP, 6));
    INTEGER(pel_pred_dim)[0] = out->nout_m;
    INTEGER(pel_pred_dim)[1] = out->nout_x;
    INTEGER(pel_pred_dim)[2] = out->nout_y;
    INTEGER(pel_pred_dim)[3] = run->no_pelagic;
    INTEGER(pel_pred_dim)[4] = run->no_pelagic;
    INTEGER(pel_pred_dim)[5] = out->nout_t;

    SEXP pel_mu_pred = PROTECT(allocVector(REALSXP, (R_xlen_t)out->nout_m * out->nout_x * out->nout_y *
                                                        run->no_pelagic * run->no_pelagic * out->nout_t));
    setAttrib(pel_mu_pred, R_DimSymbol, pel_pred_dim);
    SET_VECTOR_ELT(result, 8, pel_mu_pred);
    out->pel_mu_pred = REAL(pel_mu_pred);

    // Allocate the pelagic mortality (fishing) array in R (pel_mu_fish)
    SEXP pel_mu_fish = PROTECT(
        allocVector(REALSXP, (R_xlen_t)out->nout_m * out->nout_x * out->nout_y * run->no_pelagic * out->nout_t));
    setAttrib(pel_mu_fish, R_DimSymbol, pel_dim);
    SET_VECTOR_ELT(result, 9, pel_mu_fish);
    out->pel_mu_fish = REAL(pel_mu_fish);

    // Allocate the pelagic growth rate array in R (pel_g)
    SEXP pel_g = PROTECT(
        allocVector(REALSXP, (R_xlen_t)out->nout_m * out->nout_x * out->nout_y * run->no_pelagic * out->nout_t));
    setAttrib(pel_g, R_DimSymbol, pel_dim);
    SET_VECTOR_ELT(result, 10, pel_g);
    out->pel_g = REAL(pel_g);

    // Benthic species //
    // Allocate benthic species array dimensions (once for the group)
    SEXP ben_dim = PROTECT(allocVector(INTSXP, 5));
    INTEGER(ben_dim)[0] = out->nout_m;
    INTEGER(ben_dim)[1] = out->nout_x;
    INTEGER(ben_dim)[2] = out->nout_y;
    INTEGER(ben_dim)[3] = run->no_benthic;
    INTEGER(ben_dim)[4] = out->nout_t;

    // Allocate the benthic mass density array in R (ben_u)
    SEXP ben_u = PROTECT(
        allocVector(REALSXP, (R_xlen_t)out->nout_m * out->nout_x * out->nout_y * run->no_benthic * out->nout_t));
    setAttrib(ben_u, R_DimSymbol, ben_dim);
    SET_VECTOR_ELT(result, 11, ben_u);
    out->ben_u = REAL(ben_u);

    // Allocate the benthic mortality (predation) array in R (ben_mu_pred)
    SEXP ben_mu_pred = PROTECT(
        allocVector(REALSXP, (R_xlen_t)out->nout_m * out->nout_x * out->nout_y * run->no_benthic * out->nout_t));
    setAttrib(ben_mu_pred, R_DimSymbol, ben_dim);
    SET_VECTOR_ELT(result, 12, ben_mu_pred);
    out->ben_mu_pred = REAL(ben_mu_pred);

    // Allocate the benthic mortality (fishing) array in R (ben_mu_fish)
    SEXP ben_mu_fish = PROTECT(
        allocVector(REALSXP, (R_xlen_t)out->nout_m * out->nout_x * out->nout_y * run->no_benthic * out->nout_t));
    setAttrib(ben_mu_fish, R_DimSymbol, ben_dim);
    SET_VECTOR_ELT(result, 13, ben_mu_fish);
    out->ben_mu_fish = REAL(ben_mu_fish);

    // Allocate the benthic growth rate array in R (ben_g)
    SEXP ben_g = PROTECT(
        allocVector(REALSXP, (R_xlen_t)out->nout_m * out->nout_x * out->nout_y * run->no_benthic * out->nout_t));
    setAttrib(ben_g, R_DimSymbol, ben_dim);
    SET_VECTOR_ELT(result, 14, ben_g);
    out->ben_g = REAL(ben_g);

    // Detritus allocations
    // Allocate detritus dimensions (pooled, so no size structure)
    SEXP det_dim = PROTECT(allocVector(INTSXP, 3));
    INTEGER(det_dim)[0] = out->nout_x;
    INTEGER(det_dim)[1] = out->nout_y;
    INTEGER(det_dim)[2] = out->nout_t;

    // Allocate the detritus mass density array in R (det_w)
    SEXP det_w = PROTECT(allocVector(REALSXP, (R_xlen_t)out->nout_x * out->nout_y * out->nout_t));
    setAttrib(det_w, R_DimSymbol, det_dim);
    SET_VECTOR_ELT(result, 15, det_w);
    out->det_w = REAL(det_w);

    // Allocate the detritus growth array in R (det_g)
    SEXP det_g = PROTECT(allocVector(REALSXP, (R_xlen_t)out->nout_x * out->nout_y * out->nout_t));
    setAttrib(det_g, R_DimSymbol, det_dim);
    SET_VECTOR_ELT(result, 16, det_g);
    out->det_g = REAL(det_g);

    // Allocate the detritus out array in R (det_mu)
    SEXP det_mu = PROTECT(allocVector(REALSXP, (R_xlen_t)out->nout_x * out->nout_y * out->nout_t));
    setAttrib(det_mu, R_DimSymbol, det_dim);
    SET_VECTOR_ELT(result, 17, det_mu);
    out->det_mu = REAL(det_mu);

    // Allocate the output name fields in R so that each of the outputs have names
    // This chunk uses C indexing (i.e. above^ + 1).
    SEXP nms = PROTECT(allocVector(STRSXP, 18));
    const char *n[18] = {"time",
                         "m",
                         "x",
                         "y",
                         "pelagic_species",
                         "benthic_species",
                         "plankton_u",
                         "pelagic_u",
                         "pelagic_mu_pred",
                         "pelagic_mu_fish",
                         "pelagic_g",
                         "benthic_u",
                         "benthic_mu_pred",
                         "benthic_mu_fish",
                         "benthic_g",
                         "detritus_w",
                         "detritus_g",
                         "detritus_mu"};
    for (i = 0; i < 18; i++) {SET_STRING_ELT(nms, i, mkChar(n[i]));}
    setAttrib(result, R_NamesSymbol, nms);

    for (i = 6; i < 18; i++) {fill_na(VECTOR_ELT(result, i));}

    UNPROTECT(24);
}

// Calculate whether the current timestep is an output step
static int is_output_timestep(GRID *grid, int ts_idx)
{
    if (ts_idx == 0)
        return grid->toutmin == 0;
    return (grid->t_values[ts_idx] >= grid->toutmin && // Check that the current timestep is greater than toutmin
            grid->t_values[ts_idx] <= grid->toutmax && // Check that the current step is lower than toutmax
            ts_idx % ((int)rint(grid->toutstep / grid->tstep)) == 0 // Check that the current step is an output step
    );
}

// Calculate the interval of model steps within an outstep. e.g. if timestep outstep
// is 1 and the model timestep is 1/48 then the interval is 48. This function ensures
// that outstep is a multiple of step and at least 1 model step.
static int output_interval(double outstep, double step)
{
    int interval = (int)rint(outstep / step);
    if (interval < 1)
        error("an output step must be at least one model step");
    return (interval);
}

static void record_output(RUN *run, GRID *grid, COMMUNITY *c, OUTPUT *out, int ts_idx)
{
    int sc_idx, x_idx, y_idx;
    int pred_idx, pel_species_idx, ben_species_idx;
    if (out->next_time >= out->nout_t)
        error("more output steps than allocated");

    out->t_out[out->next_time] = grid->t_values[ts_idx];

    // Record output for plankton
    for (y_idx = 0; y_idx < out->nout_y; y_idx++)
    {
        for (x_idx = 0; x_idx < out->nout_x; x_idx++)
        {
            SKIP_IF_LAND(grid, out->x_sel[x_idx], out->y_sel[y_idx]);

            for (sc_idx = 0; sc_idx < out->nout_m; sc_idx++)
            {
                out->pla_u[sc_idx + (R_xlen_t)out->nout_m *
                                        (x_idx + out->nout_x * (y_idx + (R_xlen_t)out->nout_y * out->next_time))] =
                    c->plankton->u_values[out->m_sel[sc_idx]][out->x_sel[x_idx]][out->y_sel[y_idx]];
            }
        }
    }

    // Record outputs for pelagic species
    for (y_idx = 0; y_idx < out->nout_y; y_idx++)
    {
        for (x_idx = 0; x_idx < out->nout_x; x_idx++)
        {
            SKIP_IF_LAND(grid, out->x_sel[x_idx], out->y_sel[y_idx]);

            for (sc_idx = 0; sc_idx < out->nout_m; sc_idx++)
            {
                for (pel_species_idx = 0; pel_species_idx < run->no_pelagic; pel_species_idx++)
                {
                    out->pel_u[sc_idx + (R_xlen_t)out->nout_m *
                                            (x_idx + out->nout_x * (y_idx + (R_xlen_t)out->nout_y *
                                                                                (pel_species_idx +
                                                                                 run->no_pelagic * out->next_time)))] =
                        c->pelagic[pel_species_idx].u_values[out->m_sel[sc_idx]][out->x_sel[x_idx]][out->y_sel[y_idx]];

                    for (pred_idx = 0; pred_idx < run->no_pelagic; pred_idx++)
                    {
                        out->pel_mu_pred[sc_idx +
                                         (R_xlen_t)out->nout_m *
                                             (x_idx +
                                              out->nout_x * (y_idx + (R_xlen_t)out->nout_y *
                                                                         (pel_species_idx +
                                                                          (R_xlen_t)run->no_pelagic *
                                                                              (pred_idx + (R_xlen_t)run->no_pelagic *
                                                                                              out->next_time))))] =
                            c->pelagic[pel_species_idx]
                                .mu_pred_values[pred_idx][out->m_sel[sc_idx]][out->x_sel[x_idx]][out->y_sel[y_idx]];
                    }

                    out->pel_mu_fish[sc_idx +
                                     (R_xlen_t)out->nout_m *
                                         (x_idx + out->nout_x * (y_idx + (R_xlen_t)out->nout_y *
                                                                             (pel_species_idx +
                                                                              run->no_pelagic * out->next_time)))] =
                        c->pelagic[pel_species_idx]
                            .mu_fish_values[out->m_sel[sc_idx]][out->x_sel[x_idx]][out->y_sel[y_idx]];

                    out->pel_g[sc_idx + (R_xlen_t)out->nout_m *
                                            (x_idx + out->nout_x * (y_idx + (R_xlen_t)out->nout_y *
                                                                                (pel_species_idx +
                                                                                 run->no_pelagic * out->next_time)))] =
                        c->pelagic[pel_species_idx].g_values[out->m_sel[sc_idx]][out->x_sel[x_idx]][out->y_sel[y_idx]];
                }
            }
        }
    }

    // Record outputs for benthic species
    for (y_idx = 0; y_idx < out->nout_y; y_idx++)
    {
        for (x_idx = 0; x_idx < out->nout_x; x_idx++)
        {
            SKIP_IF_LAND(grid, out->x_sel[x_idx], out->y_sel[y_idx]);

            for (sc_idx = 0; sc_idx < out->nout_m; sc_idx++)
            {
                for (ben_species_idx = 0; ben_species_idx < run->no_benthic; ben_species_idx++)
                {
                    out->ben_u[sc_idx + (R_xlen_t)out->nout_m *
                                            (x_idx + out->nout_x * (y_idx + (R_xlen_t)out->nout_y *
                                                                                (ben_species_idx +
                                                                                 run->no_benthic * out->next_time)))] =
                        c->benthic[ben_species_idx].u_values[out->m_sel[sc_idx]][out->x_sel[x_idx]][out->y_sel[y_idx]];

                    out->ben_mu_pred[sc_idx +
                                     (R_xlen_t)out->nout_m *
                                         (x_idx + out->nout_x * (y_idx + (R_xlen_t)out->nout_y *
                                                                             (ben_species_idx +
                                                                              run->no_benthic * out->next_time)))] =
                        c->benthic[ben_species_idx]
                            .mu_pred_values[out->m_sel[sc_idx]][out->x_sel[x_idx]][out->y_sel[y_idx]];

                    out->ben_mu_fish[sc_idx +
                                     (R_xlen_t)out->nout_m *
                                         (x_idx + out->nout_x * (y_idx + (R_xlen_t)out->nout_y *
                                                                             (ben_species_idx +
                                                                              run->no_benthic * out->next_time)))] =
                        c->benthic[ben_species_idx]
                            .mu_fish_values[out->m_sel[sc_idx]][out->x_sel[x_idx]][out->y_sel[y_idx]];

                    out->ben_g[sc_idx + (R_xlen_t)out->nout_m *
                                            (x_idx + out->nout_x * (y_idx + (R_xlen_t)out->nout_y *
                                                                                (ben_species_idx +
                                                                                 run->no_benthic * out->next_time)))] =
                        c->benthic[ben_species_idx].g_values[out->m_sel[sc_idx]][out->x_sel[x_idx]][out->y_sel[y_idx]];
                }
            }
        }
    }

    if (run->no_benthic != 0) {
        // Record outputs for detritus component
        for (y_idx = 0; y_idx < out->nout_y; y_idx++)
        {
            for (x_idx = 0; x_idx < out->nout_x; x_idx++)
            {
                SKIP_IF_LAND(grid, out->x_sel[x_idx], out->y_sel[y_idx]);

                out->det_w[

                    x_idx + out->nout_x * (y_idx + (R_xlen_t)out->nout_y * out->next_time)

                ] = c->detritus->w_values[out->x_sel[x_idx]][out->y_sel[y_idx]];

                out->det_g[

                    x_idx + out->nout_x * (y_idx + (R_xlen_t)out->nout_y * out->next_time)

                ] = c->detritus->g_values[out->x_sel[x_idx]][out->y_sel[y_idx]];

                out->det_mu[

                    x_idx + out->nout_x * (y_idx + (R_xlen_t)out->nout_y * out->next_time)

                ] = c->detritus->mu_values[out->x_sel[x_idx]][out->y_sel[y_idx]];
            }
        }
    }

    out->next_time++;
}