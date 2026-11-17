#include "memory_outputs.h"


/*Allocate one array of len values per species as element `pos` of an R list (which must be protected),
  and return C pointers to their data*/
static double **species_arrays(SEXP parent, int pos, int nspecies, R_xlen_t len)
{
        int s;
        SEXP l = allocVector(VECSXP, nspecies);
        SET_VECTOR_ELT(parent, pos, l);
        double **ptr = (double **) R_alloc(nspecies, sizeof(double *));
        for (s = 0; s < nspecies; s++) {
                SEXP v = allocVector(REALSXP, len);
                SET_VECTOR_ELT(l, s, v);
                ptr[s] = REAL(v);
        }
        return ptr;
}

/*Allocate the arrays returned to R. Sizes follow the text output: every (xoutstep/xstep)-th and
  (youtstep/ystep)-th cell, every (moutstep/mstep)-th size class plus the last one, and every output time
  in [toutmin,toutmax]. The list is left PROTECTed; the .Call wrapper unprotects it when returning it.*/
void setup_output_arrays(RUN *run, GRID *grid, OUTPUT_ARRAYS *out)
{
        int i, j, k, s, n = run->no_pelagic, c = run->no_benthic;
        int m = grid->mnum, x = grid->xnum, y = grid->ynum, t = grid->tnum;
        int mstride = (int)rint(grid->moutstep / grid->mstep);
        int ostep = (int)rint(grid->toutstep / grid->tstep);
        out->xstride = (int)rint(grid->xoutstep / grid->xstep);
        out->ystride = (int)rint(grid->youtstep / grid->ystep);
        if (mstride < 1) mstride = 1;
        if (ostep < 1) ostep = 1;
        if (out->xstride < 1) out->xstride = 1;
        if (out->ystride < 1) out->ystride = 1;
        out->nxo = (x - 1) / out->xstride + 1;
        out->nyo = (y - 1) / out->ystride + 1;
        out->nmo = 0;
        for (i = 0; i < (m - 1); i += mstride) out->nmo++;
        out->nmo++;                                            /*the last size class is always output*/
        out->mi = (int *) R_alloc(out->nmo, sizeof(int));
        k = 0;
        for (i = 0; i < (m - 1); i += mstride) out->mi[k++] = i;
        out->mi[k] = m - 1;
        out->nto = (grid->toutmin == 0) ? 1 : 0;
        for (j = 1; j < t; j++) {
                if (grid->t_values[j] >= grid->toutmin && grid->t_values[j] <= grid->toutmax && j % ostep == 0) out->nto++;
        }
        out->it = 0;

        R_xlen_t ncell = (R_xlen_t)out->nxo * out->nyo;
        R_xlen_t nfull = ncell * out->nmo * out->nto;
        const char *names[] = {"t", "x", "y", "m", "plankton", "pelagic", "benthic", "detritus",
                               "pelagic_growth", "pelagic_mortality", "pelagic_predation", "pelagic_fishing",
                               "benthic_growth", "benthic_mortality", "benthic_predation", "benthic_fishing",
                               "pelagic_catch", "benthic_catch", ""};
        SEXP list = PROTECT(mkNamed(VECSXP, names));
        out->list = list;
        SEXP v;
        v = allocVector(REALSXP, out->nto); SET_VECTOR_ELT(list, 0, v); out->t = REAL(v);
        v = allocVector(REALSXP, out->nxo); SET_VECTOR_ELT(list, 1, v); out->x = REAL(v);
        v = allocVector(REALSXP, out->nyo); SET_VECTOR_ELT(list, 2, v); out->y = REAL(v);
        v = allocVector(REALSXP, out->nmo); SET_VECTOR_ELT(list, 3, v); out->m = REAL(v);
        for (k = 0; k < out->nxo; k++) out->x[k] = grid->x_values[k * out->xstride];
        for (k = 0; k < out->nyo; k++) out->y[k] = grid->y_values[k * out->ystride];
        for (k = 0; k < out->nmo; k++) out->m[k] = grid->m_values[out->mi[k]];
        for (k = 0; k < out->nto; k++) out->t[k] = NA_REAL;
        v = allocVector(REALSXP, nfull); SET_VECTOR_ELT(list, 4, v); out->pla = REAL(v);
        out->pel = species_arrays(list, 5, n, nfull);
        out->pel_g    = out->want_g    ? species_arrays(list, 8, n, nfull)  : NULL;
        out->pel_mu   = out->want_mu   ? species_arrays(list, 9, n, nfull)  : NULL;
        out->pel_fish = out->want_fish ? species_arrays(list, 11, n, nfull) : NULL;
        out->pel_catch = out->want_catch ? species_arrays(list, 16, n, nfull) : NULL;
        out->pel_pred = NULL;
        if (out->want_pred) {
                /*one list per prey species, holding one array per pelagic predator species*/
                SEXP pred = allocVector(VECSXP, n); SET_VECTOR_ELT(list, 10, pred);
                out->pel_pred = (double ***) R_alloc(n, sizeof(double **));
                for (s = 0; s < n; s++) out->pel_pred[s] = species_arrays(pred, s, n, nfull);
        }
        out->ben = NULL; out->det = NULL;
        out->ben_g = NULL; out->ben_mu = NULL; out->ben_pred = NULL; out->ben_fish = NULL; out->ben_catch = NULL;
        if (c != 0) {
                out->ben = species_arrays(list, 6, c, nfull);
                v = allocVector(REALSXP, ncell * out->nto); SET_VECTOR_ELT(list, 7, v); out->det = REAL(v);
                if (out->want_g)    out->ben_g    = species_arrays(list, 12, c, nfull);
                if (out->want_mu)   out->ben_mu   = species_arrays(list, 13, c, nfull);
                if (out->want_pred) out->ben_pred = species_arrays(list, 14, c, nfull);
                if (out->want_fish) out->ben_fish = species_arrays(list, 15, c, nfull);
                if (out->want_catch) out->ben_catch = species_arrays(list, 17, c, nfull);
        }
}

/*Copy the current state into the next output time slot of the arrays returned to R*/
void capture_output(RUN *run, GRID *grid, COMMUNITY *community, int time)
{
        OUTPUT_ARRAYS *o = run->out;
        int ii, kk, ll, s, p;
        if (o == NULL || o->mode == 0 || o->it >= o->nto) return;
        o->t[o->it] = grid->t_values[time];
        for (kk = 0; kk < o->nxo; kk++) {
                int k = kk * o->xstride;
                for (ll = 0; ll < o->nyo; ll++) {
                        int l = ll * o->ystride;
                        for (ii = 0; ii < o->nmo; ii++) {
                                int i = o->mi[ii];
                                R_xlen_t idx = (R_xlen_t)kk + (R_xlen_t)o->nxo * ((R_xlen_t)ll + (R_xlen_t)o->nyo * ((R_xlen_t)ii + (R_xlen_t)o->nmo * o->it));
                                o->pla[idx] = community->plankton->u_values[i][k][l];
                                for (s = 0; s < run->no_pelagic; s++) {
                                        PELAGIC *pe = &community->pelagic[s];
                                        o->pel[s][idx] = pe->u_values[i][k][l];
                                        if (o->pel_g != NULL)    o->pel_g[s][idx]    = pe->g_values[i][k][l];
                                        if (o->pel_mu != NULL)   o->pel_mu[s][idx]   = pe->mu_values[i][k][l];
                                        if (o->pel_fish != NULL) o->pel_fish[s][idx] = pe->mu_fish_values[i][k][l];
                                        if (o->pel_catch != NULL) o->pel_catch[s][idx] = pe->mu_fish_values[i][k][l] * pe->u_values[i][k][l] * grid->exp_m_values[i];
                                        if (o->pel_pred != NULL) {
                                                for (p = 0; p < run->no_pelagic; p++) o->pel_pred[s][p][idx] = pe->mu_pred_values[p][i][k][l];
                                        }
                                }
                                for (s = 0; s < run->no_benthic; s++) {
                                        BENTHIC *be = &community->benthic[s];
                                        o->ben[s][idx] = be->u_values[i][k][l];
                                        if (o->ben_g != NULL)    o->ben_g[s][idx]    = be->g_values[i][k][l];
                                        if (o->ben_mu != NULL)   o->ben_mu[s][idx]   = be->mu_values[i][k][l];
                                        if (o->ben_pred != NULL) o->ben_pred[s][idx] = be->mu_pred_values[i][k][l];
                                        if (o->ben_fish != NULL) o->ben_fish[s][idx] = be->mu_fish_values[i][k][l];
                                        if (o->ben_catch != NULL) o->ben_catch[s][idx] = be->mu_fish_values[i][k][l] * be->u_values[i][k][l] * grid->exp_m_values[i];
                                }
                        }
                        if (o->det != NULL) {
                                o->det[(R_xlen_t)kk + (R_xlen_t)o->nxo * ((R_xlen_t)ll + (R_xlen_t)o->nyo * o->it)] = community->detritus->w_values[k][l];
                        }
                }
        }
        o->it++;
}
