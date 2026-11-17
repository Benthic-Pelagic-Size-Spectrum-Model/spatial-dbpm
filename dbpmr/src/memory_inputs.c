#include "memory_inputs.h"

/*Names of the forcing channels in the list passed from R, in the order of the FRC_ enum*/
const char *FRC_NAMES[N_FRC] = { "forcing_pel_tempeff", "forcing_ben_tempeff", "forcing_sinking_rate",
       "forcing_depth", "forcing_residence_time", "forcing_export_attn", "forcing_export_gamma",
       "forcing_export_magg", "forcing_gravity" };

/*Fill the plankton spectrum of every cell from the R arrays for the next time block, then advance the block.
  Past the last block the final one is reused. Size classes outside the plankton range stay at zero.*/
void plankton_from_arrays(GRID *grid, PLANKTON *plankton)
{
        const ARRAY_INPUTS *arr = plankton->arr;
        int i, k, l;
        int x = grid->xnum, y = grid->ynum;
        int block = plankton->ts_block < arr->pla_nblock ? plankton->ts_block : arr->pla_nblock - 1;
        R_xlen_t layer = arr->pla_layer[block];
        for (k = 0; k < x; k++) {
                for (l = 0; l < y; l++) {
                        R_xlen_t idx = (R_xlen_t)k + (R_xlen_t)x * ((R_xlen_t)l + (R_xlen_t)y * layer);
                        double a = arr->pla_int[idx], b = arr->pla_slp[idx];
                        for (i = plankton->iplamin; i < (plankton->iplamax + 1); i++) {
                                plankton->u_values[i][k][l] = a * exp(b * grid->m_values[i]);
                        }
                }
        }
        plankton->ts_block++;
}

/*Forcing values of cell (k,l) for the current timestep from the R arrays, written into vals[] at the
  channel positions of FRC_NAMES (absent channels are left untouched). Returns the number of columns.*/
int forcing_from_arrays(RUN *run, GRID *grid, int k, int l, double *vals)
{
        const ARRAY_INPUTS *arr = run->arr;
        int c;
        int x = grid->xnum, y = grid->ynum;
        int tp = run->frc_step + 1;
        if (tp > arr->frc_nblock - 1) tp = arr->frc_nblock - 1;
        R_xlen_t layer = arr->frc_layer[tp];
        for (c = 0; c < N_FRC; c++) {
                if (arr->frc[c] == NULL) continue;
                R_xlen_t lay = arr->frc_nlayer[c] == 1 ? 0 : layer;
                vals[c] = arr->frc[c][(R_xlen_t)k + (R_xlen_t)x * ((R_xlen_t)l + (R_xlen_t)y * lay)];
        }
        return N_FRC;
}

/*A forced series held in memory must belong to a group with ts_flag set and hold x*y*ns*layers values
  (ns = number of size classes of the group; 1 for detritus)*/
void series_check(const SPECIES_INPUT *in, const char *name, int ts_flag, int x, int y, int ns)
{
        if (ts_flag != 1) error("a forced series was supplied for %s but its ts_flag is not set\n", name);
        if (in->series_len != (R_xlen_t)x * (R_xlen_t)y * (R_xlen_t)ns * (R_xlen_t)in->series_nlayer) {
                error("the forced series of %s holds %ld values; expected x*y*size classes*layers = %d*%d*%d*%d\n",
                      name, (long)in->series_len, x, y, ns, in->series_nlayer);
        }
}

/*Abundance of one group from its forced series: the layer of model time point tp (past the last time point the
  final one is reused). The series is an R array [x][y][size class of the group][layer]; imin is the index of the
  group's smallest size class on the model grid and ns its number of size classes.*/
void series_to_u(const SPECIES_INPUT *in, int tp, int x, int y, int imin, int ns, double ***u)
{
        int i, k, l;
        R_xlen_t layer = in->series_layer[tp < in->series_nblock ? tp : in->series_nblock - 1];
        for (i = 0; i < ns; i++) {
                for (k = 0; k < x; k++) {
                        for (l = 0; l < y; l++) {
                                u[imin + i][k][l] = in->series[(R_xlen_t)k + (R_xlen_t)x * ((R_xlen_t)l + (R_xlen_t)y * ((R_xlen_t)i + (R_xlen_t)ns * layer))];
                        }
                }
        }
}

/*Detritus biomass from its forced series, an R array [x][y][layer]*/
void series_to_w(const SPECIES_INPUT *in, int tp, int x, int y, double **w)
{
        int k, l;
        R_xlen_t layer = in->series_layer[tp < in->series_nblock ? tp : in->series_nblock - 1];
        for (k = 0; k < x; k++) {
                for (l = 0; l < y; l++) {
                        w[k][l] = in->series[(R_xlen_t)k + (R_xlen_t)x * ((R_xlen_t)l + (R_xlen_t)y * layer)];
                }
        }
}
