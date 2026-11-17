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
