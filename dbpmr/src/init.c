#include <stdlib.h>
#include <R.h>
#include <R_ext/Rdynload.h>

/* Forward declaration of the simulation entry point in SizeSpectra.c */
extern void SizeSpectrum(int *run_params, double *grid_params, double *pla_params,
                         double *pel_params, double *ben_params, double *det_params,
                         char **names_params, int *flags_params);

static const R_CMethodDef CEntries[] = {
    {"SizeSpectrum", (DL_FUNC) &SizeSpectrum, 8},
    {NULL, NULL, 0}
};

void R_init_dbpmr(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
