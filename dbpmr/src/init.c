#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Forward declaration of the simulation entry point in SizeSpectra.c */
extern SEXP SizeSpectrum(
        SEXP run_params_r, 
        SEXP grid_params_r, 
        SEXP plankton_params_r, 
        SEXP pelagic_params_r, 
        SEXP benthic_params_r, 
        SEXP detritus_params_r, 
        SEXP names_params_r, 
        SEXP flags_params_r
);

static const R_CallMethodDef CallEntries[] = {
    {"SizeSpectrum", (DL_FUNC) &SizeSpectrum, 8},
    {NULL, NULL, 0}
};

void R_init_dbpmr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
