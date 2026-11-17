// memory_inputs.h and memory_inputs.c contain the functions that take model
// inputs (plankton series, forcing, forced abundance series) from arrays passed from R instead of text files.

#ifndef MEMORY_INPUTS_H
#define MEMORY_INPUTS_H

#include "common.h"

void plankton_from_arrays(GRID *grid, PLANKTON *plankton);
int forcing_from_arrays(RUN *run, GRID *grid, int k, int l, double *vals);
void series_check(const SPECIES_INPUT *in, const char *name, int ts_flag, int x, int y, int ns);
void series_to_u(const SPECIES_INPUT *in, int tp, int x, int y, int imin, int ns, double ***u);
void series_to_w(const SPECIES_INPUT *in, int tp, int x, int y, double **w);

#endif
