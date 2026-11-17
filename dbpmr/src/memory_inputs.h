// memory_inputs.h and memory_inputs.c contain the functions that take model
// inputs (plankton series, forcing) from arrays passed from R instead of text files.

#ifndef MEMORY_INPUTS_H
#define MEMORY_INPUTS_H

#include "common.h"

void plankton_from_arrays(GRID *grid, PLANKTON *plankton);
int forcing_from_arrays(RUN *run, GRID *grid, int k, int l, double *vals);

#endif
