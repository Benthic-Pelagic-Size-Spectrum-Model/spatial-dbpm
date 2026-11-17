// SizeSpectra.h and SizeSpectra.c contain the main model functions that step the
// size spectrum model through time. These files contain the function (SizeSpectrum)
// that is called through R using .Call()

#ifndef SIZESPECTRA_H
#define SIZESPECTRA_H

#include "common.h"
#include "setup.h"
#include "memory_inputs.h"
#include "memory_outputs.h"
#include "file_outputs.h"

#include "biology.h"
#include "movement.h"
#include "fishing.h"

/*Protype of the internal size-spectrum function (no interface with R)*/
static void SizeSpectrum_internals(int *, double *, double *, double *, double *, double *, char **, int *, const ARRAY_INPUTS *, OUTPUT_ARRAYS *);

/*Timestepping and Output Function*/
void calculate_results(RUN *, GRID *, COMMUNITY *, MATRIX *, MATRIX*, MATRIX*, MATRIX*);

/*differencing Scheme Solver Functions*/
void mass_solver(RUN *, GRID *,COMMUNITY *, MATRIX *, MATRIX *, double ***, double ****, double ****, double **);   //Calculates u values for next time step

#endif
