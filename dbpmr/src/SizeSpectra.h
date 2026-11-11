// SizeSpectra.h and SizeSpectra.c contain the main model functions that step the
// size spectrum model through time. These files contain the function (SizeSpectrum)
// that is called through R using .Call()

#ifndef SIZESPECTRA_H
#define SIZESPECTRA_H

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "setup.h"
#include "memory_outputs.h"
#include "file_outputs.h"

#include "biology.h"
#include "movement.h"
#include "fishing.h"

/*Protype of the internal size-spectrum function (no interface with R)*/
static SEXP SizeSpectrum_internals(int *, double *, double *, double *, double *, double *, char **, int *, int *, int *, R_xlen_t);

/*Timestepping and Output Function*/
void calculate_results(RUN *, GRID *, COMMUNITY *, OUTPUT *, MATRIX *, MATRIX *, MATRIX *, MATRIX *);

/*differencing Scheme Solver Functions*/
void mass_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, MATRIX *, double ***, double ****, double ****,
                 double **);                                               // Calculates u values for next time step

#endif