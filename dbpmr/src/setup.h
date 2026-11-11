// setup.h and setup.c contains functions used to perform initial setup of model
// structs

#ifndef SETUP_H
#define SETUP_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>
#include "common.h"

//-----------------------//
// Function Declarations //
//-----------------------//

/*Setup Functions*/
void setup_run(RUN *, int *, char *, int *);
void setup_grid(GRID *, double *, int *, R_xlen_t);

void setup_plankton(RUN *, GRID *, PLANKTON *, double *, char *, int, int);         // sets up plankton values
void setup_pelagic(RUN *, GRID *, PELAGIC *, double *, char *, int, int, int, int); // sets up pelagic values
void setup_benthic(RUN *, GRID *, BENTHIC *, double *, char *, int, int, int, int); // sets up benthic values
void setup_detritus(RUN *, GRID *, DETRITUS *, double *, char *, int, int);         // sets up detritus values

double **setup_pelagic_params(int, double *); // Converts 1D vector of pelagic params to 2D array of pelagic params
double **setup_benthic_params(int, double *); // Converts 1D vector of benthic params to 2D array of benthic params
void setup_matrix(int, MATRIX *);
void set_land_to_zero(RUN *, GRID *, COMMUNITY *);

double x_start(double, double);
double y_start(double, double);

#endif