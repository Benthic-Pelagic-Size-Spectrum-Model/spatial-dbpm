// memory_outputs.h and memory_outputs.c contain functions relevant to setting up
// and writing memory outputs to return from C to R.

#ifndef MEMORY_OUTPUTS_H
#define MEMORY_OUTPUTS_H

#include <R.h>
#include <Rinternals.h>
#include "common.h"

typedef struct output_info
{
    int nout_m, nout_t, nout_x, nout_y; // Numbers of output based on model and output steps
    int *m_sel, *x_sel, *y_sel; // Holders of all indices for output indices (not every index on grid, just output ones)
    int next_time;              // The next timestep index (is iterated to)

    double *t_out; // Time values that are completed through the run

    double *pla_u; // Holder of plankton density output values (size of nout_m, nout_x, nout_y, nout_t)

    double *pel_u; // Holder of pelagic density output values (size of nout_m, nout_x, nout_y, n_pel_species, nout_t)
    double *pel_mu_pred; // Holder of pelagic mortality (from predation) values (size of nout_m, nout_x, nout_y,
                         // n_pel_species, nout_t)
    double *pel_mu_fish; // Holder of pelagic mortality (from fishing) values (size of nout_m, nout_x, nout_y,
                         // n_pel_species, nout_t)
    double *pel_g;       // Holder of pelagic growth rate values (size of nout_m, nout_x, nout_y, n_pel_species, nout_t)

    double *ben_u; // Holder of pelagic density output values (size of nout_m, nout_x, nout_y, n_ben_species, nout_t)
    double *ben_mu_pred; // Holder of benthic mortality (from predation) values (size of nout_m, nout_x, nout_y,
                         // n_ben_species, nout_t)
    double *ben_mu_fish; // Holder of benthic mortality (from fishing) values (size of nout_m, nout_x, nout_y,
                         // n_ben_species, nout_t)
    double *ben_g;       // Holder of pelagic growth rate values (size of nout_m, nout_x, nout_y, n_ben_species, nout_t)

    double *det_w;  // Holder of detritus mass density output values (size of nout_x, nout_y, nout_t)
    double *det_g;  // Holder of detritus gain values (detritus - in) (size of nout_x, nout_y, nout_t)
    double *det_mu; // Holder of detritus loss values (detritus - out) (size of nout_x, nout_y, nout_t)

} OUTPUT;

void fill_na(SEXP);
void setup_output(RUN *, GRID *, OUTPUT *, SEXP); // sets up the output containers for the run
int is_output_timestep(GRID *, int);
int output_interval(double, double);
void record_output(RUN *, GRID *, COMMUNITY *, OUTPUT *, int);

#endif