#ifndef FILE_OUTPUTS_H
#define FILE_OUTPUTS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "common.h"

void print_run(RUN *, FILE *);
void print_grid(GRID *, FILE *);
void print_plankton(PLANKTON *, FILE *);
void print_pelagic(PELAGIC *, FILE *);
void print_benthic(BENTHIC *, FILE *);
void print_detritus(DETRITUS *, FILE *);

void print_mass_header(GRID *, FILE *);
void print_timestep_plankton(GRID *, PLANKTON *, int);
void print_timestep_pelagic(RUN *, GRID *, PELAGIC *, int);
void print_timestep_benthic(GRID *, BENTHIC *, int);

void print_plankton_header(FILE *);
void print_plankton_summary(GRID *, PLANKTON *, int);
void print_pelagic_header(FILE *);
void print_pelagic_summary(GRID *, PELAGIC *, int);
void print_benthic_header(FILE *);
void print_benthic_summary(GRID *, BENTHIC *, int);
void print_detritus_header(FILE *);
void print_detritus_summary(GRID *, DETRITUS *, int);

void print_detailed_header(GRID *, FILE *);
void print_detailed_plankton(GRID *, PLANKTON *, int);
void print_detailed_pelagic(RUN *, GRID *, PELAGIC *, int);
void print_detailed_benthic(GRID *, BENTHIC *, int);

#endif