// movement.h and movement.c contains function relating to biological movement processes.

#ifndef MOVEMENT_H
#define MOVEMENT_H

#include <math.h>
#include "common.h"

double Cfun(double, double, GRID *, PELAGIC *);
double Dfun(double, double, GRID *, PELAGIC *);
double Diffun(double, double, GRID *, PELAGIC *);
int horizontal_face_valid(GRID *, int, int); // Identifies faces that are open for fluxes (both water cells)
int vertical_face_valid(GRID *, int, int); // " "

void xmove_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, double ****, int); // Calculates x movement in 1D space
void ymove_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, double ****, int); // Calculates y movement in 2D space

#endif