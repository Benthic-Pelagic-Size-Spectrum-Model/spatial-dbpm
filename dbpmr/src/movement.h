// movement.h and movement.c contains function relating to biological movement processes.

#ifndef MOVEMENT_H
#define MOVEMENT_H

#include "common.h"

/*Spatial Movement Functions*/
double Cfun(double, double, GRID *, PELAGIC *);
double Dfun(double, double, GRID *, PELAGIC *);
double Diffun(double, double, GRID *, PELAGIC *);
void xmove_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, double ****, int); //Calculates x movement in 1D space
void ymove_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, double ****, int); //Calculates y movement in 2D space

#endif
