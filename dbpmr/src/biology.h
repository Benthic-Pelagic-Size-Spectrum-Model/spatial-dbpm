// biology.h and biology.c files contain items relating to size spectrum model
// processes and base biological and ecological functions

#ifndef BIOLOGY_H
#define BIOLOGY_H

#include <math.h>
#include "common.h"

/*Predation, Growth and Renewal Functions*/
void calculate_g_and_mu(RUN *run, GRID *, COMMUNITY *); // Calculates g and mu values given u values

double pla_biomass(int, int, int, int, GRID *, COMMUNITY *); // Calculates the biomass eaten from plankton spectra
double pel_biomass(int, int, int, int, RUN *, GRID *, COMMUNITY *); // Calculates the biomass eaten from pelagic spectra
double ben_biomass(int, int, int, int, RUN *, GRID *, COMMUNITY *); // Calculates the biomass eaten from benthic spectra
double det_biomass(int, int, int, int, RUN *, GRID *, COMMUNITY *); // Calculates the biomass eaten from detritus system

double g_pel(int, int, int, int, RUN *, GRID *,
             COMMUNITY *); // Calculates g values for pelagic system from biomass values
double g_ben(int, int, int, int, RUN *, GRID *,
             COMMUNITY *);                          // Calculates g values for benthic system from biomass values
double g_det(int, int, RUN *, GRID *, COMMUNITY *); // Calculates g values for detrital system from biomass values

double mu_pel_pred(int, int, int, int, int, RUN *, GRID *, COMMUNITY *);
double mu_pel_fish(int, int, int, int, RUN *, GRID *, COMMUNITY *);
double mu_ben_pred(int, int, int, int, RUN *, GRID *, COMMUNITY *);
double mu_ben_fish(int, int, int, int, RUN *, GRID *, COMMUNITY *);

double mu_pel(int, int, int, int, RUN *, GRID *, COMMUNITY *); // Calculates mu values for pelagic system
double mu_ben(int, int, int, int, RUN *, GRID *, COMMUNITY *); // Calculates mu values for benthic system
double mu_det(int, int, RUN *, GRID *, COMMUNITY *); // Calculates mu values for detrital system from biomass values

void calculate_reproduction(RUN *, GRID *, COMMUNITY *); // Calculates reproduction for any system

#endif