#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//------------------------//
// Structure Declarations //
//------------------------//

typedef struct run_info
{

    char *filename; // root filename for runs
    int no_pelagic;
    int no_benthic;
    int spatial_dim;
    int coupled_flag; // flag for coupled info; 0-uncoupled, 1-coupled
    int diff_method;  // flag for integration method; 0-fast (poss. unstable), 1-slow (more stable)

    char fname_summ[40];
    FILE *fptr_summ;

    /* --- Environmental forcing plugin (auto-detected <run>/Input/forcing_ts.txt) ---
       A header line names the channels (any subset/order), one data row per timestep.
       Recognised channels: pel_tempeff, ben_tempeff (surface/seafloor Boltzmann-Arrhenius
       multipliers on feeding A + background mortality mu_0) and sinking_rate (fraction of
       surface-origin detritus reaching the seafloor). Add a forcing = add a column name +
       a col_* index + an apply hook. Missing channels fall back to the neutral defaults. */
    int forcing_flag;    // time-varying forcing driver active
    FILE *fptr_forcing;  // forcing time-series file ptr
    int forcing_ncol;    // number of columns in each data row
    int col_pel_tempeff; // column index of each channel (-1 if absent)
    int col_ben_tempeff;
    int col_sinking_rate;
    int col_depth;
    int col_residence_time;
    int col_export_attn;  // size-dependent detritus export (opt-in): particle of ln-mass m reaches
    int col_export_gamma; // the seafloor with fraction exp(-export_attn*exp(-export_gamma*max(m,export_magg)))
    int col_export_magg;  // aggregation floor: production below export_magg (ln-mass) sinks as marine snow of that
                          // effective size
    int col_gravity; // DBPM.md gravity fishing (opt-in): if present, calculate_fishing() scales F on each group by that
                     // group's CURRENT areal fishable-biomass share (F_g = base * B_g/sum B), recomputed per cell per
                     // timestep
    double pel_tempeff;      // current surface-temperature multiplier (default 1)
    double ben_tempeff;      // current seafloor-temperature multiplier (default 1)
    double sinking_rate;     // current detritus sinking/export fraction (default 1)
    int forcing_percell;     // per-cell forcing: forcing_ts.txt has x*y rows per timestep (env FORCING_PERCELL=1,
                             // spatial_dim>=1). 0 = one global row/timestep (default, bit-identical)
    int temp_on_mort;        // cached TEMP_ON_MORT (1 = temperature scales background mortality mu_0)
    double **pel_tempeff_c;  // per-cell pelagic temperature multiplier (used when forcing_percell; applied in
                             // feeding/mortality)
    double **ben_tempeff_c;  // per-cell benthic temperature multiplier
    double **sinking_rate_c; // per-cell detritus sinking/export fraction
    double **depth_c; // per-cell water-column depth (value used only in the gravity-fishing pelagic:benthic split)
    double **export_attn_c;  // per-cell size-dependent export attenuation number
    double **export_gamma_c; // per-cell size-dependent export sinking-velocity exponent
    double **export_magg_c;  // per-cell size-dependent export aggregation floor (ln-mass)
    double depth;            // water-column depth (m); read+available. Supplying the channel toggles Dunne burial ON
                  // (col_depth>=0); the value is NOT used to rescale the burial flux (sizemodel applies burial to the
                  // per-volume flux directly)
    double residence_time; // detritus residence time tau (yr) for the alternative first-order pool loss W/tau
                           // (dimensionally clean; supplied via the residence_time channel). 0/absent = off
    double export_attn;    // size-dependent export: attenuation number A = k(T)*z/w0 (per-timestep, from export_attn
                           // channel; larger = more remineralisation = less export)
    double export_gamma;   // size-dependent export: allometric sinking-velocity exponent (w_s ~ mass^gamma); larger
                           // particles sink faster -> higher export
    double
        export_magg; // aggregation floor (ln-mass): small production coagulates into marine snow of this effective size
                     // before sinking; export uses max(m, export_magg). Default -1e30 = no floor (pure per-particle)
    double gravity;  // gravity channel value (unused; presence via col_gravity toggles gravity fishing)

} RUN;

/*Contains all information regarding grid discretisation*/
typedef struct grid_info
{

    double mmin, mmax; // ln(mass)
    double mstep;      // step length
    int mnum;          // number of steps
    double moutstep;   // output to print every m
    double *m_values;  // actual log(m) values

    double t1, tmax; // time discretisation
    double tstep;
    int j1, tnum;
    double toutmin, toutmax;
    double toutstep;
    double *t_values;

    double xmin, xmax; // x-space discretisation (if required)
    double xstep;
    int xnum;
    double xoutstep;
    double *x_values;

    double ymin, ymax; // y-space discretisation (if required)
    double ystep;
    int ynum;
    double youtstep;
    double *y_values;

} GRID;

/*Contains all info regarding the plankton dynamics*/
typedef struct plankton_info
{

    char *filename;

    double plamin, plamax;
    int iplamin, iplamax;

    double mu_0;
    double beta;
    double u_0;
    double lambda;

    int initial_flag;
    int ts_flag;

    double ***u_values;
    double ***g_values; // for futureproofing

    // output files
    char fname_r[40];
    char fname_summ[40];

    FILE *fptr_r;
    FILE *fptr_summ;

    // input files
    char fname_ts[40];

    FILE *fptr_ts;

} PLANKTON;

/*Contains all info regarding an individual species' constants and abundances*/
typedef struct pelagic_info
{

    char *filename; // filename for each species

    double pelmin, pelmat, pelmax;
    int ipelmin, ipelmat, ipelmax;

    double A;       // volume search coefficient
    double alpha;   // volume search exponent
    double mu_0;    // intrinsic (juvenile) mortality coefficient
    double beta;    // intrinsic (juvenile) mortality exponent
    double mu_s;    // senesence mortality coefficient
    double epsilon; // senesence mortality constant

    double A_base;    // un-temperature-scaled search coefficient (temp driver)
    double mu_0_base; // un-temperature-scaled background mortality (temp driver)
    double tempeff;   // current-timestep temperature multiplier (1 if no temp driver)

    double u_0;    // initial intercept and abundance
    double lambda; // primary production slope

    double K_pla;  // growth efficiency (proportion of plankton food to growth)
    double R_pla;  // reproduction efficiency (proportion of plankton food to reproduction)
    double Ex_pla; // defecation efficiency (proportion of plankton food to defecation)
    double K_pel;  // growth efficiency (proportion of pelagic food to growth)
    double R_pel;  // reproduction efficiency (proportion of pelagic food to reproduction)
    double Ex_pel; // defecation efficiency (proportion of pelagic food to defecation)
    double K_ben;  // growth efficiency (proportion of benthic food to growth)
    double R_ben;  // reproduction efficiency (proportion of benthic food to reproduction)
    double Ex_ben; // defecation efficiency (proportion of benthic food to defecation)

    double pref_pla; // attack rate for plankton
    double pref_pel; // attack rate for pelagic
    double pref_ben; // attack rate for benthic

    double q_0;         // optimum predator prey mass ratio
    double sig;         // standard deviation of feeding kernel
    double trunc;       // number of standard deviations at which to truncate feeding kernel
    double *phi_values; // phi values

    double prey;       // prey seeking coefficient
    double pred;       // predator avoiding coefficient
    double comp;       // competition coefficient
    double gamma_prey; // prey seeking exponent
    double gamma_pred; // predator avoiding exponent
    double gamma_comp; // compeititon exponent

    int rep_method;  // reproduction method; 0-fixed amount , 1-biomass dependent
    double imm_frac; // connectivity/immigration recruitment floor as a fraction of unfished recruitment (0=off)
    double rec0; // running-peak (unfished) recruitment density at the smallest size (0-D / legacy scalar; spatial runs
                 // use rec0_c)
    int rec0_frozen; // freeze rec0 at first bind so the subsidy cannot inflate its own reference (0-D / legacy scalar)
    double **rec0_c; // PER-CELL running-peak recruitment (spatial-safe connectivity floor; each cell keeps its own
                     // reference so cells don't couple)
    int **rec0_frozen_c;   // PER-CELL freeze flag for the connectivity floor
    int refuge_on;         // reef prey-refuge mechanism active (refuge.txt supplied); 0 = off (a(w)=1, standard dbpmr)
    double refuge_min_a;   // availability floor min_A (REFUGE_MIN_A env, default 0); a(w) clamped to [min_a,1]
    double *refuge_values; // per prey-size refuge-space density (count m^-2); NULL/all-zero when off
    double ***a_refuge; // per (size,x,y) availability multiplier a(w)=max(min_a,1-refuge(w)/competitor(w)); recomputed
                        // each step
    int initial_flag;
    int ts_flag;
    int fishing_flag; // flag for fishing; 0-off, 1-on

    double ***u_values;        // number densities for each time step (outputted as file)
    double ***g_values;        // growth rates for each time step (outputted as file)
    double ***mu_values;       // mortality rates for each time step (outputted as file)
    double ****mu_pred_values; // mortality rates due to predation (outputted as file)
    double ***mu_fish_values;  // mortality rates due to fishing (outputted as file)

    double ***pla_bio; // plankton biomass eaten for each time step (used only in calculations)
    double ***pel_bio; // pelagic biomass eaten for each time step (used only in calculations)
    double ***ben_bio; // benthic biomass eaten for each time step (used only in calculations)

    double **pla_total; // total plankton biomass eaten at a spatial point (ouputted as rate in summary file)
    double **pel_total; // total pelagic biomass eaten at a spatial point (ouputted as rate in summary file)
    double **ben_total; // total benthic biomass eaten at a spatial point (ouputted as rate in summary file)

    double **pred_total;   // total biomass lost to predation at a spatial point (outputted as a rate in sumamry file)
    double **fish_total;   // total biomass lost to fishing at a spatial point (ouputted as rate in summary file)
    double **reproduction; // reproduction number density at a spatial point (ouputted as rate in summary file)

    // output files
    char fname_r[40];
    char fname_g[40];
    char fname_m[40];
    char **fname_pred;
    char fname_fish[40];
    char fname_summ[40];

    FILE *fptr_r;     // results output file ptr
    FILE *fptr_g;     // growth output file ptr
    FILE *fptr_m;     // mortality output file ptr
    FILE **fptr_pred; // predation mortality file ptr
    FILE *fptr_fish;  // fishing mortality file ptr
    FILE *fptr_summ;  // summary file ptr

    // input files
    char fname_ts[40];
    char fname_fish_ts[40];
    char fname_rep_ts[40];

    FILE *fptr_ts;      // time series file ptr
    FILE *fptr_fish_ts; // time series fishing file ptr
    FILE *fptr_rep_ts;  // time series eggs file ptr

} PELAGIC;

/*Contains all info regarding an individual species' constants and abundances*/
typedef struct benthic_info
{

    char *filename; // filename for each species

    double benmin, benmat, benmax;
    int ibenmin, ibenmat, ibenmax;

    double alpha;   // volume search exponent
    double A;       // volume search coefficient
    double beta;    // intrinsic (juvenile) mortality exponent
    double mu_0;    // intrinsic (juvenile) mortality coefficient
    double epsilon; // senesence mortality constant
    double mu_s;    // senesence mortality coefficient

    double A_base;    // un-temperature-scaled search coefficient (temp driver)
    double mu_0_base; // un-temperature-scaled background mortality (temp driver)
    double tempeff;   // current-timestep temperature multiplier (1 if no temp driver)

    double lambda; // primary production slope
    double u_0;    // initial intercept and abundance

    double K_det;  // growth efficiency (proportion of plankton food to growth)
    double R_det;  // reproduction efficiency (proportion of plankton food to reproduction)
    double Ex_det; // defecation efficiency (proportion of plankton food to defecation)

    double pref_det; // attack rate for detritus

    int rep_method;   // reproduction method; 0-fixed amount , 1-biomass dependent
    int initial_flag; // initial distribution supplied; 0-no , 1-yes
    int ts_flag;      // timeseries of distributions supplied; 0-no , 1-yes
    int fishing_flag; // flag for fishing; 0-off, 1-on

    double ***u_values;       // number densities for each time step (outputted as file)
    double ***g_values;       // growth rates for each time step (outputted as file)
    double ***mu_values;      // mortality rates for each time step (outputted as file)
    double ***mu_pred_values; // mortality rates due to predation (outputted as file)
    double ***mu_fish_values; // mortality rates due to fishing (outputted as file)

    double ***det_bio; // detritus biomass eaten for each time step (used in calculations)

    double **det_total; // total detritus biomass eaten at a spatial point (ouputted as rate in summary file)

    double **fish_total;   // total biomass lost to fishing at a spatial point (ouputted as rate in summary file)
    double **pred_total;   // total biomass lost to predation at a spatial point (ouputted as rate in summary file)
    double **reproduction; // reproduction number density at a spatial point (ouputted as rate in summary file)

    // output files
    char fname_r[40];
    char fname_g[40];
    char fname_m[40];
    char fname_pred[40];
    char fname_fish[40];
    char fname_summ[40];

    FILE *fptr_r;    // results output file ptr
    FILE *fptr_g;    // growth output file ptr
    FILE *fptr_m;    // mortality output file ptr
    FILE *fptr_pred; // predation mortality file ptr
    FILE *fptr_fish; // fishing mortality file ptr
    FILE *fptr_summ;

    // input files
    char fname_ts[40];
    char fname_fish_ts[40];
    char fname_rep_ts[40];

    FILE *fptr_ts;      // time series file ptr
    FILE *fptr_fish_ts; // time series fishing file ptr
    FILE *fptr_rep_ts;  // time series rep file ptr

} BENTHIC;

/*Contains all info regarding the plankton dynamics*/
typedef struct detritus_info
{

    char *filename;

    double w_0;

    int initial_flag;
    int ts_flag;

    double **w_values;
    double **g_values;
    double **mu_values;

    // output files
    char fname_summ[40];

    FILE *fptr_summ;

    // input files
    char fname_ts[40];

    FILE *fptr_ts; // time series file ptrr

} DETRITUS;

/*Contains all species data*/
typedef struct species_answer_info
{

    PLANKTON *plankton;
    PELAGIC *pelagic;
    BENTHIC *benthic;
    DETRITUS *detritus;

} COMMUNITY;

/*Structure for use in implicit upwind method*/
typedef struct vectors
{

    double *a;
    double *b;
    double *c;
    double *r;
    double *u;
    int size;

} MATRIX;

//-----------------------//
// Function Declarations //
//-----------------------//

/*Setup Functions*/
void setup_run(RUN *, int *, char *);
void setup_grid(GRID *, double *);
void setup_plankton(RUN *, GRID *, PLANKTON *, double *, char *, int, int);         // sets up plankton values
void setup_pelagic(RUN *, GRID *, PELAGIC *, double *, char *, int, int, int, int); // sets up pelagic values
void setup_benthic(RUN *, GRID *, BENTHIC *, double *, char *, int, int, int, int); // sets up benthic values
void setup_detritus(RUN *, GRID *, DETRITUS *, double *, char *, int, int);         // sets up detritus values

double phi(double, double, double, double); // Calculates feeding kernel for each pelagic species

double **setup_pelagic_params(int, double *); // Converts 1D vector of pelagic params to 2D array of pelagic params
double **setup_benthic_params(int, double *); // Converts 1D vector of benthic params to 2D array of benthic params
void setup_matrix(int, MATRIX *);

/*Timestepping and Output Function*/
void calculate_results(RUN *, GRID *, COMMUNITY *, MATRIX *, MATRIX *, MATRIX *, MATRIX *);

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

/*differencing Scheme Solver Functions*/
void mass_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, MATRIX *, double ***, double ****, double ****,
                 double **);                                               // Calculates u values for next time step
void xmove_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, double ****, int); // Calculates x movement in 1D space
void ymove_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, double ****, int); // Calculates y movement in 2D space
void tridag(MATRIX *);                                                     // Inverts a tridiagonal matrix
void trimul(MATRIX *);

/*Predation, Growth and Renewal Functions*/
void calculate_g_and_mu(RUN *run, GRID *, COMMUNITY *); // Calculates g and mu values given u values

double pla_biomass(int, int, int, int, GRID *, COMMUNITY *); // Calculates the biomass eaten from plankton spectra
double pel_biomass(int, int, int, int, RUN *, GRID *, COMMUNITY *); // Calculates the biomass eaten from pelagic spectra
double ben_biomass(int, int, int, int, RUN *, GRID *, COMMUNITY *); // Calculates the biomass eaten from benthic spectra
double det_biomass(int, int, int, int, RUN *, GRID *, COMMUNITY *); // Calculates the biomass eaten from detritus system

double g_pel(int, int, int, int, RUN *, GRID *,
             COMMUNITY *); // Calculates g values for pelagic system from biomass values
double g_ben(int, int, int, int, RUN *, GRID *,
             COMMUNITY *);                          // Calculates g values for benthic systsm from biomass values
double g_det(int, int, RUN *, GRID *, COMMUNITY *); // Calculates g values for detrital system from biomass values

double mu_pel_pred(int, int, int, int, int, RUN *, GRID *, COMMUNITY *);
double mu_pel_fish(int, int, int, int, RUN *, GRID *, COMMUNITY *);
double mu_ben_pred(int, int, int, int, RUN *, GRID *, COMMUNITY *);
double mu_ben_fish(int, int, int, int, RUN *, GRID *, COMMUNITY *);

double mu_pel(int, int, int, int, RUN *, GRID *, COMMUNITY *); // Calculates mu values for pelagic system
double mu_ben(int, int, int, int, RUN *, GRID *, COMMUNITY *); // Calculates mu values for benthic system
double mu_det(int, int, RUN *, GRID *, COMMUNITY *); // Calculates mu values for detrital system from biomass values

void calculate_reproduction(RUN *, GRID *, COMMUNITY *); // Calculates reproduction for any system
void calculate_fishing(RUN *, GRID *, COMMUNITY *);      // Calculates fishing mortality for any system

/*Spatial Movement Functions*/
double Cfun(double, double, GRID *, PELAGIC *);
double Dfun(double, double, GRID *, PELAGIC *);
double Diffun(double, double, GRID *, PELAGIC *);
double x_start(double, double);
double y_start(double, double);

/*Memory Management Functions*/
void free_mem(RUN *, GRID *, COMMUNITY *, MATRIX *, MATRIX *, MATRIX *, MATRIX *, double **,
              double **);                      // Frees all allocated memory
void *safe_malloc(size_t, int, char *);        // Wrappered function for memory alloaction
FILE *safe_fopen(char *, char *, int, char *); // Wrappered function for file opening

/*Protype of the internal size-spectrum function (no interface with R)*/
static void SizeSpectrum_internals(int *, double *, double *, double *, double *, double *, char **, int *, int *);

// This is the new R interface. SizeSpectrum() is now a wrapper function that
// takes the objects from R (.Call()). SizeSpectrum() checks the SEXP objects from
// R are the correct SXPtypes, sets C pointers to the R data, and then runs the
// internal function SizeSpectrum_internals() which performs the calculations and
// writes outputs. SizeSpectrum() returns NULL to R.
// TODO: return to R will be changed to be optional so that users can toggle
// file output, memory output or both.
SEXP SizeSpectrum(SEXP run_params_r, SEXP grid_params_r, SEXP plankton_params_r, SEXP pelagic_params_r,
                  SEXP benthic_params_r, SEXP detritus_params_r, SEXP names_params_r, SEXP flags_params_r,
                  SEXP output_params_r)
{
    // Ensure that all of the SEXP objects are of the correct type after
    // passing in from R.
    if (!isInteger(run_params_r) || !isInteger(flags_params_r) || !isInteger(output_params_r))
    {
        error("run, flags and output *_params_r R objects must be integer vectors");
    }
    if (!isReal(grid_params_r) || !isReal(plankton_params_r) || !isReal(pelagic_params_r) ||
        !isReal(benthic_params_r) || !isReal(detritus_params_r))
    {
        error("grid, plankton, pelagic, benthic and detritus params must be numeric vectors");
    }
    if (!isString(names_params_r))
    {
        error("names_params_r R object must be character");
    }

    // Checking the length of the parameter vectors that were passed from R.
    // The values here are taken from hardcoded indices throughout the code.
    if (LENGTH(run_params_r) != 5)
    {
        error("run_params_r must have 5 elements");
    }
    if (LENGTH(grid_params_r) != 18)
    {
        error("grid_params_r must have 18 elements");
    }
    if (LENGTH(plankton_params_r) != 6)
    {
        error("plankton_params_r must have 6 elements");
    }
    if (LENGTH(output_params_r) != 2)
    {
        error("output_params_r must have 2 elements");
    }

    int npel = INTEGER(run_params_r)[0];
    if (LENGTH(pelagic_params_r) != 32 * npel)
    {
        error("pelagic_params_r must have 32 * n_pel_species elements");
    }

    // It is possible to include no benthic species, in which case benthic_params_r
    // and detritus_params_r do not need to be checked for length.
    int nben = INTEGER(run_params_r)[1];
    if (nben > 0)
    {
        if (LENGTH(benthic_params_r) != 15 * nben)
        {
            error("benthic_params_r must have 15 * n_ben_species elements");
        }
        if (LENGTH(detritus_params_r) != 1)
        {
            error("detritus_params_r must have 1 element");
        }
        if (LENGTH(names_params_r) != (npel + nben + 1) + 2)
        {
            error("names_params_r not correct length");
        }
    }
    else if (LENGTH(names_params_r) != npel + 2)
    {
        error("names_params_r not correct length");
    }

    int nflags = 2 + 4 * (npel + nben) + (nben > 0 ? 2 : 0);
    // nflags is 2 for plankton, 4 for each pelagic / benthic species and 2
    // for detritus if there is at least 1 benthic species
    if (LENGTH(flags_params_r) != nflags)
    {
        error("flags_params_r not correct length");
    }

    // Assign pointers for the SEXP data to correct names for function use
    int *run_params = INTEGER(run_params_r);
    int *flags_params = INTEGER(flags_params_r);
    int *output_params = INTEGER(output_params_r);
    double *grid_params = REAL(grid_params_r);
    double *pla_params = REAL(plankton_params_r);
    double *pel_params = REAL(pelagic_params_r);
    double *ben_params = REAL(benthic_params_r);
    double *det_params = REAL(detritus_params_r);

    char **names_params = (char **)R_alloc(LENGTH(names_params_r), sizeof(char *));
    int name_idx;
    for (name_idx = 0; name_idx < LENGTH(names_params_r); name_idx++)
    {
        names_params[name_idx] = (char *)CHAR(STRING_ELT(names_params_r, name_idx));
    }

    SizeSpectrum_internals(run_params, grid_params, pla_params, pel_params, ben_params, det_params, names_params,
                           flags_params, output_params);
    return (R_NilValue);
}

//-----------------------------------------------------------------------------//
// This is the internal function that does the work of the program called by R //
//-----------------------------------------------------------------------------//

static void SizeSpectrum_internals(int *run_params, double *grid_params, double *pla_params, double *pel_params,
                                   double *ben_params, double *det_params, char **names_params, int *flags_params,
                                   int *output_params)
// run_params is an array of run parameters
// grid_params is an array of grid discretisation values
// pla_params is plankton params
// pel_params is pelagic params
// ben_params is benthic parameters
// det_params is detritus parameters
// names_params is a list of filenames
// flags_params is a list of selection 'flag' parameters
{
    int pel_species_idx, ben_species_idx;
    RUN run;
    GRID grid;
    COMMUNITY community;
    MATRIX *pelmatrix;
    MATRIX *benmatrix = NULL;
    MATRIX xmatrix, ymatrix;
    double **temp_pel_params;
    double **temp_ben_params;

    /*Setup run*/
    setup_run(&run, run_params, names_params[0]);

    /*Setup grid*/
    setup_grid(&grid, grid_params);

    /*Setup community*/
    community.plankton = (PLANKTON *)safe_malloc(sizeof(PLANKTON), __LINE__, __FILE__);
    community.pelagic = (PELAGIC *)safe_malloc(run.no_pelagic * sizeof(PELAGIC), __LINE__, __FILE__);
    if (run.no_benthic != 0)
    {
        community.benthic = (BENTHIC *)safe_malloc(run.no_benthic * sizeof(BENTHIC), __LINE__, __FILE__);
        community.detritus = (DETRITUS *)safe_malloc(sizeof(DETRITUS), __LINE__, __FILE__);
    }

    /*Setup plankton*/
    setup_plankton(&run, &grid, community.plankton, pla_params, names_params[1], flags_params[0], flags_params[1]);

    /*Setup pelagic*/
    temp_pel_params = setup_pelagic_params(run.no_pelagic, pel_params);
    for (pel_species_idx = 0; pel_species_idx < run.no_pelagic; pel_species_idx++)
    {
        setup_pelagic(&run, &grid, &(community.pelagic[pel_species_idx]), temp_pel_params[pel_species_idx],
                      names_params[pel_species_idx + 2], flags_params[4 * pel_species_idx + 2],
                      flags_params[4 * pel_species_idx + 3], flags_params[4 * pel_species_idx + 4],
                      flags_params[4 * pel_species_idx + 5]);
    }

    if (run.no_benthic != 0)
    {
        /*Setup benthic*/
        temp_ben_params = setup_benthic_params(run.no_benthic, ben_params);
        for (ben_species_idx = 0; ben_species_idx < run.no_benthic; ben_species_idx++)
        {
            setup_benthic(&run, &grid, &(community.benthic[ben_species_idx]), temp_ben_params[ben_species_idx],
                          names_params[ben_species_idx + run.no_pelagic + 2],
                          flags_params[4 * ben_species_idx + (2 + (4 * run.no_pelagic))],
                          flags_params[4 * ben_species_idx + (3 + (4 * run.no_pelagic))],
                          flags_params[4 * ben_species_idx + (4 + (4 * run.no_pelagic))],
                          flags_params[4 * ben_species_idx + (5 + (4 * run.no_pelagic))]);
        }

        /*Setup detritus*/
        setup_detritus(&run, &grid, community.detritus, det_params, names_params[run.no_pelagic + run.no_benthic + 2],
                       flags_params[2 + 4 * (run.no_pelagic + run.no_benthic)],
                       flags_params[3 + 4 * (run.no_pelagic + run.no_benthic)]);
    }

    /*Setup differencing matrices*/
    pelmatrix = (MATRIX *)safe_malloc(run.no_pelagic * sizeof(MATRIX), __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < run.no_pelagic; pel_species_idx++)
    {
        setup_matrix((community.pelagic[pel_species_idx].ipelmax - community.pelagic[pel_species_idx].ipelmin + 1),
                     &(pelmatrix[pel_species_idx]));
    }

    if (run.no_benthic != 0)
    {
        benmatrix = (MATRIX *)safe_malloc(run.no_benthic * sizeof(MATRIX), __LINE__, __FILE__);
        for (ben_species_idx = 0; ben_species_idx < run.no_benthic; ben_species_idx++)
        {
            setup_matrix((community.benthic[ben_species_idx].ibenmax - community.benthic[ben_species_idx].ibenmin + 1),
                         &(benmatrix[ben_species_idx]));
        }
    }

    setup_matrix(grid.xnum, &xmatrix);
    setup_matrix(grid.ynum, &ymatrix);

    /*Perform all finite differencing*/
    calculate_results(&run, &grid, &community, pelmatrix, benmatrix, &xmatrix, &ymatrix);

    /*Free all allocated memory*/
    free_mem(&run, &grid, &community, pelmatrix, benmatrix, &xmatrix, &ymatrix, temp_pel_params, temp_ben_params);
}

//-----------------//
// Setup Functions //
//-----------------//

void setup_run(RUN *run, int *run_params, char *filename)
{
    run->filename = filename;
    run->no_pelagic = run_params[0];
    run->no_benthic = run_params[1];
    run->spatial_dim = run_params[2];
    run->coupled_flag = run_params[3];
    run->diff_method = run_params[4];

    /*Create summary output filenames*/
    snprintf(run->fname_summ, sizeof(run->fname_summ), "%s/%s.txt", run->filename, "parameters");

    /*Open output files*/
    run->fptr_summ = safe_fopen(run->fname_summ, "w", __LINE__, __FILE__);

    fprintf(run->fptr_summ, "filetype:      %s\n", "parameters");
    fprintf(run->fptr_summ, "speciestype:   %s\n", "NA");
    fprintf(run->fptr_summ, "\n");

    /*Close ouput files*/
    fclose(run->fptr_summ);
}

void setup_grid(GRID *grid, double *grid_input)
{
    /*Declare Mass Discretisation*/
    grid->mmin = grid_input[0];
    grid->mmax = grid_input[1];
    grid->mstep = grid_input[2];
    // m step length
    grid->mnum = (int)rint(((grid->mmax - grid->mmin) / grid->mstep) + 1); // number of m values
    grid->moutstep = grid_input[3];

    /*Declare Time Discretisation*/
    grid->t1 = grid_input[4];
    grid->tmax = grid_input[5];                             // maximum t value
    grid->tstep = grid_input[6];                            // t step length
    grid->tnum = (int)rint((grid->tmax / grid->tstep) + 1); // number of t values
    grid->j1 = (int)rint(grid->t1 / grid->tstep);
    grid->toutmin = grid_input[7];
    grid->toutmax = grid_input[8];
    grid->toutstep = grid_input[9];

    /*Declare X-Space Discretisation*/
    grid->xmin = grid_input[10];
    grid->xmax = grid_input[11];
    grid->xstep = grid_input[12];
    grid->xnum = (int)rint(((grid->xmax - grid->xmin) / grid->xstep) + 1);
    grid->xoutstep = grid_input[13];

    /*Declare Y-Space Discretisation*/
    grid->ymin = grid_input[14];
    grid->ymax = grid_input[15];
    grid->ystep = grid_input[16];
    grid->ynum = (int)rint(((grid->ymax - grid->ymin) / grid->ystep) + 1);
    grid->youtstep = grid_input[17];

    /*Create arrays*/
    int sc_idx, ts_idx, x_idx, y_idx;
    int n_m, n_t, n_x, n_y;

    n_m = grid->mnum;
    n_t = grid->tnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    grid->m_values = (double *)safe_malloc(n_m * sizeof(double), __LINE__, __FILE__);
    grid->t_values = (double *)safe_malloc(n_t * sizeof(double), __LINE__, __FILE__);
    grid->x_values = (double *)safe_malloc(n_x * sizeof(double), __LINE__, __FILE__);
    grid->y_values = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);

    /*Setup values for m*/
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        grid->m_values[sc_idx] = grid->mmin + (sc_idx * grid->mstep);
    }
    /*Setup values for t*/
    for (ts_idx = 0; ts_idx < n_t; ts_idx++)
    {
        grid->t_values[ts_idx] = (ts_idx * grid->tstep);
    }
    /*Setup values for x*/
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        grid->x_values[x_idx] = grid->xmin + (x_idx * grid->xstep);
    }
    /*Setup values for y*/
    for (y_idx = 0; y_idx < n_y; y_idx++)
    {
        grid->y_values[y_idx] = grid->ymin + (y_idx * grid->ystep);
    }
}

double **setup_pelagic_params(int no_pelagic, double *pel_params)
{
    int param_idx, pel_species_idx;
    int no_params = 32;
    double **temp_params;
    temp_params = (double **)safe_malloc(no_pelagic * sizeof(double *), __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < no_pelagic; pel_species_idx++)
    {
        temp_params[pel_species_idx] = (double *)safe_malloc(no_params * sizeof(double), __LINE__, __FILE__);
        for (param_idx = 0; param_idx < no_params; param_idx++)
        {
            temp_params[pel_species_idx][param_idx] = pel_params[no_params * pel_species_idx + param_idx];
        }
    }
    return (temp_params);
}

double **setup_benthic_params(int no_benthic, double *ben_params)
{
    int param_idx, ben_species_idx;
    int no_params = 15;
    double **temp_params;
    temp_params = (double **)safe_malloc(no_benthic * sizeof(double *), __LINE__, __FILE__);
    for (ben_species_idx = 0; ben_species_idx < no_benthic; ben_species_idx++)
    {
        temp_params[ben_species_idx] = (double *)safe_malloc(no_params * sizeof(double), __LINE__, __FILE__);
        for (param_idx = 0; param_idx < no_params; param_idx++)
        {
            temp_params[ben_species_idx][param_idx] = ben_params[no_params * ben_species_idx + param_idx];
        }
    }
    return (temp_params);
}

void setup_plankton(RUN *run, GRID *grid, PLANKTON *plankton, double *pla_params, char *filename, int initial_flag,
                    int ts_flag)
{
    /*Declare plankton name*/
    plankton->filename = filename;

    /*Setup Plankton Sizes*/
    plankton->plamin = pla_params[0];
    plankton->plamax = pla_params[1];
    plankton->iplamin = (int)rint((plankton->plamin - grid->mmin) / grid->mstep);
    plankton->iplamax = (int)rint((plankton->plamax - grid->mmin) / grid->mstep);

    /*Setup Plankton Constants*/
    plankton->mu_0 = pla_params[2];
    plankton->beta = pla_params[3];
    plankton->u_0 = pla_params[4];
    plankton->lambda = pla_params[5];

    /*Setup Plankton Flags*/
    plankton->initial_flag = initial_flag;
    plankton->ts_flag = ts_flag;

    /*Create plankton output filenames*/
    snprintf(plankton->fname_r, sizeof(plankton->fname_r), "%s/%s/%s.txt", run->filename, plankton->filename,
             "results");
    snprintf(plankton->fname_summ, sizeof(plankton->fname_summ), "%s/%s/%s.txt", run->filename, plankton->filename,
             "summary");

    /*Open output files*/
    plankton->fptr_r = safe_fopen(plankton->fname_r, "w", __LINE__, __FILE__);
    plankton->fptr_summ = safe_fopen(plankton->fname_summ, "w", __LINE__, __FILE__);

    /*Initialise output files*/
    fprintf(plankton->fptr_r, "filetype:      %s\n", "results");
    fprintf(plankton->fptr_r, "speciestype:   %s\n", "plankton");
    fprintf(plankton->fptr_r, "\n");
    print_mass_header(grid, plankton->fptr_r);

    fprintf(plankton->fptr_summ, "filetype:      %s\n", "summary");
    fprintf(plankton->fptr_summ, "speciestype:   %s\n", "plankton");
    fprintf(plankton->fptr_summ, "\n");
    print_plankton_header(plankton->fptr_summ);

    /*Close ouput files*/
    fclose(plankton->fptr_r);
    fclose(plankton->fptr_summ);

    /*Setup And Initialise Plankton Abundance Arrays*/
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    plankton->u_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    plankton->g_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        plankton->u_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        plankton->g_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            plankton->u_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            plankton->g_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                plankton->u_values[sc_idx][x_idx][y_idx] = 0;
                plankton->g_values[sc_idx][x_idx][y_idx] = 0;
            }
        }
    }

    /*Setup initial values for u(m,0,x,y)*/
    /*Use user defined values*/
    if (plankton->initial_flag == 1 || plankton->ts_flag == 1)
    {
        /*Create plankton input filename*/
        snprintf(plankton->fname_ts, sizeof(plankton->fname_ts), "%s/Input/%s_%s.txt", run->filename,
                 plankton->filename, "ts");
        plankton->fptr_ts = safe_fopen(plankton->fname_ts, "r", __LINE__, __FILE__);

        int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
        char *input;                   // input string to store line as string
        input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
        char sep[] = ",";    // identifying the ',' char as seperating values
        char *result = NULL; // declaring a pointer for storing each individual value as a char

        /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                fgets(input, MAX_LEN, plankton->fptr_ts);
                input[strlen(input) - 1] =
                    '\0'; // replacing the newline \n char of the line of text with a null \0 char

                result = (char *)strtok(input, sep); // this is the 1st value
                sc_idx = 0;
                while (result != NULL)
                {
                    plankton->u_values[sc_idx][x_idx][y_idx] = atof(result);
                    result = (char *)strtok(NULL, sep);
                    sc_idx++;
                }
            }
        }
        free(input);

        /*If this is only an intial value then close the file now*/
        if (plankton->ts_flag == 0)
        {
            fclose(plankton->fptr_ts);
        }
    }
    /*Use standard values*/
    else
    {
        for (sc_idx = 0; sc_idx < (plankton->iplamax + 1); sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    plankton->u_values[sc_idx][x_idx][y_idx] =
                        plankton->u_0 * exp(plankton->lambda * grid->m_values[sc_idx]);
                    if (run->spatial_dim == 1 || run->spatial_dim == 2)
                    {
                        plankton->u_values[sc_idx][x_idx][y_idx] *= x_start(grid->x_values[x_idx], 0);
                        if (run->spatial_dim == 2)
                        {
                            plankton->u_values[sc_idx][x_idx][y_idx] *= y_start(grid->y_values[y_idx], 0);
                        }
                    }
                }
            }
        }
    }
}

void setup_pelagic(RUN *run, GRID *grid, PELAGIC *pelagic, double *pel_params, char *filename, int rep_method,
                   int initial_flag, int ts_flag, int fishing_flag)
{
    /*Pelagic Species Filename*/
    pelagic->filename = filename;

    /*Pelagic Species Size Range*/
    pelagic->pelmin = pel_params[0];
    pelagic->pelmat = pel_params[1];
    pelagic->pelmax = pel_params[2];
    pelagic->ipelmin = (int)rint((pelagic->pelmin - grid->mmin) / grid->mstep);
    pelagic->ipelmat = (int)rint((pelagic->pelmat - grid->mmin) / grid->mstep);
    pelagic->ipelmax = (int)rint((pelagic->pelmax - grid->mmin) / grid->mstep);

    /*Pelagic Species Constants*/
    pelagic->A = pel_params[3];         // Volume searched by unit weight
    pelagic->alpha = pel_params[4];     // Volume search allometry
    pelagic->mu_0 = pel_params[5];      // natural mortality rate
    pelagic->beta = pel_params[6];      // natural mortality allometry
    pelagic->A_base = pelagic->A;       // temp driver: base (unscaled) search rate
    pelagic->mu_0_base = pelagic->mu_0; // temp driver: base (unscaled) background mortality
    pelagic->tempeff = 1.0;             // temp driver: neutral until a temperature series is read
    pelagic->mu_s = pel_params[7];      // senescence mortality rate
    pelagic->epsilon = pel_params[8];   // senescece mortality constant
    pelagic->u_0 = pel_params[9];       // initial intercept
    pelagic->lambda = pel_params[10];   // primary production slope

    /*Pelagic Species Conversion Efficiencies*/
    pelagic->K_pla = pel_params[11]; // growth efficiency
    pelagic->R_pla = pel_params[12];
    pelagic->Ex_pla = pel_params[13];
    pelagic->K_pel = pel_params[14];
    pelagic->R_pel = pel_params[15];
    pelagic->Ex_pel = pel_params[16];
    pelagic->K_ben = pel_params[17];
    pelagic->R_ben = pel_params[18];
    pelagic->Ex_ben = pel_params[19];

    /*Pelagic Species Attack Rates*/
    pelagic->pref_pla = pel_params[20];
    pelagic->pref_pel = pel_params[21];
    pelagic->pref_ben = pel_params[22];

    /*Pelagic Species Feeding Kernel*/
    pelagic->q_0 = pel_params[23]; // modal ratio of predation sizes
    pelagic->sig = pel_params[24]; // predator-prey speciality
    pelagic->trunc = pel_params[25];

    /*Pelagic Species Movement Constants*/
    pelagic->prey = 0;
    pelagic->pred = 0;
    pelagic->comp = 0;
    pelagic->gamma_prey = 0;
    pelagic->gamma_pred = 0;
    pelagic->gamma_comp = 0;

    if (run->spatial_dim == 1 || run->spatial_dim == 2)
    {
        pelagic->prey = pel_params[26];
        pelagic->pred = pel_params[27];
        pelagic->comp = pel_params[28];
        pelagic->gamma_prey = pel_params[29];
        pelagic->gamma_pred = pel_params[30];
        pelagic->gamma_comp = pel_params[31];
    }

    /*Pelagic Species Flags and Methods*/
    pelagic->rep_method = rep_method;
    {
        char *e = getenv("PEL_IMM_FRAC");
        pelagic->imm_frac = e ? atof(e) : 0.0;
        pelagic->rec0 = 0.0;
        pelagic->rec0_frozen = 0;
    }
    pelagic->initial_flag = initial_flag;
    pelagic->ts_flag = ts_flag;
    pelagic->fishing_flag = fishing_flag;

    int pel_species_idx, n_pel_species;
    n_pel_species = run->no_pelagic;

    /*Create species output filenames*/
    snprintf(pelagic->fname_r, sizeof(pelagic->fname_r), "%s/%s/%s.txt", run->filename, pelagic->filename, "results");
    snprintf(pelagic->fname_g, sizeof(pelagic->fname_g), "%s/%s/%s.txt", run->filename, pelagic->filename, "growth");
    snprintf(pelagic->fname_m, sizeof(pelagic->fname_m), "%s/%s/%s.txt", run->filename, pelagic->filename, "mortality");
    pelagic->fname_pred = (char **)safe_malloc(n_pel_species * sizeof(char *), __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        pelagic->fname_pred[pel_species_idx] = (char *)safe_malloc(40 * sizeof(char), __LINE__, __FILE__);
        snprintf(pelagic->fname_pred[pel_species_idx], 40, "%s/%s/%s%d.txt", run->filename, pelagic->filename,
                 "predation", pel_species_idx);
    }
    snprintf(pelagic->fname_fish, sizeof(pelagic->fname_fish), "%s/%s/%s.txt", run->filename, pelagic->filename,
             "fishing");
    snprintf(pelagic->fname_summ, sizeof(pelagic->fname_summ), "%s/%s/%s.txt", run->filename, pelagic->filename,
             "summary");

    /*Open output files*/
    pelagic->fptr_r = safe_fopen(pelagic->fname_r, "w", __LINE__, __FILE__);
    pelagic->fptr_g = safe_fopen(pelagic->fname_g, "w", __LINE__, __FILE__);
    pelagic->fptr_m = safe_fopen(pelagic->fname_m, "w", __LINE__, __FILE__);
    pelagic->fptr_pred = (FILE **)safe_malloc(n_pel_species * sizeof(FILE *), __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        pelagic->fptr_pred[pel_species_idx] = safe_fopen(pelagic->fname_pred[pel_species_idx], "w", __LINE__, __FILE__);
    }
    pelagic->fptr_fish = safe_fopen(pelagic->fname_fish, "w", __LINE__, __FILE__);
    pelagic->fptr_summ = safe_fopen(pelagic->fname_summ, "w", __LINE__, __FILE__);

    /*Initialise output files*/
    fprintf(pelagic->fptr_r, "filetype:      %s\n", "results");
    fprintf(pelagic->fptr_r, "speciestype:   %s\n", "pelagic");
    fprintf(pelagic->fptr_r, "\n");
    print_mass_header(grid, pelagic->fptr_r);

    fprintf(pelagic->fptr_g, "filetype:      %s\n", "growth");
    fprintf(pelagic->fptr_g, "speciestype:   %s\n", "pelagic");
    fprintf(pelagic->fptr_g, "\n");
    print_mass_header(grid, pelagic->fptr_g);

    fprintf(pelagic->fptr_m, "filetype:      %s\n", "mortality");
    fprintf(pelagic->fptr_m, "speciestype:   %s\n", "pelagic");
    fprintf(pelagic->fptr_m, "\n");
    print_mass_header(grid, pelagic->fptr_m);

    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        fprintf(pelagic->fptr_pred[pel_species_idx], "filetype:      %s\n", "predation");
        fprintf(pelagic->fptr_pred[pel_species_idx], "speciestype:   %s\n", "pelagic");
        fprintf(pelagic->fptr_pred[pel_species_idx], "\n");
        print_mass_header(grid, pelagic->fptr_pred[pel_species_idx]);
    }

    fprintf(pelagic->fptr_fish, "filetype:      %s\n", "fishing");
    fprintf(pelagic->fptr_fish, "speciestype:   %s\n", "pelagic");
    fprintf(pelagic->fptr_fish, "\n");
    print_mass_header(grid, pelagic->fptr_fish);

    fprintf(pelagic->fptr_summ, "filetype:      %s\n", "summary");
    fprintf(pelagic->fptr_summ, "speciestype:   %s\n", "pelagic");
    fprintf(pelagic->fptr_summ, "\n");
    print_pelagic_header(pelagic->fptr_summ);

    /*Close output files*/
    fclose(pelagic->fptr_r);
    fclose(pelagic->fptr_g);
    fclose(pelagic->fptr_m);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        fclose(pelagic->fptr_pred[pel_species_idx]);
    }
    fclose(pelagic->fptr_fish);
    fclose(pelagic->fptr_summ);

    /*Setup Species Abundance Arrays*/
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Setup 3D size-space arrays*/
    pelagic->mu_pred_values = (double ****)safe_malloc(n_pel_species * sizeof(double ***), __LINE__, __FILE__);
    pelagic->u_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    pelagic->g_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    pelagic->mu_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        pelagic->mu_pred_values[pel_species_idx] = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    }
    pelagic->mu_fish_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);

    pelagic->pla_bio = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    pelagic->pel_bio = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    pelagic->ben_bio = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    pelagic->a_refuge =
        (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__); /*reef refuge availability a(w)*/

    pelagic->phi_values = (double *)safe_malloc(n_m * sizeof(double), __LINE__, __FILE__);
    pelagic->refuge_values =
        (double *)safe_malloc(n_m * sizeof(double), __LINE__, __FILE__); /*per prey-size refuge density*/
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        pelagic->u_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        pelagic->g_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        pelagic->mu_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            pelagic->mu_pred_values[pel_species_idx][sc_idx] =
                (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        }
        pelagic->mu_fish_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

        pelagic->pla_bio[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        pelagic->pel_bio[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        pelagic->ben_bio[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        pelagic->a_refuge[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

        pelagic->phi_values[sc_idx] = 0;
        pelagic->refuge_values[sc_idx] = 0;
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            pelagic->u_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            pelagic->g_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            pelagic->mu_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                pelagic->mu_pred_values[pel_species_idx][sc_idx][x_idx] =
                    (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            }
            pelagic->mu_fish_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);

            pelagic->pla_bio[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            pelagic->pel_bio[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            pelagic->ben_bio[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            pelagic->a_refuge[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                pelagic->u_values[sc_idx][x_idx][y_idx] = 0;
                pelagic->g_values[sc_idx][x_idx][y_idx] = 0;
                pelagic->mu_values[sc_idx][x_idx][y_idx] = 0;
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    pelagic->mu_pred_values[pel_species_idx][sc_idx][x_idx][y_idx] = 0;
                }
                pelagic->mu_fish_values[sc_idx][x_idx][y_idx] = 0;

                pelagic->pla_bio[sc_idx][x_idx][y_idx] = 0;
                pelagic->pel_bio[sc_idx][x_idx][y_idx] = 0;
                pelagic->ben_bio[sc_idx][x_idx][y_idx] = 0;
                pelagic->a_refuge[sc_idx][x_idx][y_idx] = 1.0; /*default fully exposed (a=1) -> standard dbpmr*/
            }
        }
    }

    /*Setup 2D space arrays*/
    pelagic->pla_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->pel_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->ben_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

    pelagic->fish_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->pred_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->reproduction = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->rec0_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->rec0_frozen_c = (int **)safe_malloc(n_x * sizeof(int *), __LINE__, __FILE__);
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        pelagic->pla_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->pel_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->ben_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);

        pelagic->fish_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->pred_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->reproduction[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->rec0_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->rec0_frozen_c[x_idx] = (int *)safe_malloc(n_y * sizeof(int), __LINE__, __FILE__);
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            pelagic->pla_total[x_idx][y_idx] = 0;
            pelagic->pel_total[x_idx][y_idx] = 0;
            pelagic->ben_total[x_idx][y_idx] = 0;

            pelagic->fish_total[x_idx][y_idx] = 0;
            pelagic->pred_total[x_idx][y_idx] = 0;
            pelagic->reproduction[x_idx][y_idx] = 0;
            pelagic->rec0_c[x_idx][y_idx] = 0.0;
            pelagic->rec0_frozen_c[x_idx][y_idx] = 0;
        }
    }

    /*Setup initial values for u(m,0,x,y)*/
    if (pelagic->initial_flag == 1)
    {
        int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
        char *input;                   // input string to store line as string
        input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
        char sep[] = ",";    // identifying the ',' char as seperating values
        char *result = NULL; // declaring a pointer for storing each individual value as a char

        /*Create plankton input filename*/
        snprintf(pelagic->fname_ts, sizeof(pelagic->fname_ts), "%s/Input/%s_%s.txt", run->filename, pelagic->filename,
                 "ts");
        pelagic->fptr_ts = safe_fopen(pelagic->fname_ts, "r", __LINE__, __FILE__);

        /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                fgets(input, MAX_LEN, pelagic->fptr_ts);
                input[strlen(input) - 1] =
                    '\0'; // replacing the newline \n char of the line of text with a null \0 char

                result = (char *)strtok(input, sep); // this is the 1st value
                sc_idx = 0;
                while (result != NULL)
                {
                    pelagic->u_values[sc_idx][x_idx][y_idx] = atof(result);
                    result = (char *)strtok(NULL, sep);
                    sc_idx++;
                }
            }
        }
        free(input);
        /*If this is only an intial value then close the file now*/
        if (pelagic->ts_flag == 0)
        {
            fclose(pelagic->fptr_ts);
        }
    }
    else
    {
        for (sc_idx = pelagic->ipelmin; sc_idx < (pelagic->ipelmax + 1); sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    pelagic->u_values[sc_idx][x_idx][y_idx] =
                        pelagic->u_0 * exp(pelagic->lambda * grid->m_values[sc_idx]);
                    if (run->spatial_dim == 1 || run->spatial_dim == 2)
                    {
                        pelagic->u_values[sc_idx][x_idx][y_idx] *= x_start(grid->x_values[x_idx], 0);
                        if (run->spatial_dim == 2)
                        {
                            pelagic->u_values[sc_idx][x_idx][y_idx] *= y_start(grid->y_values[y_idx], 0);
                        }
                    }
                }
            }
        }
    }

    /*Open pelagic fishing file*/
    if (pelagic->fishing_flag == 1)
    {
        /*Create pelagic fishing input filename*/
        snprintf(pelagic->fname_fish_ts, sizeof(pelagic->fname_fish_ts), "%s/Input/%s_%s.txt", run->filename,
                 pelagic->filename, "fishing_ts");
        pelagic->fptr_fish_ts = safe_fopen(pelagic->fname_fish_ts, "r", __LINE__, __FILE__);
    }

    /*Setup pelagic reproduction if needed*/
    if (pelagic->ts_flag == 0)
    {
        /*set eggs equal to intial values*/
        if (pelagic->rep_method == 0)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    pelagic->reproduction[x_idx][y_idx] = pelagic->u_values[pelagic->ipelmin][x_idx][y_idx];
                }
            }
        }
        if (pelagic->rep_method == 1)
        {
            /*Create pelagic fishing input filename*/
            snprintf(pelagic->fname_rep_ts, sizeof(pelagic->fname_rep_ts), "%s/Input/%s_%s.txt", run->filename,
                     pelagic->filename, "rep_ts");
            pelagic->fptr_rep_ts = safe_fopen(pelagic->fname_rep_ts, "r", __LINE__, __FILE__);
        }
    }

    /*Setup values for predator-prey fedding kernel*/
    for (sc_idx = 0; sc_idx < (pelagic->ipelmax + 1); sc_idx++)
    {
        pelagic->phi_values[sc_idx] =
            phi((grid->m_values[sc_idx] - grid->m_values[0]), pelagic->q_0, pelagic->sig, pelagic->trunc);
    }

    /* --- Reef prey-refuge plugin (opt-in; Rogers, Blanchard & Mumby 2014) ---
       Availability a(w) = max(min_A, 1 - refuge(w)/competitor(w)) scales BOTH predation mortality on,
       and predator intake of, prey of body size w (mass balance). Reads an optional per-size refuge
       density from <run>/Input/<name>_refuge.txt (one value per grid size class). Absent file (or all
       zeros) -> refuge_on=0 -> a(w)=1 everywhere -> bit-identical to standard dbpmr. min_A (floor) from
       the REFUGE_MIN_A env var (default 0). */
    pelagic->refuge_on = 0;
    pelagic->refuge_min_a = 0.0;
    {
        char *e = getenv("REFUGE_MIN_A");
        if (e)
            pelagic->refuge_min_a = atof(e);
    }
    {
        char rfn[256];
        snprintf(rfn, sizeof(rfn), "%s/Input/%s_refuge.txt", run->filename, pelagic->filename);
        FILE *rf = fopen(rfn, "r");
        if (rf)
        {
            double v;
            int cnt = 0;
            for (sc_idx = 0; sc_idx < n_m && fscanf(rf, "%lf", &v) == 1; sc_idx++)
            {
                pelagic->refuge_values[sc_idx] = v;
                if (v > 0)
                    cnt++;
            }
            fclose(rf);
            if (cnt > 0)
                pelagic->refuge_on = 1; /*only activate when a positive refuge is supplied*/
        }
    }
}

void setup_benthic(RUN *run, GRID *grid, BENTHIC *benthic, double *ben_params, char *filename, int rep_method,
                   int initial_flag, int ts_flag, int fishing_flag)
{
    /*Declare species name*/
    benthic->filename = filename;

    /*Benthic Species Size Range*/
    benthic->benmin = ben_params[0];
    benthic->benmat = ben_params[1];
    benthic->benmax = ben_params[2];
    benthic->ibenmin = (int)rint((benthic->benmin - grid->mmin) / grid->mstep);
    benthic->ibenmat = (int)rint((benthic->benmat - grid->mmin) / grid->mstep);
    benthic->ibenmax = (int)rint((benthic->benmax - grid->mmin) / grid->mstep);

    /*Benthic Species Constants*/
    benthic->A = ben_params[3];         // Volume searched by unit weight
    benthic->alpha = ben_params[4];     // Volume search allometry
    benthic->mu_0 = ben_params[5];      // natural mortality rate
    benthic->beta = ben_params[6];      // natural mortality allometry
    benthic->A_base = benthic->A;       // temp driver: base (unscaled) search rate
    benthic->mu_0_base = benthic->mu_0; // temp driver: base (unscaled) background mortality
    benthic->tempeff = 1.0;             // temp driver: neutral until a temperature series is read
    benthic->mu_s = ben_params[7];      // senescence mortality rate
    benthic->epsilon = ben_params[8];   // senescece mortality constant
    benthic->u_0 = ben_params[9];       // initial intercept
    benthic->lambda = ben_params[10];   // primary production slope

    /*Benthic Species Conversion Efficiencies*/
    benthic->K_det = ben_params[11]; // growth efficiency
    benthic->R_det = ben_params[12];
    benthic->Ex_det = ben_params[13];

    benthic->pref_det = ben_params[14];

    /*Benthic Species Flags and Methods*/
    benthic->rep_method = rep_method;
    benthic->initial_flag = initial_flag;
    benthic->ts_flag = ts_flag;
    benthic->fishing_flag = fishing_flag;

    /*Create species output filenames*/
    snprintf(benthic->fname_r, sizeof(benthic->fname_r), "%s/%s/%s.txt", run->filename, benthic->filename, "results");
    snprintf(benthic->fname_g, sizeof(benthic->fname_g), "%s/%s/%s.txt", run->filename, benthic->filename, "growth");
    snprintf(benthic->fname_m, sizeof(benthic->fname_m), "%s/%s/%s.txt", run->filename, benthic->filename, "mortality");
    snprintf(benthic->fname_pred, sizeof(benthic->fname_pred), "%s/%s/%s.txt", run->filename, benthic->filename,
             "predation");
    snprintf(benthic->fname_fish, sizeof(benthic->fname_fish), "%s/%s/%s.txt", run->filename, benthic->filename,
             "fishing");
    snprintf(benthic->fname_summ, sizeof(benthic->fname_summ), "%s/%s/%s.txt", run->filename, benthic->filename,
             "summary");

    /*Open output files*/
    benthic->fptr_r = safe_fopen(benthic->fname_r, "w", __LINE__, __FILE__);
    benthic->fptr_g = safe_fopen(benthic->fname_g, "w", __LINE__, __FILE__);
    benthic->fptr_m = safe_fopen(benthic->fname_m, "w", __LINE__, __FILE__);
    benthic->fptr_pred = safe_fopen(benthic->fname_pred, "w", __LINE__, __FILE__);
    benthic->fptr_fish = safe_fopen(benthic->fname_fish, "w", __LINE__, __FILE__);
    benthic->fptr_summ = safe_fopen(benthic->fname_summ, "w", __LINE__, __FILE__);

    /*Initialise output files*/
    fprintf(benthic->fptr_r, "filetype:      %s\n", "results");
    fprintf(benthic->fptr_r, "speciestype:   %s\n", "benthic");
    fprintf(benthic->fptr_r, "\n");
    print_mass_header(grid, benthic->fptr_r);

    fprintf(benthic->fptr_g, "filetype:      %s\n", "growth");
    fprintf(benthic->fptr_g, "speciestype:   %s\n", "benthic");
    fprintf(benthic->fptr_g, "\n");
    print_mass_header(grid, benthic->fptr_g);

    fprintf(benthic->fptr_m, "filetype:      %s\n", "mortality");
    fprintf(benthic->fptr_m, "speciestype:   %s\n", "benthic");
    fprintf(benthic->fptr_m, "\n");
    print_mass_header(grid, benthic->fptr_m);

    fprintf(benthic->fptr_pred, "filetype:      %s\n", "predation");
    fprintf(benthic->fptr_pred, "speciestype:   %s\n", "benthic");
    fprintf(benthic->fptr_pred, "\n");
    print_mass_header(grid, benthic->fptr_pred);

    fprintf(benthic->fptr_fish, "filetype:      %s\n", "fishing");
    fprintf(benthic->fptr_fish, "speciestype:   %s\n", "benthic");
    fprintf(benthic->fptr_fish, "\n");
    print_mass_header(grid, benthic->fptr_fish);

    fprintf(benthic->fptr_summ, "filetype:      %s\n", "summary");
    fprintf(benthic->fptr_summ, "speciestype:   %s\n", "benthic");
    fprintf(benthic->fptr_summ, "\n");
    print_benthic_header(benthic->fptr_summ);

    /*Close output files*/
    fclose(benthic->fptr_r);
    fclose(benthic->fptr_g);
    fclose(benthic->fptr_m);
    fclose(benthic->fptr_pred);
    fclose(benthic->fptr_fish);
    fclose(benthic->fptr_summ);

    /*Setup Species Abundance Arrays*/
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    benthic->u_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    benthic->g_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    benthic->mu_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    benthic->mu_pred_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    benthic->mu_fish_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);

    benthic->det_bio = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        benthic->u_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        benthic->g_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        benthic->mu_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        benthic->mu_pred_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        benthic->mu_fish_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

        benthic->det_bio[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            benthic->u_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            benthic->g_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            benthic->mu_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            benthic->mu_pred_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            benthic->mu_fish_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);

            benthic->det_bio[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                benthic->u_values[sc_idx][x_idx][y_idx] = 0;
                benthic->g_values[sc_idx][x_idx][y_idx] = 0;
                benthic->mu_values[sc_idx][x_idx][y_idx] = 0;
                benthic->mu_pred_values[sc_idx][x_idx][y_idx] = 0;
                benthic->mu_fish_values[sc_idx][x_idx][y_idx] = 0;

                benthic->det_bio[sc_idx][x_idx][y_idx] = 0;
            }
        }
    }
    benthic->det_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

    benthic->fish_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    benthic->pred_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    benthic->reproduction = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        benthic->det_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);

        benthic->fish_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        benthic->pred_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        benthic->reproduction[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            benthic->det_total[x_idx][y_idx] = 0;

            benthic->fish_total[x_idx][y_idx] = 0;
            benthic->pred_total[x_idx][y_idx] = 0;
            benthic->reproduction[x_idx][y_idx] = 0;
        }
    }

    /*Setup initial values for u(m,0,x,y)*/
    if (benthic->initial_flag == 1)
    {
        int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
        char *input;                   // input string to store line as string
        input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
        char sep[] = ",";    // identifying the ',' char as seperating values
        char *result = NULL; // declaring a pointer for storing each individual value as a char

        /*Create plankton input filename*/
        snprintf(benthic->fname_ts, sizeof(benthic->fname_ts), "%s/Input/%s_%s.txt", run->filename, benthic->filename,
                 "ts");
        benthic->fptr_ts = safe_fopen(benthic->fname_ts, "r", __LINE__, __FILE__);

        /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                fgets(input, MAX_LEN, benthic->fptr_ts);
                input[strlen(input) - 1] =
                    '\0'; // replacing the newline \n char of the line of text with a null \0 char

                result = (char *)strtok(input, sep); // this is the 1st value
                sc_idx = 0;
                while (result != NULL)
                {
                    benthic->u_values[sc_idx][x_idx][y_idx] = atof(result);
                    result = (char *)strtok(NULL, sep);
                    sc_idx++;
                }
            }
        }
        free(input);
        /*If this is only an intial value then close the file now*/
        if (benthic->ts_flag == 0)
        {
            fclose(benthic->fptr_ts);
        }
    }
    else
    {
        for (sc_idx = benthic->ibenmin; sc_idx < (benthic->ibenmax + 1); sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    benthic->u_values[sc_idx][x_idx][y_idx] =
                        benthic->u_0 * exp(benthic->lambda * grid->m_values[sc_idx]);
                    if (run->spatial_dim == 1 || run->spatial_dim == 2)
                    {
                        benthic->u_values[sc_idx][x_idx][y_idx] *= x_start(grid->x_values[x_idx], 0);
                        if (run->spatial_dim == 2)
                        {
                            benthic->u_values[sc_idx][x_idx][y_idx] *= y_start(grid->y_values[y_idx], 0);
                        }
                    }
                }
            }
        }
    }

    /*Setup initial values for fishing mortality*/
    if (benthic->fishing_flag == 1)
    {
        /*Create benthic fishing input filename*/
        snprintf(benthic->fname_fish_ts, sizeof(benthic->fname_fish_ts), "%s/Input/%s_%s.txt", run->filename,
                 benthic->filename, "fishing_ts");
        benthic->fptr_fish_ts = safe_fopen(benthic->fname_fish_ts, "r", __LINE__, __FILE__);
    }

    /*Setup benthic reproduction if needed*/
    if (benthic->ts_flag == 0)
    {
        /*set eggs equal to intial values*/
        if (benthic->rep_method == 0)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    benthic->reproduction[x_idx][y_idx] = benthic->u_values[benthic->ibenmin][x_idx][y_idx];
                }
            }
        }
        if (benthic->rep_method == 1)
        {
            /*Create benthic fishing input filename*/
            snprintf(benthic->fname_rep_ts, sizeof(benthic->fname_rep_ts), "%s/Input/%s_%s.txt", run->filename,
                     benthic->filename, "rep_ts");
            benthic->fptr_rep_ts = safe_fopen(benthic->fname_rep_ts, "r", __LINE__, __FILE__);
        }
    }
}

void setup_detritus(RUN *run, GRID *grid, DETRITUS *detritus, double *det_params, char *filename, int initial_flag,
                    int ts_flag)
{

    /*Setup Detritus Constants*/
    detritus->w_0 = det_params[0];

    /*Declare detritus name*/
    detritus->filename = filename;

    /*Setup Plankton Flags*/
    detritus->initial_flag = initial_flag;
    detritus->ts_flag = ts_flag;

    /*Create detritus output filenames*/
    snprintf(detritus->fname_summ, sizeof(detritus->fname_summ), "%s/%s/%s.txt", run->filename, detritus->filename,
             "summary");

    /*Open output files*/
    detritus->fptr_summ = safe_fopen(detritus->fname_summ, "w", __LINE__, __FILE__);

    fprintf(detritus->fptr_summ, "filetype:      %s\n", "summary");
    fprintf(detritus->fptr_summ, "speciestype:   %s\n", "detritus");
    fprintf(detritus->fptr_summ, "\n");
    print_detritus_header(detritus->fptr_summ);

    /*Close output files*/
    fclose(detritus->fptr_summ);

    /*Setup And Initialise detritus Abundance Arrays*/
    int x_idx, y_idx;
    int n_x, n_y;
    n_x = grid->xnum;
    n_y = grid->ynum;

    detritus->w_values = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    detritus->g_values = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    detritus->mu_values = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        detritus->w_values[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        detritus->g_values[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        detritus->mu_values[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            detritus->w_values[x_idx][y_idx] = 0;
            detritus->g_values[x_idx][y_idx] = 0;
            detritus->mu_values[x_idx][y_idx] = 0;
        }
    }

    /*Setup initial values for w(0,x,y)*/
    /*Use user defined values*/
    if (detritus->initial_flag == 1)
    {
        int MAX_LEN = 15; // this is the upper bound for the row length in characters
        char *input;      // input string to store line as string
        input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
        char sep[] = ",";    // identifying the ',' char as seperating values
        char *result = NULL; // declaring a pointer for storing each individual value as a char

        /*Create detritus input filename*/
        snprintf(detritus->fname_ts, sizeof(detritus->fname_ts), "%s/Input/%s_%s.txt", run->filename,
                 detritus->filename, "ts");
        detritus->fptr_ts = safe_fopen(detritus->fname_ts, "r", __LINE__, __FILE__);

        /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                fgets(input, MAX_LEN, detritus->fptr_ts);
                input[strlen(input) - 1] =
                    '\0'; // replacing the newline \n char of the line of text with a null \0 char

                result = (char *)strtok(input, sep); // this is the 1st value
                while (result != NULL)
                {
                    detritus->w_values[x_idx][y_idx] = atof(result);
                    result = (char *)strtok(NULL, sep);
                }
            }
        }
        free(input);
        /*If this is only an intial value then close the file now*/
        if (detritus->ts_flag == 0)
        {
            fclose(detritus->fptr_ts);
        }
    }
    /*Use standard values*/
    else
    {
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                detritus->w_values[x_idx][y_idx] = detritus->w_0;
            }
        }
    }
}

void setup_matrix(int dim, MATRIX *matrix)
/*Creates a matrix for the consumer size species*/
{
    matrix->size = dim;
    matrix->a = (double *)safe_malloc(dim * sizeof(double), __LINE__, __FILE__);
    matrix->b = (double *)safe_malloc(dim * sizeof(double), __LINE__, __FILE__);
    matrix->c = (double *)safe_malloc(dim * sizeof(double), __LINE__, __FILE__);
    matrix->r = (double *)safe_malloc(dim * sizeof(double), __LINE__, __FILE__);
    matrix->u = (double *)safe_malloc(dim * sizeof(double), __LINE__, __FILE__);
}

//-----------------------------------//
// Output Function. Not to be edited //
//-----------------------------------//

void calculate_results(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *pelmatrix, MATRIX *benmatrix,
                       MATRIX *xmatrix, MATRIX *ymatrix)
/* Sets up output txt files and implements time stepping process*/
{
    // Indices for each variable with readable names
    int sc_idx, ts_idx, x_idx, y_idx, pel_species_idx, ben_species_idx;
    // Counts for each variable with readable names
    int n_m, n_t, n_x, n_y, n_pel_species, n_ben_species;
    double time1 = 0, time2 = 0, timed = 0; // variables for countdown timer

    n_m = grid->mnum;
    n_t = grid->tnum;
    n_x = grid->xnum;
    n_y = grid->ynum;
    n_pel_species = run->no_pelagic;
    n_ben_species = run->no_benthic;

    /*Optional environmental-forcing plugin: auto-detect <run>/Input/forcing_ts.txt, whose
      first line names the channels present (comma-separated) and each subsequent line gives
      one timestep's values. We map the recognised channel names to column indices here; the
      per-timestep read + apply happens in mass_solver(). Absent file/channels -> neutral.*/
    run->forcing_flag = 0;
    run->forcing_ncol = 0;
    run->col_pel_tempeff = -1;
    run->col_ben_tempeff = -1;
    run->col_sinking_rate = -1;
    run->col_depth = -1;
    run->col_residence_time = -1;
    run->col_export_attn = -1;
    run->col_export_gamma = -1;
    run->col_export_magg = -1;
    run->col_gravity = -1;
    run->pel_tempeff = 1.0;
    run->ben_tempeff = 1.0;
    run->sinking_rate = 1.0;
    run->depth = 1.0;
    run->residence_time = 0.0;
    run->export_attn = 0.0;
    run->export_gamma = 0.0;
    run->export_magg = -1e30;
    run->gravity = 0.0;
    run->fptr_forcing = NULL;
    run->forcing_percell = 0;
    run->pel_tempeff_c = NULL;
    run->ben_tempeff_c = NULL;
    run->sinking_rate_c = NULL;
    run->depth_c = NULL;
    run->export_attn_c = NULL;
    run->export_gamma_c = NULL;
    run->export_magg_c = NULL;
    {
        char *e = getenv("TEMP_ON_MORT");
        run->temp_on_mort = e ? atoi(e) : 1;
    }
    {
        char fname_forcing[256];
        snprintf(fname_forcing, sizeof(fname_forcing), "%s/Input/forcing_ts.txt", run->filename);
        run->fptr_forcing = fopen(fname_forcing, "r");
        if (run->fptr_forcing != NULL)
        {
            run->forcing_flag = 1;
            char hdr[512];
            if (fgets(hdr, sizeof(hdr), run->fptr_forcing) != NULL)
            {
                /*strip trailing newline*/
                size_t hl = strlen(hdr);
                while (hl > 0 && (hdr[hl - 1] == '\n' || hdr[hl - 1] == '\r' || hdr[hl - 1] == ' '))
                {
                    hdr[hl - 1] = '\0';
                    hl--;
                }
                int col = 0;
                char *name = strtok(hdr, ",");
                while (name != NULL)
                {
                    /*trim leading spaces*/
                    while (*name == ' ')
                        name++;
                    if (strcmp(name, "pel_tempeff") == 0)
                        run->col_pel_tempeff = col;
                    else if (strcmp(name, "ben_tempeff") == 0)
                        run->col_ben_tempeff = col;
                    else if (strcmp(name, "sinking_rate") == 0)
                        run->col_sinking_rate = col;
                    else if (strcmp(name, "depth") == 0)
                        run->col_depth = col;
                    else if (strcmp(name, "residence_time") == 0)
                        run->col_residence_time = col;
                    else if (strcmp(name, "export_attn") == 0)
                        run->col_export_attn = col;
                    else if (strcmp(name, "export_gamma") == 0)
                        run->col_export_gamma = col;
                    else if (strcmp(name, "export_magg") == 0)
                        run->col_export_magg = col;
                    else if (strcmp(name, "gravity") == 0)
                        run->col_gravity = col;
                    col++;
                    name = strtok(NULL, ",");
                }
                run->forcing_ncol = col;
            }
        }
    }
    /* Per-cell forcing (opt-in): FORCING_PERCELL=1 + spatial_dim>=1 => forcing_ts.txt has one row per
       CELL per timestep (x*y rows, ordered x-outer y-inner, like the plankton _ts file), so temperature
       and export can vary across space. Allocate the per-cell multiplier arrays (init neutral=1). When
       off (default), the single global row is used and these arrays are unused -> behaviour unchanged. */
    if (run->forcing_flag == 1)
    {
        char *e = getenv("FORCING_PERCELL");
        run->forcing_percell = (e && atoi(e) != 0 && (run->spatial_dim == 1 || run->spatial_dim == 2)) ? 1 : 0;

        run->pel_tempeff_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->ben_tempeff_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->sinking_rate_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->depth_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->export_attn_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->export_gamma_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        run->export_magg_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            run->pel_tempeff_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->ben_tempeff_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->sinking_rate_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->depth_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->export_attn_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->export_gamma_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            run->export_magg_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                run->pel_tempeff_c[x_idx][y_idx] = 1.0;
                run->ben_tempeff_c[x_idx][y_idx] = 1.0;
                run->sinking_rate_c[x_idx][y_idx] = 1.0;
                run->depth_c[x_idx][y_idx] = 1.0;
                run->export_attn_c[x_idx][y_idx] = 0.0;
                run->export_gamma_c[x_idx][y_idx] = 0.0;
                run->export_magg_c[x_idx][y_idx] = -1e30;
            }
        }
    }

    /*Create temporary arrays for timestepping process*/
    double ***temp_plankton;
    double ****temp_pelagic;
    double ****temp_benthic;
    double **temp_detritus;

    temp_plankton = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        temp_plankton[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            temp_plankton[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                temp_plankton[sc_idx][x_idx][y_idx] = 0;
            }
        }
    }
    temp_pelagic = (double ****)safe_malloc(n_pel_species * sizeof(double ***), __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        temp_pelagic[pel_species_idx] = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
        for (sc_idx = 0; sc_idx < n_m; sc_idx++)
        {
            temp_pelagic[pel_species_idx][sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                temp_pelagic[pel_species_idx][sc_idx][x_idx] =
                    (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    temp_pelagic[pel_species_idx][sc_idx][x_idx][y_idx] = 0;
                }
            }
        }
    }

    if (run->no_benthic != 0)
    {
        temp_benthic = (double ****)safe_malloc(n_ben_species * sizeof(double ***), __LINE__, __FILE__);
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            temp_benthic[ben_species_idx] = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
            for (sc_idx = 0; sc_idx < n_m; sc_idx++)
            {
                temp_benthic[ben_species_idx][sc_idx] =
                    (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    temp_benthic[ben_species_idx][sc_idx][x_idx] =
                        (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        temp_benthic[ben_species_idx][sc_idx][x_idx][y_idx] = 0;
                    }
                }
            }
        }
        temp_detritus = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            temp_detritus[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                temp_detritus[x_idx][y_idx] = 0;
            }
        }
    }

    /*Calculate initial fishing pressures*/
    calculate_fishing(run, grid, community);

    /*Calculate starting g and mu values*/
    calculate_g_and_mu(run, grid, community);

    /*Calculate starting reproduction allocation*/
    calculate_reproduction(run, grid, community);

    /*Print run info and starting output values*/
    /*plankton*/
    if (grid->toutmin == 0)
    {
        print_timestep_plankton(grid, community->plankton, 0);
        print_plankton_summary(grid, community->plankton, 0);
    }

    /*Pelagic*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (grid->toutmin == 0)
        {
            print_timestep_pelagic(run, grid, &(community->pelagic[pel_species_idx]), 0);
            print_pelagic_summary(grid, &(community->pelagic[pel_species_idx]), 0);
        }
    }

    if (run->no_benthic != 0)
    {
        /*Benthic*/
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            if (grid->toutmin == 0)
            {
                print_timestep_benthic(grid, &(community->benthic[ben_species_idx]), 0);
                print_benthic_summary(grid, &(community->benthic[ben_species_idx]), 0);
            }
        }
        /*Detritus*/
        if (grid->toutmin == 0)
        {
            print_detritus_summary(grid, community->detritus, 0);
        }
    }

    /*Parameter File*/
    run->fptr_summ = safe_fopen(run->fname_summ, "a", __LINE__, __FILE__);

    print_run(run, run->fptr_summ);
    print_grid(grid, run->fptr_summ);
    print_plankton(community->plankton, run->fptr_summ);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        print_pelagic(&(community->pelagic[pel_species_idx]), run->fptr_summ);
    }
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            print_benthic(&(community->benthic[ben_species_idx]), run->fptr_summ);
        }
        print_detritus(community->detritus, run->fptr_summ);
    }
    fclose(run->fptr_summ);

    /*Time Stepping Loop*/
    for (ts_idx = 1; ts_idx < n_t; ts_idx++)
    {

        /*Countdown timer for the program*/

        // Checks if user has terminated function in R (e.g. by pressing escape)
        R_CheckUserInterrupt();

        /*Countdown timer for the program*/
        if (ts_idx == 1)
        {
            time2 = (double)clock() / CLOCKS_PER_SEC;
        }
        else
        {
            if (ts_idx % ((int)rint(grid->toutstep / grid->tstep)) == 0)
            {
                time1 = time2;
                time2 = (double)clock() / CLOCKS_PER_SEC;
                timed = (time2 - time1) * (grid->tstep / grid->toutstep);
            }
            Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            R_FlushConsole();
            Rprintf("%6.1f   %d", timed * (grid->tnum - ts_idx), grid->tnum - ts_idx);
            R_FlushConsole();
        }
        if (ts_idx == (n_t - 1))
        {
            Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
            R_FlushConsole();
            Rprintf("\n");
            R_FlushConsole();
        }

        /*Store previous time step in temporary array for each species*/
        for (sc_idx = 0; sc_idx < n_m; sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    temp_plankton[sc_idx][x_idx][y_idx] = community->plankton->u_values[sc_idx][x_idx][y_idx];
                }
            }
        }
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            for (sc_idx = 0; sc_idx < n_m; sc_idx++)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        temp_pelagic[pel_species_idx][sc_idx][x_idx][y_idx] =
                            community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx];
                    }
                }
            }
        }

        if (run->no_benthic != 0)
        {
            for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
            {
                for (sc_idx = 0; sc_idx < n_m; sc_idx++)
                {
                    for (x_idx = 0; x_idx < n_x; x_idx++)
                    {
                        for (y_idx = 0; y_idx < n_y; y_idx++)
                        {
                            temp_benthic[ben_species_idx][sc_idx][x_idx][y_idx] =
                                community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx];
                        }
                    }
                }
            }
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    temp_detritus[x_idx][y_idx] = community->detritus->w_values[x_idx][y_idx];
                }
            }
        }

        /*Calculate u(m,t) using previous timestep values*/
        mass_solver(run, grid, community, pelmatrix, benmatrix, temp_plankton, temp_pelagic, temp_benthic,
                    temp_detritus);

        /*Calculate fishing pressures*/
        calculate_fishing(run, grid, community);

        /*Calculate g and mu using current timestep values for use in movement and next time step*/
        calculate_g_and_mu(run, grid, community);

        /*Calculate starting reproduction allocation*/
        calculate_reproduction(run, grid, community);

        /*If Spatial Model is being used (1D = x movement only; 2D = x then y).
          NOTE: the outer guard must admit spatial_dim==2 as well — the inner spatial_dim==2
          block below only runs when this outer test passes, so a `==1`-only guard left ALL
          movement dead for 2D runs (neither x nor y solver ran). */
        if (run->spatial_dim == 1 || run->spatial_dim == 2)
        {

            /*Store itermediate time step as temp_pelagic*/

            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                for (sc_idx = 0; sc_idx < n_m; sc_idx++)
                {
                    for (x_idx = 0; x_idx < n_x; x_idx++)
                    {
                        for (y_idx = 0; y_idx < n_y; y_idx++)
                        {
                            temp_pelagic[pel_species_idx][sc_idx][x_idx][y_idx] =
                                community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx];
                        }
                    }
                }
            }

            /*Calculate spatial movements using intermediate values*/
            xmove_solver(run, grid, community, xmatrix, temp_pelagic, ts_idx);

            if (run->spatial_dim == 2)
            {

                /*Store itermediate time step as temp_pelagic*/
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
                    {
                        for (x_idx = 0; x_idx < n_x; x_idx++)
                        {
                            for (y_idx = 0; y_idx < n_y; y_idx++)
                            {
                                temp_pelagic[pel_species_idx][sc_idx][x_idx][y_idx] =
                                    community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx];
                            }
                        }
                    }
                }
                /*Calculate spatial movements using intermediate values*/
                ymove_solver(run, grid, community, ymatrix, temp_pelagic, ts_idx);
            }
        }

        /*Print current timestep values to files*/
        if (grid->t_values[ts_idx] >= grid->toutmin && grid->t_values[ts_idx] <= grid->toutmax)
        {
            if (ts_idx % ((int)rint(grid->toutstep / grid->tstep)) == 0)
            {

                print_timestep_plankton(grid, community->plankton, ts_idx);
                print_plankton_summary(grid, community->plankton, ts_idx);
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    print_timestep_pelagic(run, grid, &(community->pelagic[pel_species_idx]), ts_idx);
                    print_pelagic_summary(grid, &(community->pelagic[pel_species_idx]), ts_idx);
                }
                if (run->no_benthic != 0)
                {
                    for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
                    {
                        print_timestep_benthic(grid, &(community->benthic[ben_species_idx]), ts_idx);
                        print_benthic_summary(grid, &(community->benthic[ben_species_idx]), ts_idx);
                    }
                    print_detritus_summary(grid, community->detritus, ts_idx);
                }
            }
        }
    }

    /*Print final detailed timestep*/
    print_detailed_plankton(grid, community->plankton, n_t - 1);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        print_detailed_pelagic(run, grid, &(community->pelagic[pel_species_idx]), n_t - 1);
    }
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            print_detailed_benthic(grid, &(community->benthic[ben_species_idx]), n_t - 1);
        }
    }

    /*Free pointers*/
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            free(temp_plankton[sc_idx][x_idx]);
        }
        free(temp_plankton[sc_idx]);
    }
    free(temp_plankton);

    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        for (sc_idx = 0; sc_idx < n_m; sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                free(temp_pelagic[pel_species_idx][sc_idx][x_idx]);
            }
            free(temp_pelagic[pel_species_idx][sc_idx]);
        }
        free(temp_pelagic[pel_species_idx]);
    }
    free(temp_pelagic);

    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            for (sc_idx = 0; sc_idx < n_m; sc_idx++)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    free(temp_benthic[ben_species_idx][sc_idx][x_idx]);
                }
                free(temp_benthic[ben_species_idx][sc_idx]);
            }
            free(temp_benthic[ben_species_idx]);
        }
        free(temp_benthic);

        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            free(temp_detritus[x_idx]);
        }
        free(temp_detritus);
    }

    /*Close ts files*/
    if (community->plankton->ts_flag == 1)
    {
        fclose(community->plankton->fptr_ts);
    }
    if (run->forcing_flag == 1)
    {
        fclose(run->fptr_forcing);
        if (run->pel_tempeff_c != NULL)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                free(run->pel_tempeff_c[x_idx]);
                free(run->ben_tempeff_c[x_idx]);
                free(run->sinking_rate_c[x_idx]);
                free(run->depth_c[x_idx]);
                free(run->export_attn_c[x_idx]);
                free(run->export_gamma_c[x_idx]);
                free(run->export_magg_c[x_idx]);
            }
            free(run->pel_tempeff_c);
            free(run->ben_tempeff_c);
            free(run->sinking_rate_c);
            free(run->depth_c);
            free(run->export_attn_c);
            free(run->export_gamma_c);
            free(run->export_magg_c);
        }
    }
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (community->pelagic[pel_species_idx].ts_flag == 1)
        {
            fclose(community->pelagic[pel_species_idx].fptr_ts);
        }
    }
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            if (community->benthic[ben_species_idx].ts_flag == 1)
            {
                fclose(community->benthic[ben_species_idx].fptr_ts);
            }
        }
        if (community->detritus->ts_flag == 1)
        {
            fclose(community->detritus->fptr_ts);
        }
    }

    /*Close fishing files*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (community->pelagic[pel_species_idx].fishing_flag == 1)
        {
            fclose(community->pelagic[pel_species_idx].fptr_fish_ts);
        }
    }
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            if (community->benthic[ben_species_idx].fishing_flag == 1)
            {
                fclose(community->benthic[ben_species_idx].fptr_fish_ts);
            }
        }
    }

    /*Close reproduction files*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (community->pelagic[pel_species_idx].rep_method == 1)
        {
            if (community->pelagic[pel_species_idx].ts_flag == 0)
            {
                fclose(community->pelagic[pel_species_idx].fptr_rep_ts);
            }
        }
    }
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            if (community->benthic[ben_species_idx].rep_method == 1)
            {
                if (community->benthic[ben_species_idx].ts_flag == 0)
                {
                    fclose(community->benthic[ben_species_idx].fptr_rep_ts);
                }
            }
        }
    }
}

void print_run(RUN *run, FILE *fptr)
{
    /*Output specific C code info*/
    fprintf(fptr, "Run Parameters\n");
    fprintf(fptr, "filename:      %s\n", run->filename);
    fprintf(fptr, "no_pelagic:    %d\n", run->no_pelagic);
    fprintf(fptr, "no_benthic:    %d\n", run->no_benthic);
    fprintf(fptr, "spatial_dim:   %d\n", run->spatial_dim);
    fprintf(fptr, "coupled_flag:  %d\n", run->coupled_flag);
    fprintf(fptr, "diff_method:   %d\n", run->diff_method);
    fprintf(fptr, "\n");
}

void print_grid(GRID *grid, FILE *fptr)
{
    /*Output discretisation information*/
    fprintf(fptr, "Grid Parameters\n");
    fprintf(fptr, "mmin:          %f\n", grid->mmin);
    fprintf(fptr, "mmax:          %f\n", grid->mmax);
    fprintf(fptr, "mstep:         %f\n", grid->mstep);
    fprintf(fptr, "moutstep:      %f\n", grid->moutstep);

    fprintf(fptr, "t1:            %.10f\n", grid->t1);
    fprintf(fptr, "tmax:          %.10f\n", grid->tmax);
    fprintf(fptr, "tstep:         %.10f\n", grid->tstep);
    fprintf(fptr, "toutmin:       %.10f\n", grid->toutmin);
    fprintf(fptr, "toutmax:       %.10f\n", grid->toutmax);
    fprintf(fptr, "toutstep:      %.10f\n", grid->toutstep);

    fprintf(fptr, "xmin:          %f\n", grid->xmin);
    fprintf(fptr, "xmax:          %f\n", grid->xmax);
    fprintf(fptr, "xstep:         %f\n", grid->xstep);
    fprintf(fptr, "xoutstep:      %f\n", grid->xoutstep);

    fprintf(fptr, "ymin:          %f\n", grid->ymin);
    fprintf(fptr, "ymax:          %f\n", grid->ymax);
    fprintf(fptr, "ystep:         %f\n", grid->ystep);
    fprintf(fptr, "youtstep:      %f\n", grid->youtstep);
    fprintf(fptr, "\n");
}

void print_plankton(PLANKTON *plankton, FILE *fptr)
{
    /*Ouput plankton specific information*/
    fprintf(fptr, "Plankton Parameters\n");
    fprintf(fptr, "filename:      %s\n", plankton->filename);
    fprintf(fptr, "speciestype:   %s\n", "plankton");
    fprintf(fptr, "mmin:        %f\n", plankton->plamin);
    fprintf(fptr, "mmax:        %f\n", plankton->plamax);
    fprintf(fptr, "mu_0:          %f\n", plankton->mu_0);
    fprintf(fptr, "beta:          %f\n", plankton->beta);
    fprintf(fptr, "u_0:           %f\n", plankton->u_0);
    fprintf(fptr, "lambda:        %f\n", plankton->lambda);
    fprintf(fptr, "initial_flag:  %d\n", plankton->initial_flag);
    fprintf(fptr, "ts_flag:       %d\n", plankton->ts_flag);
    fprintf(fptr, "\n");
}

void print_pelagic(PELAGIC *pelagic, FILE *fptr)
{
    /*Ouput species specific information*/
    fprintf(fptr, "Pelagic Parameters\n");
    fprintf(fptr, "filename:      %s\n", pelagic->filename);
    fprintf(fptr, "speciestype:   %s\n", "pelagic");

    fprintf(fptr, "mmin:        %f\n", pelagic->pelmin);
    fprintf(fptr, "mmat:        %f\n", pelagic->pelmat);
    fprintf(fptr, "mmax:        %f\n", pelagic->pelmax);

    fprintf(fptr, "A:             %f\n", pelagic->A);
    fprintf(fptr, "alpha:         %f\n", pelagic->alpha);
    fprintf(fptr, "mu_0:          %f\n", pelagic->mu_0);
    fprintf(fptr, "beta:          %f\n", pelagic->beta);
    fprintf(fptr, "mu_s:          %f\n", pelagic->mu_s);
    fprintf(fptr, "epsilon:       %f\n", pelagic->epsilon);
    fprintf(fptr, "u_0:           %f\n", pelagic->u_0);
    fprintf(fptr, "lambda:        %f\n", pelagic->lambda);

    fprintf(fptr, "K_pla:         %f\n", pelagic->K_pla);
    fprintf(fptr, "R_pla:         %f\n", pelagic->R_pla);
    fprintf(fptr, "Ex_pla:        %f\n", pelagic->Ex_pla);
    fprintf(fptr, "K_pel:         %f\n", pelagic->K_pel);
    fprintf(fptr, "R_pel:         %f\n", pelagic->R_pel);
    fprintf(fptr, "Ex_pel:        %f\n", pelagic->Ex_pel);
    fprintf(fptr, "K_ben:         %f\n", pelagic->K_ben);
    fprintf(fptr, "R_ben:         %f\n", pelagic->R_ben);
    fprintf(fptr, "Ex_ben:        %f\n", pelagic->Ex_ben);

    fprintf(fptr, "pref_pla:      %f\n", pelagic->pref_pla);
    fprintf(fptr, "pref_pel:      %f\n", pelagic->pref_pel);
    fprintf(fptr, "pref_ben:      %f\n", pelagic->pref_ben);

    fprintf(fptr, "q_0:           %f\n", pelagic->q_0);
    fprintf(fptr, "sig:           %f\n", pelagic->sig);
    fprintf(fptr, "trunc:         %f\n", pelagic->trunc);

    fprintf(fptr, "prey:          %f\n", pelagic->prey);
    fprintf(fptr, "pred:          %f\n", pelagic->pred);
    fprintf(fptr, "comp:          %f\n", pelagic->comp);
    fprintf(fptr, "gamma_pred:    %f\n", pelagic->gamma_pred);
    fprintf(fptr, "gamma_prey:    %f\n", pelagic->gamma_prey);
    fprintf(fptr, "gamma_comp:    %f\n", pelagic->gamma_comp);

    fprintf(fptr, "rep_method:    %d\n", pelagic->rep_method);
    fprintf(fptr, "initial_flag:  %d\n", pelagic->initial_flag);
    fprintf(fptr, "ts_flag:       %d\n", pelagic->ts_flag);
    fprintf(fptr, "fishing_flag:  %d\n", pelagic->fishing_flag);
    fprintf(fptr, "\n");
}

void print_benthic(BENTHIC *benthic, FILE *fptr)
{
    /*Ouput species specific information*/
    fprintf(fptr, "Benthic Parameters\n");
    fprintf(fptr, "filename:      %s\n", benthic->filename);
    fprintf(fptr, "speciestype:   %s\n", "benthic");

    fprintf(fptr, "mmin:        %f\n", benthic->benmin);
    fprintf(fptr, "mmat:        %f\n", benthic->benmat);
    fprintf(fptr, "mmax:        %f\n", benthic->benmax);

    fprintf(fptr, "A:             %f\n", benthic->A);
    fprintf(fptr, "alpha:         %f\n", benthic->alpha);
    fprintf(fptr, "mu_0:          %f\n", benthic->mu_0);
    fprintf(fptr, "beta:          %f\n", benthic->beta);
    fprintf(fptr, "mu_s:          %f\n", benthic->mu_s);
    fprintf(fptr, "epsilon:       %f\n", benthic->epsilon);
    fprintf(fptr, "u_0:           %f\n", benthic->u_0);
    fprintf(fptr, "lambda:        %f\n", benthic->lambda);

    fprintf(fptr, "K_det:         %f\n", benthic->K_det);
    fprintf(fptr, "R_det:         %f\n", benthic->R_det);
    fprintf(fptr, "Ex_det:        %f\n", benthic->Ex_det);

    fprintf(fptr, "pref_det:      %f\n", benthic->pref_det);

    fprintf(fptr, "rep_method:    %d\n", benthic->rep_method);
    fprintf(fptr, "initial_flag:  %d\n", benthic->initial_flag);
    fprintf(fptr, "ts_flag:       %d\n", benthic->ts_flag);
    fprintf(fptr, "fishing_flag:  %d\n", benthic->fishing_flag);
    fprintf(fptr, "\n");
}

void print_detritus(DETRITUS *detritus, FILE *fptr)
{
    /*Ouput plankton specific information*/
    fprintf(fptr, "Detritus Parameters\n");
    fprintf(fptr, "filename:      %s\n", detritus->filename);
    fprintf(fptr, "speciestype:   %s\n", "detritus");

    fprintf(fptr, "w_0:           %f\n", detritus->w_0);

    fprintf(fptr, "initial_flag:  %d\n", detritus->initial_flag);
    fprintf(fptr, "ts_flag:       %d\n", detritus->ts_flag);
    fprintf(fptr, "\n");
}

void print_mass_header(GRID *grid, FILE *fptr)
{
    /*Output Header Row*/
    int sc_idx;
    int n_m = grid->mnum;

    fprintf(fptr, "Ouput Data\n");
    fprintf(fptr, "t,x,y,");
    for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx += (int)rint(grid->moutstep / grid->mstep))
    {
        fprintf(fptr, "%f,", grid->m_values[sc_idx]);
    }
    fprintf(fptr, "%f\n", grid->m_values[n_m - 1]);
}

void print_detailed_header(GRID *grid, FILE *fptr)
{
    int sc_idx;
    int n_m;
    n_m = grid->mnum;

    /*Print Header*/
    fprintf(fptr, "\nDetailed Timestep\n");
    fprintf(fptr, "t,x,y,");
    for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx++)
    {
        fprintf(fptr, "%f,", grid->m_values[sc_idx]);
    }
    fprintf(fptr, "%f\n", grid->m_values[n_m - 1]);
}

void print_plankton_header(FILE *fptr)
{
    fprintf(fptr, "Ouput Data\n");
    fprintf(fptr, "t,x,y");
    fprintf(fptr, ",biomass");
    fprintf(fptr, "\n");
}

void print_pelagic_header(FILE *fptr)
{
    fprintf(fptr, "Ouput Data\n");
    fprintf(fptr, "t,x,y");
    fprintf(fptr, ",biomass");
    fprintf(fptr, ",plankton");
    fprintf(fptr, ",pelagic");
    fprintf(fptr, ",benthic");
    fprintf(fptr, ",predation");
    fprintf(fptr, ",fishing");
    fprintf(fptr, ",reproduction");
    fprintf(fptr, "\n");
}

void print_benthic_header(FILE *fptr)
{
    fprintf(fptr, "Ouput Data\n");
    fprintf(fptr, "t,x,y");
    fprintf(fptr, ",biomass");
    fprintf(fptr, ",detritus");
    fprintf(fptr, ",predation");
    fprintf(fptr, ",fishing");
    fprintf(fptr, ",reproduction");
    fprintf(fptr, "\n");
}

void print_detritus_header(FILE *fptr)
{
    fprintf(fptr, "Ouput Data\n");
    fprintf(fptr, "t,x,y");
    fprintf(fptr, ",biomass");
    fprintf(fptr, ",detritus_in");
    fprintf(fptr, ",detritus_out");
    fprintf(fptr, "\n");
}

void print_timestep_plankton(GRID *grid, PLANKTON *plankton, int time)
/*Outputs timestep values*/
{
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    plankton->fptr_r = safe_fopen(plankton->fname_r, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(plankton->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx += (int)rint(grid->moutstep / grid->mstep))
            {
                fprintf(plankton->fptr_r, "%.16f,", plankton->u_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(plankton->fptr_r, "%.16f\n", plankton->u_values[n_m - 1][x_idx][y_idx]);
        }
    }

    fclose(plankton->fptr_r);
}

void print_detailed_plankton(GRID *grid, PLANKTON *plankton, int time)
/*Outputs timestep values at high resolution*/
{
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    plankton->fptr_r = safe_fopen(plankton->fname_r, "a", __LINE__, __FILE__);

    /*Print Header*/
    print_detailed_header(grid, plankton->fptr_r);

    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            fprintf(plankton->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx++)
            {
                fprintf(plankton->fptr_r, "%.16f,", plankton->u_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(plankton->fptr_r, "%.16f\n", plankton->u_values[n_m - 1][x_idx][y_idx]);
        }
    }

    fclose(plankton->fptr_r);
}

void print_timestep_pelagic(RUN *run, GRID *grid, PELAGIC *pelagic, int time)
/*Outputs timestep values*/
{
    int pel_species_idx, sc_idx, x_idx, y_idx;
    int n_pel_species, n_m, n_x, n_y;
    n_pel_species = run->no_pelagic;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    pelagic->fptr_r = safe_fopen(pelagic->fname_r, "a", __LINE__, __FILE__);
    pelagic->fptr_g = safe_fopen(pelagic->fname_g, "a", __LINE__, __FILE__);
    pelagic->fptr_m = safe_fopen(pelagic->fname_m, "a", __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        pelagic->fptr_pred[pel_species_idx] = safe_fopen(pelagic->fname_pred[pel_species_idx], "a", __LINE__, __FILE__);
    }
    pelagic->fptr_fish = safe_fopen(pelagic->fname_fish, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(pelagic->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(pelagic->fptr_g, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(pelagic->fptr_m, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                fprintf(pelagic->fptr_pred[pel_species_idx], "%.10f,%f,%f,", (grid->t_values[time]),
                        grid->x_values[x_idx], grid->y_values[y_idx]);
            }
            fprintf(pelagic->fptr_fish, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx += (int)rint(grid->moutstep / grid->mstep))
            {
                fprintf(pelagic->fptr_r, "%.16f,", pelagic->u_values[sc_idx][x_idx][y_idx]);
                fprintf(pelagic->fptr_g, "%.16f,", pelagic->g_values[sc_idx][x_idx][y_idx]);
                fprintf(pelagic->fptr_m, "%.16f,", pelagic->mu_values[sc_idx][x_idx][y_idx]);
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    fprintf(pelagic->fptr_pred[pel_species_idx], "%.16f,",
                            pelagic->mu_pred_values[pel_species_idx][sc_idx][x_idx][y_idx]);
                }
                fprintf(pelagic->fptr_fish, "%.16f,", pelagic->mu_fish_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(pelagic->fptr_r, "%.16f\n", pelagic->u_values[n_m - 1][x_idx][y_idx]);
            fprintf(pelagic->fptr_g, "%.16f\n", pelagic->g_values[n_m - 1][x_idx][y_idx]);
            fprintf(pelagic->fptr_m, "%.16f\n", pelagic->mu_values[n_m - 1][x_idx][y_idx]);
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                fprintf(pelagic->fptr_pred[pel_species_idx], "%.16f\n",
                        pelagic->mu_pred_values[pel_species_idx][n_m - 1][x_idx][y_idx]);
            }
            fprintf(pelagic->fptr_fish, "%.16f\n", pelagic->mu_fish_values[n_m - 1][x_idx][y_idx]);
        }
    }

    /*Close output files*/
    fclose(pelagic->fptr_r);
    fclose(pelagic->fptr_g);
    fclose(pelagic->fptr_m);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        fclose(pelagic->fptr_pred[pel_species_idx]);
    }
    fclose(pelagic->fptr_fish);
}

void print_detailed_pelagic(RUN *run, GRID *grid, PELAGIC *pelagic, int time)
/*Outputs timestep values*/
{
    int pel_species_idx, sc_idx, x_idx, y_idx;
    int n_pel_species, n_m, n_x, n_y;
    n_pel_species = run->no_pelagic;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    pelagic->fptr_r = safe_fopen(pelagic->fname_r, "a", __LINE__, __FILE__);
    pelagic->fptr_g = safe_fopen(pelagic->fname_g, "a", __LINE__, __FILE__);
    pelagic->fptr_m = safe_fopen(pelagic->fname_m, "a", __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        pelagic->fptr_pred[pel_species_idx] = safe_fopen(pelagic->fname_pred[pel_species_idx], "a", __LINE__, __FILE__);
    }
    pelagic->fptr_fish = safe_fopen(pelagic->fname_fish, "a", __LINE__, __FILE__);

    /*Print Header*/
    print_detailed_header(grid, pelagic->fptr_r);
    print_detailed_header(grid, pelagic->fptr_g);
    print_detailed_header(grid, pelagic->fptr_m);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        print_detailed_header(grid, pelagic->fptr_pred[pel_species_idx]);
    }
    print_detailed_header(grid, pelagic->fptr_fish);

    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            fprintf(pelagic->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(pelagic->fptr_g, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(pelagic->fptr_m, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                fprintf(pelagic->fptr_pred[pel_species_idx], "%.10f,%f,%f,", (grid->t_values[time]),
                        grid->x_values[x_idx], grid->y_values[y_idx]);
            }
            fprintf(pelagic->fptr_fish, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx++)
            {
                fprintf(pelagic->fptr_r, "%.16f,", pelagic->u_values[sc_idx][x_idx][y_idx]);
                fprintf(pelagic->fptr_g, "%.16f,", pelagic->g_values[sc_idx][x_idx][y_idx]);
                fprintf(pelagic->fptr_m, "%.16f,", pelagic->mu_values[sc_idx][x_idx][y_idx]);
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    fprintf(pelagic->fptr_pred[pel_species_idx], "%.16f,",
                            pelagic->mu_pred_values[pel_species_idx][sc_idx][x_idx][y_idx]);
                }
                fprintf(pelagic->fptr_fish, "%.16f,", pelagic->mu_fish_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(pelagic->fptr_r, "%.16f\n", pelagic->u_values[n_m - 1][x_idx][y_idx]);
            fprintf(pelagic->fptr_g, "%.16f\n", pelagic->g_values[n_m - 1][x_idx][y_idx]);
            fprintf(pelagic->fptr_m, "%.16f\n", pelagic->mu_values[n_m - 1][x_idx][y_idx]);
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                fprintf(pelagic->fptr_pred[pel_species_idx], "%.16f\n",
                        pelagic->mu_pred_values[pel_species_idx][n_m - 1][x_idx][y_idx]);
            }
            fprintf(pelagic->fptr_fish, "%.16f\n", pelagic->mu_fish_values[n_m - 1][x_idx][y_idx]);
        }
    }

    /*Close output files*/
    fclose(pelagic->fptr_r);
    fclose(pelagic->fptr_g);
    fclose(pelagic->fptr_m);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        fclose(pelagic->fptr_pred[pel_species_idx]);
    }
    fclose(pelagic->fptr_fish);
}

void print_timestep_benthic(GRID *grid, BENTHIC *benthic, int time)
/*Outputs timestep values*/
{
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    benthic->fptr_r = safe_fopen(benthic->fname_r, "a", __LINE__, __FILE__);
    benthic->fptr_g = safe_fopen(benthic->fname_g, "a", __LINE__, __FILE__);
    benthic->fptr_m = safe_fopen(benthic->fname_m, "a", __LINE__, __FILE__);
    benthic->fptr_pred = safe_fopen(benthic->fname_pred, "a", __LINE__, __FILE__);
    benthic->fptr_fish = safe_fopen(benthic->fname_fish, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(benthic->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_g, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_m, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_pred, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_fish, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx += (int)rint(grid->moutstep / grid->mstep))
            {
                fprintf(benthic->fptr_r, "%.16f,", benthic->u_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_g, "%.16f,", benthic->g_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_m, "%.16f,", benthic->mu_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_pred, "%.16f,", benthic->mu_pred_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_fish, "%.16f,", benthic->mu_fish_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(benthic->fptr_r, "%.16f\n", benthic->u_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_g, "%.16f\n", benthic->g_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_m, "%.16f\n", benthic->mu_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_pred, "%.16f\n", benthic->mu_pred_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_fish, "%.16f\n", benthic->mu_fish_values[n_m - 1][x_idx][y_idx]);
        }
    }

    /*Close output files*/
    fclose(benthic->fptr_r);
    fclose(benthic->fptr_g);
    fclose(benthic->fptr_m);
    fclose(benthic->fptr_pred);
    fclose(benthic->fptr_fish);
}
void print_detailed_benthic(GRID *grid, BENTHIC *benthic, int time)
/*Outputs timestep values*/
{
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    benthic->fptr_r = safe_fopen(benthic->fname_r, "a", __LINE__, __FILE__);
    benthic->fptr_g = safe_fopen(benthic->fname_g, "a", __LINE__, __FILE__);
    benthic->fptr_m = safe_fopen(benthic->fname_m, "a", __LINE__, __FILE__);
    benthic->fptr_pred = safe_fopen(benthic->fname_pred, "a", __LINE__, __FILE__);
    benthic->fptr_fish = safe_fopen(benthic->fname_fish, "a", __LINE__, __FILE__);

    /*Print Header*/
    print_detailed_header(grid, benthic->fptr_r);
    print_detailed_header(grid, benthic->fptr_g);
    print_detailed_header(grid, benthic->fptr_m);
    print_detailed_header(grid, benthic->fptr_pred);
    print_detailed_header(grid, benthic->fptr_fish);

    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            fprintf(benthic->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_g, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_m, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_pred, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_fish, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx++)
            {
                fprintf(benthic->fptr_r, "%.16f,", benthic->u_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_g, "%.16f,", benthic->g_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_m, "%.16f,", benthic->mu_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_pred, "%.16f,", benthic->mu_pred_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_fish, "%.16f,", benthic->mu_fish_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(benthic->fptr_r, "%.16f\n", benthic->u_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_g, "%.16f\n", benthic->g_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_m, "%.16f\n", benthic->mu_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_pred, "%.16f\n", benthic->mu_pred_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_fish, "%.16f\n", benthic->mu_fish_values[n_m - 1][x_idx][y_idx]);
        }
    }

    /*Close output files*/
    fclose(benthic->fptr_r);
    fclose(benthic->fptr_g);
    fclose(benthic->fptr_m);
    fclose(benthic->fptr_pred);
    fclose(benthic->fptr_fish);
}

void print_plankton_summary(GRID *grid, PLANKTON *plankton, int time)
{

    int sc_idx, x_idx, y_idx;
    int n_x, n_y;
    n_x = grid->xnum;
    n_y = grid->ynum;

    double temp;

    /*Open output files*/
    plankton->fptr_summ = safe_fopen(plankton->fname_summ, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(plankton->fptr_summ, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);

            /*Plankton Biomass Totals*/
            temp = 0;
            for (sc_idx = plankton->iplamin; sc_idx < (plankton->iplamax + 1); sc_idx++)
            {
                temp += plankton->u_values[sc_idx][x_idx][y_idx] * exp(grid->m_values[sc_idx]) * grid->mstep;
            }
            fprintf(plankton->fptr_summ, "%.16f", temp);
            fprintf(plankton->fptr_summ, "\n");
        }
    }

    fclose(plankton->fptr_summ);
}

void print_pelagic_summary(GRID *grid, PELAGIC *pelagic, int time)
{

    int sc_idx, x_idx, y_idx;
    int n_x, n_y;
    n_x = grid->xnum;
    n_y = grid->ynum;

    double temp;

    /*Open output files*/
    pelagic->fptr_summ = safe_fopen(pelagic->fname_summ, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(pelagic->fptr_summ, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);

            /*Pelagic Biomass Totals*/
            temp = 0;
            for (sc_idx = pelagic->ipelmin; sc_idx < (pelagic->ipelmax + 1); sc_idx++)
            {
                temp += pelagic->u_values[sc_idx][x_idx][y_idx] * exp(grid->m_values[sc_idx]) * grid->mstep;
            }
            fprintf(pelagic->fptr_summ, ",%.16f", temp);

            /*Pelagic Biomass Intake Rate*/
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->pla_total[x_idx][y_idx]);
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->pel_total[x_idx][y_idx]);
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->ben_total[x_idx][y_idx]);

            /*Biomass Lost to Predation*/
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->pred_total[x_idx][y_idx]);

            /*Biomass Lost to Fishing*/
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->fish_total[x_idx][y_idx]);

            /*Reproduction*/
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->reproduction[x_idx][y_idx] / grid->tstep);
            fprintf(pelagic->fptr_summ, "\n");
        }
    }

    fclose(pelagic->fptr_summ);
}

void print_benthic_summary(GRID *grid, BENTHIC *benthic, int time)
{
    int sc_idx, x_idx, y_idx;
    int n_x, n_y;
    n_x = grid->xnum;
    n_y = grid->ynum;

    double temp;

    /*Open output files*/
    benthic->fptr_summ = safe_fopen(benthic->fname_summ, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(benthic->fptr_summ, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);

            /*Benthic Biomass Totals*/
            temp = 0;
            for (sc_idx = benthic->ibenmin; sc_idx < (benthic->ibenmax + 1); sc_idx++)
            {
                temp += benthic->u_values[sc_idx][x_idx][y_idx] * exp(grid->m_values[sc_idx]) * grid->mstep;
            }
            fprintf(benthic->fptr_summ, ",%.16f", temp);

            /*Benthic Biomass Intake*/
            fprintf(benthic->fptr_summ, ",%.16f", benthic->det_total[x_idx][y_idx]);

            /*Biomass Lost to Predation*/
            fprintf(benthic->fptr_summ, ",%.16f", benthic->pred_total[x_idx][y_idx]);

            /*Biomass Lost to Fishing*/
            fprintf(benthic->fptr_summ, ",%.16f", benthic->fish_total[x_idx][y_idx]);

            /*Reproduction*/
            fprintf(benthic->fptr_summ, ",%.16f", benthic->reproduction[x_idx][y_idx] / grid->tstep);
            fprintf(benthic->fptr_summ, "\n");
        }
    }
    fclose(benthic->fptr_summ);
}

void print_detritus_summary(GRID *grid, DETRITUS *detritus, int time)
{
    int x_idx, y_idx;
    int n_x, n_y;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    detritus->fptr_summ = safe_fopen(detritus->fname_summ, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(detritus->fptr_summ, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);

            /*Detritus Biomass Totals*/
            fprintf(detritus->fptr_summ, ",%.16f", detritus->w_values[x_idx][y_idx]);

            /*Detritus In*/
            fprintf(detritus->fptr_summ, ",%.16f", detritus->g_values[x_idx][y_idx]);

            /*Detritus Out*/
            fprintf(detritus->fptr_summ, ",%.16f", detritus->mu_values[x_idx][y_idx]);

            fprintf(detritus->fptr_summ, "\n");
        }
    }

    /*Close output files*/
    fclose(detritus->fptr_summ);
}

//-------------------------------------------------//
// Differencing Scheme Functions. Not to be edited //
//-------------------------------------------------//

void mass_solver(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *pelmatrix, MATRIX *benmatrix,
                 double ***temp_plankton, double ****temp_pelagic, double ****temp_benthic, double **temp_detritus)
/*Impliments the differencing scheme to obtain the values at the next time step*/
{
    int sc_idx, x_idx, y_idx, pel_species_idx, ben_species_idx;
    int n_x, n_y, n_pel_species, n_ben_species;
    double dt, dm;

    dm = grid->mstep;
    dt = grid->tstep;

    n_x = grid->xnum;
    n_y = grid->ynum;
    n_pel_species = run->no_pelagic;
    n_ben_species = run->no_benthic;

    /*Advance the environmental-forcing plugin one timestep (if active), in lockstep with the
      plankton re-read below. Parse this timestep's row into the recognised channels, then:
        - rescale every species' feeding (A) and background mortality (mu_0) from their stored
          base values by pel/ben_tempeff - the growth, predation-mortality and reproduction
          integrals computed in calculate_g_and_mu() (after this solve) then carry the
          temperature effect, matching sizemodel(): feeding & other-mortality are temperature-
          scaled; senescence (mu_s) and fishing are not;
        - store sinking_rate on run for use in g_det() (detritus input).*/
    if (run->forcing_flag == 1 && run->forcing_percell)
    {
        /* PER-CELL forcing: read x*y rows this timestep (one per cell, x-outer y-inner) into the
           per-cell temperature/sinking arrays. A and mu_0 are held at BASE here; the per-cell
           temperature is applied in the feeding/mortality integrals (via run->*_tempeff_c[x_idx][l]).
           Only pel_tempeff/ben_tempeff/sinking_rate are made per-cell in this mode. */
        char fline[1024];
        double vals[64];
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                if (fgets(fline, sizeof(fline), run->fptr_forcing) == NULL)
                    break;
                int nv = 0;
                char *tok = strtok(fline, ",");
                while (tok != NULL && nv < 64)
                {
                    vals[nv++] = atof(tok);
                    tok = strtok(NULL, ",");
                }
                run->pel_tempeff_c[x_idx][y_idx] =
                    (run->col_pel_tempeff >= 0 && run->col_pel_tempeff < nv) ? vals[run->col_pel_tempeff] : 1.0;
                run->ben_tempeff_c[x_idx][y_idx] =
                    (run->col_ben_tempeff >= 0 && run->col_ben_tempeff < nv) ? vals[run->col_ben_tempeff] : 1.0;
                run->sinking_rate_c[x_idx][y_idx] =
                    (run->col_sinking_rate >= 0 && run->col_sinking_rate < nv) ? vals[run->col_sinking_rate] : 1.0;
                run->depth_c[x_idx][y_idx] = (run->col_depth >= 0 && run->col_depth < nv) ? vals[run->col_depth] : 1.0;
                run->export_attn_c[x_idx][y_idx] =
                    (run->col_export_attn >= 0 && run->col_export_attn < nv) ? vals[run->col_export_attn] : 0.0;
                run->export_gamma_c[x_idx][y_idx] =
                    (run->col_export_gamma >= 0 && run->col_export_gamma < nv) ? vals[run->col_export_gamma] : 0.0;
                run->export_magg_c[x_idx][y_idx] =
                    (run->col_export_magg >= 0 && run->col_export_magg < nv) ? vals[run->col_export_magg] : -1e30;
            }
        }
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            community->pelagic[pel_species_idx].A = community->pelagic[pel_species_idx].A_base;
            community->pelagic[pel_species_idx].mu_0 = community->pelagic[pel_species_idx].mu_0_base;
        }
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            community->benthic[ben_species_idx].A = community->benthic[ben_species_idx].A_base;
            community->benthic[ben_species_idx].mu_0 = community->benthic[ben_species_idx].mu_0_base;
        }
    }
    else if (run->forcing_flag == 1)
    {
        char fline[1024];
        if (fgets(fline, sizeof(fline), run->fptr_forcing) != NULL)
        {
            double vals[64];
            int nv = 0;
            char *tok = strtok(fline, ",");
            while (tok != NULL && nv < 64)
            {
                vals[nv++] = atof(tok);
                tok = strtok(NULL, ",");
            }
            if (run->col_pel_tempeff >= 0 && run->col_pel_tempeff < nv)
                run->pel_tempeff = vals[run->col_pel_tempeff];
            if (run->col_ben_tempeff >= 0 && run->col_ben_tempeff < nv)
                run->ben_tempeff = vals[run->col_ben_tempeff];
            if (run->col_sinking_rate >= 0 && run->col_sinking_rate < nv)
                run->sinking_rate = vals[run->col_sinking_rate];
            if (run->col_depth >= 0 && run->col_depth < nv)
                run->depth = vals[run->col_depth];
            if (run->col_residence_time >= 0 && run->col_residence_time < nv)
                run->residence_time = vals[run->col_residence_time];
            if (run->col_export_attn >= 0 && run->col_export_attn < nv)
                run->export_attn = vals[run->col_export_attn];
            if (run->col_export_gamma >= 0 && run->col_export_gamma < nv)
                run->export_gamma = vals[run->col_export_gamma];
            if (run->col_export_magg >= 0 && run->col_export_magg < nv)
                run->export_magg = vals[run->col_export_magg];
            if (run->col_gravity >= 0 && run->col_gravity < nv)
                run->gravity = vals[run->col_gravity];
        }
        /* TEMP_ON_MORT (env, default 1): whether temperature scales background mortality.
           =0 -> temperature drives FEEDING only; mu_0 held at base. Rationale: the model does
           not resolve metabolic maintenance costs (respiration) or their temperature response,
           so mu_0*tempeff is a lone, unbalanced "cost of warmth"; feeding has a clean
           encounter-rate justification, mortality does not. Read once. */
        static int tom = -1;
        if (tom < 0)
        {
            char *e = getenv("TEMP_ON_MORT");
            tom = e ? atoi(e) : 1;
        }
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            community->pelagic[pel_species_idx].A = community->pelagic[pel_species_idx].A_base * run->pel_tempeff;
            community->pelagic[pel_species_idx].mu_0 =
                tom ? community->pelagic[pel_species_idx].mu_0_base * run->pel_tempeff
                    : community->pelagic[pel_species_idx].mu_0_base;
        }
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            community->benthic[ben_species_idx].A = community->benthic[ben_species_idx].A_base * run->ben_tempeff;
            community->benthic[ben_species_idx].mu_0 =
                tom ? community->benthic[ben_species_idx].mu_0_base * run->ben_tempeff
                    : community->benthic[ben_species_idx].mu_0_base;
        }
    }

    /*Used for calculating detritus updates*/
    double **temp_det_g;
    double **temp_det_mu;

    temp_det_g = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    temp_det_mu = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        temp_det_g[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        temp_det_mu[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            temp_det_g[x_idx][y_idx] = 0;
            temp_det_mu[x_idx][y_idx] = 0;
        }
    }
    /*******************/
    /* Plankton System */
    /*******************/
    if (community->plankton->ts_flag == 0)
    {
        // Calculate values for Plankton Size Spectra using internal dynamics //
        for (sc_idx = community->plankton->iplamin; sc_idx < (community->plankton->iplamax + 1); sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    community->plankton->u_values[sc_idx][x_idx][y_idx] =
                        temp_plankton[sc_idx][x_idx][y_idx]; // All plankton systems are invariant
                }
            }
        }
    }
    if (community->plankton->ts_flag == 1)
    {
        int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
        char *input;                   // input string to store line as string
        input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
        char sep[] = ",";    // identifying the ',' char as seperating values
        char *result = NULL; // declaring a pointer for storing each individual value as a char

        /*Read in a single timestep*/
        /* each line one at a time (corresponding to spectra at a spatial point)*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                fgets(input, MAX_LEN, community->plankton->fptr_ts);
                input[strlen(input) - 1] =
                    '\0'; // replacing the newline \n char of the line of text with a null \0 char

                result = (char *)strtok(input, sep); // this is the 1st value
                sc_idx = 0;
                while (result != NULL)
                {
                    community->plankton->u_values[sc_idx][x_idx][y_idx] = atof(result);
                    result = (char *)strtok(NULL, sep);
                    sc_idx++;
                }
            }
        }
        free(input);
    }

    /******************/
    /* Pelagic System */
    /******************/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (community->pelagic[pel_species_idx].ts_flag == 0)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    for (sc_idx = 0; sc_idx < pelmatrix[pel_species_idx].size; sc_idx++)
                    {
                        pelmatrix[pel_species_idx].a[sc_idx] =
                            -(dt / dm) *
                            community->pelagic[pel_species_idx]
                                .g_values[sc_idx - 1 + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx];
                        pelmatrix[pel_species_idx].b[sc_idx] =
                            (1 +
                             (dt * community->pelagic[pel_species_idx]
                                       .mu_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx]) +
                             (dt / dm) *
                                 community->pelagic[pel_species_idx]
                                     .g_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx]);
                        pelmatrix[pel_species_idx].c[sc_idx] = 0;
                        pelmatrix[pel_species_idx].r[sc_idx] =
                            temp_pelagic[pel_species_idx][sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx]
                                        [y_idx];
                    }
                    pelmatrix[pel_species_idx].a[0] = 0;
                    tridag(&(pelmatrix[pel_species_idx]));
                    for (sc_idx = 0; sc_idx < pelmatrix[pel_species_idx].size; sc_idx++)
                    {
                        community->pelagic[pel_species_idx]
                            .u_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx] =
                            pelmatrix[pel_species_idx].u[sc_idx];
                    }
                }
            }
        }
        if (community->pelagic[pel_species_idx].ts_flag == 1)
        {
            int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
            char *input;                   // input string to store line as string
            input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
            char sep[] = ",";    // identifying the ',' char as seperating values
            char *result = NULL; // declaring a pointer for storing each individual value as a char

            /*Read in a single timestep*/
            /* each line one at a time (corresponding to spectra at a spatial point)*/
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    fgets(input, MAX_LEN, community->pelagic[pel_species_idx].fptr_ts);
                    input[strlen(input) - 1] =
                        '\0'; // replacing the newline \n char of the line of text with a null \0 char

                    result = (char *)strtok(input, sep); // this is the 1st value
                    sc_idx = 0;
                    while (result != NULL)
                    {
                        community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] = atof(result);
                        result = (char *)strtok(NULL, sep);
                        sc_idx++;
                    }
                }
            }
            free(input);
        }
    }

    /******************/
    /* Benthic System */
    /******************/
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            if (community->benthic[ben_species_idx].ts_flag == 0)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        for (sc_idx = 0; sc_idx < benmatrix[ben_species_idx].size; sc_idx++)
                        {
                            benmatrix[ben_species_idx].a[sc_idx] =
                                -(dt / dm) *
                                community->benthic[ben_species_idx]
                                    .g_values[sc_idx - 1 + community->benthic[ben_species_idx].ibenmin][x_idx][y_idx];
                            benmatrix[ben_species_idx].b[sc_idx] =
                                (1 +
                                 (dt *
                                  community->benthic[ben_species_idx]
                                      .mu_values[sc_idx + community->benthic[ben_species_idx].ibenmin][x_idx][y_idx]) +
                                 (dt / dm) *
                                     community->benthic[ben_species_idx]
                                         .g_values[sc_idx + community->benthic[ben_species_idx].ibenmin][x_idx][y_idx]);
                            benmatrix[ben_species_idx].c[sc_idx] = 0;
                            benmatrix[ben_species_idx].r[sc_idx] =
                                temp_benthic[ben_species_idx][sc_idx + community->benthic[ben_species_idx].ibenmin]
                                            [x_idx][y_idx];
                        }
                        benmatrix[ben_species_idx].a[0] = 0;
                        tridag(&(benmatrix[ben_species_idx]));
                        for (sc_idx = 0; sc_idx < benmatrix[ben_species_idx].size; sc_idx++)
                        {
                            community->benthic[ben_species_idx]
                                .u_values[sc_idx + community->benthic[ben_species_idx].ibenmin][x_idx][y_idx] =
                                benmatrix[ben_species_idx].u[sc_idx];
                        }
                    }
                }
            }
            if (community->benthic[ben_species_idx].ts_flag == 1)
            {
                int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
                char *input;                   // input string to store line as string
                input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
                char sep[] = ",";    // identifying the ',' char as seperating values
                char *result = NULL; // declaring a pointer for storing each individual value as a char

                /*Read in a single timestep*/
                /* each line one at a time (corresponding to spectra at a spatial point)*/
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        fgets(input, MAX_LEN, community->benthic[ben_species_idx].fptr_ts);
                        input[strlen(input) - 1] =
                            '\0'; // replacing the newline \n char of the line of text with a null \0 char

                        result = (char *)strtok(input, sep); // this is the 1st value
                        sc_idx = 0;
                        while (result != NULL)
                        {
                            community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] = atof(result);
                            result = (char *)strtok(NULL, sep);
                            sc_idx++;
                        }
                    }
                }
                free(input);
            }
        }

        /*******************/
        /* Detritus System */
        /*******************/
        if (community->detritus->ts_flag == 0)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    /* Analytical (exponential) integrator for the detritus pool.
                       Over a step the fish/plankton state is fixed, so detritus obeys a
                       LINEAR ODE  dW/dt = I - kappa*W , with input I = g_values (constant)
                       and loss mu_values = kappa*W (detritivory is proportional to W).
                       The exact solution  W(t+dt) = W* + (W0-W*)*exp(-kappa*dt),
                       W* = I/kappa,  is UNCONDITIONALLY STABLE and provably non-negative
                       (a convex combination of W0>=0 and W*>=0) - unlike the explicit
                       Euler step W += dt*(I-kappa*W), which overshoots to negative when
                       kappa*dt > 1 and drives the NaN/collapse cascade. Removes the
                       timestep stability limit on the detritus pool. */
                    {
                        double W0 = temp_detritus[x_idx][y_idx];
                        double In = community->detritus->g_values[x_idx][y_idx]; /* gross input rate I */
                        /* Dunne et al. (2007) eq.3 burial: a fraction of the incoming detrital
                           flux is buried (permanently lost) before joining the pool, EXACTLY as
                           sizemodel does: burial = input_w*(0.013+0.53*iw^2/(7+iw)^2) applied to
                           the model's native input_w. dbpmr's densities (like sizemodel's) are
                           per-VOLUME - anchored to the depth-averaged plankton concentration - so
                           input_w is a per-volume flux and the Dunne fraction is applied to it
                           DIRECTLY, with no areal (x depth) rescale, matching sizemodel. Burial is
                           OPT-IN: only when a `depth` forcing channel is supplied (col_depth>=0),
                           so existing runs without the forcing plugin keep their original dynamics. */
                        if (run->col_depth >= 0)
                        {
                            double bfrac = 0.013 + 0.53 * In * In / ((7.0 + In) * (7.0 + In));
                            In = In * (1.0 - bfrac); /* net of burial   */
                        }
                        double loss = community->detritus->mu_values[x_idx][y_idx]; /* = kappa*W0   */
                        double kappa = (W0 > 1e-30) ? loss / W0 : 0.0;
                        /* Alternative/opt-in residence-time closure: a first-order loss W/tau on the
                           detritus POOL (remineralisation/export from the benthic layer), added to the
                           detritivory rate. Dimensionally clean in dbpmr's per-volume units - no areal
                           ambiguity - and the exact integrator below absorbs it (K = kappa + 1/tau, so
                           the equilibrium is W* = In/(kappa+1/tau); with negligible detritivory that is
                           the classic W* = In*tau). Supplied via the residence_time channel (yr); an
                           alternative to Dunne burial (Julia's commented CMIP5 variant made this loss
                           proportional to W). Off when the channel is absent or tau<=0. */
                        if (run->col_residence_time >= 0 && run->residence_time > 0.0)
                        {
                            kappa += 1.0 / run->residence_time;
                        }
                        if (kappa > 1e-30)
                        {
                            double Wstar = In / kappa;
                            community->detritus->w_values[x_idx][y_idx] = Wstar + (W0 - Wstar) * exp(-kappa * dt);
                        }
                        else
                        {
                            community->detritus->w_values[x_idx][y_idx] = W0 + In * dt; /* kappa -> 0 limit */
                        }
                        /* defensive backstop only against a non-finite upstream rate */
                        if (!isfinite(community->detritus->w_values[x_idx][y_idx]) ||
                            community->detritus->w_values[x_idx][y_idx] < 0)
                        {
                            community->detritus->w_values[x_idx][y_idx] = 0;
                        }
                    }
                    temp_det_g[x_idx][y_idx] = community->detritus->g_values[x_idx][y_idx];
                    temp_det_mu[x_idx][y_idx] = community->detritus->mu_values[x_idx][y_idx];
                }
            }
        }
        if (community->detritus->ts_flag == 1)
        {
            int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
            char *input;                   // input string to store line as string
            input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
            char sep[] = ",";    // identifying the ',' char as seperating values
            char *result = NULL; // declaring a pointer for storing each individual value as a char

            /*Read in a single timestep*/
            /* each line one at a time (corresponding to spectra at a spatial point)*/
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    fgets(input, MAX_LEN, community->detritus->fptr_ts);
                    input[strlen(input) - 1] =
                        '\0'; // replacing the newline \n char of the line of text with a null \0 char

                    result = (char *)strtok(input, sep); // this is the 1st value
                    while (result != NULL)
                    {
                        community->detritus->w_values[x_idx][y_idx] = atof(result);
                        result = (char *)strtok(NULL, sep);
                    }
                }
            }
            free(input);
        }
    }

    // Calculate fully implicit method if required
    if (run->diff_method == 1)
    {
        // Use intermediate u values to calculate new g and mu for full implicit methods
        calculate_g_and_mu(run, grid, community);

        // Use new g and mu to claculate final u values
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            if (community->pelagic[pel_species_idx].ts_flag == 0)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        for (sc_idx = 0; sc_idx < pelmatrix[pel_species_idx].size; sc_idx++)
                        {
                            pelmatrix[pel_species_idx].a[sc_idx] =
                                -(dt / dm) *
                                community->pelagic[pel_species_idx]
                                    .g_values[sc_idx - 1 + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx];
                            pelmatrix[pel_species_idx].b[sc_idx] =
                                (1 +
                                 (dt *
                                  community->pelagic[pel_species_idx]
                                      .mu_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx]) +
                                 (dt / dm) *
                                     community->pelagic[pel_species_idx]
                                         .g_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx]);
                            pelmatrix[pel_species_idx].c[sc_idx] = 0;
                            pelmatrix[pel_species_idx].r[sc_idx] =
                                temp_pelagic[pel_species_idx][sc_idx + community->pelagic[pel_species_idx].ipelmin]
                                            [x_idx][y_idx];
                        }
                        pelmatrix[pel_species_idx].a[0] = 0;
                        tridag(&(pelmatrix[pel_species_idx]));
                        for (sc_idx = 0; sc_idx < pelmatrix[pel_species_idx].size; sc_idx++)
                        {
                            community->pelagic[pel_species_idx]
                                .u_values[sc_idx + community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx] =
                                pelmatrix[pel_species_idx].u[sc_idx];
                        }
                    }
                }
            }
        }

        if (run->no_benthic != 0)
        {
            for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
            {
                if (community->benthic[ben_species_idx].ts_flag == 0)
                {
                    for (x_idx = 0; x_idx < n_x; x_idx++)
                    {
                        for (y_idx = 0; y_idx < n_y; y_idx++)
                        {
                            for (sc_idx = 0; sc_idx < benmatrix[ben_species_idx].size; sc_idx++)
                            {
                                benmatrix[ben_species_idx].a[sc_idx] =
                                    -(dt / dm) * community->benthic[ben_species_idx]
                                                     .g_values[sc_idx - 1 + community->benthic[ben_species_idx].ibenmin]
                                                              [x_idx][y_idx];
                                benmatrix[ben_species_idx].b[sc_idx] =
                                    (1 +
                                     (dt * community->benthic[ben_species_idx]
                                               .mu_values[sc_idx + community->benthic[ben_species_idx].ibenmin][x_idx]
                                                         [y_idx]) +
                                     (dt / dm) * community->benthic[ben_species_idx]
                                                     .g_values[sc_idx + community->benthic[ben_species_idx].ibenmin]
                                                              [x_idx][y_idx]);
                                benmatrix[ben_species_idx].c[sc_idx] = 0;
                                benmatrix[ben_species_idx].r[sc_idx] =
                                    temp_benthic[ben_species_idx][sc_idx + community->benthic[ben_species_idx].ibenmin]
                                                [x_idx][y_idx];
                            }
                            benmatrix[ben_species_idx].a[0] = 0;
                            tridag(&(benmatrix[ben_species_idx]));
                            for (sc_idx = 0; sc_idx < benmatrix[ben_species_idx].size; sc_idx++)
                            {
                                community->benthic[ben_species_idx]
                                    .u_values[sc_idx + community->benthic[ben_species_idx].ibenmin][x_idx][y_idx] =
                                    benmatrix[ben_species_idx].u[sc_idx];
                            }
                        }
                    }
                }
            }

            if (community->detritus->ts_flag == 0)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        community->detritus->w_values[x_idx][y_idx] =
                            temp_detritus[x_idx][y_idx] +
                            (dt / 2) * ((community->detritus->g_values[x_idx][y_idx] + temp_det_g[x_idx][y_idx]) -
                                        (community->detritus->mu_values[x_idx][y_idx] + temp_det_mu[x_idx][y_idx]));
                    }
                }
            }
        }
    }

    /***************************/
    /* Add reproduction values */
    /***************************/
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                /*If values aren't given by timeseries*/
                if (community->pelagic[pel_species_idx].ts_flag == 0)
                {
                    /*Specifying Total Eggs*/
                    if (community->pelagic[pel_species_idx].rep_method == 0 ||
                        community->pelagic[pel_species_idx].rep_method == 1)
                    {
                        community->pelagic[pel_species_idx]
                            .u_values[community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx] =
                            community->pelagic[pel_species_idx].reproduction[x_idx][y_idx];
                    }
                    /*Adding Eggs*/
                    if (community->pelagic[pel_species_idx].rep_method == 2 ||
                        community->pelagic[pel_species_idx].rep_method == 3)
                    {
                        community->pelagic[pel_species_idx]
                            .u_values[community->pelagic[pel_species_idx].ipelmin][x_idx][y_idx] +=
                            community->pelagic[pel_species_idx].reproduction[x_idx][y_idx];
                    }
                }
            }

            if (run->no_benthic != 0)
            {
                for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
                {
                    /*If values aren't given by timeseries*/
                    if (community->benthic[ben_species_idx].ts_flag == 0)
                    {
                        /*Constant Egg amount*/
                        if (community->benthic[ben_species_idx].rep_method == 0)
                        {
                            community->benthic[ben_species_idx]
                                .u_values[community->benthic[ben_species_idx].ibenmin][x_idx][y_idx] =
                                community->benthic[ben_species_idx].u_0 *
                                exp(community->benthic[ben_species_idx].lambda *
                                    community->benthic[ben_species_idx].benmin);
                        }
                        /*Eggs produced according to biomass intake*/
                        if (community->benthic[ben_species_idx].rep_method == 1 ||
                            community->benthic[ben_species_idx].rep_method == 2)
                        {
                            community->benthic[ben_species_idx]
                                .u_values[community->benthic[ben_species_idx].ibenmin][x_idx][y_idx] +=
                                community->benthic[ben_species_idx].reproduction[x_idx][y_idx];
                        }
                    }
                }
            }
        }
    }

    /*Free pointers*/
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        free(temp_det_g[x_idx]);
        free(temp_det_mu[x_idx]);
    }
    free(temp_det_g);
    free(temp_det_mu);
}

void xmove_solver(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *xmat, double ****temp_pelagic, int time)
/* Spatial movement in x, operator-split (fix for issue #41 — the old scheme accumulated biomass at the
   x-max wall regardless of source and did not conserve mass):
     (1) ADVECTION in conservative upwind FLUX form, with ZERO FLUX through the domain walls
         (mass-conserving + stable); advective velocity v = C*dg/dx - D*dmu/dx (move up the food
         gradient g, down the mortality gradient mu);
     (2) DIFFUSION by Crank-Nicolson with ZERO-FLUX (Neumann) boundaries.
   With no movement (C=D=0 and Diffun=0) both steps are the identity -> output unchanged. */
{
    int sc_idx, x_idx, y_idx, pel_species_idx;
    int n_x, n_y, n_pel_species;
    double dt, dx, C, D, r;

    dt = grid->tstep;
    dx = grid->xstep;
    n_x = grid->xnum;
    n_y = grid->ynum;
    n_pel_species = run->no_pelagic;

    double *uadv = (double *)safe_malloc(n_x * sizeof(double), __LINE__, __FILE__); /* state after advection */
    double *flux = (double *)safe_malloc(n_x * sizeof(double), __LINE__,
                                         __FILE__); /* flux[x_idx] = flux across the face between cells k and k+1 */
    double *vface = (double *)safe_malloc(n_x * sizeof(double), __LINE__,
                                          __FILE__); /* advective velocity at face k (depends on g,mu only) */

    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        for (sc_idx = community->pelagic[pel_species_idx].ipelmin;
             sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
        {
            C = Cfun(grid->m_values[sc_idx], grid->t_values[time], grid, &(community->pelagic[pel_species_idx]));
            D = Dfun(grid->m_values[sc_idx], grid->t_values[time], grid, &(community->pelagic[pel_species_idx]));
            r = Diffun(grid->m_values[sc_idx], grid->t_values[time], grid, &(community->pelagic[pel_species_idx])) *
                (dt / (dx * dx));
            int noadv = (C == 0.0 &&
                         D == 0.0); /* no advection -> identity; also avoids 0*NaN if a crashed cell left NaN in g/mu */
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                /* --- 1. advection: conservative upwind flux form; zero flux at both walls.
                   CFL SUB-STEPPING: the explicit upwind update is stable only for Courant number
                   |v|*dt/dx <= 1; a large prey/pred coefficient would otherwise flush a cell in one
                   step and (after the non-negativity clamp) silently zero the field. Split dt into
                   nsub steps so the per-sub-step Courant number stays <=1 — strong movement then
                   saturates smoothly instead of exploding. Velocities depend on g,mu only (fixed
                   across the sub-steps); only the upwind donor value is refreshed each sub-step. */
                double vmax = 0.0;
                for (x_idx = 0; x_idx < n_x - 1; x_idx++)
                {
                    vface[x_idx] = noadv
                                       ? 0.0
                                       : (C * (community->pelagic[pel_species_idx].g_values[sc_idx][x_idx + 1][y_idx] -
                                               community->pelagic[pel_species_idx].g_values[sc_idx][x_idx][y_idx]) -
                                          D * (community->pelagic[pel_species_idx].mu_values[sc_idx][x_idx + 1][y_idx] -
                                               community->pelagic[pel_species_idx].mu_values[sc_idx][x_idx][y_idx])) /
                                             dx;
                    if (fabs(vface[x_idx]) > vmax)
                        vmax = fabs(vface[x_idx]);
                }
                for (x_idx = 0; x_idx < n_x; x_idx++)
                    uadv[x_idx] = temp_pelagic[pel_species_idx][sc_idx][x_idx][y_idx];
                int nsub = (int)ceil(vmax * dt / dx);
                if (nsub < 1)
                    nsub = 1;
                double dts = dt / (double)nsub;
                for (int sub = 0; sub < nsub; sub++)
                {
                    for (x_idx = 0; x_idx < n_x - 1; x_idx++)
                    {
                        double uup = (vface[x_idx] > 0.0) ? uadv[x_idx] : uadv[x_idx + 1]; /* upwind donor */
                        flux[x_idx] = vface[x_idx] * uup;
                    }
                    for (x_idx = 0; x_idx < n_x; x_idx++)
                    {
                        double Fr = (x_idx < n_x - 1) ? flux[x_idx] : 0.0;  /* right-face flux (0 at the x-max wall) */
                        double Fl = (x_idx > 0) ? flux[x_idx - 1] : 0.0;    /* left-face  flux (0 at the x=0   wall) */
                        uadv[x_idx] = uadv[x_idx] - (dts / dx) * (Fr - Fl); /* flux[] precomputed -> in-place is safe */
                        if (uadv[x_idx] < 0.0)
                            uadv[x_idx] = 0.0; /* clamp tiny upwind undershoots */
                    }
                }
                /* --- 2. diffusion: Crank-Nicolson with zero-flux (Neumann) BCs; conserves mass --- */
                /* RHS (I + r/2 L): xmat->r = RHS * uadv via trimul */
                xmat->b[0] = 1.0 - r / 2.0;
                xmat->c[0] = r / 2.0;
                for (x_idx = 1; x_idx < n_x - 1; x_idx++)
                {
                    xmat->a[x_idx] = r / 2.0;
                    xmat->b[x_idx] = 1.0 - r;
                    xmat->c[x_idx] = r / 2.0;
                }
                xmat->a[n_x - 1] = r / 2.0;
                xmat->b[n_x - 1] = 1.0 - r / 2.0;
                for (x_idx = 0; x_idx < n_x; x_idx++)
                    xmat->u[x_idx] = uadv[x_idx];
                trimul(xmat);
                /* LHS (I - r/2 L): solve LHS * unew = r */
                xmat->b[0] = 1.0 + r / 2.0;
                xmat->c[0] = -r / 2.0;
                for (x_idx = 1; x_idx < n_x - 1; x_idx++)
                {
                    xmat->a[x_idx] = -r / 2.0;
                    xmat->b[x_idx] = 1.0 + r;
                    xmat->c[x_idx] = -r / 2.0;
                }
                xmat->a[n_x - 1] = -r / 2.0;
                xmat->b[n_x - 1] = 1.0 + r / 2.0;
                tridag(xmat);
                for (x_idx = 0; x_idx < n_x; x_idx++)
                    community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] = xmat->u[x_idx];
            }
        }
    }
    free(uadv);
    free(flux);
    free(vface);
}

void ymove_solver(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *ymat, double ****temp_pelagic, int time)
/* Spatial movement in y, operator-split (same fix as xmove_solver — see issue #41):
   conservative upwind FLUX advection + Crank-Nicolson diffusion, both with ZERO-FLUX walls. */
{
    int sc_idx, x_idx, y_idx, pel_species_idx;
    int n_x, n_y, n_pel_species;
    double dt, dy, C, D, r;

    dt = grid->tstep;
    dy = grid->ystep;
    n_x = grid->xnum;
    n_y = grid->ynum;
    n_pel_species = run->no_pelagic;

    double *uadv = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
    double *flux = (double *)safe_malloc(n_y * sizeof(double), __LINE__,
                                         __FILE__); /* flux[l] = flux across the face between cells l and l+1 */
    double *vface = (double *)safe_malloc(n_y * sizeof(double), __LINE__,
                                          __FILE__); /* advective velocity at face l (depends on g,mu only) */

    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        for (sc_idx = community->pelagic[pel_species_idx].ipelmin;
             sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
        {
            C = Cfun(grid->m_values[sc_idx], grid->t_values[time], grid, &(community->pelagic[pel_species_idx]));
            D = Dfun(grid->m_values[sc_idx], grid->t_values[time], grid, &(community->pelagic[pel_species_idx]));
            r = Diffun(grid->m_values[sc_idx], grid->t_values[time], grid, &(community->pelagic[pel_species_idx])) *
                (dt / (dy * dy));
            int noadv = (C == 0.0 && D == 0.0); /* no advection -> identity; also avoids 0*NaN if g/mu carry NaN */
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                /* --- 1. advection: conservative upwind flux form; zero flux at both walls.
                   CFL sub-stepped (see xmove_solver for the rationale). --- */
                double vmax = 0.0;
                for (y_idx = 0; y_idx < n_y - 1; y_idx++)
                {
                    vface[y_idx] = noadv
                                       ? 0.0
                                       : (C * (community->pelagic[pel_species_idx].g_values[sc_idx][x_idx][y_idx + 1] -
                                               community->pelagic[pel_species_idx].g_values[sc_idx][x_idx][y_idx]) -
                                          D * (community->pelagic[pel_species_idx].mu_values[sc_idx][x_idx][y_idx + 1] -
                                               community->pelagic[pel_species_idx].mu_values[sc_idx][x_idx][y_idx])) /
                                             dy;
                    if (fabs(vface[y_idx]) > vmax)
                        vmax = fabs(vface[y_idx]);
                }
                for (y_idx = 0; y_idx < n_y; y_idx++)
                    uadv[y_idx] = temp_pelagic[pel_species_idx][sc_idx][x_idx][y_idx];
                int nsub = (int)ceil(vmax * dt / dy);
                if (nsub < 1)
                    nsub = 1;
                double dts = dt / (double)nsub;
                for (int sub = 0; sub < nsub; sub++)
                {
                    for (y_idx = 0; y_idx < n_y - 1; y_idx++)
                    {
                        double uup = (vface[y_idx] > 0.0) ? uadv[y_idx] : uadv[y_idx + 1]; /* upwind donor */
                        flux[y_idx] = vface[y_idx] * uup;
                    }
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        double Fr = (y_idx < n_y - 1) ? flux[y_idx] : 0.0;  /* right-face flux (0 at the y-max wall) */
                        double Fl = (y_idx > 0) ? flux[y_idx - 1] : 0.0;    /* left-face  flux (0 at the y=0   wall) */
                        uadv[y_idx] = uadv[y_idx] - (dts / dy) * (Fr - Fl); /* flux[] precomputed -> in-place is safe */
                        if (uadv[y_idx] < 0.0)
                            uadv[y_idx] = 0.0;
                    }
                }
                /* --- 2. diffusion: Crank-Nicolson with zero-flux (Neumann) BCs --- */
                ymat->b[0] = 1.0 - r / 2.0;
                ymat->c[0] = r / 2.0;
                for (y_idx = 1; y_idx < n_y - 1; y_idx++)
                {
                    ymat->a[y_idx] = r / 2.0;
                    ymat->b[y_idx] = 1.0 - r;
                    ymat->c[y_idx] = r / 2.0;
                }
                ymat->a[n_y - 1] = r / 2.0;
                ymat->b[n_y - 1] = 1.0 - r / 2.0;
                for (y_idx = 0; y_idx < n_y; y_idx++)
                    ymat->u[y_idx] = uadv[y_idx];
                trimul(ymat);
                ymat->b[0] = 1.0 + r / 2.0;
                ymat->c[0] = -r / 2.0;
                for (y_idx = 1; y_idx < n_y - 1; y_idx++)
                {
                    ymat->a[y_idx] = -r / 2.0;
                    ymat->b[y_idx] = 1.0 + r;
                    ymat->c[y_idx] = -r / 2.0;
                }
                ymat->a[n_y - 1] = -r / 2.0;
                ymat->b[n_y - 1] = 1.0 + r / 2.0;
                tridag(ymat);
                for (y_idx = 0; y_idx < n_y; y_idx++)
                    community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] = ymat->u[y_idx];
            }
        }
    }
    free(uadv);
    free(flux);
    free(vface);
}

void trimul(MATRIX *matrix)
{
    int i, n;
    n = matrix->size;

    matrix->r[0] = (matrix->b[0] * matrix->u[0]) + (matrix->c[0] * matrix->u[1]);
    for (i = 1; i < n - 1; i++)
    {
        matrix->r[i] =
            (matrix->a[i] * matrix->u[i - 1]) + (matrix->b[i] * matrix->u[i]) + (matrix->c[i] * matrix->u[i + 1]);
    }
    matrix->r[n - 1] = (matrix->a[n - 1] * matrix->u[n - 2]) + (matrix->b[n - 1] * matrix->u[n - 1]);
}

void tridag(MATRIX *matrix)
/*Fast algorithm for inverting a tridiagonal matrix*/
{
    int j, n;
    n = matrix->size;
    double bet;
    double *gam;

    gam = (double *)safe_malloc(n * sizeof(double), __LINE__, __FILE__);
    if (matrix->b[0] == 0)
    {
        error("error1");
    }
    bet = matrix->b[0];
    matrix->u[0] = matrix->r[0] / bet;
    for (j = 1; j <= (n - 1); j++)
    {
        gam[j] = matrix->c[j - 1] / bet;
        bet = matrix->b[j] - (matrix->a[j] * gam[j]);
        if (matrix->b[j] == 0)
        {
            error("error2");
        }
        matrix->u[j] = (matrix->r[j] - (matrix->a[j] * matrix->u[j - 1])) / bet;
    }
    for (j = (n - 2); j >= 0; j--)
    {
        matrix->u[j] -= (gam[j + 1] * matrix->u[j + 1]);
    }
    free(gam);
}

//--------------------------------------------//
//      Growth and Mortality Functions        //
// These contain functions that may be edited //
//--------------------------------------------//

void calculate_g_and_mu(RUN *run, GRID *grid, COMMUNITY *community)
{
    int sc_idx, x_idx, y_idx, pel_species_idx, pel_sub_idx, ben_species_idx;
    int n_x, n_y, n_pel_species, n_ben_species;

    n_x = grid->xnum;
    n_y = grid->ynum;
    n_pel_species = run->no_pelagic;
    n_ben_species = run->no_benthic;

    double pla, pel, ben, det, fish, pred;

    /* Reef refuge pre-pass: recompute the per-size availability a(w) for every pelagic species from the
       CURRENT spectrum, BEFORE the feeding/mortality integrals that consume it (so predation mortality
       and predator intake see the SAME a(w) -> mass balance). a(w)=max(min_A, 1 - refuge(w)/competitor(w)),
       competitor = U abundance density. Species with refuge off keep a=1 (set at allocation) -> no effect. */
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (!community->pelagic[pel_species_idx].refuge_on)
            continue;
        double minA = community->pelagic[pel_species_idx].refuge_min_a;
        for (sc_idx = 0; sc_idx < grid->mnum; sc_idx++)
        {
            double rf = community->pelagic[pel_species_idx].refuge_values[sc_idx];
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    double comp = community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx];
                    double a = (rf > 0 && comp > 0) ? 1.0 - rf / comp : 1.0;
                    if (a < minA)
                        a = minA;
                    if (a > 1.0)
                        a = 1.0;
                    if (a < 0.0)
                        a = 0.0;
                    community->pelagic[pel_species_idx].a_refuge[sc_idx][x_idx][y_idx] = a;
                }
            }
        }
    }

    // Pelagic changes
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        for (sc_idx = community->pelagic[pel_species_idx].ipelmin;
             sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    /* per-cell pelagic temperature multiplier on FEEDING + PREDATION (1 when off). */
                    double tmP = run->forcing_percell ? run->pel_tempeff_c[x_idx][y_idx] : 1.0;
                    for (pel_sub_idx = 0; pel_sub_idx < n_pel_species; pel_sub_idx++)
                    {
                        community->pelagic[pel_species_idx].mu_pred_values[pel_sub_idx][sc_idx][x_idx][y_idx] =
                            tmP * mu_pel_pred(pel_species_idx, pel_sub_idx, sc_idx, x_idx, y_idx, run, grid, community);
                    }
                    community->pelagic[pel_species_idx].mu_values[sc_idx][x_idx][y_idx] =
                        mu_pel(pel_species_idx, sc_idx, x_idx, y_idx, run, grid, community);

                    community->pelagic[pel_species_idx].pla_bio[sc_idx][x_idx][y_idx] =
                        tmP * pla_biomass(pel_species_idx, sc_idx, x_idx, y_idx, grid, community);
                    community->pelagic[pel_species_idx].pel_bio[sc_idx][x_idx][y_idx] =
                        tmP * pel_biomass(pel_species_idx, sc_idx, x_idx, y_idx, run, grid, community);
                    community->pelagic[pel_species_idx].ben_bio[sc_idx][x_idx][y_idx] =
                        tmP * ben_biomass(pel_species_idx, sc_idx, x_idx, y_idx, run, grid, community);
                    community->pelagic[pel_species_idx].g_values[sc_idx][x_idx][y_idx] =
                        g_pel(pel_species_idx, sc_idx, x_idx, y_idx, run, grid, community);
                }
            }
        }
    }

    // Calculate pelagic biomass totals
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                pla = 0;
                pel = 0;
                ben = 0;
                fish = 0;
                pred = 0;
                for (sc_idx = community->pelagic[pel_species_idx].ipelmin;
                     sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
                {
                    pla += community->pelagic[pel_species_idx].pla_bio[sc_idx][x_idx][y_idx] *
                           community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] * grid->mstep;
                    pel += community->pelagic[pel_species_idx].pel_bio[sc_idx][x_idx][y_idx] *
                           community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] * grid->mstep;
                    ben += community->pelagic[pel_species_idx].ben_bio[sc_idx][x_idx][y_idx] *
                           community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] * grid->mstep;
                    fish += community->pelagic[pel_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] *
                            community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] *
                            exp(grid->m_values[sc_idx]) * grid->mstep;
                    for (pel_sub_idx = 0; pel_sub_idx < n_pel_species; pel_sub_idx++)
                    {
                        pred += community->pelagic[pel_species_idx].mu_pred_values[pel_sub_idx][sc_idx][x_idx][y_idx] *
                                community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] *
                                exp(grid->m_values[sc_idx]) * grid->mstep;
                    }
                }
                community->pelagic[pel_species_idx].pla_total[x_idx][y_idx] = pla;
                community->pelagic[pel_species_idx].pel_total[x_idx][y_idx] = pel;
                community->pelagic[pel_species_idx].ben_total[x_idx][y_idx] = ben;
                community->pelagic[pel_species_idx].fish_total[x_idx][y_idx] = fish;
                community->pelagic[pel_species_idx].pred_total[x_idx][y_idx] = pred;
            }
        }
    }

    // Benthic changes
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            for (sc_idx = community->benthic[ben_species_idx].ibenmin;
                 sc_idx < (community->benthic[ben_species_idx].ibenmax + 1); sc_idx++)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        /* predation ON benthos is BY pelagic predators (pel temp); benthic feeding uses benthic temp */
                        double tmPb = run->forcing_percell ? run->pel_tempeff_c[x_idx][y_idx] : 1.0;
                        double tmB = run->forcing_percell ? run->ben_tempeff_c[x_idx][y_idx] : 1.0;
                        community->benthic[ben_species_idx].mu_pred_values[sc_idx][x_idx][y_idx] =
                            tmPb * mu_ben_pred(ben_species_idx, sc_idx, x_idx, y_idx, run, grid, community);
                        community->benthic[ben_species_idx].mu_values[sc_idx][x_idx][y_idx] =
                            mu_ben(ben_species_idx, sc_idx, x_idx, y_idx, run, grid, community);

                        community->benthic[ben_species_idx].det_bio[sc_idx][x_idx][y_idx] =
                            tmB * det_biomass(ben_species_idx, sc_idx, x_idx, y_idx, run, grid, community);
                        community->benthic[ben_species_idx].g_values[sc_idx][x_idx][y_idx] =
                            g_ben(ben_species_idx, sc_idx, x_idx, y_idx, run, grid, community);
                    }
                }
            }
        }
        // Calculate benthic biomass totals
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    det = 0;
                    fish = 0;
                    pred = 0;
                    for (sc_idx = community->benthic[ben_species_idx].ibenmin;
                         sc_idx < (community->benthic[ben_species_idx].ibenmax + 1); sc_idx++)
                    {
                        det += community->benthic[ben_species_idx].det_bio[sc_idx][x_idx][y_idx] *
                               community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] * grid->mstep;
                        fish += community->benthic[ben_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] *
                                community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] *
                                exp(grid->m_values[sc_idx]) * grid->mstep;
                        pred += community->benthic[ben_species_idx].mu_pred_values[sc_idx][x_idx][y_idx] *
                                community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] *
                                exp(grid->m_values[sc_idx]) * grid->mstep;
                    }
                    community->benthic[ben_species_idx].det_total[x_idx][y_idx] = det;
                    community->benthic[ben_species_idx].fish_total[x_idx][y_idx] = fish;
                    community->benthic[ben_species_idx].pred_total[x_idx][y_idx] = pred;
                }
            }
        }
        // Detritus Changes
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                community->detritus->g_values[x_idx][y_idx] = g_det(x_idx, y_idx, run, grid, community);
                community->detritus->mu_values[x_idx][y_idx] = mu_det(x_idx, y_idx, run, grid, community);
            }
        }
    }
}

double pla_biomass(int species, int size, int xspace, int yspace, GRID *grid, COMMUNITY *community)
// Calculates biomass intake from plankton
{
    double ans = 0;
    int sc_idx;

    /* Allows feeding from Plankton Size Spectra*/
    for (sc_idx = community->plankton->iplamin; sc_idx < size; sc_idx++)
    {
        ans += exp(grid->m_values[sc_idx]) * community->pelagic[species].phi_values[size - sc_idx] *
               community->plankton->u_values[sc_idx][xspace][yspace] * grid->mstep;
    }
    ans *= community->pelagic[species].pref_pla * community->pelagic[species].A *
           exp(community->pelagic[species].alpha * grid->m_values[size]);
    return ans;
}

double pel_biomass(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
// Calculates biomass intake from Pelagic Size Spectra
{
    double ans = 0;
    int sc_idx, pel_species_idx;

    /* Allows feeding from all Size Spectra*/
    if (run->coupled_flag == 1)
    {
        for (pel_species_idx = 0; pel_species_idx < run->no_pelagic; pel_species_idx++)
        {
            for (sc_idx = community->pelagic[pel_species_idx].ipelmin; sc_idx < size; sc_idx++)
            {
                /*reef refuge: prey of size i are available for intake only with probability a(i) (=1 when off)*/
                ans += community->pelagic[pel_species_idx].a_refuge[sc_idx][xspace][yspace] *
                       exp(grid->m_values[sc_idx]) * community->pelagic[species].phi_values[size - sc_idx] *
                       community->pelagic[pel_species_idx].u_values[sc_idx][xspace][yspace] * grid->mstep;
            }
        }
    }

    /* Allows feeding from own Size Spectra Only*/
    else
    {
        for (sc_idx = community->pelagic[species].ipelmin; sc_idx < size; sc_idx++)
        {
            ans += community->pelagic[species].a_refuge[sc_idx][xspace][yspace] * exp(grid->m_values[sc_idx]) *
                   community->pelagic[species].phi_values[size - sc_idx] *
                   community->pelagic[species].u_values[sc_idx][xspace][yspace] * grid->mstep;
        }
    }

    ans *= community->pelagic[species].pref_pel * community->pelagic[species].A *
           exp(community->pelagic[species].alpha * grid->m_values[size]);
    return ans;
}

double ben_biomass(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
// Calculates biomass intake from Pelagic Size Spectra
{
    double ans = 0;
    int sc_idx, ben_species_idx;

    if (run->no_benthic != 0)
    {
        /* Allows feeding from all Size Spectra*/
        if (run->coupled_flag == 1)
        {
            for (ben_species_idx = 0; ben_species_idx < run->no_benthic; ben_species_idx++)
            {
                for (sc_idx = community->benthic[ben_species_idx].ibenmin; sc_idx < size; sc_idx++)
                {
                    ans += exp(grid->m_values[sc_idx]) * community->pelagic[species].phi_values[size - sc_idx] *
                           community->benthic[ben_species_idx].u_values[sc_idx][xspace][yspace] * grid->mstep;
                }
            }
        }
        ans *= community->pelagic[species].pref_ben * community->pelagic[species].A *
               exp(community->pelagic[species].alpha * grid->m_values[size]);
    }

    return ans;
}

double det_biomass(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
// Calculates biomass intake from Detritus
{
    double ans = 0;
    ans += community->benthic[species].pref_det * community->benthic[species].A *
           exp(community->benthic[species].alpha * grid->m_values[size]) *
           community->detritus->w_values[xspace][yspace];
    return ans;
}

double g_pel(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
// Calculates g for size class 'size' of species 'species' using biomass values
{
    double ans = 0;
    ans += community->pelagic[species].K_pla * exp(-grid->m_values[size]) *
           community->pelagic[species].pla_bio[size][xspace][yspace];
    ans += community->pelagic[species].K_pel * exp(-grid->m_values[size]) *
           community->pelagic[species].pel_bio[size][xspace][yspace];
    ans += community->pelagic[species].K_ben * exp(-grid->m_values[size]) *
           community->pelagic[species].ben_bio[size][xspace][yspace];
    return (ans);
}

double g_ben(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
// Calculates g for size class 'size' of species 'species' using biomass values
{
    double ans = 0;
    ans += community->benthic[species].K_det * exp(-grid->m_values[size]) *
           community->benthic[species].det_bio[size][xspace][yspace];
    return ans;
}

double g_det(int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
{
    int sc_idx, pel_species_idx, ben_species_idx;
    int n_pel_species, n_ben_species;

    n_pel_species = run->no_pelagic;
    n_ben_species = run->no_benthic;

    double ans_surf = 0; /*surface-origin inputs (sink to the seafloor -> scaled by sinking_rate)*/
    double ans_bed = 0;  /*seafloor-origin inputs (dead benthos, already on the bed -> no sinking)*/

    /* Size-dependent export closure (opt-in): when the export_attn channel is supplied, each
       surface particle of ln-mass m reaches the seafloor with fraction
         export(m) = exp( -export_attn * exp(-export_gamma*m) ),
       from an allometric sinking speed w_s ~ mass^export_gamma and an attenuation number
       export_attn = k(T)*z/w0 - so large carcasses/pellets sink faster and are recycled less
       than small ones. Replaces the flat sinking_rate. Default (channel absent): flat export. */
    int size_dep = (run->col_export_attn >= 0);
    /*Per-cell export parameters when FORCING_PERCELL (spatial_dim>=1); otherwise the global values.*/
    double e_attn = run->forcing_percell ? run->export_attn_c[xspace][yspace] : run->export_attn;
    double e_gamma = run->forcing_percell ? run->export_gamma_c[xspace][yspace] : run->export_gamma;
    double e_magg = run->forcing_percell ? run->export_magg_c[xspace][yspace] : run->export_magg;

    if (run->coupled_flag == 1)
    {
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            /*Rate of det biomass in from pelagic defecation (faeces follow the
              temperature-scaled feeding, as sizemodel's defbypred does). Flat export uses
              the pre-summed *_total; size-dependent export reconstructs faeces per predator
              size (pla_bio/pel_bio/ben_bio) so export(m) can weight each size class.*/
            if (!size_dep)
            {
                ans_surf += (community->pelagic[pel_species_idx].Ex_pla *
                                 community->pelagic[pel_species_idx].pla_total[xspace][yspace] +
                             community->pelagic[pel_species_idx].Ex_pel *
                                 community->pelagic[pel_species_idx].pel_total[xspace][yspace] +
                             community->pelagic[pel_species_idx].Ex_ben *
                                 community->pelagic[pel_species_idx].ben_total[xspace][yspace]);
            }
            for (sc_idx = community->pelagic[pel_species_idx].ipelmin;
                 sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
            {
                /*Dead pelagic stuff. Use the BASE (un-temperature-scaled) background
                  mortality: sizemodel removes temperature from the mortality feeding
                  the detritus pool.*/
                double dead =
                    community->pelagic[pel_species_idx].mu_0_base *
                        exp(community->pelagic[pel_species_idx].beta * grid->m_values[sc_idx]) *
                        community->pelagic[pel_species_idx].u_values[sc_idx][xspace][yspace] *
                        exp(grid->m_values[sc_idx]) * grid->mstep +
                    community->pelagic[pel_species_idx].mu_s *
                        ((log10(exp(grid->m_values[sc_idx])) - log10(exp(community->pelagic[pel_species_idx].pelmin))) /
                         ((log10(exp(community->pelagic[pel_species_idx].pelmax)) +
                           community->pelagic[pel_species_idx].epsilon) -
                          log10(exp(grid->m_values[sc_idx])))) *
                        community->pelagic[pel_species_idx].u_values[sc_idx][xspace][yspace] *
                        exp(grid->m_values[sc_idx]) * grid->mstep;
                if (size_dep)
                {
                    double faeces = (community->pelagic[pel_species_idx].Ex_pla *
                                         community->pelagic[pel_species_idx].pla_bio[sc_idx][xspace][yspace] +
                                     community->pelagic[pel_species_idx].Ex_pel *
                                         community->pelagic[pel_species_idx].pel_bio[sc_idx][xspace][yspace] +
                                     community->pelagic[pel_species_idx].Ex_ben *
                                         community->pelagic[pel_species_idx].ben_bio[sc_idx][xspace][yspace]) *
                                    community->pelagic[pel_species_idx].u_values[sc_idx][xspace][yspace] * grid->mstep;
                    double me =
                        fmax(grid->m_values[sc_idx], e_magg); /*aggregation floor: small stuff sinks as marine snow*/
                    ans_surf += exp(-e_attn * exp(-e_gamma * me)) * (faeces + dead);
                }
                else
                {
                    ans_surf += dead;
                }
            }
            /*Dead benthos (detritivore background + senescence mortality). These bodies
              are already on/in the seafloor, so they DO NOT get the sinking_rate factor -
              matching sizemodel's benthic terms outside the sinking_rate multiplier. Uses
              base (un-temperature-scaled) background mortality, as for the pelagic.*/
            for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
            {
                for (sc_idx = community->benthic[ben_species_idx].ibenmin;
                     sc_idx < (community->benthic[ben_species_idx].ibenmax + 1); sc_idx++)
                {
                    ans_bed += community->benthic[ben_species_idx].mu_0_base *
                               exp(community->benthic[ben_species_idx].beta * grid->m_values[sc_idx]) *
                               community->benthic[ben_species_idx].u_values[sc_idx][xspace][yspace] *
                               exp(grid->m_values[sc_idx]) * grid->mstep;
                    ans_bed +=
                        community->benthic[ben_species_idx].mu_s *
                        ((log10(exp(grid->m_values[sc_idx])) - log10(exp(community->benthic[ben_species_idx].benmin))) /
                         ((log10(exp(community->benthic[ben_species_idx].benmax)) +
                           community->benthic[ben_species_idx].epsilon) -
                          log10(exp(grid->m_values[sc_idx])))) *
                        community->benthic[ben_species_idx].u_values[sc_idx][xspace][yspace] *
                        exp(grid->m_values[sc_idx]) * grid->mstep;
                }
            }
        }
    }
    for (sc_idx = community->plankton->iplamin; sc_idx < (community->plankton->iplamax + 1); sc_idx++)
    {
        /*Rate of det biomass in from dead (sinking) plankton - surface origin*/
        double pdead = community->plankton->mu_0 * exp(community->plankton->beta * grid->m_values[sc_idx]) *
                       community->plankton->u_values[sc_idx][xspace][yspace] * exp(grid->m_values[sc_idx]) *
                       grid->mstep;
        if (size_dep)
            ans_surf += exp(-e_attn * exp(-e_gamma * fmax(grid->m_values[sc_idx], e_magg))) * pdead;
        else
            ans_surf += pdead;
    }

    /*Surface-origin inputs sink: scale by the time-varying sinking_rate (= export ratio),
      exactly as sizemodel multiplies its pelagic detritus inputs by sinking_rate[sc_idx]. Seafloor
      (dead-benthos) inputs are added unscaled. Default sinking_rate 1.

      REFERENCE-THICKNESS NOTE (matches sizemodel by design): DBPM densities are per-VOLUME,
      but the pelagic groups are referenced to the productive column (min(depth,200) m) and the
      benthic groups to the benthic-habitat layer (20 m) - the output areal conversions differ
      accordingly (predators x min(depth,200), detritivores x 20). Here the two per-volume
      fluxes are summed into the (20 m-referenced) detritus pool WITHOUT an inter-reference
      conversion, so the sinking pelagic input is effectively under-weighted by ~min(depth,200)/20
      (~8-10x). This is an inherent DBPM simplification, reproduced here identically to sizemodel
      (its coupled input_w mixes the per-volume terms the same way); it is partly absorbed by the
      search-volume calibration. A dimensionally-consistent alternative (areal-referenced sinking
      + areal Dunne burial) would deviate from sizemodel and require recalibration - deliberately
      not done.*/
    /* size-dependent export already applied export(m) per size class inside ans_surf, so it
       replaces the flat sinking_rate; otherwise scale the aggregated surface flux by sinking_rate */
    double sink = run->forcing_percell ? run->sinking_rate_c[xspace][yspace] : run->sinking_rate;
    return ((size_dep ? ans_surf : sink * ans_surf) + ans_bed);
}

double mu_det(int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
{
    int ben_species_idx;
    int n_ben_species;

    n_ben_species = run->no_benthic;

    double ans = 0;

    for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
    {
        /*This is the rate of det biomass out feeding the benthos*/
        ans += (community->benthic[ben_species_idx].det_total[xspace][yspace]);
    }

    return (ans);
}

double mu_pel_pred(int species, int pred, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
// Calculates mu due to predation for size class 'size' of species 'species'
{
    double ans = 0;
    int sc_idx;

    /* Mortality due to Predation by all Size Spectra. Reef refuge: scale each predator i's contribution
       by its availability a(i) -- matching Rogers' death.u, where a.u multiplies the predator (U) vector
       INSIDE the predator integral (mortality kernel mphi = q2 = -q1 is the transpose of the feeding
       kernel gphi = q1, so the same a*(U) construction lands on the PREDATOR index in mortality but the
       PREY index in feeding/pel_biomass). a(i)=1 when refuge off -> identical to standard dbpmr. */
    if (run->coupled_flag == 1)
    {
        for (sc_idx = size; sc_idx < (community->pelagic[pred].ipelmax + 1); sc_idx++)
        {
            ans += community->pelagic[pred].a_refuge[sc_idx][xspace][yspace] *
                   (exp(community->pelagic[pred].alpha * grid->m_values[sc_idx])) *
                   community->pelagic[pred].phi_values[sc_idx - size] *
                   (community->pelagic[pred].u_values[sc_idx][xspace][yspace]) * community->pelagic[pred].A *
                   community->pelagic[pred].pref_pel * grid->mstep;
        }
    }
    /*Mortality due to Own Size Spectra Predation*/
    else
    {
        if (species == pred)
        {
            for (sc_idx = size; sc_idx < (community->pelagic[species].ipelmax + 1); sc_idx++)
            {
                ans += community->pelagic[species].a_refuge[sc_idx][xspace][yspace] *
                       (exp(community->pelagic[species].alpha * grid->m_values[sc_idx])) *
                       community->pelagic[species].phi_values[sc_idx - size] *
                       (community->pelagic[species].u_values[sc_idx][xspace][yspace]) * community->pelagic[species].A *
                       community->pelagic[species].pref_pel * grid->mstep;
            }
        }
    }

    return ans;
}

double mu_pel(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
// Calculates mu for size class 'size' of species 'species'
{
    double ans = 0;
    int pel_species_idx;
    int n_pel_species;
    n_pel_species = run->no_pelagic;

    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        ans += community->pelagic[species].mu_pred_values[pel_species_idx][size][xspace][yspace];
    }
    ans += community->pelagic[species].mu_fish_values[size][xspace][yspace];
    ans += ((run->forcing_percell && run->temp_on_mort ? run->pel_tempeff_c[xspace][yspace] : 1.0) *
            community->pelagic[species].mu_0 * exp(community->pelagic[species].beta * grid->m_values[size]));
    ans += (community->pelagic[species].mu_s *
            ((log10(exp(grid->m_values[size])) - log10(exp(community->pelagic[species].pelmin))) /
             ((log10(exp(community->pelagic[species].pelmax)) + community->pelagic[species].epsilon) -
              log10(exp(grid->m_values[size])))));

    return ans;
}

double mu_ben_pred(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
// Calculates mu due to predation for benthic species
{
    double ans = 0;
    int sc_idx, pel_species_idx;

    /*Mortality due to Predation by all Size Spectra*/
    if (run->coupled_flag == 1)
    {
        for (pel_species_idx = 0; pel_species_idx < run->no_pelagic; pel_species_idx++)
        {
            for (sc_idx = size; sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
            {
                ans = ans + (exp(community->pelagic[pel_species_idx].alpha * grid->m_values[sc_idx])) *
                                community->pelagic[pel_species_idx].phi_values[sc_idx - size] *
                                (community->pelagic[pel_species_idx].u_values[sc_idx][xspace][yspace]) *
                                community->pelagic[pel_species_idx].A * community->pelagic[pel_species_idx].pref_ben *
                                grid->mstep;
            }
        }
    }

    return ans;
}

double mu_ben(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
// Calculates mu for size class 'size' of species 'species'
{
    double ans = 0;

    ans += community->benthic[species].mu_pred_values[size][xspace][yspace];
    ans += community->benthic[species].mu_fish_values[size][xspace][yspace];
    ans += ((run->forcing_percell && run->temp_on_mort ? run->ben_tempeff_c[xspace][yspace] : 1.0) *
            community->benthic[species].mu_0 * exp(community->benthic[species].beta * grid->m_values[size]));
    ans += (community->benthic[species].mu_s *
            ((log10(exp(grid->m_values[size])) - log10(exp(community->benthic[species].benmin))) /
             ((log10(exp(community->benthic[species].benmax)) + community->benthic[species].epsilon) -
              log10(exp(grid->m_values[size])))));

    return ans;
}

void calculate_fishing(RUN *run, GRID *grid, COMMUNITY *community)
{
    int sc_idx, x_idx, y_idx, pel_species_idx, ben_species_idx;
    int n_x, n_y, n_pel_species, n_ben_species;

    n_x = grid->xnum;
    n_y = grid->ynum;
    n_pel_species = run->no_pelagic;
    n_ben_species = run->no_benthic;

    /*Pelagic Fishing*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (community->pelagic[pel_species_idx].fishing_flag == 1)
        {
            int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
            char *input;                   // input string to store line as string
            input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
            char sep[] = ",";    // identifying the ',' char as seperating values
            char *result = NULL; // declaring a pointer for storing each individual value as a char

            /*Read in each line one at a time (corresponding to fishing mortality at a spatial point)*/
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    fgets(input, MAX_LEN, community->pelagic[pel_species_idx].fptr_fish_ts);
                    input[strlen(input) - 1] =
                        '\0'; // replacing the newline \n char of the line of text with a null \0 char

                    result = (char *)strtok(input, sep); // this is the 1st value
                    sc_idx = 0;
                    while (result != NULL)
                    {
                        community->pelagic[pel_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] = atof(result);
                        result = (char *)strtok(NULL, sep);
                        sc_idx++;
                    }
                }
            }
            free(input);
        }
    }

    /*Benthic Fishing*/
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            if (community->benthic[ben_species_idx].fishing_flag == 1)
            {
                int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
                char *input;                   // input string to store line as string
                input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
                char sep[] = ",";    // identifying the ',' char as seperating values
                char *result = NULL; // declaring a pointer for storing each individual value as a char

                /*Read in each line one at a time (corresponding to fishing mortality at a spatial point)*/
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        fgets(input, MAX_LEN, community->benthic[ben_species_idx].fptr_fish_ts);
                        input[strlen(input) - 1] =
                            '\0'; // replacing the newline \n char of the line of text with a null \0 char

                        result = (char *)strtok(input, sep); // this is the 1st value
                        sc_idx = 0;
                        while (result != NULL)
                        {
                            community->benthic[ben_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] = atof(result);
                            result = (char *)strtok(NULL, sep);
                            sc_idx++;
                        }
                    }
                }
                free(input);
            }
        }
    }

    /* DBPM.md gravity fishing (opt-in via the `gravity` channel): replace the static per-group
       catchability split with the true B_{t-1} feedback. For each cell, scale each group's fishing
       mortality by that group's CURRENT areal fishable-biomass share,
         F_g = F_base * A_g/(A_pel+A_ben),  A_g = sum_{fished sizes} u*exp(m)*mstep * thickness
       (pelagic x run->depth = min(depth,200) m from the depth channel; benthic x 20 m). The base F
       read above is the un-split Q*selectivity*effort, so this yields catch_g ~ B_g^2 (effort
       concentrates where biomass is). Runs per cell -> ready for gridded. */
    if (run->col_gravity >= 0)
    {
        double bhd = 20.0; /* benthic habitat depth (sizemodel default) */
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                /* pelagic column thickness (depth channel; 1 if absent) - per-cell under FORCING_PERCELL */
                double dcorr = run->forcing_percell ? run->depth_c[x_idx][y_idx] : run->depth;
                double Bpel = 0.0, Bben = 0.0;
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    if (community->pelagic[pel_species_idx].fishing_flag != 1)
                        continue;
                    for (sc_idx = community->pelagic[pel_species_idx].ipelmin;
                         sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
                        if (community->pelagic[pel_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] > 0.0)
                            Bpel += community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] *
                                    exp(grid->m_values[sc_idx]) * grid->mstep;
                }
                for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
                {
                    if (community->benthic[ben_species_idx].fishing_flag != 1)
                        continue;
                    for (sc_idx = community->benthic[ben_species_idx].ibenmin;
                         sc_idx < (community->benthic[ben_species_idx].ibenmax + 1); sc_idx++)
                        if (community->benthic[ben_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] > 0.0)
                            Bben += community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] *
                                    exp(grid->m_values[sc_idx]) * grid->mstep;
                }
                double Ap = Bpel * dcorr, Ab = Bben * bhd, tot = Ap + Ab;
                double sp = (tot > 0.0) ? Ap / tot : 0.5;
                double sb = (tot > 0.0) ? Ab / tot : 0.5;
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                    if (community->pelagic[pel_species_idx].fishing_flag == 1)
                        for (sc_idx = community->pelagic[pel_species_idx].ipelmin;
                             sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
                            community->pelagic[pel_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] *= sp;
                for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
                    if (community->benthic[ben_species_idx].fishing_flag == 1)
                        for (sc_idx = community->benthic[ben_species_idx].ibenmin;
                             sc_idx < (community->benthic[ben_species_idx].ibenmax + 1); sc_idx++)
                            community->benthic[ben_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] *= sb;
            }
        }
    }
}

void calculate_reproduction(RUN *run, GRID *grid, COMMUNITY *community)
{
    int sc_idx, x_idx, y_idx, pel_species_idx, ben_species_idx;
    int n_x, n_y, n_pel_species, n_ben_species;

    double pla, pel, ben, det;

    n_x = grid->xnum;
    n_y = grid->ynum;
    n_pel_species = run->no_pelagic;
    n_ben_species = run->no_benthic;

    /*Pelagic Reproduction*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        // read in egg values from file
        if (community->pelagic[pel_species_idx].rep_method == 1)
        {
            int MAX_LEN = 15; // this is the upper bound for the row length in characters
            char *input;      // input string to store line as string
            input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
            char sep[] = ",";    // identifying the ',' char as seperating values
            char *result = NULL; // declaring a pointer for storing each individual value as a char

            /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    fgets(input, MAX_LEN, community->pelagic[pel_species_idx].fptr_rep_ts);
                    input[strlen(input) - 1] =
                        '\0'; // replacing the newline \n char of the line of text with a null \0 char

                    result = (char *)strtok(input, sep); // this is the 1st value
                    while (result != NULL)
                    {
                        community->pelagic[pel_species_idx].reproduction[x_idx][y_idx] = atof(result);
                        result = (char *)strtok(NULL, sep);
                    }
                }
            }
            free(input);
        }
        // biomass energy allocation
        if (community->pelagic[pel_species_idx].rep_method == 2)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    pla = 0;
                    pel = 0;
                    ben = 0;
                    for (sc_idx = community->pelagic[pel_species_idx].ipelmat;
                         sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
                    {
                        pla += community->pelagic[pel_species_idx].pla_bio[sc_idx][x_idx][y_idx] *
                               community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] * grid->mstep;
                        pel += community->pelagic[pel_species_idx].pel_bio[sc_idx][x_idx][y_idx] *
                               community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] * grid->mstep;
                        ben += community->pelagic[pel_species_idx].ben_bio[sc_idx][x_idx][y_idx] *
                               community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] * grid->mstep;
                    }
                    community->pelagic[pel_species_idx].reproduction[x_idx][y_idx] =
                        ((community->pelagic[pel_species_idx].R_pla * pla +
                          community->pelagic[pel_species_idx].R_pel * pel +
                          community->pelagic[pel_species_idx].R_ben * ben) *
                         grid->tstep) /
                        (grid->mstep * exp(community->pelagic[pel_species_idx].pelmin));
                    /*Connectivity/immigration floor (0-D): recruitment cannot fall below imm_frac of the
                      unfished-peak recruitment (external larval subsidy for open oligotrophic systems). rec0
                      tracks the running peak of NATURAL recruitment, frozen the instant the floor first binds,
                      so the subsidy can never inflate its own reference (no positive feedback).*/
                    if (community->pelagic[pel_species_idx].imm_frac > 0)
                    {
                        /*Per-cell reference (rec0_c[x_idx][y_idx]): each cell keeps its OWN running-peak
                          recruitment and freeze flag, so the floor never couples cells in a
                          spatial_dim>=1 run. For spatial_dim=0 (x=y=1) this is the [0][0] scalar,
                          bit-identical to the former per-species scalar.*/
                        if (!community->pelagic[pel_species_idx].rec0_frozen_c[x_idx][y_idx] &&
                            community->pelagic[pel_species_idx].reproduction[x_idx][y_idx] >
                                community->pelagic[pel_species_idx].rec0_c[x_idx][y_idx])
                            community->pelagic[pel_species_idx].rec0_c[x_idx][y_idx] =
                                community->pelagic[pel_species_idx].reproduction[x_idx][y_idx];
                        double fl = community->pelagic[pel_species_idx].imm_frac *
                                    community->pelagic[pel_species_idx].rec0_c[x_idx][y_idx];
                        if (community->pelagic[pel_species_idx].reproduction[x_idx][y_idx] < fl)
                        {
                            community->pelagic[pel_species_idx].reproduction[x_idx][y_idx] = fl;
                            community->pelagic[pel_species_idx].rec0_frozen_c[x_idx][y_idx] = 1;
                        }
                    }
                }
            }
        }
        // fecundity calculation
        if (community->pelagic[pel_species_idx].rep_method == 3)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    pel = 0;
                    for (sc_idx = community->pelagic[pel_species_idx].ipelmat;
                         sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
                    {
                        pel += community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] *
                               exp(grid->m_values[sc_idx]) * grid->mstep;
                    }
                    community->pelagic[pel_species_idx].reproduction[x_idx][y_idx] = pel * 492 * grid->tstep;
                }
            }
        }
    }

    /*Benthic Reproduction*/
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            // read in egg values from file
            if (community->benthic[ben_species_idx].rep_method == 1)
            {
                int MAX_LEN = 15; // this is the upper bound for the row length in characters
                char *input;      // input string to store line as string
                input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
                char sep[] = ",";    // identifying the ',' char as seperating values
                char *result = NULL; // declaring a pointer for storing each individual value as a char

                /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        fgets(input, MAX_LEN, community->benthic[ben_species_idx].fptr_rep_ts);
                        input[strlen(input) - 1] =
                            '\0'; // replacing the newline \n char of the line of text with a null \0 char

                        result = (char *)strtok(input, sep); // this is the 1st value
                        while (result != NULL)
                        {
                            community->benthic[ben_species_idx].reproduction[x_idx][y_idx] = atof(result);
                            result = (char *)strtok(NULL, sep);
                        }
                    }
                }
                free(input);
            }
            // biomass energy allocation
            if (community->benthic[ben_species_idx].rep_method == 2)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        det = 0;
                        for (sc_idx = community->benthic[ben_species_idx].ibenmat;
                             sc_idx < (community->benthic[ben_species_idx].ibenmax + 1); sc_idx++)
                        {
                            det += community->benthic[ben_species_idx].det_bio[sc_idx][x_idx][y_idx] *
                                   community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] * grid->mstep;
                        }
                        community->benthic[ben_species_idx].reproduction[x_idx][y_idx] =
                            (community->benthic[ben_species_idx].R_det * det * grid->tstep) /
                            (grid->mstep * exp(community->benthic[ben_species_idx].benmin));
                    }
                }
            }
            // fecundity calculation
            if (community->benthic[ben_species_idx].rep_method == 3)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        ben = 0;
                        for (sc_idx = community->benthic[ben_species_idx].ibenmat;
                             sc_idx < (community->benthic[ben_species_idx].ibenmax + 1); sc_idx++)
                        {
                            ben += community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] *
                                   exp(grid->m_values[sc_idx]) * grid->mstep;
                        }
                        community->benthic[ben_species_idx].reproduction[x_idx][y_idx] = ben * 492 * grid->tstep;
                    }
                }
            }
        }
    }
}

double phi(double q, double q_0, double sig, double trunc)
/*function determining prey preference distribution for the given values*/
{
    double ans;

    if (q <= 0 || q < (q_0 - (trunc * sig)) || q > (q_0 + (trunc * sig)))
    {
        ans = 0;
    }
    else
    {
        ans = exp(-(q - q_0) * (q - q_0) / (2 * sig * sig));
    }

    return ans;
}

//--------------------------------------------//
//      Spatial Movement Functions            //
// These contain functions that may be edited //
//--------------------------------------------//

double Cfun(double m, double t, GRID *grid, PELAGIC *pelagic)
/*function determining prey search*/
{
    double ans;
    if (t > grid->t1)
    {
        ans = pelagic->prey * exp(pelagic->gamma_prey * m);
    }
    else
    {
        ans = 0;
    }
    return ans;
}

double Dfun(double m, double t, GRID *grid, PELAGIC *pelagic)
/*function determining predator avoidance tendancy*/
{
    double ans;
    if (t > grid->t1)
    {
        ans = pelagic->pred * exp(pelagic->gamma_pred * m);
    }
    else
    {
        ans = 0;
    }
    return ans;
}

double Diffun(double m, double t, GRID *grid, PELAGIC *pelagic)
/*function determining competition coefficient*/
{
    double ans;
    if (t > grid->t1)
    {
        ans = pelagic->comp * exp(pelagic->gamma_comp * m);
    }
    else
    {
        ans = 0;
    }
    return ans;
}

double x_start(double x, double t)
/*function determining initial plankton distribution*/
{
    double ans;
    ans = 2 * exp(-((x - 500) * (x - 500) / 20000));
    return ans;
}

double y_start(double y, double t)
/*function determining initial plankton distribution*/
{
    double ans;
    ans = 2 * exp(-((y - 500) * (y - 500) / 20000));
    return ans;
}

//---------------------------------------------------------//
// These are memory allocation functions. Not to be edited //
//---------------------------------------------------------//

void free_mem(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *pelmatrix, MATRIX *benmatrix, MATRIX *xmatrix,
              MATRIX *ymatrix, double **temp_pel_params, double **temp_ben_params)
/*Frees all memory allocated within code*/
{
    int sc_idx, x_idx, pel_species_idx, pel_sub_idx, ben_species_idx;
    int n_m, n_x, n_pel_species, n_ben_species;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_pel_species = run->no_pelagic;
    n_ben_species = run->no_benthic;

    /*Free grid values*/
    free(grid->m_values);
    free(grid->t_values);
    free(grid->x_values);
    free(grid->y_values);

    /*Free pelagic values*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        for (sc_idx = 0; sc_idx < n_m; sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                free(community->pelagic[pel_species_idx].u_values[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].g_values[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].mu_values[sc_idx][x_idx]);
                for (pel_sub_idx = 0; pel_sub_idx < n_pel_species; pel_sub_idx++)
                {
                    free(community->pelagic[pel_species_idx].mu_pred_values[pel_sub_idx][sc_idx][x_idx]);
                }
                free(community->pelagic[pel_species_idx].mu_fish_values[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].pla_bio[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].pel_bio[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].ben_bio[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].a_refuge[sc_idx][x_idx]);
            }
            free(community->pelagic[pel_species_idx].u_values[sc_idx]);
            free(community->pelagic[pel_species_idx].g_values[sc_idx]);
            free(community->pelagic[pel_species_idx].mu_values[sc_idx]);
            for (pel_sub_idx = 0; pel_sub_idx < n_pel_species; pel_sub_idx++)
            {
                free(community->pelagic[pel_species_idx].mu_pred_values[pel_sub_idx][sc_idx]);
            }
            free(community->pelagic[pel_species_idx].mu_fish_values[sc_idx]);
            free(community->pelagic[pel_species_idx].pla_bio[sc_idx]);
            free(community->pelagic[pel_species_idx].pel_bio[sc_idx]);
            free(community->pelagic[pel_species_idx].ben_bio[sc_idx]);
            free(community->pelagic[pel_species_idx].a_refuge[sc_idx]);
        }
        free(community->pelagic[pel_species_idx].u_values);
        free(community->pelagic[pel_species_idx].g_values);
        free(community->pelagic[pel_species_idx].mu_values);
        for (pel_sub_idx = 0; pel_sub_idx < n_pel_species; pel_sub_idx++)
        {
            free(community->pelagic[pel_species_idx].mu_pred_values[pel_sub_idx]);
        }
        free(community->pelagic[pel_species_idx].mu_fish_values);
        free(community->pelagic[pel_species_idx].pla_bio);
        free(community->pelagic[pel_species_idx].pel_bio);
        free(community->pelagic[pel_species_idx].ben_bio);
        free(community->pelagic[pel_species_idx].a_refuge);
        free(community->pelagic[pel_species_idx].phi_values);
        free(community->pelagic[pel_species_idx].refuge_values);
        free(temp_pel_params[pel_species_idx]);

        free(community->pelagic[pel_species_idx].fname_pred);
        free(community->pelagic[pel_species_idx].fptr_pred);
    }
    free(temp_pel_params);

    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            free(community->pelagic[pel_species_idx].pla_total[x_idx]);
            free(community->pelagic[pel_species_idx].pel_total[x_idx]);
            free(community->pelagic[pel_species_idx].ben_total[x_idx]);
            free(community->pelagic[pel_species_idx].fish_total[x_idx]);
            free(community->pelagic[pel_species_idx].pred_total[x_idx]);
            free(community->pelagic[pel_species_idx].reproduction[x_idx]);
            free(community->pelagic[pel_species_idx].rec0_c[x_idx]);
            free(community->pelagic[pel_species_idx].rec0_frozen_c[x_idx]);
        }
        free(community->pelagic[pel_species_idx].pla_total);
        free(community->pelagic[pel_species_idx].pel_total);
        free(community->pelagic[pel_species_idx].ben_total);
        free(community->pelagic[pel_species_idx].fish_total);
        free(community->pelagic[pel_species_idx].pred_total);
        free(community->pelagic[pel_species_idx].reproduction);
        free(community->pelagic[pel_species_idx].rec0_c);
        free(community->pelagic[pel_species_idx].rec0_frozen_c);
    }
    free(community->pelagic);

    /*Free Plankton pointer*/
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            free(community->plankton->u_values[sc_idx][x_idx]);
            free(community->plankton->g_values[sc_idx][x_idx]);
        }
        free(community->plankton->u_values[sc_idx]);
        free(community->plankton->g_values[sc_idx]);
    }
    free(community->plankton->u_values);
    free(community->plankton->g_values);
    free(community->plankton);

    if (run->no_benthic != 0)
    {
        /*Free Benthic pointer*/
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            for (sc_idx = 0; sc_idx < n_m; sc_idx++)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    free(community->benthic[ben_species_idx].u_values[sc_idx][x_idx]);
                    free(community->benthic[ben_species_idx].g_values[sc_idx][x_idx]);
                    free(community->benthic[ben_species_idx].mu_values[sc_idx][x_idx]);
                    free(community->benthic[ben_species_idx].mu_pred_values[sc_idx][x_idx]);
                    free(community->benthic[ben_species_idx].mu_fish_values[sc_idx][x_idx]);
                    free(community->benthic[ben_species_idx].det_bio[sc_idx][x_idx]);
                }
                free(community->benthic[ben_species_idx].u_values[sc_idx]);
                free(community->benthic[ben_species_idx].g_values[sc_idx]);
                free(community->benthic[ben_species_idx].mu_values[sc_idx]);
                free(community->benthic[ben_species_idx].mu_pred_values[sc_idx]);
                free(community->benthic[ben_species_idx].mu_fish_values[sc_idx]);
                free(community->benthic[ben_species_idx].det_bio[sc_idx]);
            }
            free(community->benthic[ben_species_idx].u_values);
            free(community->benthic[ben_species_idx].g_values);
            free(community->benthic[ben_species_idx].mu_values);
            free(community->benthic[ben_species_idx].mu_pred_values);
            free(community->benthic[ben_species_idx].mu_fish_values);
            free(community->benthic[ben_species_idx].det_bio);
            free(temp_ben_params[ben_species_idx]);
        }
        free(temp_ben_params);
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                free(community->benthic[ben_species_idx].det_total[x_idx]);
                free(community->benthic[ben_species_idx].fish_total[x_idx]);
                free(community->benthic[ben_species_idx].pred_total[x_idx]);
                free(community->benthic[ben_species_idx].reproduction[x_idx]);
            }
            free(community->benthic[ben_species_idx].det_total);
            free(community->benthic[ben_species_idx].fish_total);
            free(community->benthic[ben_species_idx].pred_total);
            free(community->benthic[ben_species_idx].reproduction);
        }
        free(community->benthic);

        /*Free Detritus*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            free(community->detritus->w_values[x_idx]);
            free(community->detritus->g_values[x_idx]);
        }
        free(community->detritus->w_values);
        free(community->detritus->g_values);
        free(community->detritus);
    }

    /*Free matrix pointers*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        free(pelmatrix[pel_species_idx].a);
        free(pelmatrix[pel_species_idx].b);
        free(pelmatrix[pel_species_idx].c);
        free(pelmatrix[pel_species_idx].r);
        free(pelmatrix[pel_species_idx].u);
    }
    free(pelmatrix);

    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            free(benmatrix[ben_species_idx].a);
            free(benmatrix[ben_species_idx].b);
            free(benmatrix[ben_species_idx].c);
            free(benmatrix[ben_species_idx].r);
            free(benmatrix[ben_species_idx].u);
        }
        free(benmatrix);
    }

    free(xmatrix->a);
    free(xmatrix->b);
    free(xmatrix->c);
    free(xmatrix->r);
    free(xmatrix->u);

    free(ymatrix->a);
    free(ymatrix->b);
    free(ymatrix->c);
    free(ymatrix->r);
    free(ymatrix->u);
}

void *safe_malloc(size_t size, int line, char *filename)
/*Safe version of malloc with built in error checking. Exits if insufficient memory is avilable*/
{
    void *ptr;
    ptr = malloc(size);
    if (ptr == NULL)
    {
        error("Out of memory at line %d file %s\n", line, filename);
    }
    return ptr;
}

FILE *safe_fopen(char *filename, char *mode, int line, char *fname)
/*Safe version of fopen with built in error checking. Exits if the file cannot be opened*/
{
    FILE *fptr;
    fptr = fopen(filename, mode);
    if (fptr == NULL)
    {
        error("Unable to open file at line %d file %s\n", line, fname);
    }
    return fptr;
}
