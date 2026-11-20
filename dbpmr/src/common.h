// `common.h` and `common.c` files contain items that are common and used across
// multiple C files in this project. These include struct definitions, matrix
// operations, the feeding kernel and memory-management functions

#ifndef COMMON_H
#define COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rinternals.h>

//------------------------//
// Structure Declarations //
//------------------------//

/*Inputs handed over from R as arrays instead of <run>/Input text files (opt-in; NULL pointers = not supplied).
  All arrays are R arrays [x][y][layer] (column-major): index k + x*(l + y*layer).
  Plankton: the spectrum in cell (k,l) for a given layer is u(m) = a*exp(b*m), with a = pla_int and b = pla_slp.
  Forcing: one array per channel of the forcing plugin, in the order of FRC_NAMES below; a channel has either
  one layer (constant in time) or the common number of layers.*/
#define N_FRC 9
enum { FRC_PEL_TEMPEFF, FRC_BEN_TEMPEFF, FRC_SINKING_RATE, FRC_DEPTH, FRC_RESIDENCE_TIME,
       FRC_EXPORT_ATTN, FRC_EXPORT_GAMMA, FRC_EXPORT_MAGG, FRC_GRAVITY };
extern const char *FRC_NAMES[N_FRC];   //names of the forcing channels in the list passed from R (defined in memory_inputs.c)
/*Fishing mortality of one species supplied from R as an array instead of its _fishing_ts.txt file:
  the rate applies to the size classes with lo <= m <= hi and is zero for the others.*/
typedef struct fishing_array{
        const double *rate;    //fishing mortality rate, an R array [x][y][layer]; NULL = read the fishing file
        const int *layer;      //0-based layer at each model time point (0, tstep, 2*tstep, ...)
        int nlayer;            //number of layers
        int nblock;            //length of layer
        R_xlen_t len;          //length of rate (checked against x*y*layers once the grid is known)
        double lo,hi;          //fished size range (natural-log body mass)
        } FISHING_ARRAY;

typedef struct array_inputs{
        const double *pla_int; //plankton intercept a
        const double *pla_slp; //plankton slope b
        const int *pla_layer;  //0-based layer used by each time block (block 0 = initial state, then one per timestep)
        int pla_nlayer;        //number of layers in the arrays
        int pla_nblock;        //length of pla_layer
        R_xlen_t pla_len;      //length of pla_int / pla_slp (checked against x*y*layers once the grid is known)

        int have_forcing;      //1 = forcing channels supplied as arrays (forcing_ts.txt is not opened)
        const double *frc[N_FRC];   //channel arrays (NULL = channel absent)
        int frc_nlayer[N_FRC];      //layers of each channel (1 = constant in time)
        R_xlen_t frc_len[N_FRC];    //length of each channel array (checked against x*y*layers)
        const int *frc_layer;  //0-based layer at each model time point (0, tstep, 2*tstep, ...); the forcing applied
                               //in timestep r (advancing to time r*tstep) is that of time point r
        int frc_nblock;        //length of frc_layer

        FISHING_ARRAY *fish_pel;   //one entry per pelagic species (NULL = none supplied for any species)
        FISHING_ARRAY *fish_ben;   //one entry per benthic species
        } ARRAY_INPUTS;

/*Model output handed back to R as arrays (opt-in). The arrays live in R memory (elements of the list
  returned by .Call) and are filled at the same times, cells and size classes as the text result files:
  abundance arrays are [x][y][size][time] (column-major): index kk + nxo*(ll + nyo*(ii + nmo*it)).*/
typedef struct output_arrays{
        int mode;              //0 = text files only (default), 1 = text files and arrays, 2 = arrays only
        int nxo,nyo,nmo,nto;   //number of output cells in x and y, output size classes and output times
        int xstride,ystride;   //every xstride-th / ystride-th cell is output (xoutstep/xstep, youtstep/ystep)
        int *mi;               //grid index of each output size class
        int it;                //next output time slot
        double *t,*x,*y,*m;    //coordinate vectors
        int want_g,want_mu,want_pred,want_fish,want_catch;   //also return growth, total mortality, predation mortality, fishing mortality, catch
        double *pla;           //plankton abundance
        double **pel;          //abundance of each pelagic species
        double **ben;          //abundance of each benthic species
        double *det;           //detritus biomass [x][y][time]
        double **pel_g,**pel_mu,**pel_fish;   //growth rate, total mortality and fishing mortality of each pelagic species (NULL if not wanted)
        double ***pel_pred;    //pel_pred[s][p]: predation mortality on pelagic species s from pelagic species p
        double **ben_g,**ben_mu,**ben_pred,**ben_fish;   //the same for each benthic species (predation from all pelagic species)
        double **pel_catch,**ben_catch;   //catch rate density of each species: fishing mortality * abundance * body mass (per unit ln mass)
        SEXP list;             //the R list holding all of the above
        } OUTPUT_ARRAYS;

typedef struct run_info{
        
        char *filename;        //root filename for runs
        int no_pelagic;
        int no_benthic;
        int spatial_dim;
        int coupled_flag;      //flag for coupled info; 0-uncoupled, 1-coupled
        int diff_method;       //flag for integration method; 0-fast (poss. unstable), 1-slow (more stable)

        char fname_summ[40];
        FILE *fptr_summ;

        /* --- Environmental forcing plugin (auto-detected <run>/Input/forcing_ts.txt) ---
           A header line names the channels (any subset/order), one data row per timestep.
           Recognised channels: pel_tempeff, ben_tempeff (surface/seafloor Boltzmann-Arrhenius
           multipliers on feeding A + background mortality mu_0) and sinking_rate (fraction of
           surface-origin detritus reaching the seafloor). Add a forcing = add a column name +
           a col_* index + an apply hook. Missing channels fall back to the neutral defaults. */
        int forcing_flag;      //time-varying forcing driver active
        FILE *fptr_forcing;    //forcing time-series file ptr
        const ARRAY_INPUTS *arr; //inputs supplied from R as arrays (NULL or have_forcing=0 -> forcing comes from the file, if any)
        int frc_step;          //timesteps of array forcing applied so far
        int fish_block;        //number of times fishing mortality has been set (indexes the time point of array fishing)
        OUTPUT_ARRAYS *out;    //output returned to R as arrays (NULL or mode 0 -> text files only)
        int forcing_ncol;      //number of columns in each data row
        int col_pel_tempeff;   //column index of each channel (-1 if absent)
        int col_ben_tempeff;
        int col_sinking_rate;
        int col_depth;
        int col_residence_time;
        int col_export_attn;   //size-dependent detritus export (opt-in): particle of ln-mass m reaches
        int col_export_gamma;  //the seafloor with fraction exp(-export_attn*exp(-export_gamma*max(m,export_magg)))
        int col_export_magg;   //aggregation floor: production below export_magg (ln-mass) sinks as marine snow of that effective size
        int col_gravity;       //DBPM.md gravity fishing (opt-in): if present, calculate_fishing() scales F on each group by that group's CURRENT areal fishable-biomass share (F_g = base * B_g/sum B), recomputed per cell per timestep
        double pel_tempeff;    //current surface-temperature multiplier (default 1)
        double ben_tempeff;    //current seafloor-temperature multiplier (default 1)
        double sinking_rate;   //current detritus sinking/export fraction (default 1)
        int forcing_percell;   //per-cell forcing: forcing_ts.txt has x*y rows per timestep (env FORCING_PERCELL=1, spatial_dim>=1). 0 = one global row/timestep (default, bit-identical)
        int temp_on_mort;      //cached TEMP_ON_MORT (1 = temperature scales background mortality mu_0)
        double **pel_tempeff_c;//per-cell pelagic temperature multiplier (used when forcing_percell; applied in feeding/mortality)
        double **ben_tempeff_c;//per-cell benthic temperature multiplier
        double **sinking_rate_c;//per-cell detritus sinking/export fraction
        double **depth_c;      //per-cell water-column depth (value used only in the gravity-fishing pelagic:benthic split)
        double **export_attn_c;//per-cell size-dependent export attenuation number
        double **export_gamma_c;//per-cell size-dependent export sinking-velocity exponent
        double **export_magg_c;//per-cell size-dependent export aggregation floor (ln-mass)
        double depth;          //water-column depth (m); read+available. Supplying the channel toggles Dunne burial ON (col_depth>=0); the value is NOT used to rescale the burial flux (sizemodel applies burial to the per-volume flux directly)
        double residence_time; //detritus residence time tau (yr) for the alternative first-order pool loss W/tau (dimensionally clean; supplied via the residence_time channel). 0/absent = off
        double export_attn;    //size-dependent export: attenuation number A = k(T)*z/w0 (per-timestep, from export_attn channel; larger = more remineralisation = less export)
        double export_gamma;   //size-dependent export: allometric sinking-velocity exponent (w_s ~ mass^gamma); larger particles sink faster -> higher export
        double export_magg;    //aggregation floor (ln-mass): small production coagulates into marine snow of this effective size before sinking; export uses max(m, export_magg). Default -1e30 = no floor (pure per-particle)
        double gravity;        //gravity channel value (unused; presence via col_gravity toggles gravity fishing)

        } RUN;

/*Contains all information regarding grid discretisation*/
typedef struct grid_info{

        double mmin,mmax;      //ln(mass)
        double mstep;          //step length
        int mnum;              //number of steps
        double moutstep;       //output to print every m
        double *m_values;      //actual log(m) values
        double *exp_m_values;  //exp(m_values), computed once
        
        double t1,tmax;        //time discretisation
        double tstep;
        int j1,tnum;
        double toutmin, toutmax;
        double toutstep;
        double *t_values;
        
        double xmin,xmax;      //x-space discretisation (if required)
        double xstep;
        int xnum;
        double xoutstep;
        double *x_values;
        
        double ymin,ymax;      //y-space discretisation (if required)
        double ystep;
        int ynum;
        double youtstep;
        double *y_values;
        
        } GRID;

/*Contains all info regarding the plankton dynamics*/
typedef struct plankton_info{
        
        char *filename;
        
        double plamin,plamax;
        int iplamin,iplamax;
       
        double mu_0;
        double beta;
        double u_0;
        double lambda;
        
        int initial_flag;
        int ts_flag;
        const ARRAY_INPUTS *arr; //plankton series supplied from R as arrays (NULL = read the _ts file)
        int ts_block;          //next time block to load from the arrays
        
        double ***u_values;
        double ***g_values;    //for futureproofing
        
        //output files
        char fname_r[40];
        char fname_summ[40];
        
        FILE *fptr_r;
        FILE *fptr_summ;
        
        //input files
        char fname_ts[40];
        
        FILE *fptr_ts;
        
        } PLANKTON;
        
/*Contains all info regarding an individual species' constants and abundances*/
typedef struct pelagic_info{

        char *filename;        //filename for each species

        double pelmin,pelmat,pelmax;
        int ipelmin,ipelmat,ipelmax;

        double A;              //volume search coefficient
        double alpha;          //volume search exponent
        double mu_0;           //intrinsic (juvenile) mortality coefficient
        double beta;           //intrinsic (juvenile) mortality exponent
        double mu_s;           //senesence mortality coefficient
        double epsilon;        //senesence mortality constant

        double A_base;         //un-temperature-scaled search coefficient (temp driver)
        double mu_0_base;      //un-temperature-scaled background mortality (temp driver)
        double tempeff;        //current-timestep temperature multiplier (1 if no temp driver)

        double u_0;            //initial intercept and abundance
        double lambda;         //primary production slope
        
        double K_pla;          //growth efficiency (proportion of plankton food to growth)
        double R_pla;          //reproduction efficiency (proportion of plankton food to reproduction)
        double Ex_pla;         //defecation efficiency (proportion of plankton food to defecation)
        double K_pel;          //growth efficiency (proportion of pelagic food to growth)
        double R_pel;          //reproduction efficiency (proportion of pelagic food to reproduction)
        double Ex_pel;         //defecation efficiency (proportion of pelagic food to defecation)        
        double K_ben;          //growth efficiency (proportion of benthic food to growth)
        double R_ben;          //reproduction efficiency (proportion of benthic food to reproduction)
        double Ex_ben;         //defecation efficiency (proportion of benthic food to defecation)
        
        double pref_pla;       //attack rate for plankton
        double pref_pel;       //attack rate for pelagic
        double pref_ben;       //attack rate for benthic
        double **pref_ben_c;   //optional PER-CELL benthic preference [x][y] (from <run>/Input/<name>_cellpar.txt); NULL = use the scalar pref_ben everywhere
        
        double q_0;            //optimum predator prey mass ratio
        double sig;            //standard deviation of feeding kernel
        double trunc;          //number of standard deviations at which to truncate feeding kernel
        double *phi_values;    //phi values
        
        double prey;           //prey seeking coefficient
        double pred;           //predator avoiding coefficient
        double comp;           //competition coefficient
        double gamma_prey;     //prey seeking exponent
        double gamma_pred;     //predator avoiding exponent
        double gamma_comp;     //compeititon exponent
        
        int rep_method;        //reproduction method; 0-fixed amount , 1-biomass dependent
        double imm_frac;       //connectivity/immigration recruitment floor as a fraction of unfished recruitment (0=off)
        double rec0;           //running-peak (unfished) recruitment density at the smallest size (0-D / legacy scalar; spatial runs use rec0_c)
        int rec0_frozen;       //freeze rec0 at first bind so the subsidy cannot inflate its own reference (0-D / legacy scalar)
        double **rec0_c;       //PER-CELL running-peak recruitment (spatial-safe connectivity floor; each cell keeps its own reference so cells don't couple)
        int **rec0_frozen_c;   //PER-CELL freeze flag for the connectivity floor
        int refuge_on;         //reef prey-refuge mechanism active (refuge.txt supplied); 0 = off (a(w)=1, standard dbpmr)
        double refuge_min_a;   //availability floor min_A (REFUGE_MIN_A env, default 0); a(w) clamped to [min_a,1]
        double *refuge_values; //per prey-size refuge-space density (count m^-2); NULL/all-zero when off
        double ***a_refuge;    //per (size,x,y) availability multiplier a(w)=max(min_a,1-refuge(w)/competitor(w)); recomputed each step
        int initial_flag;
        int ts_flag;
        int fishing_flag;      //flag for fishing; 0-off, 1-on
        const FISHING_ARRAY *fish_arr;   //fishing mortality supplied from R as an array (NULL = read the fishing file)
        
        double ***u_values;        //number densities for each time step (outputted as file)
        double ***g_values;        //growth rates for each time step (outputted as file)
        double ***mu_values;       //mortality rates for each time step (outputted as file)
        double ****mu_pred_values;  //mortality rates due to predation (outputted as file)
        double ***mu_fish_values;  //mortality rates due to fishing (outputted as file)
        
        double ***pla_bio;     //plankton biomass eaten for each time step (used only in calculations)
        double ***pel_bio;     //pelagic biomass eaten for each time step (used only in calculations)
        double ***ben_bio;     //benthic biomass eaten for each time step (used only in calculations)
        
        double **pla_total;      //total plankton biomass eaten at a spatial point (ouputted as rate in summary file)
        double **pel_total;      //total pelagic biomass eaten at a spatial point (ouputted as rate in summary file)
        double **ben_total;      //total benthic biomass eaten at a spatial point (ouputted as rate in summary file)
        
        double **pred_total;     //total biomass lost to predation at a spatial point (outputted as a rate in sumamry file)
        double **fish_total;     //total biomass lost to fishing at a spatial point (ouputted as rate in summary file)
        double **reproduction;   //reproduction number density at a spatial point (ouputted as rate in summary file)
        
        //output files
        char fname_r[40];
        char fname_g[40];
        char fname_m[40];
        char **fname_pred;
        char fname_fish[40];
        char fname_summ[40];
        
        FILE *fptr_r;          //results output file ptr
        FILE *fptr_g;          //growth output file ptr
        FILE *fptr_m;          //mortality output file ptr
        FILE **fptr_pred;       //predation mortality file ptr
        FILE *fptr_fish;       //fishing mortality file ptr
        FILE *fptr_summ;       //summary file ptr

        //input files
        char fname_ts[40];
        char fname_fish_ts[40];
        char fname_rep_ts[40];
                
        FILE *fptr_ts;         //time series file ptr
        FILE *fptr_fish_ts;    //time series fishing file ptr
        FILE *fptr_rep_ts;     //time series eggs file ptr
        
        
        } PELAGIC;

/*Contains all info regarding an individual species' constants and abundances*/
typedef struct benthic_info{

        char *filename;        //filename for each species

        double benmin,benmat,benmax;
        int ibenmin,ibenmat,ibenmax;
        
        double alpha;          //volume search exponent
        double A;              //volume search coefficient
        double beta;           //intrinsic (juvenile) mortality exponent
        double mu_0;           //intrinsic (juvenile) mortality coefficient
        double epsilon;        //senesence mortality constant
        double mu_s;           //senesence mortality coefficient

        double A_base;         //un-temperature-scaled search coefficient (temp driver)
        double mu_0_base;      //un-temperature-scaled background mortality (temp driver)
        double tempeff;        //current-timestep temperature multiplier (1 if no temp driver)

        double lambda;         //primary production slope
        double u_0;            //initial intercept and abundance
        
        double K_det;          //growth efficiency (proportion of plankton food to growth)
        double R_det;          //reproduction efficiency (proportion of plankton food to reproduction)
        double Ex_det;         //defecation efficiency (proportion of plankton food to defecation)
        
        double pref_det;       //attack rate for detritus
        
        int rep_method;        //reproduction method; 0-fixed amount , 1-biomass dependent
        int initial_flag;      //initial distribution supplied; 0-no , 1-yes
        int ts_flag;           //timeseries of distributions supplied; 0-no , 1-yes
        int fishing_flag;      //flag for fishing; 0-off, 1-on
        const FISHING_ARRAY *fish_arr;   //fishing mortality supplied from R as an array (NULL = read the fishing file)
        
        double ***u_values;        //number densities for each time step (outputted as file)
        double ***g_values;        //growth rates for each time step (outputted as file)
        double ***mu_values;       //mortality rates for each time step (outputted as file)
        double ***mu_pred_values;  //mortality rates due to predation (outputted as file)
        double ***mu_fish_values;  //mortality rates due to fishing (outputted as file)
                
        double ***det_bio;         //detritus biomass eaten for each time step (used in calculations)
        
        double **det_total;        //total detritus biomass eaten at a spatial point (ouputted as rate in summary file)

        double **fish_total;       //total biomass lost to fishing at a spatial point (ouputted as rate in summary file)
        double **pred_total;       //total biomass lost to predation at a spatial point (ouputted as rate in summary file)
        double **reproduction;     //reproduction number density at a spatial point (ouputted as rate in summary file)

        //output files
        char fname_r[40];
        char fname_g[40];
        char fname_m[40];
        char fname_pred[40];
        char fname_fish[40];
        char fname_summ[40];
        
        FILE *fptr_r;          //results output file ptr
        FILE *fptr_g;          //growth output file ptr
        FILE *fptr_m;          //mortality output file ptr
        FILE *fptr_pred;       //predation mortality file ptr
        FILE *fptr_fish;       //fishing mortality file ptr
        FILE *fptr_summ;
        
        //input files
        char fname_ts[40];
        char fname_fish_ts[40];
        char fname_rep_ts[40];
                
        FILE *fptr_ts;         //time series file ptr
        FILE *fptr_fish_ts;    //time series fishing file ptr
        FILE *fptr_rep_ts;     //time series rep file ptr
        
        } BENTHIC;

/*Contains all info regarding the plankton dynamics*/
typedef struct detritus_info{
        
        char *filename;

        double w_0;

        int initial_flag;
        int ts_flag;
        
        double **w_values;
        double **g_values;
        double **mu_values;
        
        //output files
        char fname_summ[40];
        
        FILE *fptr_summ;
        
        //input files
        char fname_ts[40];
        
        FILE *fptr_ts;         //time series file ptrr
        
        } DETRITUS;

/*Contains all species data*/
typedef struct species_answer_info{
        
        PLANKTON *plankton;
        PELAGIC *pelagic;
        BENTHIC *benthic;
        DETRITUS *detritus;
        
        } COMMUNITY;

/*Benthic feeding preference of a pelagic species in cell (k,l): the per-cell value when a
  <name>_cellpar.txt was supplied, otherwise the scalar parameter (standard dbpmr).*/
static inline double pel_pref_ben(const PELAGIC *p, int k, int l)
{
        return p->pref_ben_c ? p->pref_ben_c[k][l] : p->pref_ben;
}
        
/*Structure for use in implicit upwind method*/
typedef struct vectors{
        
        double *a;
        double *b;
        double *c;
        double *r;
        double *u;
        int size;
        
        } MATRIX;

/*Matrix operations and the feeding kernel*/
void trimul(MATRIX *);
void tridag(MATRIX *);                                          //Inverts a tridiagonal matrix
double phi(double, double, double, double);         //Calculates feeding kernel for each pelagic species

/*Memory Management Functions*/
void free_mem(RUN *, GRID *, COMMUNITY *, MATRIX *, MATRIX *, MATRIX *, MATRIX *, double **, double **);     //Frees all allocated memory
void *safe_malloc(size_t, int, char *);           //Wrappered function for memory alloaction
FILE *safe_fopen(char *, char *, int, char *);    //Wrappered function for file opening

#endif
