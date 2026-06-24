#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>

//------------------------//
// Structure Declarations //
//------------------------//

typedef struct run_info{
        
        char *filename;        //root filename for runs
        int no_pelagic;
        int no_benthic;
        int spatial_dim;
        int coupled_flag;      //flag for coupled info; 0-uncoupled, 1-coupled
        int diff_method;       //flag for integration method; 0-fast (poss. unstable), 1-slow (more stable)

        char fname_summ[40];
        FILE *fptr_summ;

        
        } RUN;

/*Contains all information regarding grid discretisation*/
typedef struct grid_info{

        double mmin,mmax;      //ln(mass)
        double mstep;          //step length
        int mnum;              //number of steps
        double moutstep;       //output to print every m
        double *m_values;      //actual log(m) values
        
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
        int initial_flag;
        int ts_flag;
        int fishing_flag;      //flag for fishing; 0-off, 1-on
        
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
        
/*Structure for use in implicit upwind method*/
typedef struct vectors{
        
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
void setup_plankton(RUN *, GRID *, PLANKTON *, double *, char *, int, int);           //sets up plankton values
void setup_pelagic(RUN *, GRID *, PELAGIC *, double *, char *, int, int, int, int);   //sets up pelagic values
void setup_benthic(RUN *, GRID *, BENTHIC *, double *, char *, int, int, int, int);   //sets up benthic values
void setup_detritus(RUN *, GRID *, DETRITUS *, double *, char *, int, int);           //sets up detritus values

double phi(double, double, double, double);         //Calculates feeding kernel for each pelagic species

double ** setup_pelagic_params(int, double *);      //Converts 1D vector of pelagic params to 2D array of pelagic params
double ** setup_benthic_params(int, double *);      //Converts 1D vector of benthic params to 2D array of benthic params
void setup_matrix(int, MATRIX *);

/*Timestepping and Output Function*/
void calculate_results(RUN *, GRID *, COMMUNITY *, MATRIX *, MATRIX*, MATRIX*, MATRIX*);

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
void mass_solver(RUN *, GRID *,COMMUNITY *, MATRIX *, MATRIX *, double ***, double ****, double ****, double **);   //Calculates u values for next time step
void xmove_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, double ****, int); //Calculates x movement in 1D space
void ymove_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, double ****, int); //Calculates y movement in 2D space
void tridag(MATRIX *);                                          //Inverts a tridiagonal matrix
void trimul(MATRIX *);

/*Predation, Growth and Renewal Functions*/
void calculate_g_and_mu(RUN *run, GRID *, COMMUNITY *);               //Calculates g and mu values given u values

double pla_biomass(int, int, int, int, GRID *, COMMUNITY *);          //Calculates the biomass eaten from plankton spectra
double pel_biomass(int, int, int, int, RUN *, GRID *, COMMUNITY *);   //Calculates the biomass eaten from pelagic spectra
double ben_biomass(int, int, int, int, RUN *, GRID *, COMMUNITY *);   //Calculates the biomass eaten from benthic spectra
double det_biomass(int, int, int, int, RUN *, GRID *, COMMUNITY *);   //Calculates the biomass eaten from detritus system

double g_pel(int, int, int, int, RUN *, GRID *, COMMUNITY *);         //Calculates g values for pelagic system from biomass values
double g_ben(int, int, int, int, RUN *, GRID *, COMMUNITY *);         //Calculates g values for benthic systsm from biomass values
double g_det(int, int, RUN *, GRID *, COMMUNITY *);                   //Calculates g values for detrital system from biomass values

double mu_pel_pred(int, int, int, int, int, RUN *, GRID *, COMMUNITY *);
double mu_pel_fish(int, int, int, int, RUN *, GRID *, COMMUNITY *);
double mu_ben_pred(int, int, int, int, RUN *, GRID *, COMMUNITY *);
double mu_ben_fish(int, int, int, int, RUN *, GRID *, COMMUNITY *);

double mu_pel(int, int, int, int, RUN *, GRID *, COMMUNITY *);        //Calculates mu values for pelagic system
double mu_ben(int, int, int, int, RUN *, GRID *, COMMUNITY *);        //Calculates mu values for benthic system
double mu_det(int, int, RUN *, GRID *, COMMUNITY *);                   //Calculates mu values for detrital system from biomass values

void calculate_reproduction(RUN *, GRID *, COMMUNITY *);              //Calculates reproduction for any system
void calculate_fishing(RUN *, GRID *, COMMUNITY *);                   //Calculates fishing mortality for any system

/*Spatial Movement Functions*/
double Cfun(double, double, GRID *, PELAGIC *);
double Dfun(double, double, GRID *, PELAGIC *);
double Diffun(double, double, GRID *, PELAGIC *);
double x_start(double, double);
double y_start(double, double);

/*Memory Management Functions*/
void free_mem(RUN *, GRID *, COMMUNITY *, MATRIX *, MATRIX *, MATRIX *, MATRIX *, double **, double **);     //Frees all allocated memory
void *safe_malloc(size_t, int, char *);           //Wrappered function for memory alloaction
FILE *safe_fopen(char *, char *, int, char *);    //Wrappered function for file opening


//--------------------------------------//
// This is the main program called by R //
//--------------------------------------//


void SizeSpectrum(int *run_params, double *grid_params, double *pla_params, double *pel_params, double *ben_params, double * det_params, char **names_params, int *flags_params)
// run_params is an array of run parameters
// grid_params is an array of grid discretisation values
// pla_params is plankton params
// pel_params is pelagic params
// ben_params is benthic parameters
// det_params is detritus parameters
// names_params is a list of filenames
// flags_params is a list of selection 'flag' parameters
{
    int s,b;
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
    community.plankton=(PLANKTON *)safe_malloc(sizeof(PLANKTON),__LINE__,__FILE__);
    community.pelagic=(PELAGIC *)safe_malloc(run.no_pelagic*sizeof(PELAGIC),__LINE__,__FILE__);
    if(run.no_benthic!=0){
            community.benthic=(BENTHIC *)safe_malloc(run.no_benthic*sizeof(BENTHIC),__LINE__,__FILE__);
            community.detritus=(DETRITUS *)safe_malloc(sizeof(DETRITUS),__LINE__,__FILE__);
    }
    
    /*Setup plankton*/
    setup_plankton(&run, &grid, community.plankton, pla_params, names_params[1], flags_params[0], flags_params[1]); 
    
    /*Setup pelagic*/
    temp_pel_params=setup_pelagic_params(run.no_pelagic, pel_params);
    for(s=0 ; s<run.no_pelagic ; s++){
            setup_pelagic(&run, &grid, &(community.pelagic[s]), temp_pel_params[s], names_params[s+2], flags_params[4*s+2], flags_params[4*s+3], flags_params[4*s+4], flags_params[4*s+5]);
    }
    
    if(run.no_benthic!=0){
            /*Setup benthic*/
            temp_ben_params=setup_benthic_params(run.no_benthic, ben_params);
            for(b=0 ; b<run.no_benthic ; b++){
                    setup_benthic(&run, &grid, &(community.benthic[b]), temp_ben_params[b], names_params[b+run.no_pelagic+2], flags_params[4*b+(2+(4*run.no_pelagic))], flags_params[4*b+(3+(4*run.no_pelagic))], flags_params[4*b+(4+(4*run.no_pelagic))], flags_params[4*b+(5+(4*run.no_pelagic))]);
            }
    
            /*Setup detritus*/
            setup_detritus(&run, &grid, community.detritus, det_params, names_params[run.no_pelagic+run.no_benthic+2], flags_params[2+4*(run.no_pelagic+run.no_benthic)], flags_params[3+4*(run.no_pelagic+run.no_benthic)]);
    }
   
    /*Setup differencing matrices*/
    pelmatrix=(MATRIX *)safe_malloc(run.no_pelagic*sizeof(MATRIX),__LINE__,__FILE__);
    for(s=0 ; s<run.no_pelagic ; s++){
            setup_matrix((community.pelagic[s].ipelmax-community.pelagic[s].ipelmin+1), &(pelmatrix[s]));
    }
    
    if(run.no_benthic!=0){
            benmatrix=(MATRIX *)safe_malloc(run.no_benthic*sizeof(MATRIX),__LINE__,__FILE__);
            for(b=0 ; b<run.no_benthic ; b++){
                    setup_matrix((community.benthic[b].ibenmax-community.benthic[b].ibenmin+1), &(benmatrix[b]));
            }
    }
    
    setup_matrix(grid.xnum, &xmatrix);
    setup_matrix(grid.ynum, &ymatrix);

    /*Perform all finite differencing*/
    calculate_results(&run, &grid, &community, pelmatrix, benmatrix, &xmatrix, &ymatrix);

    /*Free all allocated memory*/
    free_mem(&run, &grid, &community, pelmatrix, benmatrix, &xmatrix, &ymatrix, temp_pel_params, temp_ben_params);

}

/* [truncated for brevity in commit tool payload] */
