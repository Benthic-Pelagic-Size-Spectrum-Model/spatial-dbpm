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


//-----------------//
// Setup Functions //
//-----------------//

void setup_run(RUN *run, int *run_params, char *filename)
{
     run->filename=filename;
     run->no_pelagic=run_params[0];
     run->no_benthic=run_params[1];
     run->spatial_dim=run_params[2];
     run->coupled_flag=run_params[3];
     run->diff_method=run_params[4];

     /*Create summary output filenames*/
     snprintf(run->fname_summ, sizeof(run->fname_summ),"%s/%s.txt",run->filename,"parameters");
     
     /*Open output files*/
     run->fptr_summ=safe_fopen(run->fname_summ,"w",__LINE__,__FILE__);
     
     fprintf(run->fptr_summ,"filetype:      %s\n","parameters");
     fprintf(run->fptr_summ,"speciestype:   %s\n","NA");
     fprintf(run->fptr_summ,"\n");
     
     /*Close ouput files*/
     fclose(run->fptr_summ);

}

void setup_grid(GRID *grid, double *grid_input)
{
     /*Declare Mass Discretisation*/     
     grid->mmin=grid_input[0];
     grid->mmax=grid_input[1];
     grid->mstep=grid_input[2];
                                          //m step length
     grid->mnum=(int)rint(((grid->mmax-grid->mmin)/grid->mstep)+1); //number of m values
     grid->moutstep=grid_input[3];
     
     /*Declare Time Discretisation*/
     grid->t1=grid_input[4];
     grid->tmax=grid_input[5];                                      //maximum t value
     grid->tstep=grid_input[6];                                     //t step length
     grid->tnum=(int)rint((grid->tmax/grid->tstep)+1);              //number of t values
     grid->j1=(int)rint(grid->t1/grid->tstep);
     grid->toutmin=grid_input[7];
     grid->toutmax=grid_input[8];
     grid->toutstep=grid_input[9];

     /*Declare X-Space Discretisation*/
     grid->xmin=grid_input[10];
     grid->xmax=grid_input[11];
     grid->xstep=grid_input[12];
     grid->xnum=(int)rint(((grid->xmax-grid->xmin)/grid->xstep)+1);
     grid->xoutstep=grid_input[13];
     
     /*Declare Y-Space Discretisation*/
     grid->ymin=grid_input[14];
     grid->ymax=grid_input[15];
     grid->ystep=grid_input[16];
     grid->ynum=(int)rint(((grid->ymax-grid->ymin)/grid->ystep)+1);
     grid->youtstep=grid_input[17];
     
     
     /*Create arrays*/
     int i,j,k,l;
     int m,t,x,y;
     
     m=grid->mnum;
     t=grid->tnum;
     x=grid->xnum;
     y=grid->ynum;
     
     grid->m_values=(double *)safe_malloc(m*sizeof(double),__LINE__,__FILE__);
     grid->t_values=(double *)safe_malloc(t*sizeof(double),__LINE__,__FILE__);
     grid->x_values=(double *)safe_malloc(x*sizeof(double),__LINE__,__FILE__);
     grid->y_values=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
     
     /*Setup values for m*/ 
     for(i=0 ; i<m ; i++){
             grid->m_values[i]=grid->mmin+(i*grid->mstep);
     }
     /*Setup values for t*/
     for(j=0 ; j<t ; j++){
             grid->t_values[j]=(j*grid->tstep);
     }
     /*Setup values for x*/
     for(k=0 ; k<x ; k++){
             grid->x_values[k]=grid->xmin+(k*grid->xstep);
     }
     /*Setup values for y*/
     for(l=0 ; l<y ; l++){
             grid->y_values[l]=grid->ymin+(l*grid->ystep);
     }
}

double **setup_pelagic_params(int no_pelagic, double *pel_params)
{
    int i,s;
    int no_params=32;
    double **temp_params;
    temp_params=(double **)safe_malloc(no_pelagic*sizeof(double *),__LINE__,__FILE__);
    for(s=0 ; s<no_pelagic ; s++){
            temp_params[s]=(double *)safe_malloc(no_params*sizeof(double),__LINE__,__FILE__);
            for(i=0 ; i<no_params ; i++){
                    temp_params[s][i]=pel_params[no_params*s+i];
            }
    }
    return(temp_params);
}

double **setup_benthic_params(int no_benthic, double *ben_params)
{
    int i,b;
    int no_params=15;
    double **temp_params;
    temp_params=(double **)safe_malloc(no_benthic*sizeof(double *),__LINE__,__FILE__);
    for(b=0 ; b<no_benthic ; b++){
            temp_params[b]=(double *)safe_malloc(no_params*sizeof(double),__LINE__,__FILE__);
            for(i=0 ; i<no_params ; i++){
                    temp_params[b][i]=ben_params[no_params*b+i];
            }
    }
    return(temp_params);
}

void setup_plankton(RUN *run, GRID *grid, PLANKTON *plankton, double *pla_params, char *filename, int initial_flag, int ts_flag)
{
     /*Declare plankton name*/
     plankton->filename=filename;

     /*Setup Plankton Sizes*/
     plankton->plamin=pla_params[0];
     plankton->plamax=pla_params[1];
     plankton->iplamin=(int)rint((plankton->plamin-grid->mmin)/grid->mstep);
     plankton->iplamax=(int)rint((plankton->plamax-grid->mmin)/grid->mstep);

     /*Setup Plankton Constants*/
     plankton->mu_0=pla_params[2];
     plankton->beta=pla_params[3];
     plankton->u_0=pla_params[4];
     plankton->lambda=pla_params[5];
     
     /*Setup Plankton Flags*/
     plankton->initial_flag=initial_flag;
     plankton->ts_flag=ts_flag;
     
     /*Create plankton output filenames*/
     snprintf(plankton->fname_r, sizeof(plankton->fname_r),"%s/%s/%s.txt",run->filename,plankton->filename,"results");
     snprintf(plankton->fname_summ, sizeof(plankton->fname_summ),"%s/%s/%s.txt",run->filename,plankton->filename,"summary");
     
     /*Open output files*/
     plankton->fptr_r=safe_fopen(plankton->fname_r,"w",__LINE__,__FILE__);
     plankton->fptr_summ=safe_fopen(plankton->fname_summ,"w",__LINE__,__FILE__);
     
     /*Initialise output files*/
     fprintf(plankton->fptr_r,"filetype:      %s\n","results");
     fprintf(plankton->fptr_r,"speciestype:   %s\n","plankton");
     fprintf(plankton->fptr_r,"\n");
     print_mass_header(grid, plankton->fptr_r);
     
     fprintf(plankton->fptr_summ,"filetype:      %s\n","summary");
     fprintf(plankton->fptr_summ,"speciestype:   %s\n","plankton");
     fprintf(plankton->fptr_summ,"\n");
     print_plankton_header(plankton->fptr_summ);
     
     /*Close ouput files*/
     fclose(plankton->fptr_r);
     fclose(plankton->fptr_summ);
     
     /*Setup And Initialise Plankton Abundance Arrays*/
     int i,k,l;
     int m,x,y;
     m=grid->mnum;
     x=grid->xnum;
     y=grid->ynum;
     
     plankton->u_values=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     plankton->g_values=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     for(i=0 ; i<m ; i++){
             plankton->u_values[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             plankton->g_values[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             for(k=0 ; k<x ; k++){
                     plankton->u_values[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     plankton->g_values[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     for(l=0 ; l<y ; l++){
                             plankton->u_values[i][k][l]=0;
                             plankton->g_values[i][k][l]=0;
                     }
             }
     }
          
     /*Setup initial values for u(m,0,x,y)*/
     /*Use user defined values*/
     if(plankton->initial_flag==1 || plankton->ts_flag==1){
             /*Create plankton input filename*/
             snprintf(plankton->fname_ts, sizeof(plankton->fname_ts),"%s/Input/%s_%s.txt",run->filename,plankton->filename,"ts");
             plankton->fptr_ts=safe_fopen(plankton->fname_ts,"r",__LINE__,__FILE__);

             int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
             char *input;               //input string to store line as string
             input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
             char sep[]=",";            //identifying the ',' char as seperating values
             char *result = NULL;       //declaring a pointer for storing each individual value as a char
             
             /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
             for(k=0 ; k<x ; k++){
                     for(l=0 ; l<y ; l++){
                             fgets(input,MAX_LEN,plankton->fptr_ts);
                             input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                             
                             result=(char *)strtok(input,sep);    //this is the 1st value
                             i=0;
                             while(result!=NULL){
                                     plankton->u_values[i][k][l]=atof(result);
                                     result=(char *)strtok(NULL,sep);
                                     i++;
                             }
                     }
             }
             free(input);
             
             /*If this is only an intial value then close the file now*/
             if(plankton->ts_flag==0){
                     fclose(plankton->fptr_ts);
             }
     }
     /*Use standard values*/
     else{
             for(i=0 ; i<(plankton->iplamax+1) ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     plankton->u_values[i][k][l]=plankton->u_0*exp(plankton->lambda*grid->m_values[i]);
                                     if(run->spatial_dim==1 || run->spatial_dim==2){
                                              plankton->u_values[i][k][l]*=x_start(grid->x_values[k],0);
                                              if(run->spatial_dim==2){
                                                      plankton->u_values[i][k][l]*=y_start(grid->y_values[l],0);
                                              }
                                     }
                             }
                     }
             }
     }
     
}     

void setup_pelagic(RUN *run, GRID *grid, PELAGIC *pelagic, double *pel_params, char *filename, int rep_method, int initial_flag, int ts_flag, int fishing_flag)
{
     /*Pelagic Species Filename*/
     pelagic->filename=filename;
     
     /*Pelagic Species Size Range*/
     pelagic->pelmin=pel_params[0];
     pelagic->pelmat=pel_params[1];
     pelagic->pelmax=pel_params[2];
     pelagic->ipelmin=(int)rint((pelagic->pelmin-grid->mmin)/grid->mstep);
     pelagic->ipelmat=(int)rint((pelagic->pelmat-grid->mmin)/grid->mstep);
     pelagic->ipelmax=(int)rint((pelagic->pelmax-grid->mmin)/grid->mstep);

     /*Pelagic Species Constants*/
     pelagic->A=pel_params[3];           //Volume searched by unit weight
     pelagic->alpha=pel_params[4];      //Volume search allometry
     pelagic->mu_0=pel_params[5];        //natural mortality rate
     pelagic->beta=pel_params[6];      //natural mortality allometry
     pelagic->mu_s=pel_params[7];        //senescence mortality rate
     pelagic->epsilon=pel_params[8];     //senescece mortality constant
     pelagic->u_0=pel_params[9];         //initial intercept
     pelagic->lambda=pel_params[10];      //primary production slope
     
     /*Pelagic Species Conversion Efficiencies*/
     pelagic->K_pla=pel_params[11];           //growth efficiency
     pelagic->R_pla=pel_params[12];
     pelagic->Ex_pla=pel_params[13];
     pelagic->K_pel=pel_params[14];
     pelagic->R_pel=pel_params[15];
     pelagic->Ex_pel=pel_params[16];
     pelagic->K_ben=pel_params[17];
     pelagic->R_ben=pel_params[18];
     pelagic->Ex_ben=pel_params[19];
     
     /*Pelagic Species Attack Rates*/
     pelagic->pref_pla=pel_params[20];
     pelagic->pref_pel=pel_params[21];
     pelagic->pref_ben=pel_params[22];
     
     /*Pelagic Species Feeding Kernel*/
     pelagic->q_0=pel_params[23];        //modal ratio of predation sizes
     pelagic->sig=pel_params[24];         //predator-prey speciality
     pelagic->trunc=pel_params[25];
     
     /*Pelagic Species Movement Constants*/
     pelagic->prey=0;
     pelagic->pred=0;
     pelagic->comp=0;
     pelagic->gamma_prey=0;
     pelagic->gamma_pred=0;
     pelagic->gamma_comp=0;
     
     if(run->spatial_dim==1 || run->spatial_dim==2){
             pelagic->prey=pel_params[26];
             pelagic->pred=pel_params[27];
             pelagic->comp=pel_params[28];
             pelagic->gamma_prey=pel_params[29];
             pelagic->gamma_pred=pel_params[30];
             pelagic->gamma_comp=pel_params[31];
     }
     
     /*Pelagic Species Flags and Methods*/
     pelagic->rep_method=rep_method;
     pelagic->initial_flag=initial_flag;
     pelagic->ts_flag=ts_flag;
     pelagic->fishing_flag=fishing_flag;
     
     int s,n;
     n=run->no_pelagic;
     

     /*Create species output filenames*/
     snprintf(pelagic->fname_r, sizeof(pelagic->fname_r),"%s/%s/%s.txt",run->filename,pelagic->filename,"results");
     snprintf(pelagic->fname_g, sizeof(pelagic->fname_g),"%s/%s/%s.txt",run->filename,pelagic->filename,"growth");
     snprintf(pelagic->fname_m, sizeof(pelagic->fname_m),"%s/%s/%s.txt",run->filename,pelagic->filename,"mortality");
     pelagic->fname_pred=(char **)safe_malloc(n*sizeof(char *),__LINE__,__FILE__);
     for(s=0 ; s<n ; s++){
             pelagic->fname_pred[s]=(char *)safe_malloc(40*sizeof(char),__LINE__,__FILE__);
             snprintf(pelagic->fname_pred[s], 40,"%s/%s/%s%d.txt",run->filename,pelagic->filename,"predation",s);
     }
     snprintf(pelagic->fname_fish, sizeof(pelagic->fname_fish),"%s/%s/%s.txt",run->filename,pelagic->filename,"fishing");
     snprintf(pelagic->fname_summ, sizeof(pelagic->fname_summ),"%s/%s/%s.txt",run->filename,pelagic->filename,"summary");
     
     /*Open output files*/
     pelagic->fptr_r=safe_fopen(pelagic->fname_r,"w",__LINE__,__FILE__);
     pelagic->fptr_g=safe_fopen(pelagic->fname_g,"w",__LINE__,__FILE__);
     pelagic->fptr_m=safe_fopen(pelagic->fname_m,"w",__LINE__,__FILE__);
     pelagic->fptr_pred=(FILE **)safe_malloc(n*sizeof(FILE *),__LINE__,__FILE__);
     for(s=0 ; s<n ; s++){
             pelagic->fptr_pred[s]=safe_fopen(pelagic->fname_pred[s],"w",__LINE__,__FILE__);
     }
     pelagic->fptr_fish=safe_fopen(pelagic->fname_fish,"w",__LINE__,__FILE__);
     pelagic->fptr_summ=safe_fopen(pelagic->fname_summ,"w",__LINE__,__FILE__);

     /*Initialise output files*/
     fprintf(pelagic->fptr_r,"filetype:      %s\n","results");
     fprintf(pelagic->fptr_r,"speciestype:   %s\n","pelagic");
     fprintf(pelagic->fptr_r,"\n");
     print_mass_header(grid, pelagic->fptr_r);
     
     fprintf(pelagic->fptr_g,"filetype:      %s\n","growth");
     fprintf(pelagic->fptr_g,"speciestype:   %s\n","pelagic");
     fprintf(pelagic->fptr_g,"\n");
     print_mass_header(grid, pelagic->fptr_g);
     
     fprintf(pelagic->fptr_m,"filetype:      %s\n","mortality");
     fprintf(pelagic->fptr_m,"speciestype:   %s\n","pelagic");
     fprintf(pelagic->fptr_m,"\n");
     print_mass_header(grid, pelagic->fptr_m);
     
     for(s=0 ; s<n ; s++){
             fprintf(pelagic->fptr_pred[s],"filetype:      %s\n","predation");
             fprintf(pelagic->fptr_pred[s],"speciestype:   %s\n","pelagic");
             fprintf(pelagic->fptr_pred[s],"\n");
             print_mass_header(grid, pelagic->fptr_pred[s]);
     }
     
     fprintf(pelagic->fptr_fish,"filetype:      %s\n","fishing");
     fprintf(pelagic->fptr_fish,"speciestype:   %s\n","pelagic");
     fprintf(pelagic->fptr_fish,"\n");
     print_mass_header(grid, pelagic->fptr_fish);
     
     fprintf(pelagic->fptr_summ,"filetype:      %s\n","summary");
     fprintf(pelagic->fptr_summ,"speciestype:   %s\n","pelagic");
     fprintf(pelagic->fptr_summ,"\n");
     print_pelagic_header(pelagic->fptr_summ);
     
     /*Close output files*/
     fclose(pelagic->fptr_r);
     fclose(pelagic->fptr_g);
     fclose(pelagic->fptr_m);
     for(s=0 ; s<n ; s++){
             fclose(pelagic->fptr_pred[s]);
     }
     fclose(pelagic->fptr_fish);
     fclose(pelagic->fptr_summ);
     
     /*Setup Species Abundance Arrays*/
     int i,k,l;
     int m,x,y;
     m=grid->mnum;
     x=grid->xnum;
     y=grid->ynum;
     
     /*Setup 3D size-space arrays*/
     pelagic->mu_pred_values=(double****)safe_malloc(n*sizeof(double***),__LINE__,__FILE__);
     pelagic->u_values=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     pelagic->g_values=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     pelagic->mu_values=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     for(s=0 ; s<n ; s++){
             pelagic->mu_pred_values[s]=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     }
     pelagic->mu_fish_values=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     
     pelagic->pla_bio=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     pelagic->pel_bio=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     pelagic->ben_bio=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     
     pelagic->phi_values=(double *)safe_malloc(m*sizeof(double),__LINE__,__FILE__);
     for(i=0 ; i<m ; i++){
             pelagic->u_values[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             pelagic->g_values[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             pelagic->mu_values[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             for(s=0 ; s<n ; s++){
                     pelagic->mu_pred_values[s][i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             }
             pelagic->mu_fish_values[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             
             pelagic->pla_bio[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             pelagic->pel_bio[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             pelagic->ben_bio[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             
             pelagic->phi_values[i]=0;
             for(k=0 ; k<x ; k++){
                     pelagic->u_values[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     pelagic->g_values[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     pelagic->mu_values[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     for(s=0 ; s<n ; s++){
                             pelagic->mu_pred_values[s][i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     }
                     pelagic->mu_fish_values[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     
                     pelagic->pla_bio[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     pelagic->pel_bio[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     pelagic->ben_bio[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     for(l=0 ; l<y ; l++){
                             pelagic->u_values[i][k][l]=0;
                             pelagic->g_values[i][k][l]=0;
                             pelagic->mu_values[i][k][l]=0;
                             for(s=0 ; s<n ; s++){
                                     pelagic->mu_pred_values[s][i][k][l]=0;
                             }
                             pelagic->mu_fish_values[i][k][l]=0;
                             
                             pelagic->pla_bio[i][k][l]=0;
                             pelagic->pel_bio[i][k][l]=0;
                             pelagic->ben_bio[i][k][l]=0;
                     }
             }
     }
     
     /*Setup 2D space arrays*/
     pelagic->pla_total=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     pelagic->pel_total=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     pelagic->ben_total=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     
     pelagic->fish_total=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     pelagic->pred_total=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     pelagic->reproduction=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     for(k=0 ; k<x ; k++){
             pelagic->pla_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             pelagic->pel_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             pelagic->ben_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             
             pelagic->fish_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             pelagic->pred_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             pelagic->reproduction[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             for(l=0 ; l<y ; l++){
                     pelagic->pla_total[k][l]=0;
                     pelagic->pel_total[k][l]=0;
                     pelagic->ben_total[k][l]=0;
                     
                     pelagic->fish_total[k][l]=0;
                     pelagic->pred_total[k][l]=0;
                     pelagic->reproduction[k][l]=0;
             }
     }
          
     /*Setup initial values for u(m,0,x,y)*/
     if(pelagic->initial_flag==1){
             int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
             char *input;               //input string to store line as string
             input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
             char sep[]=",";            //identifying the ',' char as seperating values
             char *result = NULL;       //declaring a pointer for storing each individual value as a char
             
             /*Create plankton input filename*/
             snprintf(pelagic->fname_ts, sizeof(pelagic->fname_ts),"%s/Input/%s_%s.txt",run->filename,pelagic->filename,"ts");
             pelagic->fptr_ts=safe_fopen(pelagic->fname_ts,"r",__LINE__,__FILE__);
             
             /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
             for(k=0 ; k<x ; k++){
                     for(l=0 ; l<y ; l++){
                             fgets(input,MAX_LEN,pelagic->fptr_ts);
                             input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                             
                             result=(char *)strtok(input,sep);    //this is the 1st value
                             i=0;
                             while(result!=NULL){
                                     pelagic->u_values[i][k][l]=atof(result);
                                     result=(char *)strtok(NULL,sep);
                                     i++;
                             }
                     }
             }
             free(input);
             /*If this is only an intial value then close the file now*/
             if(pelagic->ts_flag==0){
                     fclose(pelagic->fptr_ts);
             }
     }
     else{
             for(i=pelagic->ipelmin ; i<(pelagic->ipelmax+1) ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     pelagic->u_values[i][k][l]=pelagic->u_0*exp(pelagic->lambda*grid->m_values[i]);
                                     if(run->spatial_dim==1 || run->spatial_dim==2){
                                              pelagic->u_values[i][k][l]*=x_start(grid->x_values[k],0);
                                              if(run->spatial_dim==2){
                                                      pelagic->u_values[i][k][l]*=y_start(grid->y_values[l],0);
                                              }
                                     }
                             }
                     }
             }
     }
     
     /*Open pelagic fishing file*/
     if(pelagic->fishing_flag==1){
             /*Create pelagic fishing input filename*/
             snprintf(pelagic->fname_fish_ts, sizeof(pelagic->fname_fish_ts),"%s/Input/%s_%s.txt",run->filename,pelagic->filename,"fishing_ts");
             pelagic->fptr_fish_ts=safe_fopen(pelagic->fname_fish_ts,"r",__LINE__,__FILE__);
     }
     
     /*Setup pelagic reproduction if needed*/
     if(pelagic->ts_flag==0){
             /*set eggs equal to intial values*/
             if(pelagic->rep_method==0){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     pelagic->reproduction[k][l]=pelagic->u_values[pelagic->ipelmin][k][l];
                             }
                     }
             }
             if(pelagic->rep_method==1){
                     /*Create pelagic fishing input filename*/
                     snprintf(pelagic->fname_rep_ts, sizeof(pelagic->fname_rep_ts),"%s/Input/%s_%s.txt",run->filename,pelagic->filename,"rep_ts");
                     pelagic->fptr_rep_ts=safe_fopen(pelagic->fname_rep_ts,"r",__LINE__,__FILE__);
             }
     }

     
     /*Setup values for predator-prey fedding kernel*/
     for(i=0 ; i<(pelagic->ipelmax+1) ; i++){
             pelagic->phi_values[i]=phi((grid->m_values[i]-grid->m_values[0]),pelagic->q_0, pelagic->sig, pelagic->trunc);
     }


}

void setup_benthic(RUN *run, GRID *grid, BENTHIC *benthic, double *ben_params, char *filename, int rep_method, int initial_flag, int ts_flag, int fishing_flag)
{
     /*Declare species name*/
     benthic->filename=filename;
     
     /*Benthic Species Size Range*/
     benthic->benmin=ben_params[0];
     benthic->benmat=ben_params[1];
     benthic->benmax=ben_params[2];
     benthic->ibenmin=(int)rint((benthic->benmin-grid->mmin)/grid->mstep);
     benthic->ibenmat=(int)rint((benthic->benmat-grid->mmin)/grid->mstep);
     benthic->ibenmax=(int)rint((benthic->benmax-grid->mmin)/grid->mstep);
     
     /*Benthic Species Constants*/
     benthic->A=ben_params[3];          //Volume searched by unit weight
     benthic->alpha=ben_params[4];      //Volume search allometry
     benthic->mu_0=ben_params[5];       //natural mortality rate
     benthic->beta=ben_params[6];       //natural mortality allometry
     benthic->mu_s=ben_params[7];       //senescence mortality rate
     benthic->epsilon=ben_params[8];    //senescece mortality constant
     benthic->u_0=ben_params[9];        //initial intercept
     benthic->lambda=ben_params[10];     //primary production slope
     
     /*Benthic Species Conversion Efficiencies*/
     benthic->K_det=ben_params[11];      //growth efficiency
     benthic->R_det=ben_params[12];
     benthic->Ex_det=ben_params[13];
     
     benthic->pref_det=ben_params[14];
     
     /*Benthic Species Flags and Methods*/
     benthic->rep_method=rep_method;
     benthic->initial_flag=initial_flag;
     benthic->ts_flag=ts_flag;
     benthic->fishing_flag=fishing_flag;

     /*Create species output filenames*/
     snprintf(benthic->fname_r, sizeof(benthic->fname_r),"%s/%s/%s.txt",run->filename,benthic->filename,"results");
     snprintf(benthic->fname_g, sizeof(benthic->fname_g),"%s/%s/%s.txt",run->filename,benthic->filename,"growth");
     snprintf(benthic->fname_m, sizeof(benthic->fname_m),"%s/%s/%s.txt",run->filename,benthic->filename,"mortality");
     snprintf(benthic->fname_pred, sizeof(benthic->fname_pred),"%s/%s/%s.txt",run->filename,benthic->filename,"predation");
     snprintf(benthic->fname_fish, sizeof(benthic->fname_fish),"%s/%s/%s.txt",run->filename,benthic->filename,"fishing");
     snprintf(benthic->fname_summ, sizeof(benthic->fname_summ),"%s/%s/%s.txt",run->filename,benthic->filename,"summary");
     
     /*Open output files*/
     benthic->fptr_r=safe_fopen(benthic->fname_r,"w",__LINE__,__FILE__);
     benthic->fptr_g=safe_fopen(benthic->fname_g,"w",__LINE__,__FILE__);
     benthic->fptr_m=safe_fopen(benthic->fname_m,"w",__LINE__,__FILE__);
     benthic->fptr_pred=safe_fopen(benthic->fname_pred,"w",__LINE__,__FILE__);
     benthic->fptr_fish=safe_fopen(benthic->fname_fish,"w",__LINE__,__FILE__);
     benthic->fptr_summ=safe_fopen(benthic->fname_summ,"w",__LINE__,__FILE__);

     /*Initialise output files*/
     fprintf(benthic->fptr_r,"filetype:      %s\n","results");
     fprintf(benthic->fptr_r,"speciestype:   %s\n","benthic");
     fprintf(benthic->fptr_r,"\n");
     print_mass_header(grid, benthic->fptr_r);
     
     fprintf(benthic->fptr_g,"filetype:      %s\n","growth");
     fprintf(benthic->fptr_g,"speciestype:   %s\n","benthic");
     fprintf(benthic->fptr_g,"\n");
     print_mass_header(grid, benthic->fptr_g);
     
     fprintf(benthic->fptr_m,"filetype:      %s\n","mortality");
     fprintf(benthic->fptr_m,"speciestype:   %s\n","benthic");
     fprintf(benthic->fptr_m,"\n");
     print_mass_header(grid, benthic->fptr_m);
     
     fprintf(benthic->fptr_pred,"filetype:      %s\n","predation");
     fprintf(benthic->fptr_pred,"speciestype:   %s\n","benthic");
     fprintf(benthic->fptr_pred,"\n");
     print_mass_header(grid, benthic->fptr_pred);
     
     fprintf(benthic->fptr_fish,"filetype:      %s\n","fishing");
     fprintf(benthic->fptr_fish,"speciestype:   %s\n","benthic");
     fprintf(benthic->fptr_fish,"\n");
     print_mass_header(grid, benthic->fptr_fish);
     
     fprintf(benthic->fptr_summ,"filetype:      %s\n","summary");
     fprintf(benthic->fptr_summ,"speciestype:   %s\n","benthic");
     fprintf(benthic->fptr_summ,"\n");
     print_benthic_header(benthic->fptr_summ);
     
     /*Close output files*/
     fclose(benthic->fptr_r);
     fclose(benthic->fptr_g);
     fclose(benthic->fptr_m);
     fclose(benthic->fptr_pred);
     fclose(benthic->fptr_fish);
     fclose(benthic->fptr_summ);
     
     /*Setup Species Abundance Arrays*/
     int i,k,l;
     int m,x,y;
     m=grid->mnum;
     x=grid->xnum;
     y=grid->ynum;
     
     benthic->u_values=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     benthic->g_values=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     benthic->mu_values=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     benthic->mu_pred_values=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     benthic->mu_fish_values=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     
     benthic->det_bio=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     for(i=0 ; i<m ; i++){
             benthic->u_values[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             benthic->g_values[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             benthic->mu_values[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             benthic->mu_pred_values[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             benthic->mu_fish_values[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             
             benthic->det_bio[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             for(k=0 ; k<x ; k++){
                     benthic->u_values[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     benthic->g_values[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     benthic->mu_values[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     benthic->mu_pred_values[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     benthic->mu_fish_values[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     
                     benthic->det_bio[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     for(l=0 ; l<y ; l++){
                             benthic->u_values[i][k][l]=0;
                             benthic->g_values[i][k][l]=0;
                             benthic->mu_values[i][k][l]=0;
                             benthic->mu_pred_values[i][k][l]=0;
                             benthic->mu_fish_values[i][k][l]=0;

                             benthic->det_bio[i][k][l]=0;
                     }
             }
     }
     benthic->det_total=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     
     benthic->fish_total=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     benthic->pred_total=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     benthic->reproduction=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     for(k=0 ; k<x ; k++){
             benthic->det_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             
             benthic->fish_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             benthic->pred_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             benthic->reproduction[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             for(l=0 ; l<y ; l++){
                     benthic->det_total[k][l]=0;
                     
                     benthic->fish_total[k][l]=0;
                     benthic->pred_total[k][l]=0;
                     benthic->reproduction[k][l]=0;
             }
     }
          
     /*Setup initial values for u(m,0,x,y)*/
     if(benthic->initial_flag==1){
             int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
             char *input;               //input string to store line as string
             input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
             char sep[]=",";            //identifying the ',' char as seperating values
             char *result = NULL;       //declaring a pointer for storing each individual value as a char
             
             /*Create plankton input filename*/
             snprintf(benthic->fname_ts, sizeof(benthic->fname_ts),"%s/Input/%s_%s.txt",run->filename,benthic->filename,"ts");
             benthic->fptr_ts=safe_fopen(benthic->fname_ts,"r",__LINE__,__FILE__);
             
             /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
             for(k=0 ; k<x ; k++){
                     for(l=0 ; l<y ; l++){
                             fgets(input,MAX_LEN,benthic->fptr_ts);
                             input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                             
                             result=(char *)strtok(input,sep);    //this is the 1st value
                             i=0;
                             while(result!=NULL){
                                     benthic->u_values[i][k][l]=atof(result);
                                     result=(char *)strtok(NULL,sep);
                                     i++;
                             }
                     }
             }
             free(input);    
             /*If this is only an intial value then close the file now*/
             if(benthic->ts_flag==0){
                     fclose(benthic->fptr_ts);
             }
     }
     else{
             for(i=benthic->ibenmin ; i<(benthic->ibenmax+1) ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     benthic->u_values[i][k][l]=benthic->u_0*exp(benthic->lambda*grid->m_values[i]);
                                     if(run->spatial_dim==1 || run->spatial_dim==2){
                                              benthic->u_values[i][k][l]*=x_start(grid->x_values[k],0);
                                              if(run->spatial_dim==2){
                                                      benthic->u_values[i][k][l]*=y_start(grid->y_values[l],0);
                                              }
                                     }
                             }
                     }
             }
     }
     
     /*Setup initial values for fishing mortality*/
     if(benthic->fishing_flag==1){
             /*Create benthic fishing input filename*/
             snprintf(benthic->fname_fish_ts, sizeof(benthic->fname_fish_ts),"%s/Input/%s_%s.txt",run->filename,benthic->filename,"fishing_ts");
             benthic->fptr_fish_ts=safe_fopen(benthic->fname_fish_ts,"r",__LINE__,__FILE__);
     }
     
     /*Setup benthic reproduction if needed*/
     if(benthic->ts_flag==0){
             /*set eggs equal to intial values*/
             if(benthic->rep_method==0){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     benthic->reproduction[k][l]=benthic->u_values[benthic->ibenmin][k][l];
                             }
                     }
             }
             if(benthic->rep_method==1){
                     /*Create benthic fishing input filename*/
                     snprintf(benthic->fname_rep_ts, sizeof(benthic->fname_rep_ts),"%s/Input/%s_%s.txt",run->filename,benthic->filename,"rep_ts");
                     benthic->fptr_rep_ts=safe_fopen(benthic->fname_rep_ts,"r",__LINE__,__FILE__);
             }
     }
}

void setup_detritus(RUN *run, GRID *grid, DETRITUS *detritus, double *det_params, char *filename, int initial_flag, int ts_flag)
{
         
     /*Setup Detritus Constants*/
     detritus->w_0=det_params[0];
    
     /*Declare detritus name*/
     detritus->filename=filename;
     
     /*Setup Plankton Flags*/
     detritus->initial_flag=initial_flag;
     detritus->ts_flag=ts_flag;
     
     /*Create detritus output filenames*/
     snprintf(detritus->fname_summ, sizeof(detritus->fname_summ),"%s/%s/%s.txt",run->filename,detritus->filename,"summary");
     
     /*Open output files*/
     detritus->fptr_summ=safe_fopen(detritus->fname_summ,"w",__LINE__,__FILE__);
     
     fprintf(detritus->fptr_summ,"filetype:      %s\n","summary");
     fprintf(detritus->fptr_summ,"speciestype:   %s\n","detritus");
     fprintf(detritus->fptr_summ,"\n");
     print_detritus_header(detritus->fptr_summ);
     
     /*Close output files*/
     fclose(detritus->fptr_summ);
     
     /*Setup And Initialise detritus Abundance Arrays*/
     int k,l;
     int x,y;
     x=grid->xnum;
     y=grid->ynum;
     
     detritus->w_values=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     detritus->g_values=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     detritus->mu_values=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     for(k=0 ; k<x ; k++){
             detritus->w_values[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             detritus->g_values[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             detritus->mu_values[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             for(l=0 ; l<y ; l++){
                     detritus->w_values[k][l]=0;
                     detritus->g_values[k][l]=0;
                     detritus->mu_values[k][l]=0;
             }
     }
          
     /*Setup initial values for w(0,x,y)*/
     /*Use user defined values*/
     if(detritus->initial_flag==1){
             int MAX_LEN=15; //this is the upper bound for the row length in characters
             char *input;               //input string to store line as string
             input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
             char sep[]=",";            //identifying the ',' char as seperating values
             char *result = NULL;       //declaring a pointer for storing each individual value as a char
             
             /*Create detritus input filename*/
             snprintf(detritus->fname_ts, sizeof(detritus->fname_ts),"%s/Input/%s_%s.txt",run->filename,detritus->filename,"ts");
             detritus->fptr_ts=safe_fopen(detritus->fname_ts,"r",__LINE__,__FILE__);
             
             /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
             for(k=0 ; k<x ; k++){
                     for(l=0 ; l<y ; l++){
                             fgets(input,MAX_LEN,detritus->fptr_ts);
                             input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                             
                             result=(char *)strtok(input,sep);    //this is the 1st value
                             while(result!=NULL){
                                     detritus->w_values[k][l]=atof(result);
                                     result=(char *)strtok(NULL,sep);
                             }
                     }
             }
             free(input);
             /*If this is only an intial value then close the file now*/
             if(detritus->ts_flag==0){
                     fclose(detritus->fptr_ts);
             }
     }
     /*Use standard values*/
     else{
          for(k=0 ; k<x ; k++){
                  for(l=0 ; l<y ; l++){
                          detritus->w_values[k][l]=detritus->w_0;
                  }
          }
     }
     
}
     
void setup_matrix(int dim, MATRIX *matrix)
/*Creates a matrix for the consumer size species*/
{
     matrix->size=dim;
     matrix->a=(double *)safe_malloc(dim*sizeof(double),__LINE__,__FILE__);
     matrix->b=(double *)safe_malloc(dim*sizeof(double),__LINE__,__FILE__);
     matrix->c=(double *)safe_malloc(dim*sizeof(double),__LINE__,__FILE__);
     matrix->r=(double *)safe_malloc(dim*sizeof(double),__LINE__,__FILE__);
     matrix->u=(double *)safe_malloc(dim*sizeof(double),__LINE__,__FILE__);

}


//-----------------------------------//
// Output Function. Not to be edited //
//-----------------------------------//

void calculate_results(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *pelmatrix, MATRIX *benmatrix, MATRIX *xmatrix, MATRIX *ymatrix)
/* Sets up output txt files and implements time stepping process*/
{
     int i,j,k,l,s,b;         //indices for each variable
     int m,t,x,y,n,c;         //maximum number of steps for each variable
     
     double time1=0,time2=0,timed=0;  //variables for countdown timer
     
     m=grid->mnum;
     t=grid->tnum;
     x=grid->xnum;
     y=grid->ynum;
     n=run->no_pelagic;
     c=run->no_benthic;
     
     /*Create temporary arrays for timestepping process*/
     double ***temp_plankton;
     double ****temp_pelagic;
     double ****temp_benthic;
     double **temp_detritus;
     
     temp_plankton=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     for(i=0 ; i< m ; i++){
             temp_plankton[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             for(k=0 ; k<x ; k++){
                     temp_plankton[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     for(l=0 ; l<y ; l++){
                             temp_plankton[i][k][l]=0;
                     }
             }
     }
     temp_pelagic=(double ****)safe_malloc(n*sizeof(double ***),__LINE__,__FILE__);
     for(s=0 ; s<n ; s++){
             temp_pelagic[s]=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
             for(i=0 ; i< m ; i++){
                     temp_pelagic[s][i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
                     for(k=0 ; k<x ; k++){
                             temp_pelagic[s][i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                             for(l=0 ; l<y ; l++){
                                     temp_pelagic[s][i][k][l]=0;
                             }
                     }
             }
     }
     
     if(run->no_benthic!=0){
             temp_benthic=(double ****)safe_malloc(c*sizeof(double ***),__LINE__,__FILE__);
             for(b=0 ; b<c ; b++){
                     temp_benthic[b]=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
                     for(i=0 ; i< m ; i++){
                             temp_benthic[b][i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
                             for(k=0 ; k<x ; k++){
                                     temp_benthic[b][i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                                     for(l=0 ; l<y ; l++){
                                             temp_benthic[b][i][k][l]=0;
                                     }
                             }
                     }
             }
             temp_detritus=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             for(k=0 ; k<x ; k++){
                     temp_detritus[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     for(l=0 ; l<y ; l++){
                             temp_detritus[k][l]=0;
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
     if(grid->toutmin == 0){
             print_timestep_plankton(grid, community->plankton, 0);
             print_plankton_summary(grid, community->plankton, 0);
     }

     /*Pelagic*/
     for(s=0 ; s<n ; s++){
             if(grid->toutmin == 0){
                     print_timestep_pelagic(run, grid, &(community->pelagic[s]), 0);
                     print_pelagic_summary(grid, &(community->pelagic[s]), 0);
             }
     }
     
     if(run->no_benthic!=0){
             /*Benthic*/
             for(b=0 ; b<c ; b++){
                     if(grid->toutmin == 0){
                             print_timestep_benthic(grid, &(community->benthic[b]), 0);
                             print_benthic_summary(grid, &(community->benthic[b]), 0);
                     }
             }
             /*Detritus*/
             if(grid->toutmin == 0){
                     print_detritus_summary(grid, community->detritus, 0);
             }
     }
     
     /*Parameter File*/
     run->fptr_summ=safe_fopen(run->fname_summ,"a",__LINE__,__FILE__);
     
     print_run(run, run->fptr_summ);
     print_grid(grid, run->fptr_summ);
     print_plankton(community->plankton, run->fptr_summ);
     for(s=0 ; s<n ; s++){
             print_pelagic(&(community->pelagic[s]),run->fptr_summ);
     }
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     print_benthic(&(community->benthic[b]),run->fptr_summ);
             }
             print_detritus(community->detritus, run->fptr_summ);
     }
     fclose(run->fptr_summ);

     /*Time Stepping Loop*/
     for(j=1 ; j<t ; j++){
                         
            /*Countdown timer for the program*/
             
            //Checks if user has terminated function in R (e.g. by pressing escape)
             R_CheckUserInterrupt();
             
             /*Countdown timer for the program*/
             if(j==1){
                      time2=(double)clock()/CLOCKS_PER_SEC;
             }
             else{
                      if(j%((int)rint(grid->toutstep/grid->tstep))==0){
                             time1=time2;
                             time2=(double)clock()/CLOCKS_PER_SEC;
                             timed=(time2-time1)*(grid->tstep/grid->toutstep);
                      }
                             Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                             R_FlushConsole();
                             Rprintf("%6.1f   %d",timed*(grid->tnum-j),grid->tnum-j);
                             R_FlushConsole();
             }
             if(j==(t-1)){
                      Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                      R_FlushConsole();
                      Rprintf("\n");
                      R_FlushConsole();
             }

             /*Store previous time step in temporary array for each species*/
             for(i=0 ; i<m ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     temp_plankton[i][k][l]=community->plankton->u_values[i][k][l];
                             }
                     }
             }
             for(s=0 ; s<n ; s++){
                     for(i=0 ; i<m ; i++){
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             temp_pelagic[s][i][k][l]=community->pelagic[s].u_values[i][k][l];
                                     }
                             }
                     }
             }
             
             if(run->no_benthic!=0){
                     for(b=0 ; b<c ; b++){
                             for(i=0 ; i<m ; i++){
                                     for(k=0 ; k<x ; k++){
                                             for(l=0 ; l<y ; l++){
                                                     temp_benthic[b][i][k][l]=community->benthic[b].u_values[i][k][l];
                                             }
                                     }
                             }
                     }
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     temp_detritus[k][l]=community->detritus->w_values[k][l];
                             }
                     }
             }


                          
             /*Calculate u(m,t) using previous timestep values*/
             mass_solver(run, grid, community, pelmatrix, benmatrix, temp_plankton, temp_pelagic, temp_benthic, temp_detritus);

             /*Calculate fishing pressures*/
             calculate_fishing(run, grid, community);
             
             /*Calculate g and mu using current timestep values for use in movement and next time step*/
             calculate_g_and_mu(run, grid, community);

             /*Calculate starting reproduction allocation*/
             calculate_reproduction(run, grid, community);
             


             /*If Spatial Model is being used*/
             if(run->spatial_dim==1){
             
                     /*Store itermediate time step as temp_pelagic*/

                     for(s=0 ; s<n ; s++){
                             for(i=0 ; i<m ; i++){
                                    for(k=0 ; k<x ; k++){
                                             for(l=0 ; l<y ; l++){
                                                     temp_pelagic[s][i][k][l]=community->pelagic[s].u_values[i][k][l];
                                             }
                                     }
                             }
                     }

                     
                     /*Calculate spatial movements using intermediate values*/
                     xmove_solver(run, grid, community, xmatrix, temp_pelagic, j);
             
                     
                    if(run->spatial_dim==2){
                     
                             /*Store itermediate time step as temp_pelagic*/
                             for(s=0 ; s<n ; s++){
                                     for(i=0 ; i<m ; i++){
                                             for(k=0 ; k<x ; k++){
                                                    for(l=0 ; l<y ; l++){
                                                             temp_pelagic[s][i][k][l]=community->pelagic[s].u_values[i][k][l];
                                                     }
                                             }
                                     }
                             }
                             /*Calculate spatial movements using intermediate values*/
                             ymove_solver(run, grid, community, ymatrix, temp_pelagic, j);
             
                     }
             }
             
             
             
             /*Print current timestep values to files*/
             if(grid->t_values[j]>=grid->toutmin && grid->t_values[j]<=grid->toutmax){
                     if(j%((int)rint(grid->toutstep/grid->tstep))==0){
                             
                             print_timestep_plankton(grid, community->plankton,j);
                             print_plankton_summary(grid, community->plankton, j);
                             for(s=0 ; s<n ; s++){
                                     print_timestep_pelagic(run, grid, &(community->pelagic[s]),j);
                                     print_pelagic_summary(grid, &(community->pelagic[s]),j);
                             }
                             if(run->no_benthic!=0){
                                     for(b=0 ; b<c  ;b++){
                                             print_timestep_benthic(grid, &(community->benthic[b]),j);
                                             print_benthic_summary(grid, &(community->benthic[b]),j);
                                     }
                                     print_detritus_summary(grid, community->detritus,j);
                             }
                     }
             }
     }
     
     /*Print final detailed timestep*/
     print_detailed_plankton(grid, community->plankton, t-1);
     for(s=0 ; s<n ; s++){
             print_detailed_pelagic(run, grid, &(community->pelagic[s]), t-1);
     }
     if(run->no_benthic!=0){
             for(b=0 ; b<c  ;b++){
                     print_detailed_benthic(grid, &(community->benthic[b]),t-1);
             }
     }

     /*Free pointers*/
     for(i=0 ; i<m ; i++){
             for(k=0 ; k<x ; k++){
                     free(temp_plankton[i][k]);
             }
             free(temp_plankton[i]);
     }
     free(temp_plankton);
     
     for(s=0 ; s<n ; s++){
             for(i=0 ; i<m ; i++){
                     for(k=0 ; k<x ; k++){
                             free(temp_pelagic[s][i][k]);
                     }
                     free(temp_pelagic[s][i]);
             }
             free(temp_pelagic[s]);
     }
     free(temp_pelagic);
     
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     for(i=0 ; i<m ; i++){
                             for(k=0 ; k<x ; k++){
                                     free(temp_benthic[b][i][k]);
                             }
                             free(temp_benthic[b][i]);
                     }
                     free(temp_benthic[b]);
             }
             free(temp_benthic);
             
             for(k=0 ; k<x ; k++){
                     free(temp_detritus[k]);
             }
             free(temp_detritus);
     }
     
     /*Close ts files*/
     if(community->plankton->ts_flag==1){
             fclose(community->plankton->fptr_ts);
     }
     for(s=0 ; s<n ; s++){
             if(community->pelagic[s].ts_flag==1){
                     fclose(community->pelagic[s].fptr_ts);
             }
     }
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     if(community->benthic[b].ts_flag==1){
                             fclose(community->benthic[b].fptr_ts);
                     }
             }
             if(community->detritus->ts_flag==1){
                     fclose(community->detritus->fptr_ts);
             }
     }
     
     /*Close fishing files*/
     for(s=0 ; s<n ; s++){
             if(community->pelagic[s].fishing_flag==1){
                     fclose(community->pelagic[s].fptr_fish_ts);
             }
     }
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     if(community->benthic[b].fishing_flag==1){
                             fclose(community->benthic[b].fptr_fish_ts);
                     }
             }
     }
     
     /*Close reproduction files*/
     for(s=0 ; s<n ; s++){
             if(community->pelagic[s].rep_method==1){
                     if(community->pelagic[s].ts_flag==0){
                             fclose(community->pelagic[s].fptr_rep_ts);
                     }
             }
     }
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     if(community->benthic[b].rep_method==1){
                             if(community->benthic[b].ts_flag==0){
                                     fclose(community->benthic[b].fptr_rep_ts);
                             }
                     }
             }
     }

}

void print_run(RUN *run, FILE *fptr)
{
     /*Output specific C code info*/
     fprintf(fptr,"Run Parameters\n");
     fprintf(fptr,"filename:      %s\n",run->filename);
     fprintf(fptr,"no_pelagic:    %d\n",run->no_pelagic);
     fprintf(fptr,"no_benthic:    %d\n",run->no_benthic);
     fprintf(fptr,"spatial_dim:   %d\n",run->spatial_dim);
     fprintf(fptr,"coupled_flag:  %d\n",run->coupled_flag);
     fprintf(fptr,"diff_method:   %d\n",run->diff_method);
     fprintf(fptr,"\n");
     
}     

void print_grid(GRID *grid, FILE *fptr)
{
     /*Output discretisation information*/
     fprintf(fptr,"Grid Parameters\n");
     fprintf(fptr,"mmin:          %f\n",grid->mmin);
     fprintf(fptr,"mmax:          %f\n",grid->mmax);
     fprintf(fptr,"mstep:         %f\n",grid->mstep);
     fprintf(fptr,"moutstep:      %f\n",grid->moutstep);
     
     fprintf(fptr,"t1:            %.10f\n",grid->t1);
     fprintf(fptr,"tmax:          %.10f\n",grid->tmax);
     fprintf(fptr,"tstep:         %.10f\n",grid->tstep);
     fprintf(fptr,"toutmin:       %.10f\n",grid->toutmin);
     fprintf(fptr,"toutmax:       %.10f\n",grid->toutmax);
     fprintf(fptr,"toutstep:      %.10f\n",grid->toutstep);
     
     fprintf(fptr,"xmin:          %f\n",grid->xmin);
     fprintf(fptr,"xmax:          %f\n",grid->xmax);
     fprintf(fptr,"xstep:         %f\n",grid->xstep);
     fprintf(fptr,"xoutstep:      %f\n",grid->xoutstep);
     
     fprintf(fptr,"ymin:          %f\n",grid->ymin);
     fprintf(fptr,"ymax:          %f\n",grid->ymax);
     fprintf(fptr,"ystep:         %f\n",grid->ystep);
     fprintf(fptr,"youtstep:      %f\n",grid->youtstep);
     fprintf(fptr,"\n");
     
}     

void print_plankton(PLANKTON *plankton, FILE *fptr)
{
/*Ouput plankton specific information*/
     fprintf(fptr,"Plankton Parameters\n");
     fprintf(fptr,"filename:      %s\n",plankton->filename);
     fprintf(fptr,"speciestype:   %s\n","plankton");
     fprintf(fptr,"mmin:        %f\n",plankton->plamin);
     fprintf(fptr,"mmax:        %f\n",plankton->plamax);
     fprintf(fptr,"mu_0:          %f\n",plankton->mu_0);
     fprintf(fptr,"beta:          %f\n",plankton->beta);
     fprintf(fptr,"u_0:           %f\n",plankton->u_0);
     fprintf(fptr,"lambda:        %f\n",plankton->lambda);
     fprintf(fptr,"initial_flag:  %d\n",plankton->initial_flag);
     fprintf(fptr,"ts_flag:       %d\n",plankton->ts_flag);
     fprintf(fptr,"\n");
}

void print_pelagic(PELAGIC *pelagic, FILE *fptr)
{
     /*Ouput species specific information*/
     fprintf(fptr,"Pelagic Parameters\n");
     fprintf(fptr,"filename:      %s\n",pelagic->filename);
     fprintf(fptr,"speciestype:   %s\n","pelagic");
     
     fprintf(fptr,"mmin:        %f\n",pelagic->pelmin);
     fprintf(fptr,"mmat:        %f\n",pelagic->pelmat);
     fprintf(fptr,"mmax:        %f\n",pelagic->pelmax);
     
     fprintf(fptr,"A:             %f\n",pelagic->A);
     fprintf(fptr,"alpha:         %f\n",pelagic->alpha);
     fprintf(fptr,"mu_0:          %f\n",pelagic->mu_0);
     fprintf(fptr,"beta:          %f\n",pelagic->beta);
     fprintf(fptr,"mu_s:          %f\n",pelagic->mu_s);
     fprintf(fptr,"epsilon:       %f\n",pelagic->epsilon);
     fprintf(fptr,"u_0:           %f\n",pelagic->u_0);
     fprintf(fptr,"lambda:        %f\n",pelagic->lambda);
     
     fprintf(fptr,"K_pla:         %f\n",pelagic->K_pla);
     fprintf(fptr,"R_pla:         %f\n",pelagic->R_pla);
     fprintf(fptr,"Ex_pla:        %f\n",pelagic->Ex_pla);
     fprintf(fptr,"K_pel:         %f\n",pelagic->K_pel);
     fprintf(fptr,"R_pel:         %f\n",pelagic->R_pel);
     fprintf(fptr,"Ex_pel:        %f\n",pelagic->Ex_pel);
     fprintf(fptr,"K_ben:         %f\n",pelagic->K_ben);
     fprintf(fptr,"R_ben:         %f\n",pelagic->R_ben);
     fprintf(fptr,"Ex_ben:        %f\n",pelagic->Ex_ben);
     
     fprintf(fptr,"pref_pla:      %f\n",pelagic->pref_pla);
     fprintf(fptr,"pref_pel:      %f\n",pelagic->pref_pel);
     fprintf(fptr,"pref_ben:      %f\n",pelagic->pref_ben);
     
     fprintf(fptr,"q_0:           %f\n",pelagic->q_0);
     fprintf(fptr,"sig:           %f\n",pelagic->sig);
     fprintf(fptr,"trunc:         %f\n",pelagic->trunc);
     
     fprintf(fptr,"prey:          %f\n",pelagic->prey);
     fprintf(fptr,"pred:          %f\n",pelagic->pred);
     fprintf(fptr,"comp:          %f\n",pelagic->comp);
     fprintf(fptr,"gamma_pred:    %f\n",pelagic->gamma_pred);
     fprintf(fptr,"gamma_prey:    %f\n",pelagic->gamma_prey);
     fprintf(fptr,"gamma_comp:    %f\n",pelagic->gamma_comp);
     
     fprintf(fptr,"rep_method:    %d\n",pelagic->rep_method);
     fprintf(fptr,"initial_flag:  %d\n",pelagic->initial_flag);
     fprintf(fptr,"ts_flag:       %d\n",pelagic->ts_flag);
     fprintf(fptr,"fishing_flag:  %d\n",pelagic->fishing_flag);
     fprintf(fptr,"\n");
     
}

void print_benthic(BENTHIC *benthic, FILE *fptr)
{
     /*Ouput species specific information*/
     fprintf(fptr,"Benthic Parameters\n");
     fprintf(fptr,"filename:      %s\n",benthic->filename);
     fprintf(fptr,"speciestype:   %s\n","benthic");
     
     fprintf(fptr,"mmin:        %f\n",benthic->benmin);
     fprintf(fptr,"mmat:        %f\n",benthic->benmat);
     fprintf(fptr,"mmax:        %f\n",benthic->benmax);
     
     fprintf(fptr,"A:             %f\n",benthic->A);
     fprintf(fptr,"alpha:         %f\n",benthic->alpha);
     fprintf(fptr,"mu_0:          %f\n",benthic->mu_0);
     fprintf(fptr,"beta:          %f\n",benthic->beta);
     fprintf(fptr,"mu_s:          %f\n",benthic->mu_s);
     fprintf(fptr,"epsilon:       %f\n",benthic->epsilon);
     fprintf(fptr,"u_0:           %f\n",benthic->u_0);
     fprintf(fptr,"lambda:        %f\n",benthic->lambda);
     
     fprintf(fptr,"K_det:         %f\n",benthic->K_det);
     fprintf(fptr,"R_det:         %f\n",benthic->R_det);
     fprintf(fptr,"Ex_det:        %f\n",benthic->Ex_det);
     
     fprintf(fptr,"pref_det:      %f\n",benthic->pref_det);
     
     fprintf(fptr,"rep_method:    %d\n",benthic->rep_method);
     fprintf(fptr,"initial_flag:  %d\n",benthic->initial_flag);
     fprintf(fptr,"ts_flag:       %d\n",benthic->ts_flag);
     fprintf(fptr,"fishing_flag:  %d\n",benthic->fishing_flag);
     fprintf(fptr,"\n");
}     

void print_detritus(DETRITUS *detritus, FILE *fptr)
{
     /*Ouput plankton specific information*/
     fprintf(fptr,"Detritus Parameters\n");
     fprintf(fptr,"filename:      %s\n",detritus->filename);
     fprintf(fptr,"speciestype:   %s\n","detritus");
     
     fprintf(fptr,"w_0:           %f\n",detritus->w_0);
     
     fprintf(fptr,"initial_flag:  %d\n",detritus->initial_flag);
     fprintf(fptr,"ts_flag:       %d\n",detritus->ts_flag);
     fprintf(fptr,"\n");
}

void print_mass_header(GRID *grid, FILE *fptr)
{
     /*Output Header Row*/
     int i;
     int m=grid->mnum;
     
     fprintf(fptr,"Ouput Data\n");
     fprintf(fptr,"t,x,y,");
     for(i=0 ; i<(m-1) ; i+=(int)rint(grid->moutstep/grid->mstep)){
             fprintf(fptr,"%f,",grid->m_values[i]);
     }
     fprintf(fptr,"%f\n",grid->m_values[m-1]);
}

void print_detailed_header(GRID *grid, FILE *fptr)
{
     int i;
     int m;
     m=grid->mnum;

     /*Print Header*/
     fprintf(fptr,"\nDetailed Timestep\n");
     fprintf(fptr,"t,x,y,");
     for(i=0 ; i<(m-1) ; i++){
             fprintf(fptr,"%f,",grid->m_values[i]);
     }
     fprintf(fptr,"%f\n",grid->m_values[m-1]);
}

void print_plankton_header(FILE *fptr)
{
     fprintf(fptr,"Ouput Data\n");
     fprintf(fptr,"t,x,y");
     fprintf(fptr,",biomass");
     fprintf(fptr,"\n");
}

void print_pelagic_header(FILE *fptr)
{
     fprintf(fptr,"Ouput Data\n");
     fprintf(fptr,"t,x,y");
     fprintf(fptr,",biomass");
     fprintf(fptr,",plankton");
     fprintf(fptr,",pelagic");
     fprintf(fptr,",benthic");
     fprintf(fptr,",predation");
     fprintf(fptr,",fishing");
     fprintf(fptr,",reproduction");
     fprintf(fptr,"\n");
}

void print_benthic_header(FILE *fptr)
{
     fprintf(fptr,"Ouput Data\n");
     fprintf(fptr,"t,x,y");
     fprintf(fptr,",biomass");
     fprintf(fptr,",detritus");
     fprintf(fptr,",predation");
     fprintf(fptr,",fishing");
     fprintf(fptr,",reproduction");
     fprintf(fptr,"\n");
}

void print_detritus_header(FILE *fptr)
{
     fprintf(fptr,"Ouput Data\n");
     fprintf(fptr,"t,x,y");
     fprintf(fptr,",biomass");
     fprintf(fptr,",detritus_in");
     fprintf(fptr,",detritus_out");
     fprintf(fptr,"\n");
}

void print_timestep_plankton(GRID *grid, PLANKTON *plankton, int time)
/*Outputs timestep values*/
{
     int i,k,l;
     int m,x,y;
     m=grid->mnum;
     x=grid->xnum;
     y=grid->ynum;
     
     /*Open output files*/
     plankton->fptr_r=safe_fopen(plankton->fname_r,"a",__LINE__,__FILE__);

     for(k=0 ; k<x ; k+=(int)rint(grid->xoutstep/grid->xstep)){
             for(l=0 ; l<y ; l+=(int)rint(grid->youtstep/grid->ystep)){
                     fprintf(plankton->fptr_r,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     for(i=0 ; i<(m-1) ; i+=(int)rint(grid->moutstep/grid->mstep)){
                             fprintf(plankton->fptr_r,"%.16f,",plankton->u_values[i][k][l]);
                     }
                     fprintf(plankton->fptr_r,"%.16f\n",plankton->u_values[m-1][k][l]);
             }
     }
     
     fclose(plankton->fptr_r);
}

void print_detailed_plankton(GRID *grid, PLANKTON *plankton, int time)
/*Outputs timestep values at high resolution*/
{     
     int i,k,l;
     int m,x,y;
     m=grid->mnum;
     x=grid->xnum;
     y=grid->ynum;
     
     /*Open output files*/
     plankton->fptr_r=safe_fopen(plankton->fname_r,"a",__LINE__,__FILE__);
     
     /*Print Header*/
     print_detailed_header(grid, plankton->fptr_r);

     for(k=0 ; k<x ; k++){
             for(l=0 ; l<y ; l++){
                     fprintf(plankton->fptr_r,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     for(i=0 ; i<(m-1) ; i++){
                             fprintf(plankton->fptr_r,"%.16f,",plankton->u_values[i][k][l]);
                     }
                     fprintf(plankton->fptr_r,"%.16f\n",plankton->u_values[m-1][k][l]);
             }
     }
     
     fclose(plankton->fptr_r);
}

void print_timestep_pelagic(RUN *run, GRID *grid, PELAGIC *pelagic, int time)
/*Outputs timestep values*/
{
     int s,i,k,l;
     int n,m,x,y;
     n=run->no_pelagic;
     m=grid->mnum;
     x=grid->xnum;
     y=grid->ynum;

     /*Open output files*/
     pelagic->fptr_r=safe_fopen(pelagic->fname_r,"a",__LINE__,__FILE__);
     pelagic->fptr_g=safe_fopen(pelagic->fname_g,"a",__LINE__,__FILE__);
     pelagic->fptr_m=safe_fopen(pelagic->fname_m,"a",__LINE__,__FILE__);
     for(s=0 ; s<n ; s++){
             pelagic->fptr_pred[s]=safe_fopen(pelagic->fname_pred[s],"a",__LINE__,__FILE__);
     }
     pelagic->fptr_fish=safe_fopen(pelagic->fname_fish,"a",__LINE__,__FILE__);

     for(k=0 ; k<x ; k+=(int)rint(grid->xoutstep/grid->xstep)){
             for(l=0 ; l<y ; l+=(int)rint(grid->youtstep/grid->ystep)){
                     fprintf(pelagic->fptr_r,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(pelagic->fptr_g,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(pelagic->fptr_m,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     for(s=0 ; s<n ; s++){
                             fprintf(pelagic->fptr_pred[s],"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     }
                     fprintf(pelagic->fptr_fish,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     for(i=0 ; i<(m-1) ; i+=(int)rint(grid->moutstep/grid->mstep)){
                             fprintf(pelagic->fptr_r,"%.16f,",pelagic->u_values[i][k][l]);
                             fprintf(pelagic->fptr_g,"%.16f,",pelagic->g_values[i][k][l]);
                             fprintf(pelagic->fptr_m,"%.16f,",pelagic->mu_values[i][k][l]);
                             for(s=0 ; s<n ; s++){
                                     fprintf(pelagic->fptr_pred[s],"%.16f,",pelagic->mu_pred_values[s][i][k][l]);
                             }
                             fprintf(pelagic->fptr_fish,"%.16f,",pelagic->mu_fish_values[i][k][l]);
                     }
                     fprintf(pelagic->fptr_r,"%.16f\n",pelagic->u_values[m-1][k][l]);
                     fprintf(pelagic->fptr_g,"%.16f\n",pelagic->g_values[m-1][k][l]);
                     fprintf(pelagic->fptr_m,"%.16f\n",pelagic->mu_values[m-1][k][l]);
                     for(s=0 ; s<n ; s++){
                             fprintf(pelagic->fptr_pred[s],"%.16f\n",pelagic->mu_pred_values[s][m-1][k][l]);
                     }
                     fprintf(pelagic->fptr_fish,"%.16f\n",pelagic->mu_fish_values[m-1][k][l]);
             }
     }
     
     /*Close output files*/
     fclose(pelagic->fptr_r);
     fclose(pelagic->fptr_g);
     fclose(pelagic->fptr_m);
     for(s=0 ; s<n ; s++){
             fclose(pelagic->fptr_pred[s]);
     }
     fclose(pelagic->fptr_fish);
}

void print_detailed_pelagic(RUN *run, GRID *grid, PELAGIC *pelagic, int time)
/*Outputs timestep values*/
{
     int s,i,k,l;
     int n,m,x,y;
     n=run->no_pelagic;
     m=grid->mnum;
     x=grid->xnum;
     y=grid->ynum;

     /*Open output files*/
     pelagic->fptr_r=safe_fopen(pelagic->fname_r,"a",__LINE__,__FILE__);
     pelagic->fptr_g=safe_fopen(pelagic->fname_g,"a",__LINE__,__FILE__);
     pelagic->fptr_m=safe_fopen(pelagic->fname_m,"a",__LINE__,__FILE__);
     for(s=0 ; s<n ; s++){
             pelagic->fptr_pred[s]=safe_fopen(pelagic->fname_pred[s],"a",__LINE__,__FILE__);
     }
     pelagic->fptr_fish=safe_fopen(pelagic->fname_fish,"a",__LINE__,__FILE__);

     /*Print Header*/
     print_detailed_header(grid, pelagic->fptr_r);
     print_detailed_header(grid, pelagic->fptr_g);
     print_detailed_header(grid, pelagic->fptr_m);
     for(s=0 ; s<n ; s++){
             print_detailed_header(grid, pelagic->fptr_pred[s]);
     }
     print_detailed_header(grid, pelagic->fptr_fish);


     for(k=0 ; k<x ; k++){
             for(l=0 ; l<y ; l++){
                     fprintf(pelagic->fptr_r,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(pelagic->fptr_g,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(pelagic->fptr_m,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     for(s=0 ; s<n ; s++){
                             fprintf(pelagic->fptr_pred[s],"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     }
                     fprintf(pelagic->fptr_fish,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     for(i=0 ; i<(m-1) ; i++){
                             fprintf(pelagic->fptr_r,"%.16f,",pelagic->u_values[i][k][l]);
                             fprintf(pelagic->fptr_g,"%.16f,",pelagic->g_values[i][k][l]);
                             fprintf(pelagic->fptr_m,"%.16f,",pelagic->mu_values[i][k][l]);
                             for(s=0 ; s<n ; s++){
                                     fprintf(pelagic->fptr_pred[s],"%.16f,",pelagic->mu_pred_values[s][i][k][l]);
                             }
                             fprintf(pelagic->fptr_fish,"%.16f,",pelagic->mu_fish_values[i][k][l]);
                     }
                     fprintf(pelagic->fptr_r,"%.16f\n",pelagic->u_values[m-1][k][l]);
                     fprintf(pelagic->fptr_g,"%.16f\n",pelagic->g_values[m-1][k][l]);
                     fprintf(pelagic->fptr_m,"%.16f\n",pelagic->mu_values[m-1][k][l]);
                     for(s=0 ; s<n ; s++){
                             fprintf(pelagic->fptr_pred[s],"%.16f\n",pelagic->mu_pred_values[s][m-1][k][l]);
                     }
                     fprintf(pelagic->fptr_fish,"%.16f\n",pelagic->mu_fish_values[m-1][k][l]);
             }
     }
     
     /*Close output files*/
     fclose(pelagic->fptr_r);
     fclose(pelagic->fptr_g);
     fclose(pelagic->fptr_m);
     for(s=0 ; s<n ; s++){
             fclose(pelagic->fptr_pred[s]);
     }
     fclose(pelagic->fptr_fish);
}


void print_timestep_benthic(GRID *grid, BENTHIC *benthic, int time)
/*Outputs timestep values*/
{
     int i,k,l;
     int m,x,y;
     m=grid->mnum;
     x=grid->xnum;
     y=grid->ynum;

     /*Open output files*/
     benthic->fptr_r=safe_fopen(benthic->fname_r,"a",__LINE__,__FILE__);
     benthic->fptr_g=safe_fopen(benthic->fname_g,"a",__LINE__,__FILE__);
     benthic->fptr_m=safe_fopen(benthic->fname_m,"a",__LINE__,__FILE__);
     benthic->fptr_pred=safe_fopen(benthic->fname_pred,"a",__LINE__,__FILE__);
     benthic->fptr_fish=safe_fopen(benthic->fname_fish,"a",__LINE__,__FILE__);

     for(k=0 ; k<x ; k+=(int)rint(grid->xoutstep/grid->xstep)){
             for(l=0 ; l<y ; l+=(int)rint(grid->youtstep/grid->ystep)){
                     fprintf(benthic->fptr_r,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(benthic->fptr_g,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(benthic->fptr_m,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(benthic->fptr_pred,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(benthic->fptr_fish,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     for(i=0 ; i<(m-1) ; i+=(int)rint(grid->moutstep/grid->mstep)){
                             fprintf(benthic->fptr_r,"%.16f,",benthic->u_values[i][k][l]);
                             fprintf(benthic->fptr_g,"%.16f,",benthic->g_values[i][k][l]);
                             fprintf(benthic->fptr_m,"%.16f,",benthic->mu_values[i][k][l]);
                             fprintf(benthic->fptr_pred,"%.16f,",benthic->mu_pred_values[i][k][l]);
                             fprintf(benthic->fptr_fish,"%.16f,",benthic->mu_fish_values[i][k][l]);
                     }
                     fprintf(benthic->fptr_r,"%.16f\n",benthic->u_values[m-1][k][l]);
                     fprintf(benthic->fptr_g,"%.16f\n",benthic->g_values[m-1][k][l]);
                     fprintf(benthic->fptr_m,"%.16f\n",benthic->mu_values[m-1][k][l]);
                     fprintf(benthic->fptr_pred,"%.16f\n",benthic->mu_pred_values[m-1][k][l]);
                     fprintf(benthic->fptr_fish,"%.16f\n",benthic->mu_fish_values[m-1][k][l]);
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
     int i,k,l;
     int m,x,y;
     m=grid->mnum;
     x=grid->xnum;
     y=grid->ynum;

     /*Open output files*/
     benthic->fptr_r=safe_fopen(benthic->fname_r,"a",__LINE__,__FILE__);
     benthic->fptr_g=safe_fopen(benthic->fname_g,"a",__LINE__,__FILE__);
     benthic->fptr_m=safe_fopen(benthic->fname_m,"a",__LINE__,__FILE__);
     benthic->fptr_pred=safe_fopen(benthic->fname_pred,"a",__LINE__,__FILE__);
     benthic->fptr_fish=safe_fopen(benthic->fname_fish,"a",__LINE__,__FILE__);

     /*Print Header*/
     print_detailed_header(grid, benthic->fptr_r);
     print_detailed_header(grid, benthic->fptr_g);
     print_detailed_header(grid, benthic->fptr_m);
     print_detailed_header(grid, benthic->fptr_pred);
     print_detailed_header(grid, benthic->fptr_fish);
     
     for(k=0 ; k<x ; k++){
             for(l=0 ; l<y ; l++){
                     fprintf(benthic->fptr_r,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(benthic->fptr_g,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(benthic->fptr_m,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(benthic->fptr_pred,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     fprintf(benthic->fptr_fish,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     for(i=0 ; i<(m-1) ; i++){
                             fprintf(benthic->fptr_r,"%.16f,",benthic->u_values[i][k][l]);
                             fprintf(benthic->fptr_g,"%.16f,",benthic->g_values[i][k][l]);
                             fprintf(benthic->fptr_m,"%.16f,",benthic->mu_values[i][k][l]);
                             fprintf(benthic->fptr_pred,"%.16f,",benthic->mu_pred_values[i][k][l]);
                             fprintf(benthic->fptr_fish,"%.16f,",benthic->mu_fish_values[i][k][l]);
                     }
                     fprintf(benthic->fptr_r,"%.16f\n",benthic->u_values[m-1][k][l]);
                     fprintf(benthic->fptr_g,"%.16f\n",benthic->g_values[m-1][k][l]);
                     fprintf(benthic->fptr_m,"%.16f\n",benthic->mu_values[m-1][k][l]);
                     fprintf(benthic->fptr_pred,"%.16f\n",benthic->mu_pred_values[m-1][k][l]);
                     fprintf(benthic->fptr_fish,"%.16f\n",benthic->mu_fish_values[m-1][k][l]);
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
     
     int i,k,l;
     int x,y;
     x=grid->xnum;
     y=grid->ynum;
     
     double temp;
     
     /*Open output files*/
     plankton->fptr_summ=safe_fopen(plankton->fname_summ,"a",__LINE__,__FILE__);
     
     for(k=0 ; k<x ; k+=(int)rint(grid->xoutstep/grid->xstep)){
             for(l=0 ; l<y ; l+=(int)rint(grid->youtstep/grid->ystep)){
                     fprintf(plankton->fptr_summ,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     
                     /*Plankton Biomass Totals*/
                     temp=0;
                     for(i=plankton->iplamin ; i<(plankton->iplamax+1) ; i++){
                             temp+=plankton->u_values[i][k][l]*exp(grid->m_values[i])*grid->mstep;
                     }
                     fprintf(plankton->fptr_summ,"%.16f",temp);
                     fprintf(plankton->fptr_summ,"\n");
             }
     }
     
     fclose(plankton->fptr_summ);
}

void print_pelagic_summary(GRID *grid, PELAGIC *pelagic, int time)
{
     
     int i,k,l;
     int x,y;
     x=grid->xnum;
     y=grid->ynum;
     
     double temp;
     
     /*Open output files*/
     pelagic->fptr_summ=safe_fopen(pelagic->fname_summ,"a",__LINE__,__FILE__);

     for(k=0 ; k<x ; k+=(int)rint(grid->xoutstep/grid->xstep)){
             for(l=0 ; l<y ; l+=(int)rint(grid->youtstep/grid->ystep)){
                     fprintf(pelagic->fptr_summ,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     
                     /*Pelagic Biomass Totals*/
                     temp=0;
                     for(i=pelagic->ipelmin ; i<(pelagic->ipelmax+1) ; i++){
                             temp+=pelagic->u_values[i][k][l]*exp(grid->m_values[i])*grid->mstep;
                     }
                     fprintf(pelagic->fptr_summ,",%.16f",temp);
                     
                     /*Pelagic Biomass Intake Rate*/
                     fprintf(pelagic->fptr_summ,",%.16f",pelagic->pla_total[k][l]);
                     fprintf(pelagic->fptr_summ,",%.16f",pelagic->pel_total[k][l]);
                     fprintf(pelagic->fptr_summ,",%.16f",pelagic->ben_total[k][l]);
                     
                     /*Biomass Lost to Predation*/
                     fprintf(pelagic->fptr_summ,",%.16f",pelagic->pred_total[k][l]);
                     
                     /*Biomass Lost to Fishing*/
                     fprintf(pelagic->fptr_summ,",%.16f",pelagic->fish_total[k][l]);
                     
                     /*Reproduction*/
                     fprintf(pelagic->fptr_summ,",%.16f",pelagic->reproduction[k][l]/grid->tstep);
                     fprintf(pelagic->fptr_summ,"\n");
             }
     }
     
     fclose(pelagic->fptr_summ);
}

void print_benthic_summary(GRID *grid, BENTHIC *benthic, int time)
{
     int i,k,l;
     int x,y;
     x=grid->xnum;
     y=grid->ynum;
     
     double temp;
     
     /*Open output files*/
     benthic->fptr_summ=safe_fopen(benthic->fname_summ,"a",__LINE__,__FILE__);

     for(k=0 ; k<x ; k+=(int)rint(grid->xoutstep/grid->xstep)){
             for(l=0 ; l<y ; l+=(int)rint(grid->youtstep/grid->ystep)){
                     fprintf(benthic->fptr_summ,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);
                     
                     /*Benthic Biomass Totals*/
                     temp=0;
                     for(i=benthic->ibenmin ; i<(benthic->ibenmax+1) ; i++){
                             temp+=benthic->u_values[i][k][l]*exp(grid->m_values[i])*grid->mstep;
                     }
                     fprintf(benthic->fptr_summ,",%.16f",temp);
                     
                     /*Benthic Biomass Intake*/
                     fprintf(benthic->fptr_summ,",%.16f",benthic->det_total[k][l]);
                     
                     /*Biomass Lost to Predation*/
                     fprintf(benthic->fptr_summ,",%.16f",benthic->pred_total[k][l]);
                     
                     /*Biomass Lost to Fishing*/
                     fprintf(benthic->fptr_summ,",%.16f",benthic->fish_total[k][l]);
                     
                     /*Reproduction*/
                     fprintf(benthic->fptr_summ,",%.16f",benthic->reproduction[k][l]/grid->tstep);
                     fprintf(benthic->fptr_summ,"\n");
             }
     }
     fclose(benthic->fptr_summ);
}

void print_detritus_summary(GRID *grid, DETRITUS *detritus, int time)
{
     int k,l;
     int x,y;
     x=grid->xnum;
     y=grid->ynum;
     
     /*Open output files*/
     detritus->fptr_summ=safe_fopen(detritus->fname_summ,"a",__LINE__,__FILE__);

     for(k=0 ; k<x ; k+=(int)rint(grid->xoutstep/grid->xstep)){
             for(l=0 ; l<y ; l+=(int)rint(grid->youtstep/grid->ystep)){
                     fprintf(detritus->fptr_summ,"%.10f,%f,%f,",(grid->t_values[time]),grid->x_values[k],grid->y_values[l]);                     
                     
                     /*Detritus Biomass Totals*/
                     fprintf(detritus->fptr_summ,",%.16f",detritus->w_values[k][l]);
                     
                     /*Detritus In*/
                     fprintf(detritus->fptr_summ,",%.16f",detritus->g_values[k][l]);
                     
                     /*Detritus Out*/
                     fprintf(detritus->fptr_summ,",%.16f",detritus->mu_values[k][l]);

                     fprintf(detritus->fptr_summ,"\n");
             }
     }
     
     /*Close output files*/
     fclose(detritus->fptr_summ);
}

//-------------------------------------------------//
// Differencing Scheme Functions. Not to be edited //
//-------------------------------------------------//

void mass_solver(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *pelmatrix, MATRIX *benmatrix, double ***temp_plankton, double ****temp_pelagic, double ****temp_benthic, double **temp_detritus)
/*Impliments the differencing scheme to obtain the values at the next time step*/
{
     int i,k,l,s,b;
     int x,y,n,c;
     double dt,dm;
     
     dm=grid->mstep;
     dt=grid->tstep;
     
     x=grid->xnum;
     y=grid->ynum;
     n=run->no_pelagic;
     c=run->no_benthic;
     
     /*Used for calculating detritus updates*/
     double **temp_det_g;
     double **temp_det_mu;
     
     temp_det_g=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     temp_det_mu=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     for(k=0 ; k<x ; k++){
             temp_det_g[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             temp_det_mu[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             for(l=0 ; l<y ; l++){
                     temp_det_g[k][l]=0;
                     temp_det_mu[k][l]=0;
             }
     }
     /*******************/
     /* Plankton System */
     /*******************/
     if(community->plankton->ts_flag==0){
             // Calculate values for Plankton Size Spectra using internal dynamics //
             for(i=community->plankton->iplamin ; i<(community->plankton->iplamax+1) ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     community->plankton->u_values[i][k][l]=temp_plankton[i][k][l];  //All plankton systems are invariant
                             }
                     }
             }
     }
     if(community->plankton->ts_flag==1){
             int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
             char *input;               //input string to store line as string
             input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
             char sep[]=",";            //identifying the ',' char as seperating values
             char *result = NULL;       //declaring a pointer for storing each individual value as a char
             
             /*Read in a single timestep*/
             /* each line one at a time (corresponding to spectra at a spatial point)*/
             for(k=0 ; k<x ; k++){
                     for(l=0 ; l<y ; l++){
                             fgets(input,MAX_LEN,community->plankton->fptr_ts);
                             input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                             
                             result=(char *)strtok(input,sep);    //this is the 1st value
                             i=0;
                             while(result!=NULL){
                                     community->plankton->u_values[i][k][l]=atof(result);
                                     result=(char *)strtok(NULL,sep);
                                     i++;
                             }
                     }
             }
             free(input);
     }
     
     /******************/
     /* Pelagic System */
     /******************/     
     for(s=0 ; s<n ; s++){
             if(community->pelagic[s].ts_flag==0){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     for(i=0 ; i<pelmatrix[s].size ; i++){
                                             pelmatrix[s].a[i]=-(dt/dm)*community->pelagic[s].g_values[i-1+community->pelagic[s].ipelmin][k][l];
                                             pelmatrix[s].b[i]=(1+(dt*community->pelagic[s].mu_values[i+community->pelagic[s].ipelmin][k][l])+(dt/dm)*community->pelagic[s].g_values[i+community->pelagic[s].ipelmin][k][l]);
                                             pelmatrix[s].c[i]=0;
                                             pelmatrix[s].r[i]=temp_pelagic[s][i+community->pelagic[s].ipelmin][k][l];
                                     }
                                     pelmatrix[s].a[0]=0;
                                     tridag(&(pelmatrix[s]));
                                     for(i=0 ; i<pelmatrix[s].size ; i++){
                                             community->pelagic[s].u_values[i+community->pelagic[s].ipelmin][k][l]=pelmatrix[s].u[i];
                                     }
                             }
                     }
             }
             if(community->pelagic[s].ts_flag==1){
                     int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
                     char *input;               //input string to store line as string
                     input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                     char sep[]=",";            //identifying the ',' char as seperating values
                     char *result = NULL;       //declaring a pointer for storing each individual value as a char
             
                     /*Read in a single timestep*/
                     /* each line one at a time (corresponding to spectra at a spatial point)*/
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     fgets(input,MAX_LEN,community->pelagic[s].fptr_ts);
                                     input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                                     
                                     result=(char *)strtok(input,sep);    //this is the 1st value
                                     i=0;
                                     while(result!=NULL){
                                             community->pelagic[s].u_values[i][k][l]=atof(result);
                                             result=(char *)strtok(NULL,sep);
                                             i++;
                                     }
                             }
                     }
                     free(input);
             }
     }
     
     /******************/
     /* Benthic System */
     /******************/
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     if(community->benthic[b].ts_flag==0){
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             for(i=0 ; i<benmatrix[b].size ; i++){
                                                     benmatrix[b].a[i]=-(dt/dm)*community->benthic[b].g_values[i-1+community->benthic[b].ibenmin][k][l];
                                                     benmatrix[b].b[i]=(1+(dt*community->benthic[b].mu_values[i+community->benthic[b].ibenmin][k][l])+(dt/dm)*community->benthic[b].g_values[i+community->benthic[b].ibenmin][k][l]);
                                                     benmatrix[b].c[i]=0;
                                                     benmatrix[b].r[i]=temp_benthic[b][i+community->benthic[b].ibenmin][k][l];
                                             }
                                             benmatrix[b].a[0]=0;
                                             tridag(&(benmatrix[b]));
                                             for(i=0 ; i<benmatrix[b].size ; i++){
                                                     community->benthic[b].u_values[i+community->benthic[b].ibenmin][k][l]=benmatrix[b].u[i];
                                             }
                                     }
                             }
                     }
                     if(community->benthic[b].ts_flag==1){
                             int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
                             char *input;               //input string to store line as string
                             input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                             char sep[]=",";            //identifying the ',' char as seperating values
                             char *result = NULL;       //declaring a pointer for storing each individual value as a char
                             
                             /*Read in a single timestep*/
                             /* each line one at a time (corresponding to spectra at a spatial point)*/
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             fgets(input,MAX_LEN,community->benthic[b].fptr_ts);
                                             input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                                             
                                             result=(char *)strtok(input,sep);    //this is the 1st value
                                             i=0;
                                             while(result!=NULL){
                                                     community->benthic[b].u_values[i][k][l]=atof(result);
                                                     result=(char *)strtok(NULL,sep);
                                                     i++;
                                             }
                                     }
                             }
                             free(input);
                     }
             }
             
             /*******************/
             /* Detritus System */
             /*******************/
             if(community->detritus->ts_flag==0){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
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
                                       double W0 = temp_detritus[k][l];
                                       double In = community->detritus->g_values[k][l];        /* input rate I */
                                       double loss = community->detritus->mu_values[k][l];     /* = kappa*W0   */
                                       double kappa = (W0 > 1e-30) ? loss / W0 : 0.0;
                                       if(kappa > 1e-30){
                                               double Wstar = In / kappa;
                                               community->detritus->w_values[k][l] = Wstar + (W0 - Wstar)*exp(-kappa*dt);
                                       }else{
                                               community->detritus->w_values[k][l] = W0 + In*dt;   /* kappa -> 0 limit */
                                       }
                                       /* defensive backstop only against a non-finite upstream rate */
                                       if(!isfinite(community->detritus->w_values[k][l]) || community->detritus->w_values[k][l]<0){
                                               community->detritus->w_values[k][l]=0;
                                       }
                                     }
                                     temp_det_g[k][l]=community->detritus->g_values[k][l];
                                     temp_det_mu[k][l]=community->detritus->mu_values[k][l];
                             }
                     }
             }
             if(community->detritus->ts_flag==1){
                     int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
                     char *input;               //input string to store line as string
                     input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                     char sep[]=",";            //identifying the ',' char as seperating values
                     char *result = NULL;       //declaring a pointer for storing each individual value as a char
             
                     /*Read in a single timestep*/
                     /* each line one at a time (corresponding to spectra at a spatial point)*/
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     fgets(input,MAX_LEN,community->detritus->fptr_ts);
                                     input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                                     
                                     result=(char *)strtok(input,sep);    //this is the 1st value
                                     while(result!=NULL){
                                             community->detritus->w_values[k][l]=atof(result);
                                             result=(char *)strtok(NULL,sep);
                                     }
                             }
                     }
                     free(input);
             }
     }
     
     //Calculate fully implicit method if required
     if(run->diff_method==1){
             //Use intermediate u values to calculate new g and mu for full implicit methods
             calculate_g_and_mu(run, grid, community);
     
             //Use new g and mu to claculate final u values
             for(s=0 ; s<n ; s++){
                     if(community->pelagic[s].ts_flag==0){
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             for(i=0 ; i<pelmatrix[s].size ; i++){
                                                     pelmatrix[s].a[i]=-(dt/dm)*community->pelagic[s].g_values[i-1+community->pelagic[s].ipelmin][k][l];
                                                     pelmatrix[s].b[i]=(1+(dt*community->pelagic[s].mu_values[i+community->pelagic[s].ipelmin][k][l])+(dt/dm)*community->pelagic[s].g_values[i+community->pelagic[s].ipelmin][k][l]);
                                                     pelmatrix[s].c[i]=0;
                                                     pelmatrix[s].r[i]=temp_pelagic[s][i+community->pelagic[s].ipelmin][k][l];
                                             }
                                             pelmatrix[s].a[0]=0;
                                             tridag(&(pelmatrix[s]));
                                             for(i=0 ; i<pelmatrix[s].size ; i++){
                                                     community->pelagic[s].u_values[i+community->pelagic[s].ipelmin][k][l]=pelmatrix[s].u[i];
                                             }
                                     }
                             }
                     }
             }
             
             if(run->no_benthic!=0){
                     for(b=0 ; b<c ; b++){
                             if(community->benthic[b].ts_flag==0){
                                     for(k=0 ; k<x ; k++){
                                             for(l=0 ; l<y ; l++){
                                                     for(i=0 ; i<benmatrix[b].size ; i++){
                                                             benmatrix[b].a[i]=-(dt/dm)*community->benthic[b].g_values[i-1+community->benthic[b].ibenmin][k][l];
                                                             benmatrix[b].b[i]=(1+(dt*community->benthic[b].mu_values[i+community->benthic[b].ibenmin][k][l])+(dt/dm)*community->benthic[b].g_values[i+community->benthic[b].ibenmin][k][l]);
                                                             benmatrix[b].c[i]=0;
                                                             benmatrix[b].r[i]=temp_benthic[b][i+community->benthic[b].ibenmin][k][l];
                                                     }
                                                     benmatrix[b].a[0]=0;
                                                     tridag(&(benmatrix[b]));
                                                     for(i=0 ; i<benmatrix[b].size ; i++){
                                                             community->benthic[b].u_values[i+community->benthic[b].ibenmin][k][l]=benmatrix[b].u[i];
                                                     }
                                             }
                                     }
                             }
                     }
                     
                     if(community->detritus->ts_flag==0){
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             community->detritus->w_values[k][l]=temp_detritus[k][l]+(dt/2)*((community->detritus->g_values[k][l]+temp_det_g[k][l])-(community->detritus->mu_values[k][l]+temp_det_mu[k][l]));
                                     }
                             }
                     }
             }
     }
     
     /***************************/
     /* Add reproduction values */ 
     /***************************/
     for(k=0 ; k<x ; k++){
             for(l=0 ; l<y ; l++){
                     for(s=0 ; s<n ; s++){
                             /*If values aren't given by timeseries*/
                             if(community->pelagic[s].ts_flag==0){
                                     /*Specifying Total Eggs*/
                                     if(community->pelagic[s].rep_method==0 || community->pelagic[s].rep_method==1){
                                             community->pelagic[s].u_values[community->pelagic[s].ipelmin][k][l]=community->pelagic[s].reproduction[k][l];
                                     }
                                     /*Adding Eggs*/
                                     if(community->pelagic[s].rep_method==2 || community->pelagic[s].rep_method==3){
                                             community->pelagic[s].u_values[community->pelagic[s].ipelmin][k][l]+=community->pelagic[s].reproduction[k][l];
                                     }
                             }
                     }
                     
                     if(run->no_benthic!=0){
                             for(b=0 ; b<c ; b++){
                                     /*If values aren't given by timeseries*/
                                     if(community->benthic[b].ts_flag==0){
                                             /*Constant Egg amount*/
                                             if(community->benthic[b].rep_method==0){
                                                     community->benthic[b].u_values[community->benthic[b].ibenmin][k][l]=community->benthic[b].u_0*exp(community->benthic[b].lambda*community->benthic[b].benmin);
                                             }
                                             /*Eggs produced according to biomass intake*/
                                             if(community->benthic[b].rep_method==1 || community->benthic[b].rep_method==2){
                                                     community->benthic[b].u_values[community->benthic[b].ibenmin][k][l]+=community->benthic[b].reproduction[k][l];
                                             }
                                     }
                             }
                     }
             }
     }
     
     /*Free pointers*/
     for(k=0 ; k<x ; k++){
             free(temp_det_g[k]);
             free(temp_det_mu[k]);
     }
     free(temp_det_g);
     free(temp_det_mu);
     
}

void xmove_solver(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *xmat, double ****temp_pelagic, int time)
{
     int i,k,l,s;
     int x,y,n;
     double dt,dx;
     
     double mu_ijk,mu_ij1k,mu_ijk1;
     double g_ijk,g_ij1k,g_ijk1;
     double C,D;
     
     dt=grid->tstep;
     dx=grid->xstep;
     
     x=grid->xnum;
     y=grid->ynum;
     n=run->no_pelagic;
     
     double *alpha;
     double *beta;
     double gam;
     
     alpha=(double *)safe_malloc(x*sizeof(double),__LINE__,__FILE__);
     beta=(double *)safe_malloc(x*sizeof(double),__LINE__,__FILE__);
     
     //calculate spatial movement (using CN system)
     for(s=0 ; s<n ; s++){
             for(i=community->pelagic[s].ipelmin ; i<(community->pelagic[s].ipelmax+1) ; i++){
                     C=Cfun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]));
                     D=Dfun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]));
                     gam=Diffun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]))*(dt/(dx*dx));
                     for(l=0 ; l<y ; l++){
                             for(k=1 ; k<x-1 ; k++){
                                     g_ij1k=community->pelagic[s].g_values[i][k-1][l];
                                     g_ijk=community->pelagic[s].g_values[i][k][l];
                                     g_ijk1=community->pelagic[s].g_values[i][k+1][l];
                                     mu_ij1k=community->pelagic[s].mu_values[i][k-1][l];
                                     mu_ijk=community->pelagic[s].mu_values[i][k][l];
                                     mu_ijk1=community->pelagic[s].mu_values[i][k+1][l];
                                     
                                     alpha[k]=(-1*(dt/(4*dx*dx)))*(C*(g_ijk1-g_ij1k)-D*(mu_ijk1-mu_ij1k));
                                     beta[k]=(-1*(dt/(dx*dx)))*(C*(g_ijk1+g_ij1k-2*g_ijk)-D*(mu_ijk1+mu_ij1k-2*mu_ijk));
                             }
                             
                             alpha[0]=((-1*dt*2)/(dx*dx))*(C*(community->pelagic[s].g_values[i][1][l]-community->pelagic[s].g_values[i][0][l])-D*(community->pelagic[s].mu_values[i][1][l]-community->pelagic[s].mu_values[i][0][l]));
                             alpha[x-1]=((-1*dt*2)/(dx*dx))*(C*(community->pelagic[s].g_values[i][x-2][l]-community->pelagic[s].g_values[i][x-1][l])-D*(community->pelagic[s].mu_values[i][x-2][l]-community->pelagic[s].mu_values[i][x-1][l]));
                             
                             /*Setup RHS CN Matrix*/           
                             xmat->a[0]=0;
                             xmat->b[0]=2+alpha[0]-(2*gam);
                             xmat->c[0]=(2*gam);
                             for(k=1 ; k<x-1 ; k++){
                                     xmat->a[k]=-alpha[k]+gam;
                                     xmat->b[k]=2+beta[k]-(2*gam);
                                     xmat->c[k]=alpha[k]+gam;
                             }
                             xmat->a[x-1]=(2*gam);
                             xmat->b[x-1]=2+alpha[x-1]-(2*gam);
                             xmat->c[x-1]=0;
                             
                             for(k=0 ; k<x ; k++){
                                     xmat->u[k]=temp_pelagic[s][i][k][l];
                             }
                             
                             /*Calculate xmat r[] column*/
                             trimul(xmat);
                             
                             /*Setup LHS CN Matrix*/
                             xmat->a[0]=0;
                             xmat->b[0]=2-alpha[0]+(2*gam);
                             xmat->c[0]=-(2*gam);
                             for(k=1 ; k<x-1 ; k++){
                                     xmat->a[k]=alpha[k]-gam;
                                     xmat->b[k]=2-beta[k]+(2*gam);
                                     xmat->c[k]=-alpha[k]-gam;
                             }
                             xmat->a[x-1]=-(2*gam);
                             xmat->b[x-1]=2-alpha[x-1]+(2*gam);
                             xmat->c[x-1]=0;
                             
                             /*Calculate xmat u[] column*/
                             tridag(xmat);
                             
                             for(k=0 ; k<x ; k++){
                                     community->pelagic[s].u_values[i][k][l]=xmat->u[k];
                             }
                     }
             }
     }
     free(alpha);
     free(beta);
}

void ymove_solver(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *ymat, double ****temp_pelagic, int time)
{
     int i,k,l,s;
     int x,y,n;
     double dt,dy;
     
     double mu_ijk,mu_ij1k,mu_ijk1;
     double g_ijk,g_ij1k,g_ijk1;
     double C,D;
     
     dt=grid->tstep;
     dy=grid->ystep;
     
     x=grid->xnum;
     y=grid->ynum;
     n=run->no_pelagic;
     
     double *alpha;
     double *beta;
     double gam;
     
     alpha=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
     beta=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
     
     //calculate spatial movement (using CN system)
     for(s=0 ; s<n ; s++){
             for(i=community->pelagic[s].ipelmin ; i<(community->pelagic[s].ipelmax+1) ; i++){
                     C=Cfun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]));
                     D=Dfun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]));
                     gam=Diffun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]))*(dt/(dy*dy));
                     for(k=0 ; k<x ; k++){
                             for(l=1 ; l<y-1 ; l++){
                                     g_ij1k=community->pelagic[s].g_values[i][k][l-1];
                                     g_ijk=community->pelagic[s].g_values[i][k][l];
                                     g_ijk1=community->pelagic[s].g_values[i][k][l+1];
                                     mu_ij1k=community->pelagic[s].mu_values[i][k][l-1];
                                     mu_ijk=community->pelagic[s].mu_values[i][k][l];
                                     mu_ijk1=community->pelagic[s].mu_values[i][k][l+1];
                                     
                                     alpha[l]=(-1*(dt/(4*dy*dy)))*(C*(g_ijk1-g_ij1k)-D*(mu_ijk1-mu_ij1k));
                                     beta[l]=(-1*(dt/(dy*dy)))*(C*(g_ijk1+g_ij1k-2*g_ijk)-D*(mu_ijk1+mu_ij1k-2*mu_ijk));
                             }
                             
                             alpha[0]=((-1*dt*2)/(dy*dy))*(C*(community->pelagic[s].g_values[i][k][1]-community->pelagic[s].g_values[i][k][0])-D*(community->pelagic[s].mu_values[i][k][1]-community->pelagic[s].mu_values[i][k][0]));
                             alpha[y-1]=((-1*dt*2)/(dy*dy))*(C*(community->pelagic[s].g_values[i][k][y-2]-community->pelagic[s].g_values[i][k][y-1])-D*(community->pelagic[s].mu_values[i][k][y-2]-community->pelagic[s].mu_values[i][k][y-1]));
                             
                             /*Setup RHS CN Matrix*/           
                             ymat->a[0]=0;
                             ymat->b[0]=2+alpha[0]-(2*gam);
                             ymat->c[0]=(2*gam);
                             for(l=1 ; l<y-1 ; l++){
                                     ymat->a[l]=-alpha[l]+gam;
                                     ymat->b[l]=2+beta[l]-(2*gam);
                                     ymat->c[l]=alpha[l]+gam;
                             }
                             ymat->a[y-1]=(2*gam);
                             ymat->b[y-1]=2+alpha[y-1]-(2*gam);
                             ymat->c[y-1]=0;
                             
                             for(l=0 ; l<y ; l++){
                                     ymat->u[l]=temp_pelagic[s][i][k][l];
                             }
                             
                             /*Calculate xmat r[] column*/
                             trimul(ymat);
                             
                             /*Setup LHS CN Matrix*/
                             ymat->a[0]=0;
                             ymat->b[0]=2-alpha[0]+(2*gam);
                             ymat->c[0]=-(2*gam);
                             for(l=1 ; l<y-1 ; l++){
                                     ymat->a[l]=alpha[l]-gam;
                                     ymat->b[l]=2-beta[l]+(2*gam);
                                     ymat->c[l]=-alpha[l]-gam;
                             }
                             ymat->a[y-1]=-(2*gam);
                             ymat->b[y-1]=2-alpha[y-1]+(2*gam);
                             ymat->c[y-1]=0;
                             
                             /*Calculate xmat u[] column*/
                             tridag(ymat);
                             
                             for(l=0 ; l<y ; l++){
                                     community->pelagic[s].u_values[i][k][l]=ymat->u[l];
                             }
                     }
             }
     }
     free(alpha);
     free(beta);
}

void trimul(MATRIX *matrix)
{
     int i,n;
     n=matrix->size;
     
     matrix->r[0]=(matrix->b[0]*matrix->u[0])+(matrix->c[0]*matrix->u[1]);
     for(i=1 ; i<n-1 ; i++){
             matrix->r[i]=(matrix->a[i]*matrix->u[i-1])+(matrix->b[i]*matrix->u[i])+(matrix->c[i]*matrix->u[i+1]);
     }
     matrix->r[n-1]=(matrix->a[n-1]*matrix->u[n-2])+(matrix->b[n-1]*matrix->u[n-1]);
}


void tridag(MATRIX *matrix)
/*Fast algorithm for inverting a tridiagonal matrix*/
{
     int j,n;
     n=matrix->size;
     double bet;
     double *gam;

     gam=(double *)safe_malloc(n*sizeof(double),__LINE__,__FILE__);
     if(matrix->b[0] == 0){
             error("error1");
     }
     bet=matrix->b[0];
     matrix->u[0]=matrix->r[0]/bet;
     for(j=1 ; j<=(n-1) ; j++){
              gam[j]=matrix->c[j-1]/bet;
              bet=matrix->b[j]-(matrix->a[j]*gam[j]);
              if(matrix->b[j] == 0){
                      error("error2");
              }
              matrix->u[j]=(matrix->r[j]-(matrix->a[j]*matrix->u[j-1]))/bet;
     }
     for(j=(n-2) ; j>=0 ; j--){
               matrix->u[j] -= (gam[j+1]*matrix->u[j+1]);
     }
     free(gam);
}


//--------------------------------------------//
//      Growth and Mortality Functions        //
// These contain functions that may be edited //
//--------------------------------------------//

void calculate_g_and_mu(RUN *run, GRID *grid, COMMUNITY *community)
{
     int i,k,l,s,ss,b;
     int x,y,n,c;
     
     x=grid->xnum;
     y=grid->ynum;
     n=run->no_pelagic;
     c=run->no_benthic;
     
     double pla,pel,ben,det,fish,pred;

     //Pelagic changes
     for(s=0 ; s<n ; s++){
             for(i=community->pelagic[s].ipelmin ; i<(community->pelagic[s].ipelmax+1) ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     for(ss=0 ; ss<n ; ss++){
                                              community->pelagic[s].mu_pred_values[ss][i][k][l]=mu_pel_pred(s, ss, i, k, l, run, grid, community);
                                     }
                                     community->pelagic[s].mu_values[i][k][l]=mu_pel(s, i, k, l, run, grid, community);
                                     
                                     community->pelagic[s].pla_bio[i][k][l]=pla_biomass(s, i, k, l, grid, community);
                                     community->pelagic[s].pel_bio[i][k][l]=pel_biomass(s, i, k, l, run, grid, community);
                                     community->pelagic[s].ben_bio[i][k][l]=ben_biomass(s, i, k, l, run, grid, community);
                                     community->pelagic[s].g_values[i][k][l]=g_pel(s, i, k, l, run, grid, community);
                             }
                     }
             }
     }
     
     //Calculate pelagic biomass totals
     for(s=0 ; s<n ; s++){
             for(k=0 ; k<x ; k++){
                     for(l=0 ; l<y ; l++){
                             pla=0;
                             pel=0;
                             ben=0;
                             fish=0;
                             pred=0;
                             for(i=community->pelagic[s].ipelmin ; i<(community->pelagic[s].ipelmax+1) ; i++){
                                     pla+=community->pelagic[s].pla_bio[i][k][l]*community->pelagic[s].u_values[i][k][l]*grid->mstep;
                                     pel+=community->pelagic[s].pel_bio[i][k][l]*community->pelagic[s].u_values[i][k][l]*grid->mstep;
                                     ben+=community->pelagic[s].ben_bio[i][k][l]*community->pelagic[s].u_values[i][k][l]*grid->mstep;
                                     fish+=community->pelagic[s].mu_fish_values[i][k][l]*community->pelagic[s].u_values[i][k][l]*exp(grid->m_values[i])*grid->mstep;
                                     for(ss=0 ; ss<n ; ss++){
                                              pred+=community->pelagic[s].mu_pred_values[ss][i][k][l]*community->pelagic[s].u_values[i][k][l]*exp(grid->m_values[i])*grid->mstep;
                                     }
                             }
                             community->pelagic[s].pla_total[k][l]=pla;
                             community->pelagic[s].pel_total[k][l]=pel;
                             community->pelagic[s].ben_total[k][l]=ben;
                             community->pelagic[s].fish_total[k][l]=fish;
                             community->pelagic[s].pred_total[k][l]=pred;
                     }
             }
     }

     //Benthic changes
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     for(i=community->benthic[b].ibenmin ; i<(community->benthic[b].ibenmax+1) ; i++){
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             community->benthic[b].mu_pred_values[i][k][l]=mu_ben_pred(b, i, k, l, run, grid, community);
                                             community->benthic[b].mu_values[i][k][l]=mu_ben(b, i, k, l, run, grid, community);
                                             
                                             community->benthic[b].det_bio[i][k][l]=det_biomass(b, i, k, l, run, grid, community);
                                             community->benthic[b].g_values[i][k][l]=g_ben(b, i, k, l, run, grid, community);
                                     }
                             }
                     }
             }
             //Calculate benthic biomass totals
             for(b=0 ; b<c ; b++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     det=0;
                                     fish=0;
                                     pred=0;
                                     for(i=community->benthic[b].ibenmin ; i<(community->benthic[b].ibenmax+1) ; i++){
                                             det+=community->benthic[b].det_bio[i][k][l]*community->benthic[b].u_values[i][k][l]*grid->mstep;
                                             fish+=community->benthic[b].mu_fish_values[i][k][l]*community->benthic[b].u_values[i][k][l]*exp(grid->m_values[i])*grid->mstep;
                                             pred+=community->benthic[b].mu_pred_values[i][k][l]*community->benthic[b].u_values[i][k][l]*exp(grid->m_values[i])*grid->mstep;
                                     }
                                     community->benthic[b].det_total[k][l]=det;
                                     community->benthic[b].fish_total[k][l]=fish;
                                     community->benthic[b].pred_total[k][l]=pred;
                             }
                     }
             }
             //Detritus Changes
             for(k=0 ; k<x ; k++){
                     for(l=0 ; l<y ; l++){
                             community->detritus->g_values[k][l]=g_det(k, l, run, grid, community);
                             community->detritus->mu_values[k][l]=mu_det(k, l, run, grid, community);
                     }
             }
     }
}

double pla_biomass(int species, int size, int xspace, int yspace, GRID *grid, COMMUNITY *community)
//Calculates biomass intake from plankton
{
       double ans=0;
       int i;
       
       /* Allows feeding from Plankton Size Spectra*/
       for(i=community->plankton->iplamin ; i<size ; i++){
               ans+=exp(grid->m_values[i])*community->pelagic[species].phi_values[size-i]*community->plankton->u_values[i][xspace][yspace]*grid->mstep;
       }
       ans*=community->pelagic[species].pref_pla*community->pelagic[species].A*exp(community->pelagic[species].alpha*grid->m_values[size]);
       return ans;
}

double pel_biomass(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
//Calculates biomass intake from Pelagic Size Spectra
{
       double ans=0;
       int i,s;
       
       /* Allows feeding from all Size Spectra*/
       if(run->coupled_flag==1){
               for(s=0 ; s<run->no_pelagic ; s++){
                       for(i=community->pelagic[s].ipelmin ; i<size ; i++){
                               ans+=exp(grid->m_values[i])*community->pelagic[species].phi_values[size-i]*community->pelagic[s].u_values[i][xspace][yspace]*grid->mstep;
                       }
               }
       }
       
       /* Allows feeding from own Size Spectra Only*/
       else{
              for(i=community->pelagic[species].ipelmin ; i<size ; i++){
                      ans+=exp(grid->m_values[i])*community->pelagic[species].phi_values[size-i]*community->pelagic[species].u_values[i][xspace][yspace]*grid->mstep;
              }
       }
            
       ans*=community->pelagic[species].pref_pel*community->pelagic[species].A*exp(community->pelagic[species].alpha*grid->m_values[size]);
       return ans;
}

double ben_biomass(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
//Calculates biomass intake from Pelagic Size Spectra
{
       double ans=0;
       int i,b;
       
       if(run->no_benthic!=0){
               /* Allows feeding from all Size Spectra*/
               if(run->coupled_flag==1){
                       for(b=0 ; b<run->no_benthic ; b++){
                               for(i=community->benthic[b].ibenmin ; i<size ; i++){
                                       ans+=exp(grid->m_values[i])*community->pelagic[species].phi_values[size-i]*community->benthic[b].u_values[i][xspace][yspace]*grid->mstep;
                               }
                       }
               }
               ans*=community->pelagic[species].pref_ben*community->pelagic[species].A*exp(community->pelagic[species].alpha*grid->m_values[size]);
       }
       
       return ans;
}

double det_biomass(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
//Calculates biomass intake from Detritus
{
       double ans=0;
       ans+=community->benthic[species].pref_det*community->benthic[species].A*exp(community->benthic[species].alpha*grid->m_values[size])*community->detritus->w_values[xspace][yspace];
       return ans;
}

double g_pel(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
//Calculates g for size class 'size' of species 'species' using biomass values
{
       double ans=0;
       ans+=community->pelagic[species].K_pla*exp(-grid->m_values[size])*community->pelagic[species].pla_bio[size][xspace][yspace];
       ans+=community->pelagic[species].K_pel*exp(-grid->m_values[size])*community->pelagic[species].pel_bio[size][xspace][yspace];
       ans+=community->pelagic[species].K_ben*exp(-grid->m_values[size])*community->pelagic[species].ben_bio[size][xspace][yspace];
       return(ans);
}

double g_ben(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
//Calculates g for size class 'size' of species 'species' using biomass values
{
       double ans=0;
       ans+=community->benthic[species].K_det*exp(-grid->m_values[size])*community->benthic[species].det_bio[size][xspace][yspace];
       return ans;
}

double g_det(int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
{
     int i,s;
     int n;

     n=run->no_pelagic;
     
     double ans=0;
     
     if(run->coupled_flag==1){
             for(s=0 ; s<n ; s++){     
                     /*This is the rate of det biomass in from pelagic defecation*/
                     ans+=(community->pelagic[s].Ex_pla*community->pelagic[s].pla_total[xspace][yspace]+community->pelagic[s].Ex_pel*community->pelagic[s].pel_total[xspace][yspace]+community->pelagic[s].Ex_ben*community->pelagic[s].ben_total[xspace][yspace]);
                     for(i=community->pelagic[s].ipelmin ; i<(community->pelagic[s].ipelmax+1) ; i++){
                             /*This is the rate of det biomass in from dead pelagic stuff*/
                             ans+=community->pelagic[s].mu_0*exp(community->pelagic[s].beta*grid->m_values[i])*community->pelagic[s].u_values[i][xspace][yspace]*exp(grid->m_values[i])*grid->mstep;
                             ans+=community->pelagic[s].mu_s*((log10(exp(grid->m_values[i]))-log10(exp(community->pelagic[s].pelmin)))/((log10(exp(community->pelagic[s].pelmax))+community->pelagic[s].epsilon)-log10(exp(grid->m_values[i]))))*community->pelagic[s].u_values[i][xspace][yspace]*exp(grid->m_values[i])*grid->mstep;
                     }
             }
     }
     for(i=community->plankton->iplamin ; i<(community->plankton->iplamax+1) ; i++){
             /*This is the rate of det biomass in from dead plankton*/
             ans+=community->plankton->mu_0*exp(community->plankton->beta*grid->m_values[i])*community->plankton->u_values[i][xspace][yspace]*exp(grid->m_values[i])*grid->mstep;
     }
     
     return(ans);
}

double mu_det(int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
{
     int b;
     int c;
     
     c=run->no_benthic;

     double ans=0;
     
     for(b=0 ; b<c ; b++){
             /*This is the rate of det biomass out feeding the benthos*/
             ans+=(community->benthic[b].det_total[xspace][yspace]);
     }
     
     return(ans);
}

double mu_pel_pred(int species, int pred, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
//Calculates mu due to predation for size class 'size' of species 'species'
{
       double ans=0;
       int i;
       
       /*Mortality due to Predation by all Size Spectra*/
       if(run->coupled_flag==1){
               for(i=size ; i<(community->pelagic[pred].ipelmax+1) ; i++){
                       ans+=(exp(community->pelagic[pred].alpha*grid->m_values[i]))*community->pelagic[pred].phi_values[i-size]*(community->pelagic[pred].u_values[i][xspace][yspace])*community->pelagic[pred].A*community->pelagic[pred].pref_pel*grid->mstep;
               }
       }
       /*Mortality due to Own Size Spectra Predation*/
       else{
            if(species==pred){
                    for(i=size ; i<(community->pelagic[species].ipelmax+1) ; i++){
                            ans+=(exp(community->pelagic[species].alpha*grid->m_values[i]))*community->pelagic[species].phi_values[i-size]*(community->pelagic[species].u_values[i][xspace][yspace])*community->pelagic[species].A*community->pelagic[species].pref_pel*grid->mstep;
                    }
            }
       }
       
       
       return ans;
}

double mu_pel(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
//Calculates mu for size class 'size' of species 'species'
{
       double ans=0;
       int s;
       int n;
       n=run->no_pelagic;
       
       for(s=0 ; s<n ; s++){
               ans+=community->pelagic[species].mu_pred_values[s][size][xspace][yspace];
       }
       ans+=community->pelagic[species].mu_fish_values[size][xspace][yspace];
       ans+=(community->pelagic[species].mu_0*exp(community->pelagic[species].beta*grid->m_values[size]));
       ans+=(community->pelagic[species].mu_s*((log10(exp(grid->m_values[size]))-log10(exp(community->pelagic[species].pelmin)))/((log10(exp(community->pelagic[species].pelmax))+community->pelagic[species].epsilon)-log10(exp(grid->m_values[size])))));
            
       return ans;
}

double mu_ben_pred(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
//Calculates mu due to predation for benthic species
{
       double ans=0;
       int i,s;
       
       /*Mortality due to Predation by all Size Spectra*/
       if(run->coupled_flag==1){
               for(s=0 ; s<run->no_pelagic ; s++){
                       for(i=size ; i<(community->pelagic[s].ipelmax+1) ; i++){
                                  ans=ans+(exp(community->pelagic[s].alpha*grid->m_values[i]))*community->pelagic[s].phi_values[i-size]*(community->pelagic[s].u_values[i][xspace][yspace])*community->pelagic[s].A*community->pelagic[s].pref_ben*grid->mstep;
                       }
               }
       }
       
       return ans;
}

double mu_ben(int species, int size, int xspace, int yspace, RUN *run, GRID *grid, COMMUNITY *community)
//Calculates mu for size class 'size' of species 'species'
{
       double ans=0;
       
       ans+=community->benthic[species].mu_pred_values[size][xspace][yspace];
       ans+=community->benthic[species].mu_fish_values[size][xspace][yspace];
       ans+=(community->benthic[species].mu_0*exp(community->benthic[species].beta*grid->m_values[size]));
       ans+=(community->benthic[species].mu_s*((log10(exp(grid->m_values[size]))-log10(exp(community->benthic[species].benmin)))/((log10(exp(community->benthic[species].benmax))+community->benthic[species].epsilon)-log10(exp(grid->m_values[size])))));
       
       return ans;
}

void calculate_fishing(RUN *run, GRID *grid, COMMUNITY *community)
{
     int i,k,l,s,b;
     int x,y,n,c;
     
     x=grid->xnum;
     y=grid->ynum;
     n=run->no_pelagic;
     c=run->no_benthic;
     
     /*Pelagic Fishing*/
     for(s=0 ; s<n ; s++){
             if(community->pelagic[s].fishing_flag==1){
                     int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
                     char *input;               //input string to store line as string
                     input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                     char sep[]=",";            //identifying the ',' char as seperating values
                     char *result = NULL;       //declaring a pointer for storing each individual value as a char
                     
                     /*Read in each line one at a time (corresponding to fishing mortality at a spatial point)*/
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     fgets(input,MAX_LEN,community->pelagic[s].fptr_fish_ts);
                                     input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                                     
                                     result=(char *)strtok(input,sep);    //this is the 1st value
                                     i=0;
                                     while(result!=NULL){
                                             community->pelagic[s].mu_fish_values[i][k][l]=atof(result);
                                             result=(char *)strtok(NULL,sep);
                                             i++;
                                     }
                             }
                     }
                     free(input);
             }
     }
     
     /*Benthic Fishing*/
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     if(community->benthic[b].fishing_flag==1){
                             int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
                             char *input;               //input string to store line as string
                             input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                             char sep[]=",";            //identifying the ',' char as seperating values
                             char *result = NULL;       //declaring a pointer for storing each individual value as a char
                             
                             /*Read in each line one at a time (corresponding to fishing mortality at a spatial point)*/
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             fgets(input,MAX_LEN,community->benthic[b].fptr_fish_ts);
                                             input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                                             
                                             result=(char *)strtok(input,sep);    //this is the 1st value
                                             i=0;
                                             while(result!=NULL){
                                                     community->benthic[b].mu_fish_values[i][k][l]=atof(result);
                                                     result=(char *)strtok(NULL,sep);
                                                     i++;
                                             }
                                     }
                             }
                             free(input);
                     }
             }
     }
}
                                                             
void calculate_reproduction(RUN *run, GRID *grid, COMMUNITY *community)
{
     int i,k,l,s,b;
     int x,y,n,c;
     
     double pla,pel,ben,det;
     
     x=grid->xnum;
     y=grid->ynum;
     n=run->no_pelagic;
     c=run->no_benthic;
     
     /*Pelagic Reproduction*/
     for(s=0 ; s<n ; s++){
             //read in egg values from file
             if(community->pelagic[s].rep_method==1){
                     int MAX_LEN=15; //this is the upper bound for the row length in characters
                     char *input;               //input string to store line as string
                     input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                     char sep[]=",";            //identifying the ',' char as seperating values
                     char *result = NULL;       //declaring a pointer for storing each individual value as a char
                     
                     /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     fgets(input,MAX_LEN,community->pelagic[s].fptr_rep_ts);
                                     input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                                     
                                     result=(char *)strtok(input,sep);    //this is the 1st value
                                     while(result!=NULL){
                                             community->pelagic[s].reproduction[k][l]=atof(result);
                                             result=(char *)strtok(NULL,sep);
                                     }
                             }
                     }
                     free(input);
             }
             //biomass energy allocation
             if(community->pelagic[s].rep_method==2){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     pla=0;
                                     pel=0;
                                     ben=0;
                                     for(i=community->pelagic[s].ipelmat ; i<(community->pelagic[s].ipelmax+1) ; i++){
                                             pla+=community->pelagic[s].pla_bio[i][k][l]*community->pelagic[s].u_values[i][k][l]*grid->mstep;
                                             pel+=community->pelagic[s].pel_bio[i][k][l]*community->pelagic[s].u_values[i][k][l]*grid->mstep;
                                             ben+=community->pelagic[s].ben_bio[i][k][l]*community->pelagic[s].u_values[i][k][l]*grid->mstep;
                                     }
                                     community->pelagic[s].reproduction[k][l]=((community->pelagic[s].R_pla*pla+community->pelagic[s].R_pel*pel+community->pelagic[s].R_ben*ben)*grid->tstep)/(grid->mstep*exp(community->pelagic[s].pelmin));
                             }
                     }
             }
             //fecundity calculation
             if(community->pelagic[s].rep_method==3){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     pel=0;
                                     for(i=community->pelagic[s].ipelmat ; i<(community->pelagic[s].ipelmax+1) ; i++){
                                             pel+=community->pelagic[s].u_values[i][k][l]*exp(grid->m_values[i])*grid->mstep;
                                     }
                                     community->pelagic[s].reproduction[k][l]=pel*492*grid->tstep;
                             }
                     }
             }
     }
     
     /*Benthic Reproduction*/
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     //read in egg values from file
                     if(community->benthic[b].rep_method==1){
                             int MAX_LEN=15; //this is the upper bound for the row length in characters
                             char *input;               //input string to store line as string
                             input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                             char sep[]=",";            //identifying the ',' char as seperating values
                             char *result = NULL;       //declaring a pointer for storing each individual value as a char
                             
                             /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             fgets(input,MAX_LEN,community->benthic[b].fptr_rep_ts);
                                             input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                                             
                                             result=(char *)strtok(input,sep);    //this is the 1st value
                                             while(result!=NULL){
                                                     community->benthic[b].reproduction[k][l]=atof(result);
                                                     result=(char *)strtok(NULL,sep);
                                             }
                                     }
                             }
                             free(input);
                     }
                     //biomass energy allocation
                     if(community->benthic[b].rep_method==2){
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             det=0;
                                             for(i=community->benthic[b].ibenmat ; i<(community->benthic[b].ibenmax+1) ; i++){
                                                     det+=community->benthic[b].det_bio[i][k][l]*community->benthic[b].u_values[i][k][l]*grid->mstep;
                                             }
                                             community->benthic[b].reproduction[k][l]=(community->benthic[b].R_det*det*grid->tstep)/(grid->mstep*exp(community->benthic[b].benmin));
                                     }
                             }
                     }
                     //fecundity calculation
                     if(community->benthic[b].rep_method==3){
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             ben=0;
                                             for(i=community->benthic[b].ibenmat ; i<(community->benthic[b].ibenmax+1) ; i++){
                                                     ben+=community->benthic[b].u_values[i][k][l]*exp(grid->m_values[i])*grid->mstep;
                                             }
                                             community->benthic[b].reproduction[k][l]=ben*492*grid->tstep;
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
            
       if(q<=0 || q<(q_0-(trunc*sig)) || q>(q_0+(trunc*sig))){
            ans=0;
       }
       else{
            ans=exp(-(q-q_0)*(q-q_0)/(2*sig*sig));
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
       if(t>grid->t1){
                ans=pelagic->prey*exp(pelagic->gamma_prey*m);
       }
       else{
            ans=0;
       }
       return ans;
}


double Dfun(double m, double t, GRID *grid, PELAGIC *pelagic)
/*function determining predator avoidance tendancy*/
{
       double ans;
       if(t>grid->t1){
                ans=pelagic->pred*exp(pelagic->gamma_pred*m);
       }
       else{
            ans=0;
       }
       return ans;
}


double Diffun(double m, double t, GRID *grid, PELAGIC *pelagic)
/*function determining competition coefficient*/
{
       double ans;
       if(t>grid->t1){
                ans=pelagic->comp*exp(pelagic->gamma_comp*m);
       }
       else{
            ans=0;
       }
       return ans;
}


double x_start(double x, double t)
/*function determining initial plankton distribution*/
{
       double ans;
       ans=2*exp(-((x-500)*(x-500)/20000));
       return ans;
}

double y_start(double y, double t)
/*function determining initial plankton distribution*/
{
       double ans;
       ans=2*exp(-((y-500)*(y-500)/20000));
       return ans;
}


//---------------------------------------------------------//
// These are memory allocation functions. Not to be edited // 
//---------------------------------------------------------//

void free_mem(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *pelmatrix, MATRIX *benmatrix, MATRIX *xmatrix, MATRIX *ymatrix, double **temp_pel_params, double **temp_ben_params)
/*Frees all memory allocated within code*/
{
     int i,k,s,ss,b;
     int m,x,n,c;
     m=grid->mnum;
     x=grid->xnum;
     n=run->no_pelagic;
     c=run->no_benthic;
     
     /*Free grid values*/
     free(grid->m_values);
     free(grid->t_values);
     free(grid->x_values);
     free(grid->y_values);
     
     /*Free pelagic values*/
             for(s=0 ; s<n ; s++){
                     for(i=0 ; i<m ; i++){
                             for(k=0 ; k<x ; k++){
                                     free(community->pelagic[s].u_values[i][k]);
                                     free(community->pelagic[s].g_values[i][k]);
                                     free(community->pelagic[s].mu_values[i][k]);
                                     for(ss=0 ; ss<n ; ss++){
                                              free(community->pelagic[s].mu_pred_values[ss][i][k]);
                                     }
                                     free(community->pelagic[s].mu_fish_values[i][k]);
                                     free(community->pelagic[s].pla_bio[i][k]);
                                     free(community->pelagic[s].pel_bio[i][k]);
                                     free(community->pelagic[s].ben_bio[i][k]);
                             }
                             free(community->pelagic[s].u_values[i]);
                             free(community->pelagic[s].g_values[i]);
                             free(community->pelagic[s].mu_values[i]);
                             for(ss=0 ; ss<n ; ss++){
                                      free(community->pelagic[s].mu_pred_values[ss][i]);
                             }
                             free(community->pelagic[s].mu_fish_values[i]);
                             free(community->pelagic[s].pla_bio[i]);
                             free(community->pelagic[s].pel_bio[i]);
                             free(community->pelagic[s].ben_bio[i]);
                     }
                     free(community->pelagic[s].u_values);
                     free(community->pelagic[s].g_values);
                     free(community->pelagic[s].mu_values);
                     for(ss=0 ; ss<n ; ss++){
                              free(community->pelagic[s].mu_pred_values[ss]);
                     }
                     free(community->pelagic[s].mu_fish_values);
                     free(community->pelagic[s].pla_bio);
                     free(community->pelagic[s].pel_bio);
                     free(community->pelagic[s].ben_bio);
                     free(community->pelagic[s].phi_values);
                     free(temp_pel_params[s]);
                     
                     free(community->pelagic[s].fname_pred);
                     free(community->pelagic[s].fptr_pred);
             }
             free(temp_pel_params);
     
             for(s=0 ; s<n ; s++){
                     for(k=0 ; k<x ; k++){
                             free(community->pelagic[s].pla_total[k]);
                             free(community->pelagic[s].pel_total[k]);
                             free(community->pelagic[s].ben_total[k]);
                             free(community->pelagic[s].fish_total[k]);
                             free(community->pelagic[s].pred_total[k]);
                             free(community->pelagic[s].reproduction[k]);
                     }
                     free(community->pelagic[s].pla_total);
                     free(community->pelagic[s].pel_total);
                     free(community->pelagic[s].ben_total);
                     free(community->pelagic[s].fish_total);
                     free(community->pelagic[s].pred_total);
                     free(community->pelagic[s].reproduction);
             }
     free(community->pelagic);
     
     /*Free Plankton pointer*/     
             for(i=0 ; i<m ; i++){
                     for(k=0 ; k<x ; k++){
                             free(community->plankton->u_values[i][k]);
                             free(community->plankton->g_values[i][k]);
                     }
                     free(community->plankton->u_values[i]);
                     free(community->plankton->g_values[i]);
             }
             free(community->plankton->u_values);
             free(community->plankton->g_values);
     free(community->plankton);
     
     if(run->no_benthic!=0){
             /*Free Benthic pointer*/
             for(b=0 ; b<c ; b++){
                     for(i=0 ; i<m ; i++){
                             for(k=0 ; k<x ; k++){
                                     free(community->benthic[b].u_values[i][k]);
                                     free(community->benthic[b].g_values[i][k]);
                                     free(community->benthic[b].mu_values[i][k]);
                                     free(community->benthic[b].mu_pred_values[i][k]);
                                     free(community->benthic[b].mu_fish_values[i][k]);
                                     free(community->benthic[b].det_bio[i][k]);
                             }
                             free(community->benthic[b].u_values[i]);
                             free(community->benthic[b].g_values[i]);
                             free(community->benthic[b].mu_values[i]);
                             free(community->benthic[b].mu_pred_values[i]);
                             free(community->benthic[b].mu_fish_values[i]);
                             free(community->benthic[b].det_bio[i]);
                     }
                     free(community->benthic[b].u_values);
                     free(community->benthic[b].g_values);
                     free(community->benthic[b].mu_values);
                     free(community->benthic[b].mu_pred_values);
                     free(community->benthic[b].mu_fish_values);
                     free(community->benthic[b].det_bio);
                     free(temp_ben_params[b]);
             }
             free(temp_ben_params);
             for(b=0 ; b<c ; b++){
                     for(k=0 ; k<x ; k++){
                             free(community->benthic[b].det_total[k]);
                             free(community->benthic[b].fish_total[k]);
                             free(community->benthic[b].pred_total[k]);
                             free(community->benthic[b].reproduction[k]);
                     }
                     free(community->benthic[b].det_total);
                     free(community->benthic[b].fish_total);
                     free(community->benthic[b].pred_total);
                     free(community->benthic[b].reproduction);
             }
     free(community->benthic);
             
     /*Free Detritus*/
             for(k=0 ; k<x ; k++){
                     free(community->detritus->w_values[k]);
                     free(community->detritus->g_values[k]);
             }
             free(community->detritus->w_values);
             free(community->detritus->g_values);
     free(community->detritus);
     }
     
     /*Free matrix pointers*/
     for(s=0 ; s<n ; s++){
             free(pelmatrix[s].a);
             free(pelmatrix[s].b);
             free(pelmatrix[s].c);
             free(pelmatrix[s].r);
             free(pelmatrix[s].u);
     }
     free(pelmatrix);
     
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     free(benmatrix[b].a);
                     free(benmatrix[b].b);
                     free(benmatrix[b].c);
                     free(benmatrix[b].r);
                     free(benmatrix[b].u);
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
     ptr=malloc(size);
     if(ptr==NULL){
          error("Out of memory at line %d file %s\n",line,filename);
     }
     return ptr;
}


FILE *safe_fopen(char *filename, char *mode, int line, char *fname)
/*Safe version of fopen with built in error checking. Exits if the file cannot be opened*/
{
     FILE *fptr;
     fptr=fopen(filename,mode);
     if(fptr ==NULL){
             error("Unable to open file at line %d file %s\n",line,fname);
     }
     return fptr;
}
