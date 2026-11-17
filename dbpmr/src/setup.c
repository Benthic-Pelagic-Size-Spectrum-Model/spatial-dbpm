#include "setup.h"
#include "file_outputs.h"
#include "memory_inputs.h"


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
     
     /*Text output files (not created when the output goes to R as arrays only)*/
     if(run_writes_files(run)){
             /*Open output files*/
             run->fptr_summ=safe_fopen(run->fname_summ,"w",__LINE__,__FILE__);
     
             fprintf(run->fptr_summ,"filetype:      %s\n","parameters");
             fprintf(run->fptr_summ,"speciestype:   %s\n","NA");
             fprintf(run->fptr_summ,"\n");
     
             /*Close ouput files*/
             fclose(run->fptr_summ);
     }

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
     grid->exp_m_values=(double *)safe_malloc(m*sizeof(double),__LINE__,__FILE__);
     grid->t_values=(double *)safe_malloc(t*sizeof(double),__LINE__,__FILE__);
     grid->x_values=(double *)safe_malloc(x*sizeof(double),__LINE__,__FILE__);
     grid->y_values=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
     
     /*Setup values for m*/ 
     for(i=0 ; i<m ; i++){
             grid->m_values[i]=grid->mmin+(i*grid->mstep);
             grid->exp_m_values[i]=exp(grid->m_values[i]);
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

void setup_plankton(RUN *run, GRID *grid, PLANKTON *plankton, double *pla_params, char *filename, int initial_flag, int ts_flag, const ARRAY_INPUTS *arr)
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
     
     /*Text output files (not created when the output goes to R as arrays only)*/
     if(run_writes_files(run)){
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
     }
     
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
     /*Plankton series supplied from R as intercept/slope arrays: no _ts file is opened or read*/
     plankton->arr=NULL;
     plankton->ts_block=0;
     if(arr!=NULL && arr->pla_int!=NULL){
             if(plankton->ts_flag!=1){
                     error("plankton arrays were supplied but the plankton ts_flag is not set\n");
             }
             if(arr->pla_len!=(R_xlen_t)x*(R_xlen_t)y*(R_xlen_t)arr->pla_nlayer){
                     error("plankton arrays hold %ld values; expected x*y*layers = %d*%d*%d\n",(long)arr->pla_len,x,y,arr->pla_nlayer);
             }
             plankton->arr=arr;
             plankton_from_arrays(grid, plankton);
     }
     else if(plankton->in!=NULL && plankton->in->series!=NULL){
             /*Plankton abundance forced through time from an R array [x][y][plankton size class][layer]: no _ts file is read*/
             series_check(plankton->in, plankton->filename, plankton->ts_flag, x, y, plankton->iplamax-plankton->iplamin+1);
             series_to_u(plankton->in, 0, x, y, plankton->iplamin, plankton->iplamax-plankton->iplamin+1, plankton->u_values);
     }
     else if(plankton->in!=NULL && plankton->in->initial!=NULL){
             /*Starting plankton abundance supplied from R as an array [x][y][plankton size class]: no _ts file is read*/
             int np=plankton->iplamax-plankton->iplamin+1;
             if(plankton->initial_flag!=1 || plankton->ts_flag==1){
                     error("a starting plankton array was supplied but the plankton initial_flag is not set (or ts_flag is)\n");
             }
             if(plankton->in->initial_len!=(R_xlen_t)x*(R_xlen_t)y*(R_xlen_t)np){
                     error("the starting plankton array holds %ld values; expected x*y*size classes = %d*%d*%d\n",(long)plankton->in->initial_len,x,y,np);
             }
             for(i=0 ; i<np ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     plankton->u_values[plankton->iplamin+i][k][l]=plankton->in->initial[(R_xlen_t)k+(R_xlen_t)x*((R_xlen_t)l+(R_xlen_t)y*(R_xlen_t)i)];
                             }
                     }
             }
     }
     else if(plankton->initial_flag==1 || plankton->ts_flag==1){
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
                             if(fgets(input,MAX_LEN,plankton->fptr_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
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
     pelagic->A_base=pelagic->A;         //temp driver: base (unscaled) search rate
     pelagic->mu_0_base=pelagic->mu_0;   //temp driver: base (unscaled) background mortality
     pelagic->tempeff=1.0;               //temp driver: neutral until a temperature series is read
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
     {char *e=getenv("PEL_IMM_FRAC"); pelagic->imm_frac = e ? atof(e) : 0.0; pelagic->rec0 = 0.0; pelagic->rec0_frozen = 0;}
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
     
     /*Text output files (not created when the output goes to R as arrays only)*/
     if(run_writes_files(run)){
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
     }
     
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
     pelagic->a_refuge=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);   /*reef refuge availability a(w)*/

     pelagic->phi_values=(double *)safe_malloc(m*sizeof(double),__LINE__,__FILE__);
     pelagic->refuge_values=(double *)safe_malloc(m*sizeof(double),__LINE__,__FILE__);   /*per prey-size refuge density*/
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
             pelagic->a_refuge[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);

             pelagic->phi_values[i]=0;
             pelagic->refuge_values[i]=0;
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
                     pelagic->a_refuge[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
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
                             pelagic->a_refuge[i][k][l]=1.0;   /*default fully exposed (a=1) -> standard dbpmr*/
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
     pelagic->rec0_c=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     pelagic->rec0_frozen_c=(int **)safe_malloc(x*sizeof(int *),__LINE__,__FILE__);
     for(k=0 ; k<x ; k++){
             pelagic->pla_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             pelagic->pel_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             pelagic->ben_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);

             pelagic->fish_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             pelagic->pred_total[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             pelagic->reproduction[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             pelagic->rec0_c[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             pelagic->rec0_frozen_c[k]=(int *)safe_malloc(y*sizeof(int),__LINE__,__FILE__);
             for(l=0 ; l<y ; l++){
                     pelagic->pla_total[k][l]=0;
                     pelagic->pel_total[k][l]=0;
                     pelagic->ben_total[k][l]=0;

                     pelagic->fish_total[k][l]=0;
                     pelagic->pred_total[k][l]=0;
                     pelagic->reproduction[k][l]=0;
                     pelagic->rec0_c[k][l]=0.0;
                     pelagic->rec0_frozen_c[k][l]=0;
             }
     }
          
     /*Setup initial values for u(m,0,x,y)*/
     if(pelagic->in!=NULL && pelagic->in->series!=NULL){
             /*Abundance forced through time from an R array [x][y][size class of the species][layer]: no _ts file is read*/
             series_check(pelagic->in, pelagic->filename, pelagic->ts_flag, x, y, pelagic->ipelmax-pelagic->ipelmin+1);
             series_to_u(pelagic->in, 0, x, y, pelagic->ipelmin, pelagic->ipelmax-pelagic->ipelmin+1, pelagic->u_values);
     }
     else if(pelagic->in!=NULL && pelagic->in->initial!=NULL){
             /*Starting abundance supplied from R as an array [x][y][size class of the species]: no _ts file is read*/
             int ns=pelagic->ipelmax-pelagic->ipelmin+1;
             if(pelagic->initial_flag!=1 || pelagic->ts_flag==1){
                     error("a starting array was supplied for %s but its initial_flag is not set (or ts_flag is)\n",pelagic->filename);
             }
             if(pelagic->in->initial_len!=(R_xlen_t)x*(R_xlen_t)y*(R_xlen_t)ns){
                     error("the starting array of %s holds %ld values; expected x*y*size classes = %d*%d*%d\n",pelagic->filename,(long)pelagic->in->initial_len,x,y,ns);
             }
             for(i=0 ; i<ns ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     pelagic->u_values[pelagic->ipelmin+i][k][l]=pelagic->in->initial[(R_xlen_t)k+(R_xlen_t)x*((R_xlen_t)l+(R_xlen_t)y*(R_xlen_t)i)];
                             }
                     }
             }
     }
     else if(pelagic->initial_flag==1){
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
                             if(fgets(input,MAX_LEN,pelagic->fptr_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
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
     
     /*Open pelagic fishing file (not needed when fishing mortality was supplied from R as an array)*/
     if(pelagic->fishing_flag==1 && pelagic->fish_arr!=NULL){
             if(pelagic->fish_arr->len!=(R_xlen_t)grid->xnum*(R_xlen_t)grid->ynum*(R_xlen_t)pelagic->fish_arr->nlayer){
                     error("fishing array of %s holds %ld values; expected x*y*layers = %d*%d*%d\n",pelagic->filename,(long)pelagic->fish_arr->len,grid->xnum,grid->ynum,pelagic->fish_arr->nlayer);
             }
     }
     else if(pelagic->fishing_flag==1){
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
             if(pelagic->rep_method==1 && pelagic->in!=NULL && pelagic->in->rep!=NULL){
                     /*Egg density supplied from R as an array [x][y][layer]: no _rep_ts file is opened*/
                     if(pelagic->in->rep_len!=(R_xlen_t)x*(R_xlen_t)y*(R_xlen_t)pelagic->in->rep_nlayer){
                             error("the reproduction array of %s holds %ld values; expected x*y*layers = %d*%d*%d\n",pelagic->filename,(long)pelagic->in->rep_len,x,y,pelagic->in->rep_nlayer);
                     }
             }
             else if(pelagic->rep_method==1){
                     /*Create pelagic fishing input filename*/
                     snprintf(pelagic->fname_rep_ts, sizeof(pelagic->fname_rep_ts),"%s/Input/%s_%s.txt",run->filename,pelagic->filename,"rep_ts");
                     pelagic->fptr_rep_ts=safe_fopen(pelagic->fname_rep_ts,"r",__LINE__,__FILE__);
             }
     }

     
     /*Setup values for predator-prey fedding kernel*/
     for(i=0 ; i<(pelagic->ipelmax+1) ; i++){
             pelagic->phi_values[i]=phi((grid->m_values[i]-grid->m_values[0]),pelagic->q_0, pelagic->sig, pelagic->trunc);
     }

     /* --- Reef prey-refuge plugin (opt-in; Rogers, Blanchard & Mumby 2014) ---
        Availability a(w) = max(min_A, 1 - refuge(w)/competitor(w)) scales BOTH predation mortality on,
        and predator intake of, prey of body size w (mass balance). Reads an optional per-size refuge
        density from <run>/Input/<name>_refuge.txt (one value per grid size class). Absent file (or all
        zeros) -> refuge_on=0 -> a(w)=1 everywhere -> bit-identical to standard dbpmr. min_A (floor) from
        the REFUGE_MIN_A env var (default 0). */
     pelagic->refuge_on = 0;
     pelagic->refuge_min_a = 0.0;
     { char *e = getenv("REFUGE_MIN_A"); if(e) pelagic->refuge_min_a = atof(e); }
     if(pelagic->in!=NULL && pelagic->in->refuge!=NULL){
             /*Refuge density supplied from R as an array [x][y][size class of the grid], so it can differ between
               cells: no _refuge file is read and the floor comes with the array, not from the environment*/
             R_xlen_t j, len=(R_xlen_t)x*(R_xlen_t)y*(R_xlen_t)grid->mnum;
             if(pelagic->in->refuge_len!=len){
                     error("the refuge array of %s holds %ld values; expected x*y*size classes = %d*%d*%d\n",pelagic->filename,(long)pelagic->in->refuge_len,x,y,grid->mnum);
             }
             pelagic->refuge_min_a=pelagic->in->refuge_min_a;
             for(j=0 ; j<len ; j++){
                     if(pelagic->in->refuge[j]>0){ pelagic->refuge_on = 1; break; }   /*only activate when a positive refuge is supplied*/
             }
     }
     else{
             char rfn[256];
             snprintf(rfn, sizeof(rfn), "%s/Input/%s_refuge.txt", run->filename, pelagic->filename);
             FILE *rf = fopen(rfn, "r");
             if(rf){
                     double v; int cnt=0;
                     for(i=0 ; i<grid->mnum && fscanf(rf,"%lf",&v)==1 ; i++){
                             pelagic->refuge_values[i]=v; if(v>0) cnt++;
                     }
                     fclose(rf);
                     if(cnt>0) pelagic->refuge_on = 1;   /*only activate when a positive refuge is supplied*/
             }
     }

     /* --- Per-cell model parameters (opt-in) ---
        Optional <run>/Input/<name>_cellpar.txt: a header line naming the columns, then one row per
        spatial cell (x*y rows, x-outer y-inner - the plankton _ts / per-cell forcing convention).
        Recognised column: pref_ben (benthic feeding preference, e.g. depth-dependent). The values are
        static (read once). Absent file or column -> pref_ben_c stays NULL -> the scalar pref_ben is
        used in every cell -> bit-identical to standard dbpmr. */
     pelagic->pref_ben_c = NULL;
     if(pelagic->in!=NULL && pelagic->in->pref_ben!=NULL){
             /*Per-cell benthic feeding preference supplied from R as a matrix [x][y]: no _cellpar file is read*/
             if(pelagic->in->pref_ben_len!=(R_xlen_t)x*(R_xlen_t)y){
                     error("the per-cell pref_ben of %s holds %ld values; expected x*y = %d*%d\n",pelagic->filename,(long)pelagic->in->pref_ben_len,x,y);
             }
             pelagic->pref_ben_c=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             for(k=0 ; k<x ; k++){
                     pelagic->pref_ben_c[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     for(l=0 ; l<y ; l++){
                             pelagic->pref_ben_c[k][l]=pelagic->in->pref_ben[(R_xlen_t)k+(R_xlen_t)x*(R_xlen_t)l];
                     }
             }
     }
     else{
             char cfn[256];
             snprintf(cfn, sizeof(cfn), "%s/Input/%s_cellpar.txt", run->filename, pelagic->filename);
             FILE *cf = fopen(cfn, "r");
             if(cf){
                     char cline[1024];
                     int col_pb=-1;
                     if(fgets(cline,sizeof(cline),cf)!=NULL){
                             int col=0;
                             char *name=strtok(cline,", \t\r\n");
                             while(name!=NULL){
                                     if(strcmp(name,"pref_ben")==0) col_pb=col;
                                     col++;
                                     name=strtok(NULL,", \t\r\n");
                             }
                     }
                     if(col_pb>=0){
                             int xx=grid->xnum, yy=grid->ynum, kk, ll;
                             pelagic->pref_ben_c=(double **)safe_malloc(xx*sizeof(double *),__LINE__,__FILE__);
                             for(kk=0 ; kk<xx ; kk++){
                                     pelagic->pref_ben_c[kk]=(double *)safe_malloc(yy*sizeof(double),__LINE__,__FILE__);
                             }
                             for(kk=0 ; kk<xx ; kk++){
                                     for(ll=0 ; ll<yy ; ll++){
                                             double cvals[64];
                                             int nv=0;
                                             if(fgets(cline,sizeof(cline),cf)==NULL){
                                                     fclose(cf);
                                                     error("%s: expected %d rows (one per cell, x-outer y-inner), found %d\n",cfn,xx*yy,(kk*yy)+ll);
                                             }
                                             char *tok=strtok(cline,", \t\r\n");
                                             while(tok!=NULL && nv<64){ cvals[nv++]=atof(tok); tok=strtok(NULL,", \t\r\n"); }
                                             if(col_pb>=nv){
                                                     fclose(cf);
                                                     error("%s: row %d has no value for column pref_ben\n",cfn,(kk*yy)+ll+1);
                                             }
                                             pelagic->pref_ben_c[kk][ll]=cvals[col_pb];
                                     }
                             }
                     }
                     fclose(cf);
             }
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
     benthic->A_base=benthic->A;         //temp driver: base (unscaled) search rate
     benthic->mu_0_base=benthic->mu_0;   //temp driver: base (unscaled) background mortality
     benthic->tempeff=1.0;               //temp driver: neutral until a temperature series is read
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
     
     /*Text output files (not created when the output goes to R as arrays only)*/
     if(run_writes_files(run)){
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
     }
     
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
     if(benthic->in!=NULL && benthic->in->series!=NULL){
             /*Abundance forced through time from an R array [x][y][size class of the species][layer]: no _ts file is read*/
             series_check(benthic->in, benthic->filename, benthic->ts_flag, x, y, benthic->ibenmax-benthic->ibenmin+1);
             series_to_u(benthic->in, 0, x, y, benthic->ibenmin, benthic->ibenmax-benthic->ibenmin+1, benthic->u_values);
     }
     else if(benthic->in!=NULL && benthic->in->initial!=NULL){
             /*Starting abundance supplied from R as an array [x][y][size class of the species]: no _ts file is read*/
             int ns=benthic->ibenmax-benthic->ibenmin+1;
             if(benthic->initial_flag!=1 || benthic->ts_flag==1){
                     error("a starting array was supplied for %s but its initial_flag is not set (or ts_flag is)\n",benthic->filename);
             }
             if(benthic->in->initial_len!=(R_xlen_t)x*(R_xlen_t)y*(R_xlen_t)ns){
                     error("the starting array of %s holds %ld values; expected x*y*size classes = %d*%d*%d\n",benthic->filename,(long)benthic->in->initial_len,x,y,ns);
             }
             for(i=0 ; i<ns ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     benthic->u_values[benthic->ibenmin+i][k][l]=benthic->in->initial[(R_xlen_t)k+(R_xlen_t)x*((R_xlen_t)l+(R_xlen_t)y*(R_xlen_t)i)];
                             }
                     }
             }
     }
     else if(benthic->initial_flag==1){
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
                             if(fgets(input,MAX_LEN,benthic->fptr_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
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
     
     /*Setup initial values for fishing mortality (no file when it was supplied from R as an array)*/
     if(benthic->fishing_flag==1 && benthic->fish_arr!=NULL){
             if(benthic->fish_arr->len!=(R_xlen_t)grid->xnum*(R_xlen_t)grid->ynum*(R_xlen_t)benthic->fish_arr->nlayer){
                     error("fishing array of %s holds %ld values; expected x*y*layers = %d*%d*%d\n",benthic->filename,(long)benthic->fish_arr->len,grid->xnum,grid->ynum,benthic->fish_arr->nlayer);
             }
     }
     else if(benthic->fishing_flag==1){
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
             if(benthic->rep_method==1 && benthic->in!=NULL && benthic->in->rep!=NULL){
                     /*Egg density supplied from R as an array [x][y][layer]: no _rep_ts file is opened*/
                     if(benthic->in->rep_len!=(R_xlen_t)x*(R_xlen_t)y*(R_xlen_t)benthic->in->rep_nlayer){
                             error("the reproduction array of %s holds %ld values; expected x*y*layers = %d*%d*%d\n",benthic->filename,(long)benthic->in->rep_len,x,y,benthic->in->rep_nlayer);
                     }
             }
             else if(benthic->rep_method==1){
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
     
     /*Text output files (not created when the output goes to R as arrays only)*/
     if(run_writes_files(run)){
             /*Open output files*/
             detritus->fptr_summ=safe_fopen(detritus->fname_summ,"w",__LINE__,__FILE__);
     
             fprintf(detritus->fptr_summ,"filetype:      %s\n","summary");
             fprintf(detritus->fptr_summ,"speciestype:   %s\n","detritus");
             fprintf(detritus->fptr_summ,"\n");
             print_detritus_header(detritus->fptr_summ);
     
             /*Close output files*/
             fclose(detritus->fptr_summ);
     }
     
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
     if(detritus->in!=NULL && detritus->in->series!=NULL){
             /*Detritus biomass forced through time from an R array [x][y][layer]: no _ts file is read*/
             series_check(detritus->in, detritus->filename, detritus->ts_flag, x, y, 1);
             series_to_w(detritus->in, 0, x, y, detritus->w_values);
     }
     else if(detritus->in!=NULL && detritus->in->initial!=NULL){
             /*Starting detritus biomass supplied from R as a matrix [x][y]: no _ts file is read*/
             if(detritus->initial_flag!=1 || detritus->ts_flag==1){
                     error("a starting detritus matrix was supplied but the detritus initial_flag is not set (or ts_flag is)\n");
             }
             if(detritus->in->initial_len!=(R_xlen_t)x*(R_xlen_t)y){
                     error("the starting detritus matrix holds %ld values; expected x*y = %d*%d\n",(long)detritus->in->initial_len,x,y);
             }
             for(k=0 ; k<x ; k++){
                     for(l=0 ; l<y ; l++){
                             detritus->w_values[k][l]=detritus->in->initial[(R_xlen_t)k+(R_xlen_t)x*(R_xlen_t)l];
                     }
             }
     }
     else if(detritus->initial_flag==1){
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
                             if(fgets(input,MAX_LEN,detritus->fptr_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
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
