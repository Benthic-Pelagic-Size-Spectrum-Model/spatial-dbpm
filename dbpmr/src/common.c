#include "common.h"


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
     
     /*Free the gravity multiplier on fishing mortality*/
     if(run->fish_share!=NULL){
             for(k=0 ; k<x ; k++) free(run->fish_share[k]);
             free(run->fish_share);
     }
     if(run->share_u_pel!=NULL){
             for(s=0 ; s<n ; s++) free(run->share_u_pel[s]);
             free(run->share_u_pel);
     }
     if(run->share_u_ben!=NULL){
             for(b=0 ; b<c ; b++) free(run->share_u_ben[b]);
             free(run->share_u_ben);
     }

     /*Free grid values*/
     free(grid->m_values);
     free(grid->exp_m_values);
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
                                     free(community->pelagic[s].a_refuge[i][k]);
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
                             free(community->pelagic[s].a_refuge[i]);
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
                     free(community->pelagic[s].a_refuge);
                     free(community->pelagic[s].phi_values);
                     free(community->pelagic[s].refuge_values);
                     if(community->pelagic[s].pref_ben_c){
                             for(k=0 ; k<x ; k++){
                                     free(community->pelagic[s].pref_ben_c[k]);
                             }
                             free(community->pelagic[s].pref_ben_c);
                     }
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
                             free(community->pelagic[s].rec0_c[k]);
                             free(community->pelagic[s].rec0_frozen_c[k]);
                     }
                     free(community->pelagic[s].pla_total);
                     free(community->pelagic[s].pel_total);
                     free(community->pelagic[s].ben_total);
                     free(community->pelagic[s].fish_total);
                     free(community->pelagic[s].pred_total);
                     free(community->pelagic[s].reproduction);
                     free(community->pelagic[s].rec0_c);
                     free(community->pelagic[s].rec0_frozen_c);
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
