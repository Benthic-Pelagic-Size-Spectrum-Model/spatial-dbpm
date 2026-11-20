#include "fishing.h"


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
                                     if(fgets(input,MAX_LEN,community->pelagic[s].fptr_fish_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
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
                                             if(fgets(input,MAX_LEN,community->benthic[b].fptr_fish_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
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

     /* DBPM.md gravity fishing (opt-in via the `gravity` channel): replace the static per-group
        catchability split with the true B_{t-1} feedback. For each cell, scale each group's fishing
        mortality by that group's CURRENT areal fishable-biomass share,
          F_g = F_base * A_g/(A_pel+A_ben),  A_g = sum_{fished sizes} u*exp(m)*mstep * thickness
        (pelagic x run->depth = min(depth,200) m from the depth channel; benthic x 20 m). The base F
        read above is the un-split Q*selectivity*effort, so this yields catch_g ~ B_g^2 (effort
        concentrates where biomass is). Runs per cell -> ready for gridded. */
     if(run->col_gravity>=0){
             double bhd=20.0;                 /* benthic habitat depth (sizemodel default) */
             for(k=0 ; k<x ; k++){
                     for(l=0 ; l<y ; l++){
                             /* pelagic column thickness (depth channel; 1 if absent) - per-cell under FORCING_PERCELL */
                             double dcorr=run->forcing_percell ? run->depth_c[k][l] : run->depth;
                             double Bpel=0.0, Bben=0.0;
                             for(s=0 ; s<n ; s++){
                                     if(community->pelagic[s].fishing_flag!=1) continue;
                                     for(i=community->pelagic[s].ipelmin ; i<(community->pelagic[s].ipelmax+1) ; i++)
                                             if(community->pelagic[s].mu_fish_values[i][k][l]>0.0)
                                                     Bpel += community->pelagic[s].u_values[i][k][l]*grid->exp_m_values[i]*grid->mstep;
                             }
                             for(b=0 ; b<c ; b++){
                                     if(community->benthic[b].fishing_flag!=1) continue;
                                     for(i=community->benthic[b].ibenmin ; i<(community->benthic[b].ibenmax+1) ; i++)
                                             if(community->benthic[b].mu_fish_values[i][k][l]>0.0)
                                                     Bben += community->benthic[b].u_values[i][k][l]*grid->exp_m_values[i]*grid->mstep;
                             }
                             double Ap=Bpel*dcorr, Ab=Bben*bhd, tot=Ap+Ab;
                             double sp = (tot>0.0) ? Ap/tot : 0.5;
                             double sb = (tot>0.0) ? Ab/tot : 0.5;
                             for(s=0 ; s<n ; s++)
                                     if(community->pelagic[s].fishing_flag==1)
                                             for(i=community->pelagic[s].ipelmin ; i<(community->pelagic[s].ipelmax+1) ; i++)
                                                     community->pelagic[s].mu_fish_values[i][k][l] *= sp;
                             for(b=0 ; b<c ; b++)
                                     if(community->benthic[b].fishing_flag==1)
                                             for(i=community->benthic[b].ibenmin ; i<(community->benthic[b].ibenmax+1) ; i++)
                                                     community->benthic[b].mu_fish_values[i][k][l] *= sb;
                     }
             }
     }
}
