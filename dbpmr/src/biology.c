#include "biology.h"


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

     /* Reef refuge pre-pass: recompute the per-size availability a(w) for every pelagic species from the
        CURRENT spectrum, BEFORE the feeding/mortality integrals that consume it (so predation mortality
        and predator intake see the SAME a(w) -> mass balance). a(w)=max(min_A, 1 - refuge(w)/competitor(w)),
        competitor = U abundance density. Species with refuge off keep a=1 (set at allocation) -> no effect. */
     for(s=0 ; s<n ; s++){
             if(!community->pelagic[s].refuge_on) continue;
             double minA=community->pelagic[s].refuge_min_a;
             const double *rfc=(community->pelagic[s].in!=NULL) ? community->pelagic[s].in->refuge : NULL;   /*refuge per cell, [x][y][m]*/
             for(i=0 ; i<grid->mnum ; i++){
                     double rf=community->pelagic[s].refuge_values[i];
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     if(rfc!=NULL) rf=rfc[(R_xlen_t)k+(R_xlen_t)x*((R_xlen_t)l+(R_xlen_t)y*(R_xlen_t)i)];
                                     double comp=community->pelagic[s].u_values[i][k][l];
                                     double a=(rf>0 && comp>0) ? 1.0-rf/comp : 1.0;
                                     if(a<minA) a=minA; if(a>1.0) a=1.0; if(a<0.0) a=0.0;
                                     community->pelagic[s].a_refuge[i][k][l]=a;
                             }
                     }
             }
     }

     //Pelagic changes
     for(s=0 ; s<n ; s++){
             for(i=community->pelagic[s].ipelmin ; i<(community->pelagic[s].ipelmax+1) ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     /* per-cell pelagic temperature multiplier on FEEDING + PREDATION (1 when off). */
                                     double tmP = run->forcing_percell ? run->pel_tempeff_c[k][l] : 1.0;
                                     for(ss=0 ; ss<n ; ss++){
                                              community->pelagic[s].mu_pred_values[ss][i][k][l]=tmP*mu_pel_pred(s, ss, i, k, l, run, grid, community);
                                     }
                                     community->pelagic[s].mu_values[i][k][l]=mu_pel(s, i, k, l, run, grid, community);

                                     community->pelagic[s].pla_bio[i][k][l]=tmP*pla_biomass(s, i, k, l, grid, community);
                                     community->pelagic[s].pel_bio[i][k][l]=tmP*pel_biomass(s, i, k, l, run, grid, community);
                                     community->pelagic[s].ben_bio[i][k][l]=tmP*ben_biomass(s, i, k, l, run, grid, community);
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
                                     fish+=community->pelagic[s].mu_fish_values[i][k][l]*community->pelagic[s].u_values[i][k][l]*grid->exp_m_values[i]*grid->mstep;
                                     for(ss=0 ; ss<n ; ss++){
                                              pred+=community->pelagic[s].mu_pred_values[ss][i][k][l]*community->pelagic[s].u_values[i][k][l]*grid->exp_m_values[i]*grid->mstep;
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
                                             /* predation ON benthos is BY pelagic predators (pel temp); benthic feeding uses benthic temp */
                                             double tmPb = run->forcing_percell ? run->pel_tempeff_c[k][l] : 1.0;
                                             double tmB  = run->forcing_percell ? run->ben_tempeff_c[k][l] : 1.0;
                                             community->benthic[b].mu_pred_values[i][k][l]=tmPb*mu_ben_pred(b, i, k, l, run, grid, community);
                                             community->benthic[b].mu_values[i][k][l]=mu_ben(b, i, k, l, run, grid, community);

                                             community->benthic[b].det_bio[i][k][l]=tmB*det_biomass(b, i, k, l, run, grid, community);
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
                                             fish+=community->benthic[b].mu_fish_values[i][k][l]*community->benthic[b].u_values[i][k][l]*grid->exp_m_values[i]*grid->mstep;
                                             pred+=community->benthic[b].mu_pred_values[i][k][l]*community->benthic[b].u_values[i][k][l]*grid->exp_m_values[i]*grid->mstep;
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
               ans+=grid->exp_m_values[i]*community->pelagic[species].phi_values[size-i]*community->plankton->u_values[i][xspace][yspace]*grid->mstep;
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
                               /*reef refuge: prey of size i are available for intake only with probability a(i) (=1 when off)*/
                               ans+=community->pelagic[s].a_refuge[i][xspace][yspace]*grid->exp_m_values[i]*community->pelagic[species].phi_values[size-i]*community->pelagic[s].u_values[i][xspace][yspace]*grid->mstep;
                       }
               }
       }

       /* Allows feeding from own Size Spectra Only*/
       else{
              for(i=community->pelagic[species].ipelmin ; i<size ; i++){
                      ans+=community->pelagic[species].a_refuge[i][xspace][yspace]*grid->exp_m_values[i]*community->pelagic[species].phi_values[size-i]*community->pelagic[species].u_values[i][xspace][yspace]*grid->mstep;
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
                                       ans+=grid->exp_m_values[i]*community->pelagic[species].phi_values[size-i]*community->benthic[b].u_values[i][xspace][yspace]*grid->mstep;
                               }
                       }
               }
               ans*=pel_pref_ben(&community->pelagic[species],xspace,yspace)*community->pelagic[species].A*exp(community->pelagic[species].alpha*grid->m_values[size]);
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
     int i,s,b;
     int n,c;

     n=run->no_pelagic;
     c=run->no_benthic;

     double ans_surf=0;   /*surface-origin inputs (sink to the seafloor -> scaled by sinking_rate)*/
     double ans_bed=0;    /*seafloor-origin inputs (dead benthos, already on the bed -> no sinking)*/

     /* Size-dependent export closure (opt-in): when the export_attn channel is supplied, each
        surface particle of ln-mass m reaches the seafloor with fraction
          export(m) = exp( -export_attn * exp(-export_gamma*m) ),
        from an allometric sinking speed w_s ~ mass^export_gamma and an attenuation number
        export_attn = k(T)*z/w0 - so large carcasses/pellets sink faster and are recycled less
        than small ones. Replaces the flat sinking_rate. Default (channel absent): flat export. */
     int size_dep = (run->col_export_attn>=0);
     /*Per-cell export parameters when FORCING_PERCELL (spatial_dim>=1); otherwise the global values.*/
     double e_attn  = run->forcing_percell ? run->export_attn_c[xspace][yspace]  : run->export_attn;
     double e_gamma = run->forcing_percell ? run->export_gamma_c[xspace][yspace] : run->export_gamma;
     double e_magg  = run->forcing_percell ? run->export_magg_c[xspace][yspace]  : run->export_magg;

     if(run->coupled_flag==1){
             for(s=0 ; s<n ; s++){
                     /*Rate of det biomass in from pelagic defecation (faeces follow the
                       temperature-scaled feeding, as sizemodel's defbypred does). Flat export uses
                       the pre-summed *_total; size-dependent export reconstructs faeces per predator
                       size (pla_bio/pel_bio/ben_bio) so export(m) can weight each size class.*/
                     if(!size_dep){
                             ans_surf+=(community->pelagic[s].Ex_pla*community->pelagic[s].pla_total[xspace][yspace]+community->pelagic[s].Ex_pel*community->pelagic[s].pel_total[xspace][yspace]+community->pelagic[s].Ex_ben*community->pelagic[s].ben_total[xspace][yspace]);
                     }
                     for(i=community->pelagic[s].ipelmin ; i<(community->pelagic[s].ipelmax+1) ; i++){
                             /*Dead pelagic stuff. Use the BASE (un-temperature-scaled) background
                               mortality: sizemodel removes temperature from the mortality feeding
                               the detritus pool.*/
                             double dead = community->pelagic[s].mu_0_base*exp(community->pelagic[s].beta*grid->m_values[i])*community->pelagic[s].u_values[i][xspace][yspace]*grid->exp_m_values[i]*grid->mstep
                                         + community->pelagic[s].mu_s*((log10(grid->exp_m_values[i])-log10(exp(community->pelagic[s].pelmin)))/((log10(exp(community->pelagic[s].pelmax))+community->pelagic[s].epsilon)-log10(grid->exp_m_values[i])))*community->pelagic[s].u_values[i][xspace][yspace]*grid->exp_m_values[i]*grid->mstep;
                             if(size_dep){
                                     double faeces = (community->pelagic[s].Ex_pla*community->pelagic[s].pla_bio[i][xspace][yspace]+community->pelagic[s].Ex_pel*community->pelagic[s].pel_bio[i][xspace][yspace]+community->pelagic[s].Ex_ben*community->pelagic[s].ben_bio[i][xspace][yspace])*community->pelagic[s].u_values[i][xspace][yspace]*grid->mstep;
                                     double me = fmax(grid->m_values[i], e_magg);   /*aggregation floor: small stuff sinks as marine snow*/
                                     ans_surf += exp(-e_attn*exp(-e_gamma*me))*(faeces + dead);
                             }else{
                                     ans_surf += dead;
                             }
                     }
                     /*Dead benthos (detritivore background + senescence mortality). These bodies
                       are already on/in the seafloor, so they DO NOT get the sinking_rate factor -
                       matching sizemodel's benthic terms outside the sinking_rate multiplier. Uses
                       base (un-temperature-scaled) background mortality, as for the pelagic.*/
                     for(b=0 ; b<c ; b++){
                             for(i=community->benthic[b].ibenmin ; i<(community->benthic[b].ibenmax+1) ; i++){
                                     ans_bed+=community->benthic[b].mu_0_base*exp(community->benthic[b].beta*grid->m_values[i])*community->benthic[b].u_values[i][xspace][yspace]*grid->exp_m_values[i]*grid->mstep;
                                     ans_bed+=community->benthic[b].mu_s*((log10(grid->exp_m_values[i])-log10(exp(community->benthic[b].benmin)))/((log10(exp(community->benthic[b].benmax))+community->benthic[b].epsilon)-log10(grid->exp_m_values[i])))*community->benthic[b].u_values[i][xspace][yspace]*grid->exp_m_values[i]*grid->mstep;
                             }
                     }
             }
     }
     for(i=community->plankton->iplamin ; i<(community->plankton->iplamax+1) ; i++){
             /*Rate of det biomass in from dead (sinking) plankton - surface origin*/
             double pdead = community->plankton->mu_0*exp(community->plankton->beta*grid->m_values[i])*community->plankton->u_values[i][xspace][yspace]*grid->exp_m_values[i]*grid->mstep;
             if(size_dep) ans_surf += exp(-e_attn*exp(-e_gamma*fmax(grid->m_values[i], e_magg)))*pdead;
             else         ans_surf += pdead;
     }

     /*Surface-origin inputs sink: scale by the time-varying sinking_rate (= export ratio),
       exactly as sizemodel multiplies its pelagic detritus inputs by sinking_rate[i]. Seafloor
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
     return((size_dep ? ans_surf : sink * ans_surf) + ans_bed);
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
       
       /* Mortality due to Predation by all Size Spectra. Reef refuge: scale each predator i's contribution
          by its availability a(i) -- matching Rogers' death.u, where a.u multiplies the predator (U) vector
          INSIDE the predator integral (mortality kernel mphi = q2 = -q1 is the transpose of the feeding
          kernel gphi = q1, so the same a*(U) construction lands on the PREDATOR index in mortality but the
          PREY index in feeding/pel_biomass). a(i)=1 when refuge off -> identical to standard dbpmr. */
       if(run->coupled_flag==1){
               for(i=size ; i<(community->pelagic[pred].ipelmax+1) ; i++){
                       ans+=community->pelagic[pred].a_refuge[i][xspace][yspace]*(exp(community->pelagic[pred].alpha*grid->m_values[i]))*community->pelagic[pred].phi_values[i-size]*(community->pelagic[pred].u_values[i][xspace][yspace])*community->pelagic[pred].A*community->pelagic[pred].pref_pel*grid->mstep;
               }
       }
       /*Mortality due to Own Size Spectra Predation*/
       else{
            if(species==pred){
                    for(i=size ; i<(community->pelagic[species].ipelmax+1) ; i++){
                            ans+=community->pelagic[species].a_refuge[i][xspace][yspace]*(exp(community->pelagic[species].alpha*grid->m_values[i]))*community->pelagic[species].phi_values[i-size]*(community->pelagic[species].u_values[i][xspace][yspace])*community->pelagic[species].A*community->pelagic[species].pref_pel*grid->mstep;
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
       ans+=((run->forcing_percell && run->temp_on_mort ? run->pel_tempeff_c[xspace][yspace] : 1.0)*community->pelagic[species].mu_0*exp(community->pelagic[species].beta*grid->m_values[size]));
       ans+=(community->pelagic[species].mu_s*((log10(grid->exp_m_values[size])-log10(exp(community->pelagic[species].pelmin)))/((log10(exp(community->pelagic[species].pelmax))+community->pelagic[species].epsilon)-log10(grid->exp_m_values[size]))));
            
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
                                  ans=ans+(exp(community->pelagic[s].alpha*grid->m_values[i]))*community->pelagic[s].phi_values[i-size]*(community->pelagic[s].u_values[i][xspace][yspace])*community->pelagic[s].A*pel_pref_ben(&community->pelagic[s],xspace,yspace)*grid->mstep;
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
       ans+=((run->forcing_percell && run->temp_on_mort ? run->ben_tempeff_c[xspace][yspace] : 1.0)*community->benthic[species].mu_0*exp(community->benthic[species].beta*grid->m_values[size]));
       ans+=(community->benthic[species].mu_s*((log10(grid->exp_m_values[size])-log10(exp(community->benthic[species].benmin)))/((log10(exp(community->benthic[species].benmax))+community->benthic[species].epsilon)-log10(grid->exp_m_values[size]))));
       
       return ans;
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
             //egg values supplied from R as an array [x][y][layer]: the layer of the current time point
             if(community->pelagic[s].rep_method==1 && community->pelagic[s].in!=NULL && community->pelagic[s].in->rep!=NULL){
                     const SPECIES_INPUT *in=community->pelagic[s].in;
                     R_xlen_t layer=in->rep_layer[run->rep_block < in->rep_nblock ? run->rep_block : in->rep_nblock-1];
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     community->pelagic[s].reproduction[k][l]=in->rep[(R_xlen_t)k+(R_xlen_t)x*((R_xlen_t)l+(R_xlen_t)y*layer)];
                             }
                     }
             }
             //read in egg values from file
             else if(community->pelagic[s].rep_method==1){
                     int MAX_LEN=15; //this is the upper bound for the row length in characters
                     char *input;               //input string to store line as string
                     input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                     char sep[]=",";            //identifying the ',' char as seperating values
                     char *result = NULL;       //declaring a pointer for storing each individual value as a char
                     
                     /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     if(fgets(input,MAX_LEN,community->pelagic[s].fptr_rep_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
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
                                     /*Connectivity/immigration floor (0-D): recruitment cannot fall below imm_frac of the
                                       unfished-peak recruitment (external larval subsidy for open oligotrophic systems). rec0
                                       tracks the running peak of NATURAL recruitment, frozen the instant the floor first binds,
                                       so the subsidy can never inflate its own reference (no positive feedback).*/
                                     if(community->pelagic[s].imm_frac>0){
                                             /*Per-cell reference (rec0_c[k][l]): each cell keeps its OWN running-peak
                                               recruitment and freeze flag, so the floor never couples cells in a
                                               spatial_dim>=1 run. For spatial_dim=0 (x=y=1) this is the [0][0] scalar,
                                               bit-identical to the former per-species scalar.*/
                                             if(!community->pelagic[s].rec0_frozen_c[k][l] &&
                                                community->pelagic[s].reproduction[k][l]>community->pelagic[s].rec0_c[k][l])
                                                     community->pelagic[s].rec0_c[k][l]=community->pelagic[s].reproduction[k][l];
                                             double fl=community->pelagic[s].imm_frac*community->pelagic[s].rec0_c[k][l];
                                             if(community->pelagic[s].reproduction[k][l]<fl){
                                                     community->pelagic[s].reproduction[k][l]=fl;
                                                     community->pelagic[s].rec0_frozen_c[k][l]=1;
                                             }
                                     }
                             }
                     }
             }
             //fecundity calculation
             if(community->pelagic[s].rep_method==3){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     pel=0;
                                     for(i=community->pelagic[s].ipelmat ; i<(community->pelagic[s].ipelmax+1) ; i++){
                                             pel+=community->pelagic[s].u_values[i][k][l]*grid->exp_m_values[i]*grid->mstep;
                                     }
                                     community->pelagic[s].reproduction[k][l]=pel*492*grid->tstep;
                             }
                     }
             }
     }
     
     /*Benthic Reproduction*/
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     //egg values supplied from R as an array [x][y][layer]: the layer of the current time point
                     if(community->benthic[b].rep_method==1 && community->benthic[b].in!=NULL && community->benthic[b].in->rep!=NULL){
                             const SPECIES_INPUT *in=community->benthic[b].in;
                             R_xlen_t layer=in->rep_layer[run->rep_block < in->rep_nblock ? run->rep_block : in->rep_nblock-1];
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             community->benthic[b].reproduction[k][l]=in->rep[(R_xlen_t)k+(R_xlen_t)x*((R_xlen_t)l+(R_xlen_t)y*layer)];
                                     }
                             }
                     }
                     //read in egg values from file
                     else if(community->benthic[b].rep_method==1){
                             int MAX_LEN=15; //this is the upper bound for the row length in characters
                             char *input;               //input string to store line as string
                             input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                             char sep[]=",";            //identifying the ',' char as seperating values
                             char *result = NULL;       //declaring a pointer for storing each individual value as a char
                             
                             /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             if(fgets(input,MAX_LEN,community->benthic[b].fptr_rep_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
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
                                                     ben+=community->benthic[b].u_values[i][k][l]*grid->exp_m_values[i]*grid->mstep;
                                             }
                                             community->benthic[b].reproduction[k][l]=ben*492*grid->tstep;
                                     }
                             }
                     }
             }
     }

     run->rep_block++;   /*one call per model time point (the first is the initial state)*/
}
