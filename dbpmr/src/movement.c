#include "movement.h"


void xmove_solver(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *xmat, double ****temp_pelagic, int time)
/* Spatial movement in x, operator-split (fix for issue #41 — the old scheme accumulated biomass at the
   x-max wall regardless of source and did not conserve mass):
     (1) ADVECTION in conservative upwind FLUX form, with ZERO FLUX through the domain walls
         (mass-conserving + stable); advective velocity v = C*dg/dx - D*dmu/dx (move up the food
         gradient g, down the mortality gradient mu);
     (2) DIFFUSION by Crank-Nicolson with ZERO-FLUX (Neumann) boundaries.
   With no movement (C=D=0 and Diffun=0) both steps are the identity -> output unchanged. */
{
     int i,k,l,s;
     int x,y,n;
     double dt,dx,C,D,r;

     dt=grid->tstep;
     dx=grid->xstep;
     x=grid->xnum;
     y=grid->ynum;
     n=run->no_pelagic;

     double *uadv=(double *)safe_malloc(x*sizeof(double),__LINE__,__FILE__);   /* state after advection */
     double *flux=(double *)safe_malloc(x*sizeof(double),__LINE__,__FILE__);   /* flux[k] = flux across the face between cells k and k+1 */
     double *vface=(double *)safe_malloc(x*sizeof(double),__LINE__,__FILE__);  /* advective velocity at face k (depends on g,mu only) */

     for(s=0 ; s<n ; s++){
             for(i=community->pelagic[s].ipelmin ; i<(community->pelagic[s].ipelmax+1) ; i++){
                     C=Cfun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]));
                     D=Dfun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]));
                     r=Diffun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]))*(dt/(dx*dx));
                     int noadv=(C==0.0 && D==0.0);   /* no advection -> identity; also avoids 0*NaN if a crashed cell left NaN in g/mu */
                     for(l=0 ; l<y ; l++){
                             /* --- 1. advection: conservative upwind flux form; zero flux at both walls.
                                CFL SUB-STEPPING: the explicit upwind update is stable only for Courant number
                                |v|*dt/dx <= 1; a large prey/pred coefficient would otherwise flush a cell in one
                                step and (after the non-negativity clamp) silently zero the field. Split dt into
                                nsub steps so the per-sub-step Courant number stays <=1 — strong movement then
                                saturates smoothly instead of exploding. Velocities depend on g,mu only (fixed
                                across the sub-steps); only the upwind donor value is refreshed each sub-step. */
                             double vmax=0.0;
                             for(k=0 ; k<x-1 ; k++){
                                     vface[k]= noadv ? 0.0 :
                                               (C*(community->pelagic[s].g_values[i][k+1][l]-community->pelagic[s].g_values[i][k][l])
                                               -D*(community->pelagic[s].mu_values[i][k+1][l]-community->pelagic[s].mu_values[i][k][l]))/dx;
                                     if(fabs(vface[k])>vmax) vmax=fabs(vface[k]);
                             }
                             for(k=0 ; k<x ; k++) uadv[k]=temp_pelagic[s][i][k][l];
                             int nsub=(int)ceil(vmax*dt/dx); if(nsub<1) nsub=1;
                             double dts=dt/(double)nsub;
                             for(int sub=0 ; sub<nsub ; sub++){
                                     for(k=0 ; k<x-1 ; k++){
                                             double uup=(vface[k]>0.0)? uadv[k] : uadv[k+1];   /* upwind donor */
                                             flux[k]=vface[k]*uup;
                                     }
                                     for(k=0 ; k<x ; k++){
                                             double Fr=(k<x-1)? flux[k]   : 0.0;   /* right-face flux (0 at the x-max wall) */
                                             double Fl=(k>0)  ? flux[k-1] : 0.0;   /* left-face  flux (0 at the x=0   wall) */
                                             uadv[k]=uadv[k]-(dts/dx)*(Fr-Fl);     /* flux[] precomputed -> in-place is safe */
                                             if(uadv[k]<0.0) uadv[k]=0.0;          /* clamp tiny upwind undershoots */
                                     }
                             }
                             /* --- 2. diffusion: Crank-Nicolson with zero-flux (Neumann) BCs; conserves mass --- */
                             /* RHS (I + r/2 L): xmat->r = RHS * uadv via trimul */
                             xmat->b[0]=1.0-r/2.0;  xmat->c[0]=r/2.0;
                             for(k=1 ; k<x-1 ; k++){ xmat->a[k]=r/2.0; xmat->b[k]=1.0-r; xmat->c[k]=r/2.0; }
                             xmat->a[x-1]=r/2.0;    xmat->b[x-1]=1.0-r/2.0;
                             for(k=0 ; k<x ; k++) xmat->u[k]=uadv[k];
                             trimul(xmat);
                             /* LHS (I - r/2 L): solve LHS * unew = r */
                             xmat->b[0]=1.0+r/2.0;  xmat->c[0]=-r/2.0;
                             for(k=1 ; k<x-1 ; k++){ xmat->a[k]=-r/2.0; xmat->b[k]=1.0+r; xmat->c[k]=-r/2.0; }
                             xmat->a[x-1]=-r/2.0;   xmat->b[x-1]=1.0+r/2.0;
                             tridag(xmat);
                             for(k=0 ; k<x ; k++) community->pelagic[s].u_values[i][k][l]=xmat->u[k];
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
     int i,k,l,s;
     int x,y,n;
     double dt,dy,C,D,r;

     dt=grid->tstep;
     dy=grid->ystep;
     x=grid->xnum;
     y=grid->ynum;
     n=run->no_pelagic;

     double *uadv=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
     double *flux=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);   /* flux[l] = flux across the face between cells l and l+1 */
     double *vface=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);  /* advective velocity at face l (depends on g,mu only) */

     for(s=0 ; s<n ; s++){
             for(i=community->pelagic[s].ipelmin ; i<(community->pelagic[s].ipelmax+1) ; i++){
                     C=Cfun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]));
                     D=Dfun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]));
                     r=Diffun(grid->m_values[i], grid->t_values[time], grid, &(community->pelagic[s]))*(dt/(dy*dy));
                     int noadv=(C==0.0 && D==0.0);   /* no advection -> identity; also avoids 0*NaN if g/mu carry NaN */
                     for(k=0 ; k<x ; k++){
                             /* --- 1. advection: conservative upwind flux form; zero flux at both walls.
                                CFL sub-stepped (see xmove_solver for the rationale). --- */
                             double vmax=0.0;
                             for(l=0 ; l<y-1 ; l++){
                                     vface[l]= noadv ? 0.0 :
                                               (C*(community->pelagic[s].g_values[i][k][l+1]-community->pelagic[s].g_values[i][k][l])
                                               -D*(community->pelagic[s].mu_values[i][k][l+1]-community->pelagic[s].mu_values[i][k][l]))/dy;
                                     if(fabs(vface[l])>vmax) vmax=fabs(vface[l]);
                             }
                             for(l=0 ; l<y ; l++) uadv[l]=temp_pelagic[s][i][k][l];
                             int nsub=(int)ceil(vmax*dt/dy); if(nsub<1) nsub=1;
                             double dts=dt/(double)nsub;
                             for(int sub=0 ; sub<nsub ; sub++){
                                     for(l=0 ; l<y-1 ; l++){
                                             double uup=(vface[l]>0.0)? uadv[l] : uadv[l+1];   /* upwind donor */
                                             flux[l]=vface[l]*uup;
                                     }
                                     for(l=0 ; l<y ; l++){
                                             double Fr=(l<y-1)? flux[l]   : 0.0;   /* right-face flux (0 at the y-max wall) */
                                             double Fl=(l>0)  ? flux[l-1] : 0.0;   /* left-face  flux (0 at the y=0   wall) */
                                             uadv[l]=uadv[l]-(dts/dy)*(Fr-Fl);     /* flux[] precomputed -> in-place is safe */
                                             if(uadv[l]<0.0) uadv[l]=0.0;
                                     }
                             }
                             /* --- 2. diffusion: Crank-Nicolson with zero-flux (Neumann) BCs --- */
                             ymat->b[0]=1.0-r/2.0;  ymat->c[0]=r/2.0;
                             for(l=1 ; l<y-1 ; l++){ ymat->a[l]=r/2.0; ymat->b[l]=1.0-r; ymat->c[l]=r/2.0; }
                             ymat->a[y-1]=r/2.0;    ymat->b[y-1]=1.0-r/2.0;
                             for(l=0 ; l<y ; l++) ymat->u[l]=uadv[l];
                             trimul(ymat);
                             ymat->b[0]=1.0+r/2.0;  ymat->c[0]=-r/2.0;
                             for(l=1 ; l<y-1 ; l++){ ymat->a[l]=-r/2.0; ymat->b[l]=1.0+r; ymat->c[l]=-r/2.0; }
                             ymat->a[y-1]=-r/2.0;   ymat->b[y-1]=1.0+r/2.0;
                             tridag(ymat);
                             for(l=0 ; l<y ; l++) community->pelagic[s].u_values[i][k][l]=ymat->u[l];
                     }
             }
     }
     free(uadv);
     free(flux);
     free(vface);
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
