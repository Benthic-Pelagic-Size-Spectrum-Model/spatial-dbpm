#include "movement.h"

//--------------------------------------------//
//      Spatial Movement Functions            //
// These contain functions that may be edited //
//--------------------------------------------//

static double Cfun(double m, double t, GRID *grid, PELAGIC *pelagic)
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

static double Dfun(double m, double t, GRID *grid, PELAGIC *pelagic)
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

static double Diffun(double m, double t, GRID *grid, PELAGIC *pelagic)
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

// Helper function that checks if two neighbouring (horizontal) cells are both water
// and would therefore accept movement fluxes across their boundary
static int horizontal_face_valid(GRID *grid, int x_idx, int y_idx) {
    return !is_land(grid, x_idx, y_idx) && !is_land(grid, x_idx + 1, y_idx);
}

// Helper function that checks if two neighbouring (vertical) cells are both water
// and would therefore accept movement fluxes across their boundary
static int vertical_face_valid(GRID *grid, int x_idx, int y_idx) {
    return !is_land(grid, x_idx, y_idx) && !is_land(grid, x_idx, y_idx + 1);
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
    int *face_open = (int *)safe_malloc(n_x * sizeof(int), __LINE__, __FILE__); // Holds whether faces are open for movement (no land cells)

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
                for (x_idx=0; x_idx<n_x - 1; x_idx++) {
                    face_open[x_idx] = horizontal_face_valid(grid, x_idx, y_idx);
                }
                
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
                    vface[x_idx] = noadv || !face_open[x_idx]
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
                for (x_idx=0; x_idx < n_x; x_idx++) {
                    // Use face_open [0,1] as a multiplier to disable diffusion across closed faces (faces where land cells are present).
                    double rl = (x_idx > 0) ? r * face_open[x_idx - 1] : 0.0; // Left face of cell
                    double rr = (x_idx < n_x - 1) ? r * face_open[x_idx] : 0.0; // Right face of cell

                    xmat->a[x_idx] = rl / 2.0;
                    xmat->b[x_idx] = 1.0 - (rl + rr) / 2.0;
                    xmat->c[x_idx] = rr / 2.0;
                    xmat->u[x_idx] = uadv[x_idx];
                }
                trimul(xmat);

                /* LHS (I - r/2 L): solve LHS * unew = r */
                for (x_idx=0; x_idx < n_x; x_idx++) {
                    // Use face_open [0,1] as a multiplier to disable diffusion across closed faces (faces where land cells are present).
                    double rl = (x_idx > 0) ? r * face_open[x_idx - 1] : 0.0; // Left face of cell
                    double rr = (x_idx < n_x - 1) ? r * face_open[x_idx] : 0.0; // Right face of cell

                    xmat->a[x_idx] = -rl / 2.0;
                    xmat->b[x_idx] = 1.0 + (rl + rr) / 2.0;
                    xmat->c[x_idx] = -rr / 2.0;
                }
                tridag(xmat);

                for (x_idx = 0; x_idx < n_x; x_idx++) {
                    if (!is_land(grid, x_idx, y_idx)) {
                        community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] = xmat->u[x_idx];
                    }
                }
            }
        }
    }
    free(uadv);
    free(flux);
    free(vface);
    free(face_open);
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
    int *face_open = (int *)safe_malloc(n_y * sizeof(int), __LINE__, __FILE__); // Holds whether faces are open for movement (no land cells)

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
                for (y_idx=0; y_idx < n_y - 1; y_idx++) { // Precompute where valid (water to water) faces occur between cells
                    face_open[y_idx] = vertical_face_valid(grid, x_idx, y_idx);
                }
                /* --- 1. advection: conservative upwind flux form; zero flux at both walls.
                   CFL sub-stepped (see xmove_solver for the rationale). --- */
                double vmax = 0.0;
                for (y_idx = 0; y_idx < n_y - 1; y_idx++)
                {
                    vface[y_idx] = (noadv || !face_open[y_idx])
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
                // RHS (I + r/2 L):
                for (y_idx=0; y_idx<n_y; y_idx++) {
                    double rl = (y_idx > 0) ? r * face_open[y_idx - 1] : 0.0; // face with prior cell
                    double rr = (y_idx < n_y - 1) ? r * face_open[y_idx] : 0.0; // face with next cell

                    ymat->a[y_idx] = rl / 2.0;
                    ymat->b[y_idx] = 1.0 - (rl + rr) / 2.0;
                    ymat->c[y_idx] = rr / 2.0;
                    ymat->u[y_idx] = uadv[y_idx];
                }
                trimul(ymat);

                // LHS (I - r/2 L)
                for (y_idx=0; y_idx<n_y; y_idx++) {
                    double rl = (y_idx > 0) ? r * face_open[y_idx - 1] : 0.0; // face with prior cell
                    double rr = (y_idx < n_y - 1) ? r * face_open[y_idx] : 0.0; // face with next cell

                    ymat->a[y_idx] = -rl / 2.0;
                    ymat->b[y_idx] = 1.0 + (rl + rr) / 2.0;
                    ymat->c[y_idx] = -rr / 2.0;
                }
                tridag(ymat);

                for (y_idx = 0; y_idx < n_y; y_idx++) {
                    if (!is_land(grid, x_idx, y_idx)) {
                        community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] = ymat->u[y_idx];
                    }
                }                
            }
        }
    }
    free(uadv);
    free(flux);
    free(vface);
    free(face_open);
}