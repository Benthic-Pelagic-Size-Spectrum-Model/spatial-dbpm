#include "biology.h"

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
                    SKIP_IF_LAND(grid, x_idx, y_idx);

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
                    SKIP_IF_LAND(grid, x_idx, y_idx);

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
                SKIP_IF_LAND(grid, x_idx, y_idx);

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
                            grid->exp_m_values[sc_idx] * grid->mstep;
                    for (pel_sub_idx = 0; pel_sub_idx < n_pel_species; pel_sub_idx++)
                    {
                        pred += community->pelagic[pel_species_idx].mu_pred_values[pel_sub_idx][sc_idx][x_idx][y_idx] *
                                community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] *
                                grid->exp_m_values[sc_idx] * grid->mstep;
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
                        SKIP_IF_LAND(grid, x_idx, y_idx);

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
                    SKIP_IF_LAND(grid, x_idx, y_idx);

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
                                grid->exp_m_values[sc_idx] * grid->mstep;
                        pred += community->benthic[ben_species_idx].mu_pred_values[sc_idx][x_idx][y_idx] *
                                community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] *
                                grid->exp_m_values[sc_idx] * grid->mstep;
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
                SKIP_IF_LAND(grid, x_idx, y_idx);
                
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
        ans += grid->exp_m_values[sc_idx] * community->pelagic[species].phi_values[size - sc_idx] *
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
                       grid->exp_m_values[sc_idx] * community->pelagic[species].phi_values[size - sc_idx] *
                       community->pelagic[pel_species_idx].u_values[sc_idx][xspace][yspace] * grid->mstep;
            }
        }
    }

    /* Allows feeding from own Size Spectra Only*/
    else
    {
        for (sc_idx = community->pelagic[species].ipelmin; sc_idx < size; sc_idx++)
        {
            ans += community->pelagic[species].a_refuge[sc_idx][xspace][yspace] * grid->exp_m_values[sc_idx] *
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
                    ans += grid->exp_m_values[sc_idx] * community->pelagic[species].phi_values[size - sc_idx] *
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
                        grid->exp_m_values[sc_idx] * grid->mstep +
                    community->pelagic[pel_species_idx].mu_s *
                        ((log10(grid->exp_m_values[sc_idx]) - log10(exp(community->pelagic[pel_species_idx].pelmin))) /
                         ((log10(exp(community->pelagic[pel_species_idx].pelmax)) +
                           community->pelagic[pel_species_idx].epsilon) -
                          log10(grid->exp_m_values[sc_idx]))) *
                        community->pelagic[pel_species_idx].u_values[sc_idx][xspace][yspace] *
                        grid->exp_m_values[sc_idx] * grid->mstep;
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
                               grid->exp_m_values[sc_idx] * grid->mstep;
                    ans_bed +=
                        community->benthic[ben_species_idx].mu_s *
                        ((log10(grid->exp_m_values[sc_idx]) - log10(exp(community->benthic[ben_species_idx].benmin))) /
                         ((log10(exp(community->benthic[ben_species_idx].benmax)) +
                           community->benthic[ben_species_idx].epsilon) -
                          log10(grid->exp_m_values[sc_idx]))) *
                        community->benthic[ben_species_idx].u_values[sc_idx][xspace][yspace] *
                        grid->exp_m_values[sc_idx] * grid->mstep;
                }
            }
        }
    }
    for (sc_idx = community->plankton->iplamin; sc_idx < (community->plankton->iplamax + 1); sc_idx++)
    {
        /*Rate of det biomass in from dead (sinking) plankton - surface origin*/
        double pdead = community->plankton->mu_0 * exp(community->plankton->beta * grid->m_values[sc_idx]) *
                       community->plankton->u_values[sc_idx][xspace][yspace] * grid->exp_m_values[sc_idx] * grid->mstep;
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
            ((log10(grid->exp_m_values[size]) - log10(exp(community->pelagic[species].pelmin))) /
             ((log10(exp(community->pelagic[species].pelmax)) + community->pelagic[species].epsilon) -
              log10(grid->exp_m_values[size]))));

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
            ((log10(grid->exp_m_values[size]) - log10(exp(community->benthic[species].benmin))) /
             ((log10(exp(community->benthic[species].benmax)) + community->benthic[species].epsilon) -
              log10(grid->exp_m_values[size]))));

    return ans;
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
                    SKIP_IF_LAND(grid, x_idx, y_idx);

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
                    SKIP_IF_LAND(grid, x_idx, y_idx);

                    pel = 0;
                    for (sc_idx = community->pelagic[pel_species_idx].ipelmat;
                         sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
                    {
                        pel += community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] *
                               grid->exp_m_values[sc_idx] * grid->mstep;
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
                        SKIP_IF_LAND(grid, x_idx, y_idx);

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
                        SKIP_IF_LAND(grid, x_idx, y_idx);

                        ben = 0;
                        for (sc_idx = community->benthic[ben_species_idx].ibenmat;
                             sc_idx < (community->benthic[ben_species_idx].ibenmax + 1); sc_idx++)
                        {
                            ben += community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] *
                                   grid->exp_m_values[sc_idx] * grid->mstep;
                        }
                        community->benthic[ben_species_idx].reproduction[x_idx][y_idx] = ben * 492 * grid->tstep;
                    }
                }
            }
        }
    }
}