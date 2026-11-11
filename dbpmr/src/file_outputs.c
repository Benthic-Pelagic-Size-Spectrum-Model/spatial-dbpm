#include <file_outputs.h>

void print_run(RUN *run, FILE *fptr)
{
    /*Output specific C code info*/
    fprintf(fptr, "Run Parameters\n");
    fprintf(fptr, "filename:      %s\n", run->filename);
    fprintf(fptr, "no_pelagic:    %d\n", run->no_pelagic);
    fprintf(fptr, "no_benthic:    %d\n", run->no_benthic);
    fprintf(fptr, "spatial_dim:   %d\n", run->spatial_dim);
    fprintf(fptr, "coupled_flag:  %d\n", run->coupled_flag);
    fprintf(fptr, "diff_method:   %d\n", run->diff_method);
    fprintf(fptr, "\n");
}

void print_grid(GRID *grid, FILE *fptr)
{
    /*Output discretisation information*/
    fprintf(fptr, "Grid Parameters\n");
    fprintf(fptr, "mmin:          %f\n", grid->mmin);
    fprintf(fptr, "mmax:          %f\n", grid->mmax);
    fprintf(fptr, "mstep:         %f\n", grid->mstep);
    fprintf(fptr, "moutstep:      %f\n", grid->moutstep);

    fprintf(fptr, "t1:            %.10f\n", grid->t1);
    fprintf(fptr, "tmax:          %.10f\n", grid->tmax);
    fprintf(fptr, "tstep:         %.10f\n", grid->tstep);
    fprintf(fptr, "toutmin:       %.10f\n", grid->toutmin);
    fprintf(fptr, "toutmax:       %.10f\n", grid->toutmax);
    fprintf(fptr, "toutstep:      %.10f\n", grid->toutstep);

    fprintf(fptr, "xmin:          %f\n", grid->xmin);
    fprintf(fptr, "xmax:          %f\n", grid->xmax);
    fprintf(fptr, "xstep:         %f\n", grid->xstep);
    fprintf(fptr, "xoutstep:      %f\n", grid->xoutstep);

    fprintf(fptr, "ymin:          %f\n", grid->ymin);
    fprintf(fptr, "ymax:          %f\n", grid->ymax);
    fprintf(fptr, "ystep:         %f\n", grid->ystep);
    fprintf(fptr, "youtstep:      %f\n", grid->youtstep);
    fprintf(fptr, "\n");
}

void print_plankton(PLANKTON *plankton, FILE *fptr)
{
    /*Ouput plankton specific information*/
    fprintf(fptr, "Plankton Parameters\n");
    fprintf(fptr, "filename:      %s\n", plankton->filename);
    fprintf(fptr, "speciestype:   %s\n", "plankton");
    fprintf(fptr, "mmin:        %f\n", plankton->plamin);
    fprintf(fptr, "mmax:        %f\n", plankton->plamax);
    fprintf(fptr, "mu_0:          %f\n", plankton->mu_0);
    fprintf(fptr, "beta:          %f\n", plankton->beta);
    fprintf(fptr, "u_0:           %f\n", plankton->u_0);
    fprintf(fptr, "lambda:        %f\n", plankton->lambda);
    fprintf(fptr, "initial_flag:  %d\n", plankton->initial_flag);
    fprintf(fptr, "ts_flag:       %d\n", plankton->ts_flag);
    fprintf(fptr, "\n");
}

void print_pelagic(PELAGIC *pelagic, FILE *fptr)
{
    /*Ouput species specific information*/
    fprintf(fptr, "Pelagic Parameters\n");
    fprintf(fptr, "filename:      %s\n", pelagic->filename);
    fprintf(fptr, "speciestype:   %s\n", "pelagic");

    fprintf(fptr, "mmin:        %f\n", pelagic->pelmin);
    fprintf(fptr, "mmat:        %f\n", pelagic->pelmat);
    fprintf(fptr, "mmax:        %f\n", pelagic->pelmax);

    fprintf(fptr, "A:             %f\n", pelagic->A);
    fprintf(fptr, "alpha:         %f\n", pelagic->alpha);
    fprintf(fptr, "mu_0:          %f\n", pelagic->mu_0);
    fprintf(fptr, "beta:          %f\n", pelagic->beta);
    fprintf(fptr, "mu_s:          %f\n", pelagic->mu_s);
    fprintf(fptr, "epsilon:       %f\n", pelagic->epsilon);
    fprintf(fptr, "u_0:           %f\n", pelagic->u_0);
    fprintf(fptr, "lambda:        %f\n", pelagic->lambda);

    fprintf(fptr, "K_pla:         %f\n", pelagic->K_pla);
    fprintf(fptr, "R_pla:         %f\n", pelagic->R_pla);
    fprintf(fptr, "Ex_pla:        %f\n", pelagic->Ex_pla);
    fprintf(fptr, "K_pel:         %f\n", pelagic->K_pel);
    fprintf(fptr, "R_pel:         %f\n", pelagic->R_pel);
    fprintf(fptr, "Ex_pel:        %f\n", pelagic->Ex_pel);
    fprintf(fptr, "K_ben:         %f\n", pelagic->K_ben);
    fprintf(fptr, "R_ben:         %f\n", pelagic->R_ben);
    fprintf(fptr, "Ex_ben:        %f\n", pelagic->Ex_ben);

    fprintf(fptr, "pref_pla:      %f\n", pelagic->pref_pla);
    fprintf(fptr, "pref_pel:      %f\n", pelagic->pref_pel);
    fprintf(fptr, "pref_ben:      %f\n", pelagic->pref_ben);

    fprintf(fptr, "q_0:           %f\n", pelagic->q_0);
    fprintf(fptr, "sig:           %f\n", pelagic->sig);
    fprintf(fptr, "trunc:         %f\n", pelagic->trunc);

    fprintf(fptr, "prey:          %f\n", pelagic->prey);
    fprintf(fptr, "pred:          %f\n", pelagic->pred);
    fprintf(fptr, "comp:          %f\n", pelagic->comp);
    fprintf(fptr, "gamma_pred:    %f\n", pelagic->gamma_pred);
    fprintf(fptr, "gamma_prey:    %f\n", pelagic->gamma_prey);
    fprintf(fptr, "gamma_comp:    %f\n", pelagic->gamma_comp);

    fprintf(fptr, "rep_method:    %d\n", pelagic->rep_method);
    fprintf(fptr, "initial_flag:  %d\n", pelagic->initial_flag);
    fprintf(fptr, "ts_flag:       %d\n", pelagic->ts_flag);
    fprintf(fptr, "fishing_flag:  %d\n", pelagic->fishing_flag);
    fprintf(fptr, "\n");
}

void print_benthic(BENTHIC *benthic, FILE *fptr)
{
    /*Ouput species specific information*/
    fprintf(fptr, "Benthic Parameters\n");
    fprintf(fptr, "filename:      %s\n", benthic->filename);
    fprintf(fptr, "speciestype:   %s\n", "benthic");

    fprintf(fptr, "mmin:        %f\n", benthic->benmin);
    fprintf(fptr, "mmat:        %f\n", benthic->benmat);
    fprintf(fptr, "mmax:        %f\n", benthic->benmax);

    fprintf(fptr, "A:             %f\n", benthic->A);
    fprintf(fptr, "alpha:         %f\n", benthic->alpha);
    fprintf(fptr, "mu_0:          %f\n", benthic->mu_0);
    fprintf(fptr, "beta:          %f\n", benthic->beta);
    fprintf(fptr, "mu_s:          %f\n", benthic->mu_s);
    fprintf(fptr, "epsilon:       %f\n", benthic->epsilon);
    fprintf(fptr, "u_0:           %f\n", benthic->u_0);
    fprintf(fptr, "lambda:        %f\n", benthic->lambda);

    fprintf(fptr, "K_det:         %f\n", benthic->K_det);
    fprintf(fptr, "R_det:         %f\n", benthic->R_det);
    fprintf(fptr, "Ex_det:        %f\n", benthic->Ex_det);

    fprintf(fptr, "pref_det:      %f\n", benthic->pref_det);

    fprintf(fptr, "rep_method:    %d\n", benthic->rep_method);
    fprintf(fptr, "initial_flag:  %d\n", benthic->initial_flag);
    fprintf(fptr, "ts_flag:       %d\n", benthic->ts_flag);
    fprintf(fptr, "fishing_flag:  %d\n", benthic->fishing_flag);
    fprintf(fptr, "\n");
}

void print_detritus(DETRITUS *detritus, FILE *fptr)
{
    /*Ouput plankton specific information*/
    fprintf(fptr, "Detritus Parameters\n");
    fprintf(fptr, "filename:      %s\n", detritus->filename);
    fprintf(fptr, "speciestype:   %s\n", "detritus");

    fprintf(fptr, "w_0:           %f\n", detritus->w_0);

    fprintf(fptr, "initial_flag:  %d\n", detritus->initial_flag);
    fprintf(fptr, "ts_flag:       %d\n", detritus->ts_flag);
    fprintf(fptr, "\n");
}

void print_mass_header(GRID *grid, FILE *fptr)
{
    /*Output Header Row*/
    int sc_idx;
    int n_m = grid->mnum;

    fprintf(fptr, "Ouput Data\n");
    fprintf(fptr, "t,x,y,");
    for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx += (int)rint(grid->moutstep / grid->mstep))
    {
        fprintf(fptr, "%f,", grid->m_values[sc_idx]);
    }
    fprintf(fptr, "%f\n", grid->m_values[n_m - 1]);
}

void print_detailed_header(GRID *grid, FILE *fptr)
{
    int sc_idx;
    int n_m;
    n_m = grid->mnum;

    /*Print Header*/
    fprintf(fptr, "\nDetailed Timestep\n");
    fprintf(fptr, "t,x,y,");
    for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx++)
    {
        fprintf(fptr, "%f,", grid->m_values[sc_idx]);
    }
    fprintf(fptr, "%f\n", grid->m_values[n_m - 1]);
}

void print_plankton_header(FILE *fptr)
{
    fprintf(fptr, "Ouput Data\n");
    fprintf(fptr, "t,x,y");
    fprintf(fptr, ",biomass");
    fprintf(fptr, "\n");
}

void print_pelagic_header(FILE *fptr)
{
    fprintf(fptr, "Ouput Data\n");
    fprintf(fptr, "t,x,y");
    fprintf(fptr, ",biomass");
    fprintf(fptr, ",plankton");
    fprintf(fptr, ",pelagic");
    fprintf(fptr, ",benthic");
    fprintf(fptr, ",predation");
    fprintf(fptr, ",fishing");
    fprintf(fptr, ",reproduction");
    fprintf(fptr, "\n");
}

void print_benthic_header(FILE *fptr)
{
    fprintf(fptr, "Ouput Data\n");
    fprintf(fptr, "t,x,y");
    fprintf(fptr, ",biomass");
    fprintf(fptr, ",detritus");
    fprintf(fptr, ",predation");
    fprintf(fptr, ",fishing");
    fprintf(fptr, ",reproduction");
    fprintf(fptr, "\n");
}

void print_detritus_header(FILE *fptr)
{
    fprintf(fptr, "Ouput Data\n");
    fprintf(fptr, "t,x,y");
    fprintf(fptr, ",biomass");
    fprintf(fptr, ",detritus_in");
    fprintf(fptr, ",detritus_out");
    fprintf(fptr, "\n");
}

void print_timestep_plankton(GRID *grid, PLANKTON *plankton, int time)
/*Outputs timestep values*/
{
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    plankton->fptr_r = safe_fopen(plankton->fname_r, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(plankton->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx += (int)rint(grid->moutstep / grid->mstep))
            {
                fprintf(plankton->fptr_r, "%.16f,", plankton->u_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(plankton->fptr_r, "%.16f\n", plankton->u_values[n_m - 1][x_idx][y_idx]);
        }
    }

    fclose(plankton->fptr_r);
}

void print_detailed_plankton(GRID *grid, PLANKTON *plankton, int time)
/*Outputs timestep values at high resolution*/
{
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    plankton->fptr_r = safe_fopen(plankton->fname_r, "a", __LINE__, __FILE__);

    /*Print Header*/
    print_detailed_header(grid, plankton->fptr_r);

    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            fprintf(plankton->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx++)
            {
                fprintf(plankton->fptr_r, "%.16f,", plankton->u_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(plankton->fptr_r, "%.16f\n", plankton->u_values[n_m - 1][x_idx][y_idx]);
        }
    }

    fclose(plankton->fptr_r);
}

void print_timestep_pelagic(RUN *run, GRID *grid, PELAGIC *pelagic, int time)
/*Outputs timestep values*/
{
    int pel_species_idx, sc_idx, x_idx, y_idx;
    int n_pel_species, n_m, n_x, n_y;
    n_pel_species = run->no_pelagic;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    pelagic->fptr_r = safe_fopen(pelagic->fname_r, "a", __LINE__, __FILE__);
    pelagic->fptr_g = safe_fopen(pelagic->fname_g, "a", __LINE__, __FILE__);
    pelagic->fptr_m = safe_fopen(pelagic->fname_m, "a", __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        pelagic->fptr_pred[pel_species_idx] = safe_fopen(pelagic->fname_pred[pel_species_idx], "a", __LINE__, __FILE__);
    }
    pelagic->fptr_fish = safe_fopen(pelagic->fname_fish, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(pelagic->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(pelagic->fptr_g, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(pelagic->fptr_m, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                fprintf(pelagic->fptr_pred[pel_species_idx], "%.10f,%f,%f,", (grid->t_values[time]),
                        grid->x_values[x_idx], grid->y_values[y_idx]);
            }
            fprintf(pelagic->fptr_fish, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx += (int)rint(grid->moutstep / grid->mstep))
            {
                fprintf(pelagic->fptr_r, "%.16f,", pelagic->u_values[sc_idx][x_idx][y_idx]);
                fprintf(pelagic->fptr_g, "%.16f,", pelagic->g_values[sc_idx][x_idx][y_idx]);
                fprintf(pelagic->fptr_m, "%.16f,", pelagic->mu_values[sc_idx][x_idx][y_idx]);
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    fprintf(pelagic->fptr_pred[pel_species_idx], "%.16f,",
                            pelagic->mu_pred_values[pel_species_idx][sc_idx][x_idx][y_idx]);
                }
                fprintf(pelagic->fptr_fish, "%.16f,", pelagic->mu_fish_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(pelagic->fptr_r, "%.16f\n", pelagic->u_values[n_m - 1][x_idx][y_idx]);
            fprintf(pelagic->fptr_g, "%.16f\n", pelagic->g_values[n_m - 1][x_idx][y_idx]);
            fprintf(pelagic->fptr_m, "%.16f\n", pelagic->mu_values[n_m - 1][x_idx][y_idx]);
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                fprintf(pelagic->fptr_pred[pel_species_idx], "%.16f\n",
                        pelagic->mu_pred_values[pel_species_idx][n_m - 1][x_idx][y_idx]);
            }
            fprintf(pelagic->fptr_fish, "%.16f\n", pelagic->mu_fish_values[n_m - 1][x_idx][y_idx]);
        }
    }

    /*Close output files*/
    fclose(pelagic->fptr_r);
    fclose(pelagic->fptr_g);
    fclose(pelagic->fptr_m);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        fclose(pelagic->fptr_pred[pel_species_idx]);
    }
    fclose(pelagic->fptr_fish);
}

void print_detailed_pelagic(RUN *run, GRID *grid, PELAGIC *pelagic, int time)
/*Outputs timestep values*/
{
    int pel_species_idx, sc_idx, x_idx, y_idx;
    int n_pel_species, n_m, n_x, n_y;
    n_pel_species = run->no_pelagic;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    pelagic->fptr_r = safe_fopen(pelagic->fname_r, "a", __LINE__, __FILE__);
    pelagic->fptr_g = safe_fopen(pelagic->fname_g, "a", __LINE__, __FILE__);
    pelagic->fptr_m = safe_fopen(pelagic->fname_m, "a", __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        pelagic->fptr_pred[pel_species_idx] = safe_fopen(pelagic->fname_pred[pel_species_idx], "a", __LINE__, __FILE__);
    }
    pelagic->fptr_fish = safe_fopen(pelagic->fname_fish, "a", __LINE__, __FILE__);

    /*Print Header*/
    print_detailed_header(grid, pelagic->fptr_r);
    print_detailed_header(grid, pelagic->fptr_g);
    print_detailed_header(grid, pelagic->fptr_m);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        print_detailed_header(grid, pelagic->fptr_pred[pel_species_idx]);
    }
    print_detailed_header(grid, pelagic->fptr_fish);

    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            fprintf(pelagic->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(pelagic->fptr_g, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(pelagic->fptr_m, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                fprintf(pelagic->fptr_pred[pel_species_idx], "%.10f,%f,%f,", (grid->t_values[time]),
                        grid->x_values[x_idx], grid->y_values[y_idx]);
            }
            fprintf(pelagic->fptr_fish, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx++)
            {
                fprintf(pelagic->fptr_r, "%.16f,", pelagic->u_values[sc_idx][x_idx][y_idx]);
                fprintf(pelagic->fptr_g, "%.16f,", pelagic->g_values[sc_idx][x_idx][y_idx]);
                fprintf(pelagic->fptr_m, "%.16f,", pelagic->mu_values[sc_idx][x_idx][y_idx]);
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    fprintf(pelagic->fptr_pred[pel_species_idx], "%.16f,",
                            pelagic->mu_pred_values[pel_species_idx][sc_idx][x_idx][y_idx]);
                }
                fprintf(pelagic->fptr_fish, "%.16f,", pelagic->mu_fish_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(pelagic->fptr_r, "%.16f\n", pelagic->u_values[n_m - 1][x_idx][y_idx]);
            fprintf(pelagic->fptr_g, "%.16f\n", pelagic->g_values[n_m - 1][x_idx][y_idx]);
            fprintf(pelagic->fptr_m, "%.16f\n", pelagic->mu_values[n_m - 1][x_idx][y_idx]);
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                fprintf(pelagic->fptr_pred[pel_species_idx], "%.16f\n",
                        pelagic->mu_pred_values[pel_species_idx][n_m - 1][x_idx][y_idx]);
            }
            fprintf(pelagic->fptr_fish, "%.16f\n", pelagic->mu_fish_values[n_m - 1][x_idx][y_idx]);
        }
    }

    /*Close output files*/
    fclose(pelagic->fptr_r);
    fclose(pelagic->fptr_g);
    fclose(pelagic->fptr_m);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        fclose(pelagic->fptr_pred[pel_species_idx]);
    }
    fclose(pelagic->fptr_fish);
}

void print_timestep_benthic(GRID *grid, BENTHIC *benthic, int time)
/*Outputs timestep values*/
{
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    benthic->fptr_r = safe_fopen(benthic->fname_r, "a", __LINE__, __FILE__);
    benthic->fptr_g = safe_fopen(benthic->fname_g, "a", __LINE__, __FILE__);
    benthic->fptr_m = safe_fopen(benthic->fname_m, "a", __LINE__, __FILE__);
    benthic->fptr_pred = safe_fopen(benthic->fname_pred, "a", __LINE__, __FILE__);
    benthic->fptr_fish = safe_fopen(benthic->fname_fish, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(benthic->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_g, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_m, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_pred, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_fish, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx += (int)rint(grid->moutstep / grid->mstep))
            {
                fprintf(benthic->fptr_r, "%.16f,", benthic->u_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_g, "%.16f,", benthic->g_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_m, "%.16f,", benthic->mu_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_pred, "%.16f,", benthic->mu_pred_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_fish, "%.16f,", benthic->mu_fish_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(benthic->fptr_r, "%.16f\n", benthic->u_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_g, "%.16f\n", benthic->g_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_m, "%.16f\n", benthic->mu_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_pred, "%.16f\n", benthic->mu_pred_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_fish, "%.16f\n", benthic->mu_fish_values[n_m - 1][x_idx][y_idx]);
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
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    benthic->fptr_r = safe_fopen(benthic->fname_r, "a", __LINE__, __FILE__);
    benthic->fptr_g = safe_fopen(benthic->fname_g, "a", __LINE__, __FILE__);
    benthic->fptr_m = safe_fopen(benthic->fname_m, "a", __LINE__, __FILE__);
    benthic->fptr_pred = safe_fopen(benthic->fname_pred, "a", __LINE__, __FILE__);
    benthic->fptr_fish = safe_fopen(benthic->fname_fish, "a", __LINE__, __FILE__);

    /*Print Header*/
    print_detailed_header(grid, benthic->fptr_r);
    print_detailed_header(grid, benthic->fptr_g);
    print_detailed_header(grid, benthic->fptr_m);
    print_detailed_header(grid, benthic->fptr_pred);
    print_detailed_header(grid, benthic->fptr_fish);

    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            fprintf(benthic->fptr_r, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_g, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_m, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_pred, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            fprintf(benthic->fptr_fish, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);
            for (sc_idx = 0; sc_idx < (n_m - 1); sc_idx++)
            {
                fprintf(benthic->fptr_r, "%.16f,", benthic->u_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_g, "%.16f,", benthic->g_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_m, "%.16f,", benthic->mu_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_pred, "%.16f,", benthic->mu_pred_values[sc_idx][x_idx][y_idx]);
                fprintf(benthic->fptr_fish, "%.16f,", benthic->mu_fish_values[sc_idx][x_idx][y_idx]);
            }
            fprintf(benthic->fptr_r, "%.16f\n", benthic->u_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_g, "%.16f\n", benthic->g_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_m, "%.16f\n", benthic->mu_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_pred, "%.16f\n", benthic->mu_pred_values[n_m - 1][x_idx][y_idx]);
            fprintf(benthic->fptr_fish, "%.16f\n", benthic->mu_fish_values[n_m - 1][x_idx][y_idx]);
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

    int sc_idx, x_idx, y_idx;
    int n_x, n_y;
    n_x = grid->xnum;
    n_y = grid->ynum;

    double temp;

    /*Open output files*/
    plankton->fptr_summ = safe_fopen(plankton->fname_summ, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(plankton->fptr_summ, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);

            /*Plankton Biomass Totals*/
            temp = 0;
            for (sc_idx = plankton->iplamin; sc_idx < (plankton->iplamax + 1); sc_idx++)
            {
                temp += plankton->u_values[sc_idx][x_idx][y_idx] * grid->exp_m_values[sc_idx] * grid->mstep;
            }
            fprintf(plankton->fptr_summ, "%.16f", temp);
            fprintf(plankton->fptr_summ, "\n");
        }
    }

    fclose(plankton->fptr_summ);
}

void print_pelagic_summary(GRID *grid, PELAGIC *pelagic, int time)
{

    int sc_idx, x_idx, y_idx;
    int n_x, n_y;
    n_x = grid->xnum;
    n_y = grid->ynum;

    double temp;

    /*Open output files*/
    pelagic->fptr_summ = safe_fopen(pelagic->fname_summ, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(pelagic->fptr_summ, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);

            /*Pelagic Biomass Totals*/
            temp = 0;
            for (sc_idx = pelagic->ipelmin; sc_idx < (pelagic->ipelmax + 1); sc_idx++)
            {
                temp += pelagic->u_values[sc_idx][x_idx][y_idx] * grid->exp_m_values[sc_idx] * grid->mstep;
            }
            fprintf(pelagic->fptr_summ, ",%.16f", temp);

            /*Pelagic Biomass Intake Rate*/
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->pla_total[x_idx][y_idx]);
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->pel_total[x_idx][y_idx]);
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->ben_total[x_idx][y_idx]);

            /*Biomass Lost to Predation*/
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->pred_total[x_idx][y_idx]);

            /*Biomass Lost to Fishing*/
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->fish_total[x_idx][y_idx]);

            /*Reproduction*/
            fprintf(pelagic->fptr_summ, ",%.16f", pelagic->reproduction[x_idx][y_idx] / grid->tstep);
            fprintf(pelagic->fptr_summ, "\n");
        }
    }

    fclose(pelagic->fptr_summ);
}

void print_benthic_summary(GRID *grid, BENTHIC *benthic, int time)
{
    int sc_idx, x_idx, y_idx;
    int n_x, n_y;
    n_x = grid->xnum;
    n_y = grid->ynum;

    double temp;

    /*Open output files*/
    benthic->fptr_summ = safe_fopen(benthic->fname_summ, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(benthic->fptr_summ, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);

            /*Benthic Biomass Totals*/
            temp = 0;
            for (sc_idx = benthic->ibenmin; sc_idx < (benthic->ibenmax + 1); sc_idx++)
            {
                temp += benthic->u_values[sc_idx][x_idx][y_idx] * grid->exp_m_values[sc_idx] * grid->mstep;
            }
            fprintf(benthic->fptr_summ, ",%.16f", temp);

            /*Benthic Biomass Intake*/
            fprintf(benthic->fptr_summ, ",%.16f", benthic->det_total[x_idx][y_idx]);

            /*Biomass Lost to Predation*/
            fprintf(benthic->fptr_summ, ",%.16f", benthic->pred_total[x_idx][y_idx]);

            /*Biomass Lost to Fishing*/
            fprintf(benthic->fptr_summ, ",%.16f", benthic->fish_total[x_idx][y_idx]);

            /*Reproduction*/
            fprintf(benthic->fptr_summ, ",%.16f", benthic->reproduction[x_idx][y_idx] / grid->tstep);
            fprintf(benthic->fptr_summ, "\n");
        }
    }
    fclose(benthic->fptr_summ);
}

void print_detritus_summary(GRID *grid, DETRITUS *detritus, int time)
{
    int x_idx, y_idx;
    int n_x, n_y;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Open output files*/
    detritus->fptr_summ = safe_fopen(detritus->fname_summ, "a", __LINE__, __FILE__);

    for (x_idx = 0; x_idx < n_x; x_idx += (int)rint(grid->xoutstep / grid->xstep))
    {
        for (y_idx = 0; y_idx < n_y; y_idx += (int)rint(grid->youtstep / grid->ystep))
        {
            fprintf(detritus->fptr_summ, "%.10f,%f,%f,", (grid->t_values[time]), grid->x_values[x_idx],
                    grid->y_values[y_idx]);

            /*Detritus Biomass Totals*/
            fprintf(detritus->fptr_summ, ",%.16f", detritus->w_values[x_idx][y_idx]);

            /*Detritus In*/
            fprintf(detritus->fptr_summ, ",%.16f", detritus->g_values[x_idx][y_idx]);

            /*Detritus Out*/
            fprintf(detritus->fptr_summ, ",%.16f", detritus->mu_values[x_idx][y_idx]);

            fprintf(detritus->fptr_summ, "\n");
        }
    }

    /*Close output files*/
    fclose(detritus->fptr_summ);
}