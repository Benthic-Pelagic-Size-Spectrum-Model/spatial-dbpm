#include "common.h"

int is_land(GRID *grid, int x_idx, int y_idx) {
    return (grid->land_mask[x_idx][y_idx] == 0);
}

void trimul(MATRIX *matrix)
{
    int i, n;
    n = matrix->size;

    matrix->r[0] = (matrix->b[0] * matrix->u[0]) + (matrix->c[0] * matrix->u[1]);
    for (i = 1; i < n - 1; i++)
    {
        matrix->r[i] =
            (matrix->a[i] * matrix->u[i - 1]) + (matrix->b[i] * matrix->u[i]) + (matrix->c[i] * matrix->u[i + 1]);
    }
    matrix->r[n - 1] = (matrix->a[n - 1] * matrix->u[n - 2]) + (matrix->b[n - 1] * matrix->u[n - 1]);
}

void tridag(MATRIX *matrix)
/*Fast algorithm for inverting a tridiagonal matrix*/
{
    int j, n;
    n = matrix->size;
    double bet;
    double *gam;

    gam = (double *)safe_malloc(n * sizeof(double), __LINE__, __FILE__);
    if (matrix->b[0] == 0)
    {
        error("error1");
    }
    bet = matrix->b[0];
    matrix->u[0] = matrix->r[0] / bet;
    for (j = 1; j <= (n - 1); j++)
    {
        gam[j] = matrix->c[j - 1] / bet;
        bet = matrix->b[j] - (matrix->a[j] * gam[j]);
        if (matrix->b[j] == 0)
        {
            error("error2");
        }
        matrix->u[j] = (matrix->r[j] - (matrix->a[j] * matrix->u[j - 1])) / bet;
    }
    for (j = (n - 2); j >= 0; j--)
    {
        matrix->u[j] -= (gam[j + 1] * matrix->u[j + 1]);
    }
    free(gam);
}

//---------------------------------------------------------//
// These are memory allocation functions. Not to be edited //
//---------------------------------------------------------//

void free_mem(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *pelmatrix, MATRIX *benmatrix, MATRIX *xmatrix,
              MATRIX *ymatrix, double **temp_pel_params, double **temp_ben_params)
/*Frees all memory allocated within code*/
{
    int sc_idx, x_idx, pel_species_idx, pel_sub_idx, ben_species_idx;
    int n_m, n_x, n_pel_species, n_ben_species;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_pel_species = run->no_pelagic;
    n_ben_species = run->no_benthic;

    /*Free grid values*/
    free(grid->m_values);
    free(grid->exp_m_values);
    free(grid->t_values);
    free(grid->x_values);
    free(grid->y_values);

    for (x_idx=0; x_idx<n_x; x_idx++) {
        free(grid->land_mask[x_idx]);
    }
    free(grid->land_mask);

    /*Free pelagic values*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        for (sc_idx = 0; sc_idx < n_m; sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                free(community->pelagic[pel_species_idx].u_values[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].g_values[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].mu_values[sc_idx][x_idx]);
                for (pel_sub_idx = 0; pel_sub_idx < n_pel_species; pel_sub_idx++)
                {
                    free(community->pelagic[pel_species_idx].mu_pred_values[pel_sub_idx][sc_idx][x_idx]);
                }
                free(community->pelagic[pel_species_idx].mu_fish_values[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].pla_bio[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].pel_bio[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].ben_bio[sc_idx][x_idx]);
                free(community->pelagic[pel_species_idx].a_refuge[sc_idx][x_idx]);
            }
            free(community->pelagic[pel_species_idx].u_values[sc_idx]);
            free(community->pelagic[pel_species_idx].g_values[sc_idx]);
            free(community->pelagic[pel_species_idx].mu_values[sc_idx]);
            for (pel_sub_idx = 0; pel_sub_idx < n_pel_species; pel_sub_idx++)
            {
                free(community->pelagic[pel_species_idx].mu_pred_values[pel_sub_idx][sc_idx]);
            }
            free(community->pelagic[pel_species_idx].mu_fish_values[sc_idx]);
            free(community->pelagic[pel_species_idx].pla_bio[sc_idx]);
            free(community->pelagic[pel_species_idx].pel_bio[sc_idx]);
            free(community->pelagic[pel_species_idx].ben_bio[sc_idx]);
            free(community->pelagic[pel_species_idx].a_refuge[sc_idx]);
        }
        free(community->pelagic[pel_species_idx].u_values);
        free(community->pelagic[pel_species_idx].g_values);
        free(community->pelagic[pel_species_idx].mu_values);
        for (pel_sub_idx = 0; pel_sub_idx < n_pel_species; pel_sub_idx++)
        {
            free(community->pelagic[pel_species_idx].mu_pred_values[pel_sub_idx]);
        }
        free(community->pelagic[pel_species_idx].mu_fish_values);
        free(community->pelagic[pel_species_idx].pla_bio);
        free(community->pelagic[pel_species_idx].pel_bio);
        free(community->pelagic[pel_species_idx].ben_bio);
        free(community->pelagic[pel_species_idx].a_refuge);
        free(community->pelagic[pel_species_idx].phi_values);
        free(community->pelagic[pel_species_idx].refuge_values);
        free(temp_pel_params[pel_species_idx]);

        if (run->write_files)
        {
            free(community->pelagic[pel_species_idx].fname_pred);
            free(community->pelagic[pel_species_idx].fptr_pred);
        }
    }
    free(temp_pel_params);

    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            free(community->pelagic[pel_species_idx].pla_total[x_idx]);
            free(community->pelagic[pel_species_idx].pel_total[x_idx]);
            free(community->pelagic[pel_species_idx].ben_total[x_idx]);
            free(community->pelagic[pel_species_idx].fish_total[x_idx]);
            free(community->pelagic[pel_species_idx].pred_total[x_idx]);
            free(community->pelagic[pel_species_idx].reproduction[x_idx]);
            free(community->pelagic[pel_species_idx].rec0_c[x_idx]);
            free(community->pelagic[pel_species_idx].rec0_frozen_c[x_idx]);
        }
        free(community->pelagic[pel_species_idx].pla_total);
        free(community->pelagic[pel_species_idx].pel_total);
        free(community->pelagic[pel_species_idx].ben_total);
        free(community->pelagic[pel_species_idx].fish_total);
        free(community->pelagic[pel_species_idx].pred_total);
        free(community->pelagic[pel_species_idx].reproduction);
        free(community->pelagic[pel_species_idx].rec0_c);
        free(community->pelagic[pel_species_idx].rec0_frozen_c);
    }
    free(community->pelagic);

    /*Free Plankton pointer*/
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            free(community->plankton->u_values[sc_idx][x_idx]);
            free(community->plankton->g_values[sc_idx][x_idx]);
        }
        free(community->plankton->u_values[sc_idx]);
        free(community->plankton->g_values[sc_idx]);
    }
    free(community->plankton->u_values);
    free(community->plankton->g_values);
    free(community->plankton);

    if (run->no_benthic != 0)
    {
        /*Free Benthic pointer*/
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            for (sc_idx = 0; sc_idx < n_m; sc_idx++)
            {
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    free(community->benthic[ben_species_idx].u_values[sc_idx][x_idx]);
                    free(community->benthic[ben_species_idx].g_values[sc_idx][x_idx]);
                    free(community->benthic[ben_species_idx].mu_values[sc_idx][x_idx]);
                    free(community->benthic[ben_species_idx].mu_pred_values[sc_idx][x_idx]);
                    free(community->benthic[ben_species_idx].mu_fish_values[sc_idx][x_idx]);
                    free(community->benthic[ben_species_idx].det_bio[sc_idx][x_idx]);
                }
                free(community->benthic[ben_species_idx].u_values[sc_idx]);
                free(community->benthic[ben_species_idx].g_values[sc_idx]);
                free(community->benthic[ben_species_idx].mu_values[sc_idx]);
                free(community->benthic[ben_species_idx].mu_pred_values[sc_idx]);
                free(community->benthic[ben_species_idx].mu_fish_values[sc_idx]);
                free(community->benthic[ben_species_idx].det_bio[sc_idx]);
            }
            free(community->benthic[ben_species_idx].u_values);
            free(community->benthic[ben_species_idx].g_values);
            free(community->benthic[ben_species_idx].mu_values);
            free(community->benthic[ben_species_idx].mu_pred_values);
            free(community->benthic[ben_species_idx].mu_fish_values);
            free(community->benthic[ben_species_idx].det_bio);
            free(temp_ben_params[ben_species_idx]);
        }
        free(temp_ben_params);
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                free(community->benthic[ben_species_idx].det_total[x_idx]);
                free(community->benthic[ben_species_idx].fish_total[x_idx]);
                free(community->benthic[ben_species_idx].pred_total[x_idx]);
                free(community->benthic[ben_species_idx].reproduction[x_idx]);
            }
            free(community->benthic[ben_species_idx].det_total);
            free(community->benthic[ben_species_idx].fish_total);
            free(community->benthic[ben_species_idx].pred_total);
            free(community->benthic[ben_species_idx].reproduction);
        }
        free(community->benthic);

        /*Free Detritus*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            free(community->detritus->w_values[x_idx]);
            free(community->detritus->g_values[x_idx]);
        }
        free(community->detritus->w_values);
        free(community->detritus->g_values);
        free(community->detritus);
    }

    /*Free matrix pointers*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        free(pelmatrix[pel_species_idx].a);
        free(pelmatrix[pel_species_idx].b);
        free(pelmatrix[pel_species_idx].c);
        free(pelmatrix[pel_species_idx].r);
        free(pelmatrix[pel_species_idx].u);
    }
    free(pelmatrix);

    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            free(benmatrix[ben_species_idx].a);
            free(benmatrix[ben_species_idx].b);
            free(benmatrix[ben_species_idx].c);
            free(benmatrix[ben_species_idx].r);
            free(benmatrix[ben_species_idx].u);
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
    ptr = malloc(size);
    if (ptr == NULL)
    {
        error("Out of memory at line %d file %s\n", line, filename);
    }
    return ptr;
}

FILE *safe_fopen(char *filename, char *mode, int line, char *fname)
/*Safe version of fopen with built in error checking. Exits if the file cannot be opened*/
{
    FILE *fptr;
    fptr = fopen(filename, mode);
    if (fptr == NULL)
    {
        error("Unable to open file at line %d file %s\n", line, fname);
    }
    return fptr;
}

// Biological function that is used across both biology.c and setup.c files
double phi(double q, double q_0, double sig, double trunc)
/*function determining prey preference distribution for the given values*/
{
    double ans;

    if (q <= 0 || q < (q_0 - (trunc * sig)) || q > (q_0 + (trunc * sig)))
    {
        ans = 0;
    }
    else
    {
        ans = exp(-(q - q_0) * (q - q_0) / (2 * sig * sig));
    }

    return ans;
}