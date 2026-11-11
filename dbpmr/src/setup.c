#include "setup.h"

//-----------------//
// Setup Functions //
//-----------------//

void setup_run(RUN *run, int *run_params, char *filename, int *output_params)
{
    run->filename = filename;
    run->no_pelagic = run_params[0];
    run->no_benthic = run_params[1];
    run->spatial_dim = run_params[2];
    run->coupled_flag = run_params[3];
    run->diff_method = run_params[4];
    run->write_files = output_params[0];
    run->return_memory = output_params[1];
    run->movement_flag = run_params[5];

    /*Create summary output filenames*/
    snprintf(run->fname_summ, sizeof(run->fname_summ), "%s/%s.txt", run->filename, "parameters");

    /*Open output files*/
    run->fptr_summ = safe_fopen(run->fname_summ, "w", __LINE__, __FILE__);

    fprintf(run->fptr_summ, "filetype:      %s\n", "parameters");
    fprintf(run->fptr_summ, "speciestype:   %s\n", "NA");
    fprintf(run->fptr_summ, "\n");

    /*Close ouput files*/
    fclose(run->fptr_summ);
}

void setup_grid(GRID *grid, double *grid_input, int *land_mask, R_xlen_t mask_length)
{
    /*Declare Mass Discretisation*/
    grid->mmin = grid_input[0];
    grid->mmax = grid_input[1];
    grid->mstep = grid_input[2];
    // m step length
    grid->mnum = (int)rint(((grid->mmax - grid->mmin) / grid->mstep) + 1); // number of m values
    grid->moutstep = grid_input[3];

    /*Declare Time Discretisation*/
    grid->t1 = grid_input[4];
    grid->tmax = grid_input[5];                             // maximum t value
    grid->tstep = grid_input[6];                            // t step length
    grid->tnum = (int)rint((grid->tmax / grid->tstep) + 1); // number of t values
    grid->j1 = (int)rint(grid->t1 / grid->tstep);
    grid->toutmin = grid_input[7];
    grid->toutmax = grid_input[8];
    grid->toutstep = grid_input[9];

    /*Declare X-Space Discretisation*/
    grid->xmin = grid_input[10];
    grid->xmax = grid_input[11];
    grid->xstep = grid_input[12];
    grid->xnum = (int)rint(((grid->xmax - grid->xmin) / grid->xstep) + 1);
    grid->xoutstep = grid_input[13];

    /*Declare Y-Space Discretisation*/
    grid->ymin = grid_input[14];
    grid->ymax = grid_input[15];
    grid->ystep = grid_input[16];
    grid->ynum = (int)rint(((grid->ymax - grid->ymin) / grid->ystep) + 1);
    grid->youtstep = grid_input[17];

    /*Create arrays*/
    int sc_idx, ts_idx, x_idx, y_idx;
    int n_m, n_t, n_x, n_y;

    n_m = grid->mnum;
    n_t = grid->tnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    grid->m_values = (double *)safe_malloc(n_m * sizeof(double), __LINE__, __FILE__);
    grid->exp_m_values = (double *)safe_malloc(n_m * sizeof(double), __LINE__, __FILE__);
    grid->t_values = (double *)safe_malloc(n_t * sizeof(double), __LINE__, __FILE__);
    grid->x_values = (double *)safe_malloc(n_x * sizeof(double), __LINE__, __FILE__);
    grid->y_values = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);

    /*Setup values for m*/
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        grid->m_values[sc_idx] = grid->mmin + (sc_idx * grid->mstep);
        grid->exp_m_values[sc_idx] = exp(grid->m_values[sc_idx]);
    }
    /*Setup values for t*/
    for (ts_idx = 0; ts_idx < n_t; ts_idx++)
    {
        grid->t_values[ts_idx] = (ts_idx * grid->tstep);
    }
    /*Setup values for x*/
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        grid->x_values[x_idx] = grid->xmin + (x_idx * grid->xstep);
    }
    /*Setup values for y*/
    for (y_idx = 0; y_idx < n_y; y_idx++)
    {
        grid->y_values[y_idx] = grid->ymin + (y_idx * grid->ystep);
    }

    // Land mask setup -

    // Check the length of the mask
    if (mask_length != grid->xnum * grid->ynum) {
        error("land_mask must have a length of n_x * n_y.");
    }

    // Setup cell land mask vector from R inputs (the order of cell indices matches
    // the current order in C loops - Y changes fastest so cell indices are e.g.
    // [1,1], [1,2], [2,1], [2,2]). 
    // This index value is calculated two ways:
    // In R this index would be calculated `(x_idx - 1) * n_y + y_idx` but in C
    // it is calculated `x_idx * n_y + y_idx`
    grid->land_mask = (int **)safe_malloc(n_x * sizeof(int *), __LINE__, __FILE__);

    for (x_idx=0; x_idx < n_x; x_idx++) {
        grid->land_mask[x_idx] = (int *)safe_malloc(n_y * sizeof(int), __LINE__, __FILE__);
        for (y_idx=0; y_idx < n_y; y_idx++) {
            grid->land_mask[x_idx][y_idx] = land_mask[x_idx * n_y + y_idx];
        }
    }
}

void setup_plankton(RUN *run, GRID *grid, PLANKTON *plankton, double *pla_params, char *filename, int initial_flag,
                    int ts_flag)
{
    /*Declare plankton name*/
    plankton->filename = filename;

    /*Setup Plankton Sizes*/
    plankton->plamin = pla_params[0];
    plankton->plamax = pla_params[1];
    plankton->iplamin = (int)rint((plankton->plamin - grid->mmin) / grid->mstep);
    plankton->iplamax = (int)rint((plankton->plamax - grid->mmin) / grid->mstep);

    /*Setup Plankton Constants*/
    plankton->mu_0 = pla_params[2];
    plankton->beta = pla_params[3];
    plankton->u_0 = pla_params[4];
    plankton->lambda = pla_params[5];

    /*Setup Plankton Flags*/
    plankton->initial_flag = initial_flag;
    plankton->ts_flag = ts_flag;

    if (run->write_files)
    {
        /*Create plankton output filenames*/
        snprintf(plankton->fname_r, sizeof(plankton->fname_r), "%s/%s/%s.txt", run->filename, plankton->filename,
                 "results");
        snprintf(plankton->fname_summ, sizeof(plankton->fname_summ), "%s/%s/%s.txt", run->filename, plankton->filename,
                 "summary");

        /*Open output files*/
        plankton->fptr_r = safe_fopen(plankton->fname_r, "w", __LINE__, __FILE__);
        plankton->fptr_summ = safe_fopen(plankton->fname_summ, "w", __LINE__, __FILE__);

        /*Initialise output files*/
        fprintf(plankton->fptr_r, "filetype:      %s\n", "results");
        fprintf(plankton->fptr_r, "speciestype:   %s\n", "plankton");
        fprintf(plankton->fptr_r, "\n");
        print_mass_header(grid, plankton->fptr_r);

        fprintf(plankton->fptr_summ, "filetype:      %s\n", "summary");
        fprintf(plankton->fptr_summ, "speciestype:   %s\n", "plankton");
        fprintf(plankton->fptr_summ, "\n");
        print_plankton_header(plankton->fptr_summ);

        /*Close ouput files*/
        fclose(plankton->fptr_r);
        fclose(plankton->fptr_summ);
    }

    /*Setup And Initialise Plankton Abundance Arrays*/
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    plankton->u_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    plankton->g_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        plankton->u_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        plankton->g_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            plankton->u_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            plankton->g_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                plankton->u_values[sc_idx][x_idx][y_idx] = 0;
                plankton->g_values[sc_idx][x_idx][y_idx] = 0;
            }
        }
    }

    /*Setup initial values for u(m,0,x,y)*/
    /*Use user defined values*/
    if (plankton->initial_flag == 1 || plankton->ts_flag == 1)
    {
        /*Create plankton input filename*/
        snprintf(plankton->fname_ts, sizeof(plankton->fname_ts), "%s/Input/%s_%s.txt", run->filename,
                 plankton->filename, "ts");
        plankton->fptr_ts = safe_fopen(plankton->fname_ts, "r", __LINE__, __FILE__);

        int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
        char *input;                   // input string to store line as string
        input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
        char sep[] = ",";    // identifying the ',' char as seperating values
        char *result = NULL; // declaring a pointer for storing each individual value as a char

        /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                fgets(input, MAX_LEN, plankton->fptr_ts);
                input[strlen(input) - 1] =
                    '\0'; // replacing the newline \n char of the line of text with a null \0 char

                result = (char *)strtok(input, sep); // this is the 1st value
                sc_idx = 0;
                while (result != NULL)
                {
                    plankton->u_values[sc_idx][x_idx][y_idx] = atof(result);
                    result = (char *)strtok(NULL, sep);
                    sc_idx++;
                }
            }
        }
        free(input);

        /*If this is only an intial value then close the file now*/
        if (plankton->ts_flag == 0)
        {
            fclose(plankton->fptr_ts);
        }
    }
    /*Use standard values*/
    else
    {
        for (sc_idx = 0; sc_idx < (plankton->iplamax + 1); sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    plankton->u_values[sc_idx][x_idx][y_idx] =
                        plankton->u_0 * exp(plankton->lambda * grid->m_values[sc_idx]);
                    if (run->spatial_dim == 1 || run->spatial_dim == 2)
                    {
                        plankton->u_values[sc_idx][x_idx][y_idx] *= x_start(grid->x_values[x_idx], 0);
                        if (run->spatial_dim == 2)
                        {
                            plankton->u_values[sc_idx][x_idx][y_idx] *= y_start(grid->y_values[y_idx], 0);
                        }
                    }
                }
            }
        }
    }
}

void setup_pelagic(RUN *run, GRID *grid, PELAGIC *pelagic, double *pel_params, char *filename, int rep_method,
                   int initial_flag, int ts_flag, int fishing_flag)
{
    /*Pelagic Species Filename*/
    pelagic->filename = filename;

    /*Pelagic Species Size Range*/
    pelagic->pelmin = pel_params[0];
    pelagic->pelmat = pel_params[1];
    pelagic->pelmax = pel_params[2];
    pelagic->ipelmin = (int)rint((pelagic->pelmin - grid->mmin) / grid->mstep);
    pelagic->ipelmat = (int)rint((pelagic->pelmat - grid->mmin) / grid->mstep);
    pelagic->ipelmax = (int)rint((pelagic->pelmax - grid->mmin) / grid->mstep);

    /*Pelagic Species Constants*/
    pelagic->A = pel_params[3];         // Volume searched by unit weight
    pelagic->alpha = pel_params[4];     // Volume search allometry
    pelagic->mu_0 = pel_params[5];      // natural mortality rate
    pelagic->beta = pel_params[6];      // natural mortality allometry
    pelagic->A_base = pelagic->A;       // temp driver: base (unscaled) search rate
    pelagic->mu_0_base = pelagic->mu_0; // temp driver: base (unscaled) background mortality
    pelagic->tempeff = 1.0;             // temp driver: neutral until a temperature series is read
    pelagic->mu_s = pel_params[7];      // senescence mortality rate
    pelagic->epsilon = pel_params[8];   // senescece mortality constant
    pelagic->u_0 = pel_params[9];       // initial intercept
    pelagic->lambda = pel_params[10];   // primary production slope

    /*Pelagic Species Conversion Efficiencies*/
    pelagic->K_pla = pel_params[11]; // growth efficiency
    pelagic->R_pla = pel_params[12];
    pelagic->Ex_pla = pel_params[13];
    pelagic->K_pel = pel_params[14];
    pelagic->R_pel = pel_params[15];
    pelagic->Ex_pel = pel_params[16];
    pelagic->K_ben = pel_params[17];
    pelagic->R_ben = pel_params[18];
    pelagic->Ex_ben = pel_params[19];

    /*Pelagic Species Attack Rates*/
    pelagic->pref_pla = pel_params[20];
    pelagic->pref_pel = pel_params[21];
    pelagic->pref_ben = pel_params[22];

    /*Pelagic Species Feeding Kernel*/
    pelagic->q_0 = pel_params[23]; // modal ratio of predation sizes
    pelagic->sig = pel_params[24]; // predator-prey speciality
    pelagic->trunc = pel_params[25];

    /*Pelagic Species Movement Constants*/
    pelagic->prey = 0;
    pelagic->pred = 0;
    pelagic->comp = 0;
    pelagic->gamma_prey = 0;
    pelagic->gamma_pred = 0;
    pelagic->gamma_comp = 0;

    if (run->spatial_dim == 1 || run->spatial_dim == 2)
    {
        pelagic->prey = pel_params[26];
        pelagic->pred = pel_params[27];
        pelagic->comp = pel_params[28];
        pelagic->gamma_prey = pel_params[29];
        pelagic->gamma_pred = pel_params[30];
        pelagic->gamma_comp = pel_params[31];
    }

    /*Pelagic Species Flags and Methods*/
    pelagic->rep_method = rep_method;
    {
        char *e = getenv("PEL_IMM_FRAC");
        pelagic->imm_frac = e ? atof(e) : 0.0;
        pelagic->rec0 = 0.0;
        pelagic->rec0_frozen = 0;
    }
    pelagic->initial_flag = initial_flag;
    pelagic->ts_flag = ts_flag;
    pelagic->fishing_flag = fishing_flag;

    int pel_species_idx, n_pel_species;
    n_pel_species = run->no_pelagic;

    if (run->write_files)
    {
        /*Create species output filenames*/
        snprintf(pelagic->fname_r, sizeof(pelagic->fname_r), "%s/%s/%s.txt", run->filename, pelagic->filename,
                 "results");
        snprintf(pelagic->fname_g, sizeof(pelagic->fname_g), "%s/%s/%s.txt", run->filename, pelagic->filename,
                 "growth");
        snprintf(pelagic->fname_m, sizeof(pelagic->fname_m), "%s/%s/%s.txt", run->filename, pelagic->filename,
                 "mortality");
        pelagic->fname_pred = (char **)safe_malloc(n_pel_species * sizeof(char *), __LINE__, __FILE__);
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            pelagic->fname_pred[pel_species_idx] = (char *)safe_malloc(40 * sizeof(char), __LINE__, __FILE__);
            snprintf(pelagic->fname_pred[pel_species_idx], 40, "%s/%s/%s%d.txt", run->filename, pelagic->filename,
                     "predation", pel_species_idx);
        }
        snprintf(pelagic->fname_fish, sizeof(pelagic->fname_fish), "%s/%s/%s.txt", run->filename, pelagic->filename,
                 "fishing");
        snprintf(pelagic->fname_summ, sizeof(pelagic->fname_summ), "%s/%s/%s.txt", run->filename, pelagic->filename,
                 "summary");

        /*Open output files*/
        pelagic->fptr_r = safe_fopen(pelagic->fname_r, "w", __LINE__, __FILE__);
        pelagic->fptr_g = safe_fopen(pelagic->fname_g, "w", __LINE__, __FILE__);
        pelagic->fptr_m = safe_fopen(pelagic->fname_m, "w", __LINE__, __FILE__);
        pelagic->fptr_pred = (FILE **)safe_malloc(n_pel_species * sizeof(FILE *), __LINE__, __FILE__);
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            pelagic->fptr_pred[pel_species_idx] =
                safe_fopen(pelagic->fname_pred[pel_species_idx], "w", __LINE__, __FILE__);
        }
        pelagic->fptr_fish = safe_fopen(pelagic->fname_fish, "w", __LINE__, __FILE__);
        pelagic->fptr_summ = safe_fopen(pelagic->fname_summ, "w", __LINE__, __FILE__);

        /*Initialise output files*/
        fprintf(pelagic->fptr_r, "filetype:      %s\n", "results");
        fprintf(pelagic->fptr_r, "speciestype:   %s\n", "pelagic");
        fprintf(pelagic->fptr_r, "\n");
        print_mass_header(grid, pelagic->fptr_r);

        fprintf(pelagic->fptr_g, "filetype:      %s\n", "growth");
        fprintf(pelagic->fptr_g, "speciestype:   %s\n", "pelagic");
        fprintf(pelagic->fptr_g, "\n");
        print_mass_header(grid, pelagic->fptr_g);

        fprintf(pelagic->fptr_m, "filetype:      %s\n", "mortality");
        fprintf(pelagic->fptr_m, "speciestype:   %s\n", "pelagic");
        fprintf(pelagic->fptr_m, "\n");
        print_mass_header(grid, pelagic->fptr_m);

        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            fprintf(pelagic->fptr_pred[pel_species_idx], "filetype:      %s\n", "predation");
            fprintf(pelagic->fptr_pred[pel_species_idx], "speciestype:   %s\n", "pelagic");
            fprintf(pelagic->fptr_pred[pel_species_idx], "\n");
            print_mass_header(grid, pelagic->fptr_pred[pel_species_idx]);
        }

        fprintf(pelagic->fptr_fish, "filetype:      %s\n", "fishing");
        fprintf(pelagic->fptr_fish, "speciestype:   %s\n", "pelagic");
        fprintf(pelagic->fptr_fish, "\n");
        print_mass_header(grid, pelagic->fptr_fish);

        fprintf(pelagic->fptr_summ, "filetype:      %s\n", "summary");
        fprintf(pelagic->fptr_summ, "speciestype:   %s\n", "pelagic");
        fprintf(pelagic->fptr_summ, "\n");
        print_pelagic_header(pelagic->fptr_summ);

        /*Close output files*/
        fclose(pelagic->fptr_r);
        fclose(pelagic->fptr_g);
        fclose(pelagic->fptr_m);
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            fclose(pelagic->fptr_pred[pel_species_idx]);
        }
        fclose(pelagic->fptr_fish);
        fclose(pelagic->fptr_summ);
    }

    /*Setup Species Abundance Arrays*/
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    /*Setup 3D size-space arrays*/
    pelagic->mu_pred_values = (double ****)safe_malloc(n_pel_species * sizeof(double ***), __LINE__, __FILE__);
    pelagic->u_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    pelagic->g_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    pelagic->mu_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        pelagic->mu_pred_values[pel_species_idx] = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    }
    pelagic->mu_fish_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);

    pelagic->pla_bio = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    pelagic->pel_bio = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    pelagic->ben_bio = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    pelagic->a_refuge =
        (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__); /*reef refuge availability a(w)*/

    pelagic->phi_values = (double *)safe_malloc(n_m * sizeof(double), __LINE__, __FILE__);
    pelagic->refuge_values =
        (double *)safe_malloc(n_m * sizeof(double), __LINE__, __FILE__); /*per prey-size refuge density*/
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        pelagic->u_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        pelagic->g_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        pelagic->mu_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
        {
            pelagic->mu_pred_values[pel_species_idx][sc_idx] =
                (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        }
        pelagic->mu_fish_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

        pelagic->pla_bio[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        pelagic->pel_bio[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        pelagic->ben_bio[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        pelagic->a_refuge[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

        pelagic->phi_values[sc_idx] = 0;
        pelagic->refuge_values[sc_idx] = 0;
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            pelagic->u_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            pelagic->g_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            pelagic->mu_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
            {
                pelagic->mu_pred_values[pel_species_idx][sc_idx][x_idx] =
                    (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            }
            pelagic->mu_fish_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);

            pelagic->pla_bio[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            pelagic->pel_bio[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            pelagic->ben_bio[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            pelagic->a_refuge[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                pelagic->u_values[sc_idx][x_idx][y_idx] = 0;
                pelagic->g_values[sc_idx][x_idx][y_idx] = 0;
                pelagic->mu_values[sc_idx][x_idx][y_idx] = 0;
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    pelagic->mu_pred_values[pel_species_idx][sc_idx][x_idx][y_idx] = 0;
                }
                pelagic->mu_fish_values[sc_idx][x_idx][y_idx] = 0;

                pelagic->pla_bio[sc_idx][x_idx][y_idx] = 0;
                pelagic->pel_bio[sc_idx][x_idx][y_idx] = 0;
                pelagic->ben_bio[sc_idx][x_idx][y_idx] = 0;
                pelagic->a_refuge[sc_idx][x_idx][y_idx] = 1.0; /*default fully exposed (a=1) -> standard dbpmr*/
            }
        }
    }

    /*Setup 2D space arrays*/
    pelagic->pla_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->pel_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->ben_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

    pelagic->fish_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->pred_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->reproduction = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->rec0_c = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    pelagic->rec0_frozen_c = (int **)safe_malloc(n_x * sizeof(int *), __LINE__, __FILE__);
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        pelagic->pla_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->pel_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->ben_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);

        pelagic->fish_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->pred_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->reproduction[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->rec0_c[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        pelagic->rec0_frozen_c[x_idx] = (int *)safe_malloc(n_y * sizeof(int), __LINE__, __FILE__);
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            pelagic->pla_total[x_idx][y_idx] = 0;
            pelagic->pel_total[x_idx][y_idx] = 0;
            pelagic->ben_total[x_idx][y_idx] = 0;

            pelagic->fish_total[x_idx][y_idx] = 0;
            pelagic->pred_total[x_idx][y_idx] = 0;
            pelagic->reproduction[x_idx][y_idx] = 0;
            pelagic->rec0_c[x_idx][y_idx] = 0.0;
            pelagic->rec0_frozen_c[x_idx][y_idx] = 0;
        }
    }

    /*Setup initial values for u(m,0,x,y)*/
    if (pelagic->initial_flag == 1)
    {
        int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
        char *input;                   // input string to store line as string
        input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
        char sep[] = ",";    // identifying the ',' char as seperating values
        char *result = NULL; // declaring a pointer for storing each individual value as a char

        /*Create plankton input filename*/
        snprintf(pelagic->fname_ts, sizeof(pelagic->fname_ts), "%s/Input/%s_%s.txt", run->filename, pelagic->filename,
                 "ts");
        pelagic->fptr_ts = safe_fopen(pelagic->fname_ts, "r", __LINE__, __FILE__);

        /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                fgets(input, MAX_LEN, pelagic->fptr_ts);
                input[strlen(input) - 1] =
                    '\0'; // replacing the newline \n char of the line of text with a null \0 char

                result = (char *)strtok(input, sep); // this is the 1st value
                sc_idx = 0;
                while (result != NULL)
                {
                    pelagic->u_values[sc_idx][x_idx][y_idx] = atof(result);
                    result = (char *)strtok(NULL, sep);
                    sc_idx++;
                }
            }
        }
        free(input);
        /*If this is only an intial value then close the file now*/
        if (pelagic->ts_flag == 0)
        {
            fclose(pelagic->fptr_ts);
        }
    }
    else
    {
        for (sc_idx = pelagic->ipelmin; sc_idx < (pelagic->ipelmax + 1); sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    pelagic->u_values[sc_idx][x_idx][y_idx] =
                        pelagic->u_0 * exp(pelagic->lambda * grid->m_values[sc_idx]);
                    if (run->spatial_dim == 1 || run->spatial_dim == 2)
                    {
                        pelagic->u_values[sc_idx][x_idx][y_idx] *= x_start(grid->x_values[x_idx], 0);
                        if (run->spatial_dim == 2)
                        {
                            pelagic->u_values[sc_idx][x_idx][y_idx] *= y_start(grid->y_values[y_idx], 0);
                        }
                    }
                }
            }
        }
    }

    /*Open pelagic fishing file*/
    if (pelagic->fishing_flag == 1)
    {
        /*Create pelagic fishing input filename*/
        snprintf(pelagic->fname_fish_ts, sizeof(pelagic->fname_fish_ts), "%s/Input/%s_%s.txt", run->filename,
                 pelagic->filename, "fishing_ts");
        pelagic->fptr_fish_ts = safe_fopen(pelagic->fname_fish_ts, "r", __LINE__, __FILE__);
    }

    /*Setup pelagic reproduction if needed*/
    if (pelagic->ts_flag == 0)
    {
        /*set eggs equal to intial values*/
        if (pelagic->rep_method == 0)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    pelagic->reproduction[x_idx][y_idx] = pelagic->u_values[pelagic->ipelmin][x_idx][y_idx];
                }
            }
        }
        if (pelagic->rep_method == 1)
        {
            /*Create pelagic fishing input filename*/
            snprintf(pelagic->fname_rep_ts, sizeof(pelagic->fname_rep_ts), "%s/Input/%s_%s.txt", run->filename,
                     pelagic->filename, "rep_ts");
            pelagic->fptr_rep_ts = safe_fopen(pelagic->fname_rep_ts, "r", __LINE__, __FILE__);
        }
    }

    /*Setup values for predator-prey fedding kernel*/
    for (sc_idx = 0; sc_idx < (pelagic->ipelmax + 1); sc_idx++)
    {
        pelagic->phi_values[sc_idx] =
            phi((grid->m_values[sc_idx] - grid->m_values[0]), pelagic->q_0, pelagic->sig, pelagic->trunc);
    }

    /* --- Reef prey-refuge plugin (opt-in; Rogers, Blanchard & Mumby 2014) ---
       Availability a(w) = max(min_A, 1 - refuge(w)/competitor(w)) scales BOTH predation mortality on,
       and predator intake of, prey of body size w (mass balance). Reads an optional per-size refuge
       density from <run>/Input/<name>_refuge.txt (one value per grid size class). Absent file (or all
       zeros) -> refuge_on=0 -> a(w)=1 everywhere -> bit-identical to standard dbpmr. min_A (floor) from
       the REFUGE_MIN_A env var (default 0). */
    pelagic->refuge_on = 0;
    pelagic->refuge_min_a = 0.0;
    {
        char *e = getenv("REFUGE_MIN_A");
        if (e)
            pelagic->refuge_min_a = atof(e);
    }
    {
        char rfn[256];
        snprintf(rfn, sizeof(rfn), "%s/Input/%s_refuge.txt", run->filename, pelagic->filename);
        FILE *rf = fopen(rfn, "r");
        if (rf)
        {
            double v;
            int cnt = 0;
            for (sc_idx = 0; sc_idx < n_m && fscanf(rf, "%lf", &v) == 1; sc_idx++)
            {
                pelagic->refuge_values[sc_idx] = v;
                if (v > 0)
                    cnt++;
            }
            fclose(rf);
            if (cnt > 0)
                pelagic->refuge_on = 1; /*only activate when a positive refuge is supplied*/
        }
    }
}

void setup_benthic(RUN *run, GRID *grid, BENTHIC *benthic, double *ben_params, char *filename, int rep_method,
                   int initial_flag, int ts_flag, int fishing_flag)
{
    /*Declare species name*/
    benthic->filename = filename;

    /*Benthic Species Size Range*/
    benthic->benmin = ben_params[0];
    benthic->benmat = ben_params[1];
    benthic->benmax = ben_params[2];
    benthic->ibenmin = (int)rint((benthic->benmin - grid->mmin) / grid->mstep);
    benthic->ibenmat = (int)rint((benthic->benmat - grid->mmin) / grid->mstep);
    benthic->ibenmax = (int)rint((benthic->benmax - grid->mmin) / grid->mstep);

    /*Benthic Species Constants*/
    benthic->A = ben_params[3];         // Volume searched by unit weight
    benthic->alpha = ben_params[4];     // Volume search allometry
    benthic->mu_0 = ben_params[5];      // natural mortality rate
    benthic->beta = ben_params[6];      // natural mortality allometry
    benthic->A_base = benthic->A;       // temp driver: base (unscaled) search rate
    benthic->mu_0_base = benthic->mu_0; // temp driver: base (unscaled) background mortality
    benthic->tempeff = 1.0;             // temp driver: neutral until a temperature series is read
    benthic->mu_s = ben_params[7];      // senescence mortality rate
    benthic->epsilon = ben_params[8];   // senescece mortality constant
    benthic->u_0 = ben_params[9];       // initial intercept
    benthic->lambda = ben_params[10];   // primary production slope

    /*Benthic Species Conversion Efficiencies*/
    benthic->K_det = ben_params[11]; // growth efficiency
    benthic->R_det = ben_params[12];
    benthic->Ex_det = ben_params[13];

    benthic->pref_det = ben_params[14];

    /*Benthic Species Flags and Methods*/
    benthic->rep_method = rep_method;
    benthic->initial_flag = initial_flag;
    benthic->ts_flag = ts_flag;
    benthic->fishing_flag = fishing_flag;

    if (run->write_files)
    {
        /*Create species output filenames*/
        snprintf(benthic->fname_r, sizeof(benthic->fname_r), "%s/%s/%s.txt", run->filename, benthic->filename,
                 "results");
        snprintf(benthic->fname_g, sizeof(benthic->fname_g), "%s/%s/%s.txt", run->filename, benthic->filename,
                 "growth");
        snprintf(benthic->fname_m, sizeof(benthic->fname_m), "%s/%s/%s.txt", run->filename, benthic->filename,
                 "mortality");
        snprintf(benthic->fname_pred, sizeof(benthic->fname_pred), "%s/%s/%s.txt", run->filename, benthic->filename,
                 "predation");
        snprintf(benthic->fname_fish, sizeof(benthic->fname_fish), "%s/%s/%s.txt", run->filename, benthic->filename,
                 "fishing");
        snprintf(benthic->fname_summ, sizeof(benthic->fname_summ), "%s/%s/%s.txt", run->filename, benthic->filename,
                 "summary");

        /*Open output files*/
        benthic->fptr_r = safe_fopen(benthic->fname_r, "w", __LINE__, __FILE__);
        benthic->fptr_g = safe_fopen(benthic->fname_g, "w", __LINE__, __FILE__);
        benthic->fptr_m = safe_fopen(benthic->fname_m, "w", __LINE__, __FILE__);
        benthic->fptr_pred = safe_fopen(benthic->fname_pred, "w", __LINE__, __FILE__);
        benthic->fptr_fish = safe_fopen(benthic->fname_fish, "w", __LINE__, __FILE__);
        benthic->fptr_summ = safe_fopen(benthic->fname_summ, "w", __LINE__, __FILE__);

        /*Initialise output files*/
        fprintf(benthic->fptr_r, "filetype:      %s\n", "results");
        fprintf(benthic->fptr_r, "speciestype:   %s\n", "benthic");
        fprintf(benthic->fptr_r, "\n");
        print_mass_header(grid, benthic->fptr_r);

        fprintf(benthic->fptr_g, "filetype:      %s\n", "growth");
        fprintf(benthic->fptr_g, "speciestype:   %s\n", "benthic");
        fprintf(benthic->fptr_g, "\n");
        print_mass_header(grid, benthic->fptr_g);

        fprintf(benthic->fptr_m, "filetype:      %s\n", "mortality");
        fprintf(benthic->fptr_m, "speciestype:   %s\n", "benthic");
        fprintf(benthic->fptr_m, "\n");
        print_mass_header(grid, benthic->fptr_m);

        fprintf(benthic->fptr_pred, "filetype:      %s\n", "predation");
        fprintf(benthic->fptr_pred, "speciestype:   %s\n", "benthic");
        fprintf(benthic->fptr_pred, "\n");
        print_mass_header(grid, benthic->fptr_pred);

        fprintf(benthic->fptr_fish, "filetype:      %s\n", "fishing");
        fprintf(benthic->fptr_fish, "speciestype:   %s\n", "benthic");
        fprintf(benthic->fptr_fish, "\n");
        print_mass_header(grid, benthic->fptr_fish);

        fprintf(benthic->fptr_summ, "filetype:      %s\n", "summary");
        fprintf(benthic->fptr_summ, "speciestype:   %s\n", "benthic");
        fprintf(benthic->fptr_summ, "\n");
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
    int sc_idx, x_idx, y_idx;
    int n_m, n_x, n_y;
    n_m = grid->mnum;
    n_x = grid->xnum;
    n_y = grid->ynum;

    benthic->u_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    benthic->g_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    benthic->mu_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    benthic->mu_pred_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    benthic->mu_fish_values = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);

    benthic->det_bio = (double ***)safe_malloc(n_m * sizeof(double **), __LINE__, __FILE__);
    for (sc_idx = 0; sc_idx < n_m; sc_idx++)
    {
        benthic->u_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        benthic->g_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        benthic->mu_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        benthic->mu_pred_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        benthic->mu_fish_values[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

        benthic->det_bio[sc_idx] = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            benthic->u_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            benthic->g_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            benthic->mu_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            benthic->mu_pred_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            benthic->mu_fish_values[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);

            benthic->det_bio[sc_idx][x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                benthic->u_values[sc_idx][x_idx][y_idx] = 0;
                benthic->g_values[sc_idx][x_idx][y_idx] = 0;
                benthic->mu_values[sc_idx][x_idx][y_idx] = 0;
                benthic->mu_pred_values[sc_idx][x_idx][y_idx] = 0;
                benthic->mu_fish_values[sc_idx][x_idx][y_idx] = 0;

                benthic->det_bio[sc_idx][x_idx][y_idx] = 0;
            }
        }
    }
    benthic->det_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);

    benthic->fish_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    benthic->pred_total = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    benthic->reproduction = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        benthic->det_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);

        benthic->fish_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        benthic->pred_total[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        benthic->reproduction[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            benthic->det_total[x_idx][y_idx] = 0;

            benthic->fish_total[x_idx][y_idx] = 0;
            benthic->pred_total[x_idx][y_idx] = 0;
            benthic->reproduction[x_idx][y_idx] = 0;
        }
    }

    /*Setup initial values for u(m,0,x,y)*/
    if (benthic->initial_flag == 1)
    {
        int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
        char *input;                   // input string to store line as string
        input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
        char sep[] = ",";    // identifying the ',' char as seperating values
        char *result = NULL; // declaring a pointer for storing each individual value as a char

        /*Create plankton input filename*/
        snprintf(benthic->fname_ts, sizeof(benthic->fname_ts), "%s/Input/%s_%s.txt", run->filename, benthic->filename,
                 "ts");
        benthic->fptr_ts = safe_fopen(benthic->fname_ts, "r", __LINE__, __FILE__);

        /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                fgets(input, MAX_LEN, benthic->fptr_ts);
                input[strlen(input) - 1] =
                    '\0'; // replacing the newline \n char of the line of text with a null \0 char

                result = (char *)strtok(input, sep); // this is the 1st value
                sc_idx = 0;
                while (result != NULL)
                {
                    benthic->u_values[sc_idx][x_idx][y_idx] = atof(result);
                    result = (char *)strtok(NULL, sep);
                    sc_idx++;
                }
            }
        }
        free(input);
        /*If this is only an intial value then close the file now*/
        if (benthic->ts_flag == 0)
        {
            fclose(benthic->fptr_ts);
        }
    }
    else
    {
        for (sc_idx = benthic->ibenmin; sc_idx < (benthic->ibenmax + 1); sc_idx++)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    benthic->u_values[sc_idx][x_idx][y_idx] =
                        benthic->u_0 * exp(benthic->lambda * grid->m_values[sc_idx]);
                    if (run->spatial_dim == 1 || run->spatial_dim == 2)
                    {
                        benthic->u_values[sc_idx][x_idx][y_idx] *= x_start(grid->x_values[x_idx], 0);
                        if (run->spatial_dim == 2)
                        {
                            benthic->u_values[sc_idx][x_idx][y_idx] *= y_start(grid->y_values[y_idx], 0);
                        }
                    }
                }
            }
        }
    }

    /*Setup initial values for fishing mortality*/
    if (benthic->fishing_flag == 1)
    {
        /*Create benthic fishing input filename*/
        snprintf(benthic->fname_fish_ts, sizeof(benthic->fname_fish_ts), "%s/Input/%s_%s.txt", run->filename,
                 benthic->filename, "fishing_ts");
        benthic->fptr_fish_ts = safe_fopen(benthic->fname_fish_ts, "r", __LINE__, __FILE__);
    }

    /*Setup benthic reproduction if needed*/
    if (benthic->ts_flag == 0)
    {
        /*set eggs equal to intial values*/
        if (benthic->rep_method == 0)
        {
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    benthic->reproduction[x_idx][y_idx] = benthic->u_values[benthic->ibenmin][x_idx][y_idx];
                }
            }
        }
        if (benthic->rep_method == 1)
        {
            /*Create benthic fishing input filename*/
            snprintf(benthic->fname_rep_ts, sizeof(benthic->fname_rep_ts), "%s/Input/%s_%s.txt", run->filename,
                     benthic->filename, "rep_ts");
            benthic->fptr_rep_ts = safe_fopen(benthic->fname_rep_ts, "r", __LINE__, __FILE__);
        }
    }
}

void setup_detritus(RUN *run, GRID *grid, DETRITUS *detritus, double *det_params, char *filename, int initial_flag,
                    int ts_flag)
{

    /*Setup Detritus Constants*/
    detritus->w_0 = det_params[0];

    /*Declare detritus name*/
    detritus->filename = filename;

    /*Setup Plankton Flags*/
    detritus->initial_flag = initial_flag;
    detritus->ts_flag = ts_flag;

    if (run->write_files)
    {
        /*Create detritus output filenames*/
        snprintf(detritus->fname_summ, sizeof(detritus->fname_summ), "%s/%s/%s.txt", run->filename, detritus->filename,
                 "summary");

        /*Open output files*/
        detritus->fptr_summ = safe_fopen(detritus->fname_summ, "w", __LINE__, __FILE__);

        fprintf(detritus->fptr_summ, "filetype:      %s\n", "summary");
        fprintf(detritus->fptr_summ, "speciestype:   %s\n", "detritus");
        fprintf(detritus->fptr_summ, "\n");
        print_detritus_header(detritus->fptr_summ);

        /*Close output files*/
        fclose(detritus->fptr_summ);
    }

    /*Setup And Initialise detritus Abundance Arrays*/
    int x_idx, y_idx;
    int n_x, n_y;
    n_x = grid->xnum;
    n_y = grid->ynum;

    detritus->w_values = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    detritus->g_values = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    detritus->mu_values = (double **)safe_malloc(n_x * sizeof(double *), __LINE__, __FILE__);
    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        detritus->w_values[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        detritus->g_values[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        detritus->mu_values[x_idx] = (double *)safe_malloc(n_y * sizeof(double), __LINE__, __FILE__);
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            detritus->w_values[x_idx][y_idx] = 0;
            detritus->g_values[x_idx][y_idx] = 0;
            detritus->mu_values[x_idx][y_idx] = 0;
        }
    }

    /*Setup initial values for w(0,x,y)*/
    /*Use user defined values*/
    if (detritus->initial_flag == 1)
    {
        int MAX_LEN = 15; // this is the upper bound for the row length in characters
        char *input;      // input string to store line as string
        input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
        char sep[] = ",";    // identifying the ',' char as seperating values
        char *result = NULL; // declaring a pointer for storing each individual value as a char

        /*Create detritus input filename*/
        snprintf(detritus->fname_ts, sizeof(detritus->fname_ts), "%s/Input/%s_%s.txt", run->filename,
                 detritus->filename, "ts");
        detritus->fptr_ts = safe_fopen(detritus->fname_ts, "r", __LINE__, __FILE__);

        /*Read in each line one at a time (corresponding to spectra at a spatial point)*/
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                fgets(input, MAX_LEN, detritus->fptr_ts);
                input[strlen(input) - 1] =
                    '\0'; // replacing the newline \n char of the line of text with a null \0 char

                result = (char *)strtok(input, sep); // this is the 1st value
                while (result != NULL)
                {
                    detritus->w_values[x_idx][y_idx] = atof(result);
                    result = (char *)strtok(NULL, sep);
                }
            }
        }
        free(input);
        /*If this is only an intial value then close the file now*/
        if (detritus->ts_flag == 0)
        {
            fclose(detritus->fptr_ts);
        }
    }
    /*Use standard values*/
    else
    {
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                detritus->w_values[x_idx][y_idx] = detritus->w_0;
            }
        }
    }
}

double **setup_pelagic_params(int no_pelagic, double *pel_params)
{
    int param_idx, pel_species_idx;
    int no_params = 32;
    double **temp_params;
    temp_params = (double **)safe_malloc(no_pelagic * sizeof(double *), __LINE__, __FILE__);
    for (pel_species_idx = 0; pel_species_idx < no_pelagic; pel_species_idx++)
    {
        temp_params[pel_species_idx] = (double *)safe_malloc(no_params * sizeof(double), __LINE__, __FILE__);
        for (param_idx = 0; param_idx < no_params; param_idx++)
        {
            temp_params[pel_species_idx][param_idx] = pel_params[no_params * pel_species_idx + param_idx];
        }
    }
    return (temp_params);
}

double **setup_benthic_params(int no_benthic, double *ben_params)
{
    int param_idx, ben_species_idx;
    int no_params = 15;
    double **temp_params;
    temp_params = (double **)safe_malloc(no_benthic * sizeof(double *), __LINE__, __FILE__);
    for (ben_species_idx = 0; ben_species_idx < no_benthic; ben_species_idx++)
    {
        temp_params[ben_species_idx] = (double *)safe_malloc(no_params * sizeof(double), __LINE__, __FILE__);
        for (param_idx = 0; param_idx < no_params; param_idx++)
        {
            temp_params[ben_species_idx][param_idx] = ben_params[no_params * ben_species_idx + param_idx];
        }
    }
    return (temp_params);
}

void setup_matrix(int dim, MATRIX *matrix)
/*Creates a matrix for the consumer size species*/
{
    matrix->size = dim;
    matrix->a = (double *)safe_malloc(dim * sizeof(double), __LINE__, __FILE__);
    matrix->b = (double *)safe_malloc(dim * sizeof(double), __LINE__, __FILE__);
    matrix->c = (double *)safe_malloc(dim * sizeof(double), __LINE__, __FILE__);
    matrix->r = (double *)safe_malloc(dim * sizeof(double), __LINE__, __FILE__);
    matrix->u = (double *)safe_malloc(dim * sizeof(double), __LINE__, __FILE__);
}


// Set land cells to zero after setup. If land cells are non-zero after reading
// in input values etc. this could create confusion. So, after setting up the community,
// we set land cell values to 0 for clarity.
void set_land_to_zero(RUN *run, GRID *grid, COMMUNITY *comm)
{
    int sc_idx, x_idx, y_idx, pel_species_idx, ben_species_idx;
    int n_m = grid->mnum, n_x = grid->xnum, n_y = grid->ynum;
    int n_pel_species = run->no_pelagic, n_ben_species = run->no_benthic;

    for (x_idx = 0; x_idx < n_x; x_idx++)
    {
        for (y_idx = 0; y_idx < n_y; y_idx++)
        {
            if (!is_land(grid, x_idx, y_idx)) {continue;} // Only zero land cells

            for (sc_idx = 0; sc_idx < n_m; sc_idx++)
            {
                comm->plankton->u_values[sc_idx][x_idx][y_idx] = 0;

                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    comm->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] = 0;
                }

                if (n_ben_species != 0)
                {
                    for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
                    {
                        comm->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] = 0;
                    }
                }
            }

            if (n_ben_species != 0)
            {
                comm->detritus->w_values[x_idx][y_idx] = 0;
            }
        }
    }
}

double x_start(double x, double t)
/*function determining initial plankton distribution*/
{
    double ans;
    ans = 2 * exp(-((x - 500) * (x - 500) / 20000));
    return ans;
}

double y_start(double y, double t)
/*function determining initial plankton distribution*/
{
    double ans;
    ans = 2 * exp(-((y - 500) * (y - 500) / 20000));
    return ans;
}