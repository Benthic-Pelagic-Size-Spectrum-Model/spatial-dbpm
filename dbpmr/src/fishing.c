#include "fishing.h"

void calculate_fishing(RUN *run, GRID *grid, COMMUNITY *community)
{
    int sc_idx, x_idx, y_idx, pel_species_idx, ben_species_idx;
    int n_x, n_y, n_pel_species, n_ben_species;

    n_x = grid->xnum;
    n_y = grid->ynum;
    n_pel_species = run->no_pelagic;
    n_ben_species = run->no_benthic;

    /*Pelagic Fishing*/
    for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
    {
        if (community->pelagic[pel_species_idx].fishing_flag == 1)
        {
            int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
            char *input;                   // input string to store line as string
            input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
            char sep[] = ",";    // identifying the ',' char as seperating values
            char *result = NULL; // declaring a pointer for storing each individual value as a char

            /*Read in each line one at a time (corresponding to fishing mortality at a spatial point)*/
            for (x_idx = 0; x_idx < n_x; x_idx++)
            {
                for (y_idx = 0; y_idx < n_y; y_idx++)
                {
                    fgets(input, MAX_LEN, community->pelagic[pel_species_idx].fptr_fish_ts);
                    input[strlen(input) - 1] =
                        '\0'; // replacing the newline \n char of the line of text with a null \0 char

                    result = (char *)strtok(input, sep); // this is the 1st value
                    sc_idx = 0;
                    SKIP_IF_LAND(grid, x_idx, y_idx);

                    while (result != NULL)
                    {
                        community->pelagic[pel_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] = atof(result);
                        result = (char *)strtok(NULL, sep);
                        sc_idx++;
                    }
                }
            }
            free(input);
        }
    }

    /*Benthic Fishing*/
    if (run->no_benthic != 0)
    {
        for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
        {
            if (community->benthic[ben_species_idx].fishing_flag == 1)
            {
                int MAX_LEN = grid->mnum * 15; // this is the upper bound for the row length in characters
                char *input;                   // input string to store line as string
                input = (char *)safe_malloc(MAX_LEN * sizeof(char), __LINE__, __FILE__);
                char sep[] = ",";    // identifying the ',' char as seperating values
                char *result = NULL; // declaring a pointer for storing each individual value as a char

                /*Read in each line one at a time (corresponding to fishing mortality at a spatial point)*/
                for (x_idx = 0; x_idx < n_x; x_idx++)
                {
                    for (y_idx = 0; y_idx < n_y; y_idx++)
                    {
                        fgets(input, MAX_LEN, community->benthic[ben_species_idx].fptr_fish_ts);
                        input[strlen(input) - 1] =
                            '\0'; // replacing the newline \n char of the line of text with a null \0 char

                        result = (char *)strtok(input, sep); // this is the 1st value
                        sc_idx = 0;
                        SKIP_IF_LAND(grid, x_idx, y_idx);

                        while (result != NULL)
                        {
                            community->benthic[ben_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] = atof(result);
                            result = (char *)strtok(NULL, sep);
                            sc_idx++;
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
    if (run->col_gravity >= 0)
    {
        double bhd = 20.0; /* benthic habitat depth (sizemodel default) */
        for (x_idx = 0; x_idx < n_x; x_idx++)
        {
            for (y_idx = 0; y_idx < n_y; y_idx++)
            {
                SKIP_IF_LAND(grid, x_idx, y_idx);

                /* pelagic column thickness (depth channel; 1 if absent) - per-cell under FORCING_PERCELL */
                double dcorr = run->forcing_percell ? run->depth_c[x_idx][y_idx] : run->depth;
                double Bpel = 0.0, Bben = 0.0;
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                {
                    if (community->pelagic[pel_species_idx].fishing_flag != 1)
                        continue;
                    for (sc_idx = community->pelagic[pel_species_idx].ipelmin;
                         sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
                        if (community->pelagic[pel_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] > 0.0)
                            Bpel += community->pelagic[pel_species_idx].u_values[sc_idx][x_idx][y_idx] *
                                    grid->exp_m_values[sc_idx] * grid->mstep;
                }
                for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
                {
                    if (community->benthic[ben_species_idx].fishing_flag != 1)
                        continue;
                    for (sc_idx = community->benthic[ben_species_idx].ibenmin;
                         sc_idx < (community->benthic[ben_species_idx].ibenmax + 1); sc_idx++)
                        if (community->benthic[ben_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] > 0.0)
                            Bben += community->benthic[ben_species_idx].u_values[sc_idx][x_idx][y_idx] *
                                    grid->exp_m_values[sc_idx] * grid->mstep;
                }
                double Ap = Bpel * dcorr, Ab = Bben * bhd, tot = Ap + Ab;
                double sp = (tot > 0.0) ? Ap / tot : 0.5;
                double sb = (tot > 0.0) ? Ab / tot : 0.5;
                for (pel_species_idx = 0; pel_species_idx < n_pel_species; pel_species_idx++)
                    if (community->pelagic[pel_species_idx].fishing_flag == 1)
                        for (sc_idx = community->pelagic[pel_species_idx].ipelmin;
                             sc_idx < (community->pelagic[pel_species_idx].ipelmax + 1); sc_idx++)
                            community->pelagic[pel_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] *= sp;
                for (ben_species_idx = 0; ben_species_idx < n_ben_species; ben_species_idx++)
                    if (community->benthic[ben_species_idx].fishing_flag == 1)
                        for (sc_idx = community->benthic[ben_species_idx].ibenmin;
                             sc_idx < (community->benthic[ben_species_idx].ibenmax + 1); sc_idx++)
                            community->benthic[ben_species_idx].mu_fish_values[sc_idx][x_idx][y_idx] *= sb;
            }
        }
    }
}