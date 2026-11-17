#include "SizeSpectra.h"


/*Fishing mortality arrays for the nsp species of one type, from element `name` of the inputs list: a list
  with one entry per species, each NULL (fishing read from file, or not fished) or a list holding rate
  (numeric [x, y, layer]), layer (0-based layer at each time point), nlayer and range (lo, hi).*/
static FISHING_ARRAY *read_fishing_arrays(SEXP inputs, const char *name, int nsp);

/*Element of an R list by name (R_NilValue if absent)*/
static SEXP get_list_element(SEXP list, const char *name)
{
        SEXP names = getAttrib(list, R_NamesSymbol);
        int i;
        if (names == R_NilValue) return R_NilValue;
        for (i = 0; i < LENGTH(list); i++) {
                if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) return VECTOR_ELT(list, i);
        }
        return R_NilValue;
}

static FISHING_ARRAY *read_fishing_arrays(SEXP inputs, const char *name, int nsp)
{
        SEXP lst = get_list_element(inputs, name);
        int s, j;
        if (lst == R_NilValue || nsp < 1) return NULL;
        if (!isNewList(lst) || LENGTH(lst) != nsp) error("%s must be a list with one entry per species", name);
        FISHING_ARRAY *fa = (FISHING_ARRAY *) R_alloc(nsp, sizeof(FISHING_ARRAY));
        for (s = 0; s < nsp; s++) {
                SEXP el = VECTOR_ELT(lst, s);
                fa[s].rate = NULL; fa[s].layer = NULL; fa[s].nlayer = 0; fa[s].nblock = 0; fa[s].len = 0;
                fa[s].lo = R_NegInf; fa[s].hi = R_PosInf;
                if (el == R_NilValue) continue;
                SEXP rate = get_list_element(el, "rate"), layer = get_list_element(el, "layer");
                SEXP nl = get_list_element(el, "nlayer"), range = get_list_element(el, "range");
                if (!isReal(rate) || !isInteger(layer) || !isInteger(nl) || LENGTH(nl) != 1 || !isReal(range) || LENGTH(range) != 2) {
                        error("%s[[%d]]: rate and range must be numeric, layer and nlayer integer", name, s + 1);
                }
                int nlayer = INTEGER(nl)[0];
                if (nlayer < 1 || LENGTH(layer) < 1 || XLENGTH(rate) < nlayer) error("%s[[%d]]: inconsistent lengths", name, s + 1);
                for (j = 0; j < LENGTH(layer); j++) {
                        if (INTEGER(layer)[j] < 0 || INTEGER(layer)[j] >= nlayer) error("%s[[%d]]: layer index out of range at time point %d", name, s + 1, j + 1);
                }
                fa[s].rate = REAL(rate); fa[s].layer = INTEGER(layer); fa[s].nlayer = nlayer; fa[s].nblock = LENGTH(layer);
                fa[s].len = XLENGTH(rate); fa[s].lo = REAL(range)[0]; fa[s].hi = REAL(range)[1];
        }
        return fa;
}

// This is the new R interface. SizeSpectrum() is now a wrapper function that
// takes the objects from R (.Call()). SizeSpectrum() checks the SEXP objects from
// R are the correct SXPtypes, sets C pointers to the R data, and then runs the
// internal function SizeSpectrum_internals() which performs the calculations and
// writes outputs. SizeSpectrum() returns NULL to R.
// TODO: return to R will be changed to be optional so that users can toggle
// file output, memory output or both. 
SEXP SizeSpectrum(
        SEXP run_params_r, 
        SEXP grid_params_r, 
        SEXP plankton_params_r, 
        SEXP pelagic_params_r, 
        SEXP benthic_params_r, 
        SEXP detritus_params_r, 
        SEXP names_params_r, 
        SEXP flags_params_r,
        SEXP inputs_r,
        SEXP output_r
) {
        // Ensure that all of the SEXP objects are of the correct type after
        // passing in from R.
        if (!isInteger(run_params_r) || !isInteger(flags_params_r)) {
                error("run_params_r and flags_params_r R objects must be integer vectors");
        }
        if (
                !isReal(grid_params_r) || !isReal(plankton_params_r) ||
                !isReal(pelagic_params_r) || !isReal(benthic_params_r) ||
                !isReal(detritus_params_r)
        ) {error("grid, plankton, pelagic, benthic and detritus params must be numeric vectors");}
        if (!isString(names_params_r)) {
                error("names_params_r R object must be character");
        }

        // Checking the length of the parameter vectors that were passed from R.
        // The values here are taken from hardcoded indices throughout the code.
        if (LENGTH(run_params_r) != 5) {
                error("run_params_r must have 5 elements");
        }
        if (LENGTH(grid_params_r) != 18) {
                error("grid_params_r must have 18 elements");
        }
        if (LENGTH(plankton_params_r) != 6) {
                error("plankton_params_r must have 6 elements");
        }

        int npel = INTEGER(run_params_r)[0];
        if (LENGTH(pelagic_params_r) != 32 * npel) {
                error("pelagic_params_r must have 32 * n_pel_species elements");
        }
        
        // It is possible to include no benthic species, in which case benthic_params_r
        // and detritus_params_r do not need to be checked for length.
        int nben = INTEGER(run_params_r)[1];
        if (nben > 0) {
                if (LENGTH(benthic_params_r) != 15 * nben) {
                        error("benthic_params_r must have 15 * n_ben_species elements");
                }
                if (LENGTH(detritus_params_r) != 1) {
                        error("detritus_params_r must have 1 element");
                }
                if (LENGTH(names_params_r) != (npel + nben + 1) + 2) {
                        error("names_params_r not correct length");
                }
        } else if (LENGTH(names_params_r) != npel + 2) {
                error("names_params_r not correct length");
        }
        
        int nflags = 2 + 4 * (npel + nben) + (nben > 0 ? 2 : 0);
        // nflags is 2 for plankton, 4 for each pelagic / benthic species and 2
        // for detritus if there is at least 1 benthic species
        if (LENGTH(flags_params_r) != nflags) {
                error("flags_params_r not correct length");
        }

        // Assign pointers for the SEXP data to correct names for function use
        int *run_params = INTEGER(run_params_r);
        int *flags_params = INTEGER(flags_params_r);
        double *grid_params = REAL(grid_params_r);
        double *pla_params = REAL(plankton_params_r);
        double *pel_params = REAL(pelagic_params_r);
        double *ben_params = REAL(benthic_params_r);
        double *det_params = REAL(detritus_params_r);

        char **names_params = (char **) R_alloc(
                LENGTH(names_params_r), 
                sizeof(char *)
        );
        int i;
        for (i = 0; i < LENGTH(names_params_r); i++) {
                names_params[i] = (char *) CHAR(STRING_ELT(names_params_r, i));
        }
        
        // Output: c(mode, growth, mortality, predation, fishing, catch). mode: 0 = text files only (returns
        // NULL), 1 = text files and arrays, 2 = arrays only. The other five are 0/1 flags for the optional
        // variables returned alongside abundance.
        if (!isInteger(output_r) || LENGTH(output_r) != 6 || INTEGER(output_r)[0] < 0 || INTEGER(output_r)[0] > 2) {
                error("output_r must be 6 integers: mode (0 files, 1 files and arrays, 2 arrays) and five 0/1 variable flags");
        }
        OUTPUT_ARRAYS out;
        out.mode = INTEGER(output_r)[0];
        out.want_g = INTEGER(output_r)[1] != 0;
        out.want_mu = INTEGER(output_r)[2] != 0;
        out.want_pred = INTEGER(output_r)[3] != 0;
        out.want_fish = INTEGER(output_r)[4] != 0;
        out.want_catch = INTEGER(output_r)[5] != 0;
        out.list = R_NilValue;

        // Optional inputs passed as R arrays instead of text files (inputs_r is NULL or a named list).
        // Currently: the plankton series as intercept/slope arrays [x, y, layer] plus, for each time
        // block, the (0-based) layer it uses. The arrays stay owned by R for the duration of the call.
        ARRAY_INPUTS arr;
        arr.pla_int = NULL; arr.pla_slp = NULL; arr.pla_layer = NULL;
        arr.pla_nlayer = 0; arr.pla_nblock = 0; arr.pla_len = 0;
        arr.have_forcing = 0; arr.frc_layer = NULL; arr.frc_nblock = 0;
        arr.fish_pel = NULL; arr.fish_ben = NULL;
        arr.gravity = NULL;
        for (i = 0; i < N_FRC; i++) { arr.frc[i] = NULL; arr.frc_nlayer[i] = 0; arr.frc_len[i] = 0; }
        if (inputs_r != R_NilValue) {
                if (!isNewList(inputs_r)) {
                        error("inputs_r must be NULL or a list");
                }
                SEXP pint = get_list_element(inputs_r, "plankton_intercept");
                SEXP pslp = get_list_element(inputs_r, "plankton_slope");
                SEXP play = get_list_element(inputs_r, "plankton_layer");
                SEXP pnl  = get_list_element(inputs_r, "plankton_nlayer");
                if (pint != R_NilValue || pslp != R_NilValue || play != R_NilValue || pnl != R_NilValue) {
                        if (!isReal(pint) || !isReal(pslp) || !isInteger(play) || !isInteger(pnl) || LENGTH(pnl) != 1) {
                                error("plankton arrays: intercept and slope must be numeric, layer and nlayer integer");
                        }
                        int nlayer = INTEGER(pnl)[0];
                        if (nlayer < 1 || LENGTH(play) < 1 || XLENGTH(pint) != XLENGTH(pslp) || XLENGTH(pint) < nlayer) {
                                error("plankton arrays: inconsistent lengths");
                        }
                        for (i = 0; i < LENGTH(play); i++) {
                                if (INTEGER(play)[i] < 0 || INTEGER(play)[i] >= nlayer) {
                                        error("plankton arrays: layer index out of range at time block %d", i + 1);
                                }
                        }
                        arr.pla_int = REAL(pint); arr.pla_slp = REAL(pslp); arr.pla_layer = INTEGER(play);
                        arr.pla_nlayer = nlayer; arr.pla_nblock = LENGTH(play); arr.pla_len = XLENGTH(pint);
                }

                // Forcing channels as arrays: forcing_layer (0-based layer at each time point),
                // forcing_nlayer (layers of each channel, in FRC_NAMES order) and one array per channel present.
                SEXP flay = get_list_element(inputs_r, "forcing_layer");
                SEXP fnl  = get_list_element(inputs_r, "forcing_nlayer");
                if (flay != R_NilValue || fnl != R_NilValue) {
                        if (!isInteger(flay) || !isInteger(fnl) || LENGTH(fnl) != N_FRC || LENGTH(flay) < 1) {
                                error("forcing arrays: forcing_layer and forcing_nlayer must be integer, the latter of length %d", N_FRC);
                        }
                        int maxlayer = 1, nch = 0, c;
                        for (c = 0; c < N_FRC; c++) {
                                SEXP ch = get_list_element(inputs_r, FRC_NAMES[c]);
                                if (ch == R_NilValue) continue;
                                int nl = INTEGER(fnl)[c];
                                if (!isReal(ch) || nl < 1 || XLENGTH(ch) < nl) {
                                        error("forcing arrays: %s must be numeric with at least one layer", FRC_NAMES[c]);
                                }
                                arr.frc[c] = REAL(ch); arr.frc_nlayer[c] = nl; arr.frc_len[c] = XLENGTH(ch);
                                if (nl > maxlayer) maxlayer = nl;
                                nch++;
                        }
                        for (c = 0; c < N_FRC; c++) {
                                if (arr.frc[c] != NULL && arr.frc_nlayer[c] != 1 && arr.frc_nlayer[c] != maxlayer) {
                                        error("forcing arrays: channels must have 1 layer or the same number of layers");
                                }
                        }
                        for (i = 0; i < LENGTH(flay); i++) {
                                if (INTEGER(flay)[i] < 0 || INTEGER(flay)[i] >= maxlayer) {
                                        error("forcing arrays: layer index out of range at time point %d", i + 1);
                                }
                        }
                        if (nch > 0) {
                                arr.have_forcing = 1; arr.frc_layer = INTEGER(flay); arr.frc_nblock = LENGTH(flay);
                        }
                }

                // Fishing mortality as arrays, per species.
                arr.fish_pel = read_fishing_arrays(inputs_r, "fishing_pelagic", npel);
                arr.fish_ben = read_fishing_arrays(inputs_r, "fishing_benthic", nben);

                // Effort allocated among cells by gravity.
                SEXP gv_r = get_list_element(inputs_r, "gravity");
                if (gv_r != R_NilValue) {
                        SEXP ga = get_list_element(gv_r, "access"), gp = get_list_element(gv_r, "thick_pel");
                        SEXP gb = get_list_element(gv_r, "thick_ben"), gn = get_list_element(gv_r, "nupdate");
                        if (!isReal(ga) || !isReal(gp) || !isReal(gb) || !isInteger(gn) || LENGTH(gn) != 1 ||
                            XLENGTH(gp) != XLENGTH(ga) || XLENGTH(gb) != XLENGTH(ga) || INTEGER(gn)[0] < 1) {
                                error("gravity: access, thick_pel and thick_ben must be numeric of equal length and nupdate a positive integer");
                        }
                        GRAVITY_INPUT *gv = (GRAVITY_INPUT *) R_alloc(1, sizeof(GRAVITY_INPUT));
                        gv->access = REAL(ga); gv->thick_pel = REAL(gp); gv->thick_ben = REAL(gb);
                        gv->len = XLENGTH(ga); gv->nupdate = INTEGER(gn)[0];
                        arr.gravity = gv;
                }
        }

        SizeSpectrum_internals(
                run_params,
                grid_params,
                pla_params,
                pel_params,
                ben_params,
                det_params,
                names_params,
                flags_params,
                &arr,
                &out
        );
        if (out.mode != 0) {
                UNPROTECT(1);   // the list protected in setup_output_arrays()
                return(out.list);
        }
        return(R_NilValue);
}

//-----------------------------------------------------------------------------//
// This is the internal function that does the work of the program called by R //
//-----------------------------------------------------------------------------//

static void SizeSpectrum_internals(int *run_params, double *grid_params, double *pla_params, double *pel_params, double *ben_params, double * det_params, char **names_params, int *flags_params, const ARRAY_INPUTS *arr, OUTPUT_ARRAYS *out)
// run_params is an array of run parameters
// grid_params is an array of grid discretisation values
// pla_params is plankton params
// pel_params is pelagic params
// ben_params is benthic parameters
// det_params is detritus parameters
// names_params is a list of filenames
// flags_params is a list of selection 'flag' parameters
{
    int s,b;
    RUN run;    
    GRID grid;
    COMMUNITY community;
    MATRIX *pelmatrix;
    MATRIX *benmatrix = NULL;
    MATRIX xmatrix, ymatrix;
    double **temp_pel_params;
    double **temp_ben_params;
        
    /*Setup run*/
    setup_run(&run, run_params, names_params[0]);
    run.arr=arr;          //inputs supplied from R as arrays (may hold none)
    run.frc_step=0;
    run.fish_block=0;
    run.out=out;          //output returned to R as arrays (mode 0 = none)
    
    /*Setup grid*/
    setup_grid(&grid, grid_params);

    /*Effort allocated among cells by gravity: the multiplier on fishing mortality in each cell*/
    run.fish_share=NULL;
    run.share_u_pel=NULL;
    run.share_u_ben=NULL;
    if(arr!=NULL && arr->gravity!=NULL){
            int k,l;
            if(arr->gravity->len!=(R_xlen_t)grid.xnum*(R_xlen_t)grid.ynum) error("gravity: access and thickness must have one value per cell (x*y = %d)", grid.xnum*grid.ynum);
            run.fish_share=(double **)safe_malloc(grid.xnum*sizeof(double *),__LINE__,__FILE__);
            for(k=0 ; k<grid.xnum ; k++){
                    run.fish_share[k]=(double *)safe_malloc(grid.ynum*sizeof(double),__LINE__,__FILE__);
                    for(l=0 ; l<grid.ynum ; l++) run.fish_share[k][l]=1.0;
            }
            if(arr->gravity->nupdate>1){
                    size_t len=(size_t)grid.mnum*(size_t)grid.xnum*(size_t)grid.ynum, j;
                    run.share_u_pel=(double **)safe_malloc(run.no_pelagic*sizeof(double *),__LINE__,__FILE__);
                    for(s=0 ; s<run.no_pelagic ; s++){
                            run.share_u_pel[s]=(double *)safe_malloc(len*sizeof(double),__LINE__,__FILE__);
                            for(j=0 ; j<len ; j++) run.share_u_pel[s][j]=0.0;
                    }
                    if(run.no_benthic!=0){
                            run.share_u_ben=(double **)safe_malloc(run.no_benthic*sizeof(double *),__LINE__,__FILE__);
                            for(b=0 ; b<run.no_benthic ; b++){
                                    run.share_u_ben[b]=(double *)safe_malloc(len*sizeof(double),__LINE__,__FILE__);
                                    for(j=0 ; j<len ; j++) run.share_u_ben[b][j]=0.0;
                            }
                    }
            }
    }

    /*Setup community*/
    community.plankton=(PLANKTON *)safe_malloc(sizeof(PLANKTON),__LINE__,__FILE__);
    community.pelagic=(PELAGIC *)safe_malloc(run.no_pelagic*sizeof(PELAGIC),__LINE__,__FILE__);
    if(run.no_benthic!=0){
            community.benthic=(BENTHIC *)safe_malloc(run.no_benthic*sizeof(BENTHIC),__LINE__,__FILE__);
            community.detritus=(DETRITUS *)safe_malloc(sizeof(DETRITUS),__LINE__,__FILE__);
    }
    
    /*Setup plankton*/
    setup_plankton(&run, &grid, community.plankton, pla_params, names_params[1], flags_params[0], flags_params[1], arr); 
    
    /*Setup pelagic*/
    temp_pel_params=setup_pelagic_params(run.no_pelagic, pel_params);
    for(s=0 ; s<run.no_pelagic ; s++){
            community.pelagic[s].fish_arr=(arr!=NULL && arr->fish_pel!=NULL && arr->fish_pel[s].rate!=NULL) ? &arr->fish_pel[s] : NULL;
            setup_pelagic(&run, &grid, &(community.pelagic[s]), temp_pel_params[s], names_params[s+2], flags_params[4*s+2], flags_params[4*s+3], flags_params[4*s+4], flags_params[4*s+5]);
    }
    
    if(run.no_benthic!=0){
            /*Setup benthic*/
            temp_ben_params=setup_benthic_params(run.no_benthic, ben_params);
            for(b=0 ; b<run.no_benthic ; b++){
                    community.benthic[b].fish_arr=(arr!=NULL && arr->fish_ben!=NULL && arr->fish_ben[b].rate!=NULL) ? &arr->fish_ben[b] : NULL;
                    setup_benthic(&run, &grid, &(community.benthic[b]), temp_ben_params[b], names_params[b+run.no_pelagic+2], flags_params[4*b+(2+(4*run.no_pelagic))], flags_params[4*b+(3+(4*run.no_pelagic))], flags_params[4*b+(4+(4*run.no_pelagic))], flags_params[4*b+(5+(4*run.no_pelagic))]);
            }
    
            /*Setup detritus*/
            setup_detritus(&run, &grid, community.detritus, det_params, names_params[run.no_pelagic+run.no_benthic+2], flags_params[2+4*(run.no_pelagic+run.no_benthic)], flags_params[3+4*(run.no_pelagic+run.no_benthic)]);
    }
   
    /*Setup differencing matrices*/
    pelmatrix=(MATRIX *)safe_malloc(run.no_pelagic*sizeof(MATRIX),__LINE__,__FILE__);
    for(s=0 ; s<run.no_pelagic ; s++){
            setup_matrix((community.pelagic[s].ipelmax-community.pelagic[s].ipelmin+1), &(pelmatrix[s]));
    }
    
    if(run.no_benthic!=0){
            benmatrix=(MATRIX *)safe_malloc(run.no_benthic*sizeof(MATRIX),__LINE__,__FILE__);
            for(b=0 ; b<run.no_benthic ; b++){
                    setup_matrix((community.benthic[b].ibenmax-community.benthic[b].ibenmin+1), &(benmatrix[b]));
            }
    }
    
    setup_matrix(grid.xnum, &xmatrix);
    setup_matrix(grid.ynum, &ymatrix);

    /*Arrays returned to R (allocated now that the grid and the number of species are known)*/
    if(out!=NULL && out->mode!=0){
            setup_output_arrays(&run, &grid, out);
    }

    /*Perform all finite differencing*/
    calculate_results(&run, &grid, &community, pelmatrix, benmatrix, &xmatrix, &ymatrix);

    /*Free all allocated memory*/
    free_mem(&run, &grid, &community, pelmatrix, benmatrix, &xmatrix, &ymatrix, temp_pel_params, temp_ben_params);

}

//-----------------------------------//
// Output Function. Not to be edited //
//-----------------------------------//

void calculate_results(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *pelmatrix, MATRIX *benmatrix, MATRIX *xmatrix, MATRIX *ymatrix)
/* Sets up output txt files and implements time stepping process*/
{
     int i,j,k,l,s,b;         //indices for each variable
     int m,t,x,y,n,c;         //maximum number of steps for each variable
     
     double time1=0,time2=0,timed=0;  //variables for countdown timer
     
     m=grid->mnum;
     t=grid->tnum;
     x=grid->xnum;
     y=grid->ynum;
     n=run->no_pelagic;
     c=run->no_benthic;

     /*Optional environmental-forcing plugin: auto-detect <run>/Input/forcing_ts.txt, whose
       first line names the channels present (comma-separated) and each subsequent line gives
       one timestep's values. We map the recognised channel names to column indices here; the
       per-timestep read + apply happens in mass_solver(). Absent file/channels -> neutral.*/
     run->forcing_flag=0;
     run->forcing_ncol=0;
     run->col_pel_tempeff=-1;
     run->col_ben_tempeff=-1;
     run->col_sinking_rate=-1;
     run->col_depth=-1;
     run->col_residence_time=-1;
     run->col_export_attn=-1;
     run->col_export_gamma=-1;
     run->col_export_magg=-1;
     run->col_gravity=-1;
     run->pel_tempeff=1.0;
     run->ben_tempeff=1.0;
     run->sinking_rate=1.0;
     run->depth=1.0;
     run->residence_time=0.0;
     run->export_attn=0.0;
     run->export_gamma=0.0;
     run->export_magg=-1e30;
     run->gravity=0.0;
     run->fptr_forcing=NULL;
     run->forcing_percell=0;
     run->pel_tempeff_c=NULL; run->ben_tempeff_c=NULL; run->sinking_rate_c=NULL;
     run->depth_c=NULL; run->export_attn_c=NULL; run->export_gamma_c=NULL; run->export_magg_c=NULL;
     { char *e=getenv("TEMP_ON_MORT"); run->temp_on_mort = e ? atoi(e) : 1; }
     if(run->arr!=NULL && run->arr->have_forcing){
             /*Forcing channels supplied from R as arrays: forcing_ts.txt is not opened. A channel that is
               present gets its own position as column index, so everything keyed on col_* works unchanged.*/
             int ch;
             for(ch=0 ; ch<N_FRC ; ch++){
                     if(run->arr->frc[ch]!=NULL && run->arr->frc_len[ch]!=(R_xlen_t)grid->xnum*(R_xlen_t)grid->ynum*(R_xlen_t)run->arr->frc_nlayer[ch]){
                             error("%s holds %ld values; expected x*y*layers = %d*%d*%d\n",FRC_NAMES[ch],(long)run->arr->frc_len[ch],grid->xnum,grid->ynum,run->arr->frc_nlayer[ch]);
                     }
             }
             run->forcing_flag=1;
             run->forcing_ncol=N_FRC;
             run->col_pel_tempeff    = run->arr->frc[FRC_PEL_TEMPEFF]!=NULL    ? FRC_PEL_TEMPEFF    : -1;
             run->col_ben_tempeff    = run->arr->frc[FRC_BEN_TEMPEFF]!=NULL    ? FRC_BEN_TEMPEFF    : -1;
             run->col_sinking_rate   = run->arr->frc[FRC_SINKING_RATE]!=NULL   ? FRC_SINKING_RATE   : -1;
             run->col_depth          = run->arr->frc[FRC_DEPTH]!=NULL          ? FRC_DEPTH          : -1;
             run->col_residence_time = run->arr->frc[FRC_RESIDENCE_TIME]!=NULL ? FRC_RESIDENCE_TIME : -1;
             run->col_export_attn    = run->arr->frc[FRC_EXPORT_ATTN]!=NULL    ? FRC_EXPORT_ATTN    : -1;
             run->col_export_gamma   = run->arr->frc[FRC_EXPORT_GAMMA]!=NULL   ? FRC_EXPORT_GAMMA   : -1;
             run->col_export_magg    = run->arr->frc[FRC_EXPORT_MAGG]!=NULL    ? FRC_EXPORT_MAGG    : -1;
             run->col_gravity        = run->arr->frc[FRC_GRAVITY]!=NULL        ? FRC_GRAVITY        : -1;
     }
     else{
             char fname_forcing[256];
             snprintf(fname_forcing, sizeof(fname_forcing), "%s/Input/forcing_ts.txt", run->filename);
             run->fptr_forcing=fopen(fname_forcing,"r");
             if(run->fptr_forcing!=NULL){
                     run->forcing_flag=1;
                     char hdr[512];
                     if(fgets(hdr,sizeof(hdr),run->fptr_forcing)!=NULL){
                             /*strip trailing newline*/
                             size_t hl=strlen(hdr);
                             while(hl>0 && (hdr[hl-1]=='\n'||hdr[hl-1]=='\r'||hdr[hl-1]==' ')){ hdr[hl-1]='\0'; hl--; }
                             int col=0;
                             char *name=strtok(hdr,",");
                             while(name!=NULL){
                                     /*trim leading spaces*/
                                     while(*name==' ') name++;
                                     if(strcmp(name,"pel_tempeff")==0)  run->col_pel_tempeff=col;
                                     else if(strcmp(name,"ben_tempeff")==0)  run->col_ben_tempeff=col;
                                     else if(strcmp(name,"sinking_rate")==0)  run->col_sinking_rate=col;
                                     else if(strcmp(name,"depth")==0)  run->col_depth=col;
                                     else if(strcmp(name,"residence_time")==0)  run->col_residence_time=col;
                                     else if(strcmp(name,"export_attn")==0)  run->col_export_attn=col;
                                     else if(strcmp(name,"export_gamma")==0)  run->col_export_gamma=col;
                                     else if(strcmp(name,"export_magg")==0)  run->col_export_magg=col;
                                     else if(strcmp(name,"gravity")==0)  run->col_gravity=col;
                                     col++;
                                     name=strtok(NULL,",");
                             }
                             run->forcing_ncol=col;
                     }
             }
     }
     /* Per-cell forcing (opt-in): FORCING_PERCELL=1 + spatial_dim>=1 => forcing_ts.txt has one row per
        CELL per timestep (x*y rows, ordered x-outer y-inner, like the plankton _ts file), so temperature
        and export can vary across space. Allocate the per-cell multiplier arrays (init neutral=1). When
        off (default), the single global row is used and these arrays are unused -> behaviour unchanged. */
     if(run->forcing_flag==1){
             char *e=getenv("FORCING_PERCELL");
             int from_arrays = (run->arr!=NULL && run->arr->have_forcing);   /*array forcing is always per cell in a spatial run*/
             run->forcing_percell = ((from_arrays || (e && atoi(e)!=0)) && (run->spatial_dim==1 || run->spatial_dim==2)) ? 1 : 0;
             int xx=grid->xnum, yy=grid->ynum, kk, ll;
             run->pel_tempeff_c=(double**)safe_malloc(xx*sizeof(double*),__LINE__,__FILE__);
             run->ben_tempeff_c=(double**)safe_malloc(xx*sizeof(double*),__LINE__,__FILE__);
             run->sinking_rate_c=(double**)safe_malloc(xx*sizeof(double*),__LINE__,__FILE__);
             run->depth_c=(double**)safe_malloc(xx*sizeof(double*),__LINE__,__FILE__);
             run->export_attn_c=(double**)safe_malloc(xx*sizeof(double*),__LINE__,__FILE__);
             run->export_gamma_c=(double**)safe_malloc(xx*sizeof(double*),__LINE__,__FILE__);
             run->export_magg_c=(double**)safe_malloc(xx*sizeof(double*),__LINE__,__FILE__);
             for(kk=0;kk<xx;kk++){
                     run->pel_tempeff_c[kk]=(double*)safe_malloc(yy*sizeof(double),__LINE__,__FILE__);
                     run->ben_tempeff_c[kk]=(double*)safe_malloc(yy*sizeof(double),__LINE__,__FILE__);
                     run->sinking_rate_c[kk]=(double*)safe_malloc(yy*sizeof(double),__LINE__,__FILE__);
                     run->depth_c[kk]=(double*)safe_malloc(yy*sizeof(double),__LINE__,__FILE__);
                     run->export_attn_c[kk]=(double*)safe_malloc(yy*sizeof(double),__LINE__,__FILE__);
                     run->export_gamma_c[kk]=(double*)safe_malloc(yy*sizeof(double),__LINE__,__FILE__);
                     run->export_magg_c[kk]=(double*)safe_malloc(yy*sizeof(double),__LINE__,__FILE__);
                     for(ll=0;ll<yy;ll++){ run->pel_tempeff_c[kk][ll]=1.0; run->ben_tempeff_c[kk][ll]=1.0; run->sinking_rate_c[kk][ll]=1.0;
                             run->depth_c[kk][ll]=1.0; run->export_attn_c[kk][ll]=0.0; run->export_gamma_c[kk][ll]=0.0; run->export_magg_c[kk][ll]=-1e30; }
             }
     }

     /*Create temporary arrays for timestepping process*/
     double ***temp_plankton;
     double ****temp_pelagic;
     double ****temp_benthic;
     double **temp_detritus;
     
     temp_plankton=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
     for(i=0 ; i< m ; i++){
             temp_plankton[i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             for(k=0 ; k<x ; k++){
                     temp_plankton[i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     for(l=0 ; l<y ; l++){
                             temp_plankton[i][k][l]=0;
                     }
             }
     }
     temp_pelagic=(double ****)safe_malloc(n*sizeof(double ***),__LINE__,__FILE__);
     for(s=0 ; s<n ; s++){
             temp_pelagic[s]=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
             for(i=0 ; i< m ; i++){
                     temp_pelagic[s][i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
                     for(k=0 ; k<x ; k++){
                             temp_pelagic[s][i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                             for(l=0 ; l<y ; l++){
                                     temp_pelagic[s][i][k][l]=0;
                             }
                     }
             }
     }
     
     if(run->no_benthic!=0){
             temp_benthic=(double ****)safe_malloc(c*sizeof(double ***),__LINE__,__FILE__);
             for(b=0 ; b<c ; b++){
                     temp_benthic[b]=(double ***)safe_malloc(m*sizeof(double **),__LINE__,__FILE__);
                     for(i=0 ; i< m ; i++){
                             temp_benthic[b][i]=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
                             for(k=0 ; k<x ; k++){
                                     temp_benthic[b][i][k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                                     for(l=0 ; l<y ; l++){
                                             temp_benthic[b][i][k][l]=0;
                                     }
                             }
                     }
             }
             temp_detritus=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
             for(k=0 ; k<x ; k++){
                     temp_detritus[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
                     for(l=0 ; l<y ; l++){
                             temp_detritus[k][l]=0;
                     }
             }
     }
     
     /*Calculate initial fishing pressures*/
     calculate_fishing(run, grid, community);

     /*Calculate starting g and mu values*/
     calculate_g_and_mu(run, grid, community);
     
     /*Calculate starting reproduction allocation*/
     calculate_reproduction(run, grid, community);

     /*Print run info and starting output values*/
     int write_files=(run->out==NULL || run->out->mode!=2);   /*arrays-only output skips the per-timestep text files*/
     if(grid->toutmin == 0){
             capture_output(run, grid, community, 0);
     }
     /*plankton*/
     if(write_files && grid->toutmin == 0){
             print_timestep_plankton(grid, community->plankton, 0);
             print_plankton_summary(grid, community->plankton, 0);
     }

     /*Pelagic*/
     for(s=0 ; s<n ; s++){
             if(write_files && grid->toutmin == 0){
                     print_timestep_pelagic(run, grid, &(community->pelagic[s]), 0);
                     print_pelagic_summary(grid, &(community->pelagic[s]), 0);
             }
     }
     
     if(run->no_benthic!=0){
             /*Benthic*/
             for(b=0 ; b<c ; b++){
                     if(write_files && grid->toutmin == 0){
                             print_timestep_benthic(grid, &(community->benthic[b]), 0);
                             print_benthic_summary(grid, &(community->benthic[b]), 0);
                     }
             }
             /*Detritus*/
             if(write_files && grid->toutmin == 0){
                     print_detritus_summary(grid, community->detritus, 0);
             }
     }
     
     /*Parameter File*/
     run->fptr_summ=safe_fopen(run->fname_summ,"a",__LINE__,__FILE__);
     
     print_run(run, run->fptr_summ);
     print_grid(grid, run->fptr_summ);
     print_plankton(community->plankton, run->fptr_summ);
     for(s=0 ; s<n ; s++){
             print_pelagic(&(community->pelagic[s]),run->fptr_summ);
     }
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     print_benthic(&(community->benthic[b]),run->fptr_summ);
             }
             print_detritus(community->detritus, run->fptr_summ);
     }
     fclose(run->fptr_summ);

     /*Time Stepping Loop*/
     for(j=1 ; j<t ; j++){
                         
            /*Countdown timer for the program*/
             
            //Checks if user has terminated function in R (e.g. by pressing escape)
             R_CheckUserInterrupt();
             
             /*Countdown timer for the program*/
             if(j==1){
                      time2=(double)clock()/CLOCKS_PER_SEC;
             }
             else{
                      if(j%((int)rint(grid->toutstep/grid->tstep))==0){
                             time1=time2;
                             time2=(double)clock()/CLOCKS_PER_SEC;
                             timed=(time2-time1)*(grid->tstep/grid->toutstep);
                      }
                             Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                             R_FlushConsole();
                             Rprintf("%6.1f   %d",timed*(grid->tnum-j),grid->tnum-j);
                             R_FlushConsole();
             }
             if(j==(t-1)){
                      Rprintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b");
                      R_FlushConsole();
                      Rprintf("\n");
                      R_FlushConsole();
             }

             /*Store previous time step in temporary array for each species*/
             for(i=0 ; i<m ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     temp_plankton[i][k][l]=community->plankton->u_values[i][k][l];
                             }
                     }
             }
             for(s=0 ; s<n ; s++){
                     for(i=0 ; i<m ; i++){
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             temp_pelagic[s][i][k][l]=community->pelagic[s].u_values[i][k][l];
                                     }
                             }
                     }
             }
             
             if(run->no_benthic!=0){
                     for(b=0 ; b<c ; b++){
                             for(i=0 ; i<m ; i++){
                                     for(k=0 ; k<x ; k++){
                                             for(l=0 ; l<y ; l++){
                                                     temp_benthic[b][i][k][l]=community->benthic[b].u_values[i][k][l];
                                             }
                                     }
                             }
                     }
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     temp_detritus[k][l]=community->detritus->w_values[k][l];
                             }
                     }
             }


                          
             /*Calculate u(m,t) using previous timestep values*/
             mass_solver(run, grid, community, pelmatrix, benmatrix, temp_plankton, temp_pelagic, temp_benthic, temp_detritus);

             /*Calculate fishing pressures*/
             calculate_fishing(run, grid, community);
             
             /*Calculate g and mu using current timestep values for use in movement and next time step*/
             calculate_g_and_mu(run, grid, community);

             /*Calculate starting reproduction allocation*/
             calculate_reproduction(run, grid, community);
             


             /*If Spatial Model is being used (1D = x movement only; 2D = x then y).
               NOTE: the outer guard must admit spatial_dim==2 as well — the inner spatial_dim==2
               block below only runs when this outer test passes, so a `==1`-only guard left ALL
               movement dead for 2D runs (neither x nor y solver ran). */
             if(run->spatial_dim==1 || run->spatial_dim==2){

                     /*Store itermediate time step as temp_pelagic*/

                     for(s=0 ; s<n ; s++){
                             for(i=0 ; i<m ; i++){
                                    for(k=0 ; k<x ; k++){
                                             for(l=0 ; l<y ; l++){
                                                     temp_pelagic[s][i][k][l]=community->pelagic[s].u_values[i][k][l];
                                             }
                                     }
                             }
                     }


                     /*Calculate spatial movements using intermediate values*/
                     xmove_solver(run, grid, community, xmatrix, temp_pelagic, j);
             
                     
                    if(run->spatial_dim==2){
                     
                             /*Store itermediate time step as temp_pelagic*/
                             for(s=0 ; s<n ; s++){
                                     for(i=0 ; i<m ; i++){
                                             for(k=0 ; k<x ; k++){
                                                    for(l=0 ; l<y ; l++){
                                                             temp_pelagic[s][i][k][l]=community->pelagic[s].u_values[i][k][l];
                                                     }
                                             }
                                     }
                             }
                             /*Calculate spatial movements using intermediate values*/
                             ymove_solver(run, grid, community, ymatrix, temp_pelagic, j);
             
                     }
             }
             
             
             
             /*Print current timestep values to files*/
             if(grid->t_values[j]>=grid->toutmin && grid->t_values[j]<=grid->toutmax){
                     if(j%((int)rint(grid->toutstep/grid->tstep))==0){
                             capture_output(run, grid, community, j);
                     }
                     if(write_files && j%((int)rint(grid->toutstep/grid->tstep))==0){
                             
                             print_timestep_plankton(grid, community->plankton,j);
                             print_plankton_summary(grid, community->plankton, j);
                             for(s=0 ; s<n ; s++){
                                     print_timestep_pelagic(run, grid, &(community->pelagic[s]),j);
                                     print_pelagic_summary(grid, &(community->pelagic[s]),j);
                             }
                             if(run->no_benthic!=0){
                                     for(b=0 ; b<c  ;b++){
                                             print_timestep_benthic(grid, &(community->benthic[b]),j);
                                             print_benthic_summary(grid, &(community->benthic[b]),j);
                                     }
                                     print_detritus_summary(grid, community->detritus,j);
                             }
                     }
             }
     }
     
     /*Print final detailed timestep*/
     if(write_files){
             print_detailed_plankton(grid, community->plankton, t-1);
             for(s=0 ; s<n ; s++){
                     print_detailed_pelagic(run, grid, &(community->pelagic[s]), t-1);
             }
             if(run->no_benthic!=0){
                     for(b=0 ; b<c  ;b++){
                             print_detailed_benthic(grid, &(community->benthic[b]),t-1);
                     }
             }
     }

     /*Free pointers*/
     for(i=0 ; i<m ; i++){
             for(k=0 ; k<x ; k++){
                     free(temp_plankton[i][k]);
             }
             free(temp_plankton[i]);
     }
     free(temp_plankton);
     
     for(s=0 ; s<n ; s++){
             for(i=0 ; i<m ; i++){
                     for(k=0 ; k<x ; k++){
                             free(temp_pelagic[s][i][k]);
                     }
                     free(temp_pelagic[s][i]);
             }
             free(temp_pelagic[s]);
     }
     free(temp_pelagic);
     
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     for(i=0 ; i<m ; i++){
                             for(k=0 ; k<x ; k++){
                                     free(temp_benthic[b][i][k]);
                             }
                             free(temp_benthic[b][i]);
                     }
                     free(temp_benthic[b]);
             }
             free(temp_benthic);
             
             for(k=0 ; k<x ; k++){
                     free(temp_detritus[k]);
             }
             free(temp_detritus);
     }
     
     /*Close ts files*/
     if(community->plankton->ts_flag==1 && community->plankton->arr==NULL){
             fclose(community->plankton->fptr_ts);
     }
     if(run->forcing_flag==1 && run->fptr_forcing!=NULL){
             fclose(run->fptr_forcing);
             if(run->pel_tempeff_c!=NULL){
                     for(k=0;k<x;k++){ free(run->pel_tempeff_c[k]); free(run->ben_tempeff_c[k]); free(run->sinking_rate_c[k]);
                             free(run->depth_c[k]); free(run->export_attn_c[k]); free(run->export_gamma_c[k]); free(run->export_magg_c[k]); }
                     free(run->pel_tempeff_c); free(run->ben_tempeff_c); free(run->sinking_rate_c);
                     free(run->depth_c); free(run->export_attn_c); free(run->export_gamma_c); free(run->export_magg_c);
             }
     }
     for(s=0 ; s<n ; s++){
             if(community->pelagic[s].ts_flag==1){
                     fclose(community->pelagic[s].fptr_ts);
             }
     }
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     if(community->benthic[b].ts_flag==1){
                             fclose(community->benthic[b].fptr_ts);
                     }
             }
             if(community->detritus->ts_flag==1){
                     fclose(community->detritus->fptr_ts);
             }
     }
     
     /*Close fishing files*/
     for(s=0 ; s<n ; s++){
             if(community->pelagic[s].fishing_flag==1 && community->pelagic[s].fish_arr==NULL){
                     fclose(community->pelagic[s].fptr_fish_ts);
             }
     }
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     if(community->benthic[b].fishing_flag==1 && community->benthic[b].fish_arr==NULL){
                             fclose(community->benthic[b].fptr_fish_ts);
                     }
             }
     }
     
     /*Close reproduction files*/
     for(s=0 ; s<n ; s++){
             if(community->pelagic[s].rep_method==1){
                     if(community->pelagic[s].ts_flag==0){
                             fclose(community->pelagic[s].fptr_rep_ts);
                     }
             }
     }
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     if(community->benthic[b].rep_method==1){
                             if(community->benthic[b].ts_flag==0){
                                     fclose(community->benthic[b].fptr_rep_ts);
                             }
                     }
             }
     }

}

//-------------------------------------------------//
// Differencing Scheme Functions. Not to be edited //
//-------------------------------------------------//

void mass_solver(RUN *run, GRID *grid, COMMUNITY *community, MATRIX *pelmatrix, MATRIX *benmatrix, double ***temp_plankton, double ****temp_pelagic, double ****temp_benthic, double **temp_detritus)
/*Impliments the differencing scheme to obtain the values at the next time step*/
{
     int i,k,l,s,b;
     int x,y,n,c;
     double dt,dm;
     
     dm=grid->mstep;
     dt=grid->tstep;
     
     x=grid->xnum;
     y=grid->ynum;
     n=run->no_pelagic;
     c=run->no_benthic;

     /*Advance the environmental-forcing plugin one timestep (if active), in lockstep with the
       plankton re-read below. Parse this timestep's row into the recognised channels, then:
         - rescale every species' feeding (A) and background mortality (mu_0) from their stored
           base values by pel/ben_tempeff - the growth, predation-mortality and reproduction
           integrals computed in calculate_g_and_mu() (after this solve) then carry the
           temperature effect, matching sizemodel(): feeding & other-mortality are temperature-
           scaled; senescence (mu_s) and fishing are not;
         - store sinking_rate on run for use in g_det() (detritus input).*/
     if(run->forcing_flag==1 && run->forcing_percell){
             /* PER-CELL forcing: read x*y rows this timestep (one per cell, x-outer y-inner) into the
                per-cell temperature/sinking arrays. A and mu_0 are held at BASE here; the per-cell
                temperature is applied in the feeding/mortality integrals (via run->*_tempeff_c[k][l]).
                Only pel_tempeff/ben_tempeff/sinking_rate are made per-cell in this mode. */
             char fline[1024]; double vals[64];
             for(k=0 ; k<x ; k++){
                     for(l=0 ; l<y ; l++){
                             int nv=0;
                             if(run->arr!=NULL && run->arr->have_forcing){
                                     nv=forcing_from_arrays(run, grid, k, l, vals);
                             }
                             else{
                                     if(fgets(fline,sizeof(fline),run->fptr_forcing)==NULL) break;
                                     char *tok=strtok(fline,",");
                                     while(tok!=NULL && nv<64){ vals[nv++]=atof(tok); tok=strtok(NULL,","); }
                             }
                             run->pel_tempeff_c[k][l]  = (run->col_pel_tempeff>=0 && run->col_pel_tempeff<nv) ? vals[run->col_pel_tempeff] : 1.0;
                             run->ben_tempeff_c[k][l]  = (run->col_ben_tempeff>=0 && run->col_ben_tempeff<nv) ? vals[run->col_ben_tempeff] : 1.0;
                             run->sinking_rate_c[k][l] = (run->col_sinking_rate>=0 && run->col_sinking_rate<nv) ? vals[run->col_sinking_rate] : 1.0;
                             run->depth_c[k][l]        = (run->col_depth>=0 && run->col_depth<nv) ? vals[run->col_depth] : 1.0;
                             run->export_attn_c[k][l]  = (run->col_export_attn>=0 && run->col_export_attn<nv) ? vals[run->col_export_attn] : 0.0;
                             run->export_gamma_c[k][l] = (run->col_export_gamma>=0 && run->col_export_gamma<nv) ? vals[run->col_export_gamma] : 0.0;
                             run->export_magg_c[k][l]  = (run->col_export_magg>=0 && run->col_export_magg<nv) ? vals[run->col_export_magg] : -1e30;
                     }
             }
             if(run->arr!=NULL && run->arr->have_forcing){
                     /*channels with no per-cell storage are taken from the first cell (R checks they are uniform in space)*/
                     forcing_from_arrays(run, grid, 0, 0, vals);
                     if(run->col_residence_time>=0) run->residence_time=vals[run->col_residence_time];
                     if(run->col_gravity>=0) run->gravity=vals[run->col_gravity];
                     run->frc_step++;
             }
             for(s=0 ; s<n ; s++){ community->pelagic[s].A=community->pelagic[s].A_base; community->pelagic[s].mu_0=community->pelagic[s].mu_0_base; }
             for(b=0 ; b<c ; b++){ community->benthic[b].A=community->benthic[b].A_base; community->benthic[b].mu_0=community->benthic[b].mu_0_base; }
     }
     else if(run->forcing_flag==1){
             char fline[1024];
             double vals[64];
             int nv=0;
             int got=0;
             if(run->arr!=NULL && run->arr->have_forcing){
                     nv=forcing_from_arrays(run, grid, 0, 0, vals);
                     run->frc_step++;
                     got=1;
             }
             else if(fgets(fline,sizeof(fline),run->fptr_forcing)!=NULL){
                     char *tok=strtok(fline,",");
                     while(tok!=NULL && nv<64){ vals[nv++]=atof(tok); tok=strtok(NULL,","); }
                     got=1;
             }
             if(got){
                     if(run->col_pel_tempeff>=0 && run->col_pel_tempeff<nv) run->pel_tempeff=vals[run->col_pel_tempeff];
                     if(run->col_ben_tempeff>=0 && run->col_ben_tempeff<nv) run->ben_tempeff=vals[run->col_ben_tempeff];
                     if(run->col_sinking_rate>=0 && run->col_sinking_rate<nv) run->sinking_rate=vals[run->col_sinking_rate];
                     if(run->col_depth>=0 && run->col_depth<nv) run->depth=vals[run->col_depth];
                     if(run->col_residence_time>=0 && run->col_residence_time<nv) run->residence_time=vals[run->col_residence_time];
                     if(run->col_export_attn>=0 && run->col_export_attn<nv) run->export_attn=vals[run->col_export_attn];
                     if(run->col_export_gamma>=0 && run->col_export_gamma<nv) run->export_gamma=vals[run->col_export_gamma];
                     if(run->col_export_magg>=0 && run->col_export_magg<nv) run->export_magg=vals[run->col_export_magg];
                     if(run->col_gravity>=0 && run->col_gravity<nv) run->gravity=vals[run->col_gravity];
             }
             /* TEMP_ON_MORT (env, default 1): whether temperature scales background mortality.
                =0 -> temperature drives FEEDING only; mu_0 held at base. Rationale: the model does
                not resolve metabolic maintenance costs (respiration) or their temperature response,
                so mu_0*tempeff is a lone, unbalanced "cost of warmth"; feeding has a clean
                encounter-rate justification, mortality does not. Read once. */
             static int tom = -1;
             if(tom < 0){ char *e = getenv("TEMP_ON_MORT"); tom = e ? atoi(e) : 1; }
             for(s=0 ; s<n ; s++){
                     community->pelagic[s].A    = community->pelagic[s].A_base    * run->pel_tempeff;
                     community->pelagic[s].mu_0 = tom ? community->pelagic[s].mu_0_base * run->pel_tempeff
                                                      : community->pelagic[s].mu_0_base;
             }
             for(b=0 ; b<c ; b++){
                     community->benthic[b].A    = community->benthic[b].A_base    * run->ben_tempeff;
                     community->benthic[b].mu_0 = tom ? community->benthic[b].mu_0_base * run->ben_tempeff
                                                      : community->benthic[b].mu_0_base;
             }
     }

     /*Used for calculating detritus updates*/
     double **temp_det_g;
     double **temp_det_mu;
     
     temp_det_g=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     temp_det_mu=(double **)safe_malloc(x*sizeof(double *),__LINE__,__FILE__);
     for(k=0 ; k<x ; k++){
             temp_det_g[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             temp_det_mu[k]=(double *)safe_malloc(y*sizeof(double),__LINE__,__FILE__);
             for(l=0 ; l<y ; l++){
                     temp_det_g[k][l]=0;
                     temp_det_mu[k][l]=0;
             }
     }
     /*******************/
     /* Plankton System */
     /*******************/
     if(community->plankton->ts_flag==0){
             // Calculate values for Plankton Size Spectra using internal dynamics //
             for(i=community->plankton->iplamin ; i<(community->plankton->iplamax+1) ; i++){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     community->plankton->u_values[i][k][l]=temp_plankton[i][k][l];  //All plankton systems are invariant
                             }
                     }
             }
     }
     if(community->plankton->ts_flag==1 && community->plankton->arr!=NULL){
             /*Next time block from the R arrays*/
             plankton_from_arrays(grid, community->plankton);
     }
     else if(community->plankton->ts_flag==1){
             int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
             char *input;               //input string to store line as string
             input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
             char sep[]=",";            //identifying the ',' char as seperating values
             char *result = NULL;       //declaring a pointer for storing each individual value as a char
             
             /*Read in a single timestep*/
             /* each line one at a time (corresponding to spectra at a spatial point)*/
             for(k=0 ; k<x ; k++){
                     for(l=0 ; l<y ; l++){
                             if(fgets(input,MAX_LEN,community->plankton->fptr_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
                             input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                             
                             result=(char *)strtok(input,sep);    //this is the 1st value
                             i=0;
                             while(result!=NULL){
                                     community->plankton->u_values[i][k][l]=atof(result);
                                     result=(char *)strtok(NULL,sep);
                                     i++;
                             }
                     }
             }
             free(input);
     }
     
     /******************/
     /* Pelagic System */
     /******************/     
     for(s=0 ; s<n ; s++){
             if(community->pelagic[s].ts_flag==0){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     for(i=0 ; i<pelmatrix[s].size ; i++){
                                             pelmatrix[s].a[i]=-(dt/dm)*community->pelagic[s].g_values[i-1+community->pelagic[s].ipelmin][k][l];
                                             pelmatrix[s].b[i]=(1+(dt*community->pelagic[s].mu_values[i+community->pelagic[s].ipelmin][k][l])+(dt/dm)*community->pelagic[s].g_values[i+community->pelagic[s].ipelmin][k][l]);
                                             pelmatrix[s].c[i]=0;
                                             pelmatrix[s].r[i]=temp_pelagic[s][i+community->pelagic[s].ipelmin][k][l];
                                     }
                                     pelmatrix[s].a[0]=0;
                                     tridag(&(pelmatrix[s]));
                                     for(i=0 ; i<pelmatrix[s].size ; i++){
                                             community->pelagic[s].u_values[i+community->pelagic[s].ipelmin][k][l]=pelmatrix[s].u[i];
                                     }
                             }
                     }
             }
             if(community->pelagic[s].ts_flag==1){
                     int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
                     char *input;               //input string to store line as string
                     input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                     char sep[]=",";            //identifying the ',' char as seperating values
                     char *result = NULL;       //declaring a pointer for storing each individual value as a char
             
                     /*Read in a single timestep*/
                     /* each line one at a time (corresponding to spectra at a spatial point)*/
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     if(fgets(input,MAX_LEN,community->pelagic[s].fptr_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
                                     input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                                     
                                     result=(char *)strtok(input,sep);    //this is the 1st value
                                     i=0;
                                     while(result!=NULL){
                                             community->pelagic[s].u_values[i][k][l]=atof(result);
                                             result=(char *)strtok(NULL,sep);
                                             i++;
                                     }
                             }
                     }
                     free(input);
             }
     }
     
     /******************/
     /* Benthic System */
     /******************/
     if(run->no_benthic!=0){
             for(b=0 ; b<c ; b++){
                     if(community->benthic[b].ts_flag==0){
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             for(i=0 ; i<benmatrix[b].size ; i++){
                                                     benmatrix[b].a[i]=-(dt/dm)*community->benthic[b].g_values[i-1+community->benthic[b].ibenmin][k][l];
                                                     benmatrix[b].b[i]=(1+(dt*community->benthic[b].mu_values[i+community->benthic[b].ibenmin][k][l])+(dt/dm)*community->benthic[b].g_values[i+community->benthic[b].ibenmin][k][l]);
                                                     benmatrix[b].c[i]=0;
                                                     benmatrix[b].r[i]=temp_benthic[b][i+community->benthic[b].ibenmin][k][l];
                                             }
                                             benmatrix[b].a[0]=0;
                                             tridag(&(benmatrix[b]));
                                             for(i=0 ; i<benmatrix[b].size ; i++){
                                                     community->benthic[b].u_values[i+community->benthic[b].ibenmin][k][l]=benmatrix[b].u[i];
                                             }
                                     }
                             }
                     }
                     if(community->benthic[b].ts_flag==1){
                             int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
                             char *input;               //input string to store line as string
                             input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                             char sep[]=",";            //identifying the ',' char as seperating values
                             char *result = NULL;       //declaring a pointer for storing each individual value as a char
                             
                             /*Read in a single timestep*/
                             /* each line one at a time (corresponding to spectra at a spatial point)*/
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             if(fgets(input,MAX_LEN,community->benthic[b].fptr_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
                                             input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                                             
                                             result=(char *)strtok(input,sep);    //this is the 1st value
                                             i=0;
                                             while(result!=NULL){
                                                     community->benthic[b].u_values[i][k][l]=atof(result);
                                                     result=(char *)strtok(NULL,sep);
                                                     i++;
                                             }
                                     }
                             }
                             free(input);
                     }
             }
             
             /*******************/
             /* Detritus System */
             /*******************/
             if(community->detritus->ts_flag==0){
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     /* Analytical (exponential) integrator for the detritus pool.
                                        Over a step the fish/plankton state is fixed, so detritus obeys a
                                        LINEAR ODE  dW/dt = I - kappa*W , with input I = g_values (constant)
                                        and loss mu_values = kappa*W (detritivory is proportional to W).
                                        The exact solution  W(t+dt) = W* + (W0-W*)*exp(-kappa*dt),
                                        W* = I/kappa,  is UNCONDITIONALLY STABLE and provably non-negative
                                        (a convex combination of W0>=0 and W*>=0) - unlike the explicit
                                        Euler step W += dt*(I-kappa*W), which overshoots to negative when
                                        kappa*dt > 1 and drives the NaN/collapse cascade. Removes the
                                        timestep stability limit on the detritus pool. */
                                     {
                                       double W0 = temp_detritus[k][l];
                                       double In = community->detritus->g_values[k][l];        /* gross input rate I */
                                       /* Dunne et al. (2007) eq.3 burial: a fraction of the incoming detrital
                                          flux is buried (permanently lost) before joining the pool, EXACTLY as
                                          sizemodel does: burial = input_w*(0.013+0.53*iw^2/(7+iw)^2) applied to
                                          the model's native input_w. dbpmr's densities (like sizemodel's) are
                                          per-VOLUME - anchored to the depth-averaged plankton concentration - so
                                          input_w is a per-volume flux and the Dunne fraction is applied to it
                                          DIRECTLY, with no areal (x depth) rescale, matching sizemodel. Burial is
                                          OPT-IN: only when a `depth` forcing channel is supplied (col_depth>=0),
                                          so existing runs without the forcing plugin keep their original dynamics. */
                                       if(run->col_depth>=0){
                                               double bfrac = 0.013 + 0.53*In*In/((7.0+In)*(7.0+In));
                                               In = In * (1.0 - bfrac);                         /* net of burial   */
                                       }
                                       double loss = community->detritus->mu_values[k][l];     /* = kappa*W0   */
                                       double kappa = (W0 > 1e-30) ? loss / W0 : 0.0;
                                       /* Alternative/opt-in residence-time closure: a first-order loss W/tau on the
                                          detritus POOL (remineralisation/export from the benthic layer), added to the
                                          detritivory rate. Dimensionally clean in dbpmr's per-volume units - no areal
                                          ambiguity - and the exact integrator below absorbs it (K = kappa + 1/tau, so
                                          the equilibrium is W* = In/(kappa+1/tau); with negligible detritivory that is
                                          the classic W* = In*tau). Supplied via the residence_time channel (yr); an
                                          alternative to Dunne burial (Julia's commented CMIP5 variant made this loss
                                          proportional to W). Off when the channel is absent or tau<=0. */
                                       if(run->col_residence_time>=0 && run->residence_time>0.0){
                                               kappa += 1.0/run->residence_time;
                                       }
                                       if(kappa > 1e-30){
                                               double Wstar = In / kappa;
                                               community->detritus->w_values[k][l] = Wstar + (W0 - Wstar)*exp(-kappa*dt);
                                       }else{
                                               community->detritus->w_values[k][l] = W0 + In*dt;   /* kappa -> 0 limit */
                                       }
                                       /* defensive backstop only against a non-finite upstream rate */
                                       if(!isfinite(community->detritus->w_values[k][l]) || community->detritus->w_values[k][l]<0){
                                               community->detritus->w_values[k][l]=0;
                                       }
                                     }
                                     temp_det_g[k][l]=community->detritus->g_values[k][l];
                                     temp_det_mu[k][l]=community->detritus->mu_values[k][l];
                             }
                     }
             }
             if(community->detritus->ts_flag==1){
                     int MAX_LEN=grid->mnum*15; //this is the upper bound for the row length in characters
                     char *input;               //input string to store line as string
                     input=(char *)safe_malloc(MAX_LEN*sizeof(char),__LINE__,__FILE__); 
                     char sep[]=",";            //identifying the ',' char as seperating values
                     char *result = NULL;       //declaring a pointer for storing each individual value as a char
             
                     /*Read in a single timestep*/
                     /* each line one at a time (corresponding to spectra at a spatial point)*/
                     for(k=0 ; k<x ; k++){
                             for(l=0 ; l<y ; l++){
                                     if(fgets(input,MAX_LEN,community->detritus->fptr_ts)==NULL){/* unexpected EOF: leave buffer as-is (pre-existing behaviour) */}
                                     input[strlen(input)-1]='\0'; //replacing the newline \n char of the line of text with a null \0 char
                                     
                                     result=(char *)strtok(input,sep);    //this is the 1st value
                                     while(result!=NULL){
                                             community->detritus->w_values[k][l]=atof(result);
                                             result=(char *)strtok(NULL,sep);
                                     }
                             }
                     }
                     free(input);
             }
     }
     
     //Calculate fully implicit method if required
     if(run->diff_method==1){
             //Use intermediate u values to calculate new g and mu for full implicit methods
             calculate_g_and_mu(run, grid, community);
     
             //Use new g and mu to claculate final u values
             for(s=0 ; s<n ; s++){
                     if(community->pelagic[s].ts_flag==0){
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             for(i=0 ; i<pelmatrix[s].size ; i++){
                                                     pelmatrix[s].a[i]=-(dt/dm)*community->pelagic[s].g_values[i-1+community->pelagic[s].ipelmin][k][l];
                                                     pelmatrix[s].b[i]=(1+(dt*community->pelagic[s].mu_values[i+community->pelagic[s].ipelmin][k][l])+(dt/dm)*community->pelagic[s].g_values[i+community->pelagic[s].ipelmin][k][l]);
                                                     pelmatrix[s].c[i]=0;
                                                     pelmatrix[s].r[i]=temp_pelagic[s][i+community->pelagic[s].ipelmin][k][l];
                                             }
                                             pelmatrix[s].a[0]=0;
                                             tridag(&(pelmatrix[s]));
                                             for(i=0 ; i<pelmatrix[s].size ; i++){
                                                     community->pelagic[s].u_values[i+community->pelagic[s].ipelmin][k][l]=pelmatrix[s].u[i];
                                             }
                                     }
                             }
                     }
             }
             
             if(run->no_benthic!=0){
                     for(b=0 ; b<c ; b++){
                             if(community->benthic[b].ts_flag==0){
                                     for(k=0 ; k<x ; k++){
                                             for(l=0 ; l<y ; l++){
                                                     for(i=0 ; i<benmatrix[b].size ; i++){
                                                             benmatrix[b].a[i]=-(dt/dm)*community->benthic[b].g_values[i-1+community->benthic[b].ibenmin][k][l];
                                                             benmatrix[b].b[i]=(1+(dt*community->benthic[b].mu_values[i+community->benthic[b].ibenmin][k][l])+(dt/dm)*community->benthic[b].g_values[i+community->benthic[b].ibenmin][k][l]);
                                                             benmatrix[b].c[i]=0;
                                                             benmatrix[b].r[i]=temp_benthic[b][i+community->benthic[b].ibenmin][k][l];
                                                     }
                                                     benmatrix[b].a[0]=0;
                                                     tridag(&(benmatrix[b]));
                                                     for(i=0 ; i<benmatrix[b].size ; i++){
                                                             community->benthic[b].u_values[i+community->benthic[b].ibenmin][k][l]=benmatrix[b].u[i];
                                                     }
                                             }
                                     }
                             }
                     }
                     
                     if(community->detritus->ts_flag==0){
                             for(k=0 ; k<x ; k++){
                                     for(l=0 ; l<y ; l++){
                                             community->detritus->w_values[k][l]=temp_detritus[k][l]+(dt/2)*((community->detritus->g_values[k][l]+temp_det_g[k][l])-(community->detritus->mu_values[k][l]+temp_det_mu[k][l]));
                                     }
                             }
                     }
             }
     }
     
     /***************************/
     /* Add reproduction values */ 
     /***************************/
     for(k=0 ; k<x ; k++){
             for(l=0 ; l<y ; l++){
                     for(s=0 ; s<n ; s++){
                             /*If values aren't given by timeseries*/
                             if(community->pelagic[s].ts_flag==0){
                                     /*Specifying Total Eggs*/
                                     if(community->pelagic[s].rep_method==0 || community->pelagic[s].rep_method==1){
                                             community->pelagic[s].u_values[community->pelagic[s].ipelmin][k][l]=community->pelagic[s].reproduction[k][l];
                                     }
                                     /*Adding Eggs*/
                                     if(community->pelagic[s].rep_method==2 || community->pelagic[s].rep_method==3){
                                             community->pelagic[s].u_values[community->pelagic[s].ipelmin][k][l]+=community->pelagic[s].reproduction[k][l];
                                     }
                             }
                     }
                     
                     if(run->no_benthic!=0){
                             for(b=0 ; b<c ; b++){
                                     /*If values aren't given by timeseries*/
                                     if(community->benthic[b].ts_flag==0){
                                             /*Constant Egg amount*/
                                             if(community->benthic[b].rep_method==0){
                                                     community->benthic[b].u_values[community->benthic[b].ibenmin][k][l]=community->benthic[b].u_0*exp(community->benthic[b].lambda*community->benthic[b].benmin);
                                             }
                                             /*Eggs produced according to biomass intake*/
                                             if(community->benthic[b].rep_method==1 || community->benthic[b].rep_method==2){
                                                     community->benthic[b].u_values[community->benthic[b].ibenmin][k][l]+=community->benthic[b].reproduction[k][l];
                                             }
                                     }
                             }
                     }
             }
     }
     
     /*Free pointers*/
     for(k=0 ; k<x ; k++){
             free(temp_det_g[k]);
             free(temp_det_mu[k]);
     }
     free(temp_det_g);
     free(temp_det_mu);
     
}
