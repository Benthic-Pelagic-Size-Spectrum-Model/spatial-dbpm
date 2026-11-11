static double Cfun(double, double, GRID *, PELAGIC *);
static double Dfun(double, double, GRID *, PELAGIC *);
static double Diffun(double, double, GRID *, PELAGIC *);
static int horizontal_face_valid(GRID *, int, int); // Identifies faces that are open for fluxes (both water cells)
static int vertical_face_valid(GRID *, int, int); // " "

// Global matrix integration functions (also used in SizeSpectra.c)
void tridag(MATRIX *);                                                     // Inverts a tridiagonal matrix
void trimul(MATRIX *);

void xmove_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, double ****, int); // Calculates x movement in 1D space
void ymove_solver(RUN *, GRID *, COMMUNITY *, MATRIX *, double ****, int); // Calculates y movement in 2D space
