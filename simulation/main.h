#define N_PART 512 // number of particles in simulation box
#define FRICTION_0 1.0 // friction coefficient for the coordinates
#define K 1.0 // Boltzman constant
#define T_0 1.0 // temperature for coordinates evaluation
#define SIGMA_0 1.0 // mean for gaussian distribution of sigmas
#define SIGMA_STD 0.2 // standart deviation for gaussian distribution of sigmas
#define N_PROD 1e7 // number of steps in simulation (production run) //1e4
#define N_EQ 1e6 // number of steps for equilibration run
#define EPS 500.0
#define TIMESTEP 1e-4
#define FREQOUT 1000 // out frequency for production run
#define FREQOUTEQ 10000 // out frequency for equilibration

#define ALPHA_T 1.0 // constant that determines T_x temperature // 0.1 0.2 0.5 1 2 5 10
#define DENSITY 0.95 // number density // 0.02   //0.20 //0.60 // 0.95 // 1.30
#define LAMBDA 0.2 //   constant that determines T_sigma temperature // 0.1 0.2 0.5 1 2 5 10

#define NUM_VAL_G 100 // number of divided segments of rij values to calc RDF
#define NUM_VAL_N 100


double x[N_PART], y[N_PART], z[N_PART]; // x, y, z coordinates of each particle
double fpair_x[N_PART], fpair_y[N_PART], fpair_z[N_PART]; // forces to update coordinates
double fpair_sigma[N_PART], fsingle_sigma[N_PART]; // forces to update sigma
double sigma[N_PART]; // sizes of each particle
double rij[N_PART][N_PART]; // distances between particles
double box_size; // size of simulation box
double energy_sys; //energy of the system

int sigma_vanish; // count times when sigma value becomes less than zero
int num_steps; //number of steps saved

double g[NUM_VAL_G]; // RDF
double g_perstep[NUM_VAL_G];
double delta_r;
double r_values[NUM_VAL_G];

double N_sigma_perstep[NUM_VAL_N];
double N_sigma[NUM_VAL_N];
double delta_sigma;
double sigma_val[NUM_VAL_N];
