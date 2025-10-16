#include <stdio.h>
#include <math.h>
#include "initialize.h"
#include "main.h"
#include "mt19937ar.h"

// set initial coordinates of particles
void initialize_Cubic_BD(){
// cubic root of number of particles, determine number of particles in 1 string or row
  int cbrt_N = cbrt(N_PART);
  int i, j, k;
  int index; // particle index
  double lattice_step = box_size / (1.0 * cbrt_N); // spacing between particles

  for (i=0; i<cbrt_N ; i++){
    for (j=0; j<cbrt_N; j++){
      for (k=0; k<cbrt_N; k++){
        index = cbrt_N * cbrt_N * k + cbrt_N * j + i;
        // set initial particle coordinates and forces
        x[index] = lattice_step * (0.5 + i);
        y[index] = lattice_step * (0.5 + j);
        z[index] = lattice_step * (0.5 + k);
        fpair_x[index] = 0.0;
        fpair_y[index] = 0.0;
        fpair_z[index] = 0.0;
      }
    }
  }
}

// set initial sizes sigmas of particles
void initialize_sigmas(){

  int i;
  for (i=0; i<N_PART; i++){
    fpair_sigma[i] = 0.0;
    fsingle_sigma[i] = 0.0;

// set sigma values randomly according to standart distribution with mean SIGMA_0 and 
// standart deviation SIGMA_STD
    sigma[i] = SIGMA_0 + SIGMA_STD * GaussrandMarsaglia();
    // Get rid of too small and too large sigma values
    if (sigma[i] < 0.3 && sigma[i] > 1.7){
        sigma[i] = SIGMA_0;
    }
  }
}

double GaussrandMarsaglia()
{
	/* https://en.wikipedia.org/wiki/Marsaglia_polar_method
	 * adapted from (after tests) + alternative (Box-Muller / CLT approaches) in
	 * http://c-faq.com/lib/gaussian.html */
	/* returns Gaussian random numbers with mean 0 and standard dev 1.0 */

    static double V1, V2, S;
    static unsigned long long phase = 0;
    double X;

    if(phase == 0) {
        do {
			/* The choice of uniform random number generator doesn't matter with
			 * the Marsaglia polar method. Still, keeping MT for possible future use */
			/* good uniform random nums using Mersenne-Twister */
            double U1 = genrand_real1();
            double U2 = genrand_real1();
			/* standard rand() */
			//double U1 = (double)rand() / RAND_MAX;
			//double U2 = (double)rand() / RAND_MAX;

            V1 = 2.0 * U1 - 1.0;
            V2 = 2.0 * U2 - 1.0;
            S = V1 * V1 + V2 * V2;
            } while(S >= 1.0 || S == 0.0);

        X = V1 * sqrt(-2.0 * log(S) / S);
    }
    else
        X = V2 * sqrt(-2.0 * log(S) / S);

    phase = 1.0 - phase;
    return X;
}
