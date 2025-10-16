#include <stdio.h>
#include <math.h>
#include "initialize.h"
#include "main.h"
#include "evolve.h"
#include "mt19937ar.h"
#include "force.h"

// coordinate step
void update_coord(){
  // printf("\n New step\n\n");
  int i;
  double friction_sigma[N_PART]; //friction coefficient for each particle
  double rand_coef; // coefficient for random force with mean 0 and variance sqrt(2 D dt)

  // 1. calc forces fpair for coord and sigma
  force_pair();

  for(i=0; i<N_PART; i++){

    friction_sigma[i] = FRICTION_0 * sigma[i] / SIGMA_0;
    rand_coef = sqrt(2.0 * K * T_0 * ALPHA_T * TIMESTEP / friction_sigma[i]);

    // 2. update coordinates:
    x[i] += (fpair_x[i] / friction_sigma[i]) * TIMESTEP + rand_coef * GaussrandMarsaglia();
    y[i] += (fpair_y[i] / friction_sigma[i]) * TIMESTEP + rand_coef * GaussrandMarsaglia();
    z[i] += (fpair_z[i] / friction_sigma[i]) * TIMESTEP + rand_coef * GaussrandMarsaglia();

  }
}

// property sigma step
void update_sigma(){
  int i;
  double force_sigma[N_PART];
  double friction_sigma[N_PART];
  double rand_coef;

  // 1. Calc fsingle_sigma
  force_property();
  double old_sigma;

  for (i=0; i<N_PART; i++){
    old_sigma = sigma[i]; // save old sigma to return to it in case of possibility of negative sigma

	// Changes UB : Unlike with translation, there is no physical reasoning why friction_sigma
	// should show Stokes' behavior. So we are using constant friction_sigma in the property
	// evolution. If we decide to change this later, we should discuss and identify a physically
	// motivated form.
    friction_sigma[i] = FRICTION_0;

    // Calc common force for sigma
    force_sigma[i] = fpair_sigma[i] + fsingle_sigma[i];
    rand_coef = sqrt(2.0 * K * T_0 * LAMBDA * TIMESTEP / friction_sigma[i]); // coef for random force
    // 2. update sigma
    sigma[i] += (force_sigma[i] / friction_sigma[i]) * TIMESTEP + rand_coef * GaussrandMarsaglia();

    if (sigma[i] <= 0.05){
      sigma[i] = old_sigma; // in case of too small sigma not to get negative value return to the old value of sigma
      sigma_vanish++;
      //printf("Warning: Sigma value is less than zero!\n");
	  // Comment UB : it would be good to print an error message here and keep track of what
	  // fraction of times this hard boundary is encountered. If too much, the simulation results
	  // will be wrong.
    }

  }

}
