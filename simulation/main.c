#include <stdio.h>
#include <math.h>
#include "main.h"
#include "initialize.h"
#include "evolve.h"
#include "mt19937ar.h"
#include "force.h"
#include "output.h"
#include "PropCalc.h"

int main()
 {
  int n, l, j, i;
  double fraction_vanish; // fraction of time when perticle vanishes (sigma<0)
  box_size = cbrt(N_PART / DENSITY);
  char TrajFileEqWr[600], TrajFileEqUnwr[600], TrajFileProdWr[600], TrajFileProdUnwr[600], rdf[600], SigmaDistrib[600], EnergyValues[600];

  num_steps = N_PROD / FREQOUT;
  int seed = 2777;
  init_genrand(seed);

  // distribute particles on a 3D lattice and set sigmas
  initialize_Cubic_BD();
  initialize_sigmas();
  initialize_RDF_N(); // initialize values that will be used to calculate RDF and N(sigma)

// with index 2 it's when we don't have sigma-tenperature in psi potential. Potential scales with k_B T_0
  sprintf(TrajFileEqWr,"./Eq_traj_wrapped2_density_%.2f_Npart_%d_lambda_%.1f_Freqout_%d.dump", DENSITY, N_PART, LAMBDA, FREQOUT);
  sprintf(TrajFileEqUnwr,"./Eq_traj_unwrapped2_density_%.2f_Npart_%d_lambda_%.1f_Freqout_%d.dump", DENSITY, N_PART, LAMBDA, FREQOUT);
  int percent;

  // Equilibration run:
  sigma_vanish = 0;
  printf("Equilibration run:\n");
  for (j=0; j<N_EQ; j++){
    update_coord();
    update_sigma();

    if (j%FREQOUTEQ == 0){
      percent = 100 * j / N_EQ;
      printf("Step: %d - %d percent\n", j, percent);
      print_wrapped_coords(j, TrajFileEqWr);
      print_unwrapped_coords(j, TrajFileEqUnwr);
    }
  }
  fraction_vanish = sigma_vanish / N_EQ;
  printf("Fraction of time when sigma becomes less than zero: %g\n", fraction_vanish);

  // Production run
  sprintf(TrajFileProdWr,"./Prod_traj_wrapped2_density_%.2f_Npart_%d_lambda_%.1f_Freqout_%d.dump", DENSITY, N_PART, LAMBDA, FREQOUT);
  sprintf(TrajFileProdUnwr,"./Prod_traj_unwrapped2_density_%.2f_Npart_%d_lambda_%.1f_Freqout_%d.dump", DENSITY, N_PART, LAMBDA, FREQOUT);

  sprintf(rdf,"./RDF2_density_%.2f_Npart_%d_lambda_%.1f_Freqout_%d.txt", DENSITY, N_PART, LAMBDA, FREQOUT);
  sprintf(SigmaDistrib,"./N_sigma2_density_%.2f_Npart_%d_lambda_%.1f_Freqout_%d.txt", DENSITY, N_PART, LAMBDA, FREQOUT);
  sprintf(EnergyValues,"./Energy2_density_%.2f_Npart_%d_lambda_%.1f_Freqout_%d.txt", DENSITY, N_PART, LAMBDA, FREQOUT);
  sigma_vanish = 0;

  printf("\nProduction run:\n");
   for (i=0; i<N_PROD; i++){
     update_coord();
     update_sigma();

     if (i%FREQOUT == 0){
       percent = 100 * i / N_PROD;
       printf("Step: %d - %d percent\n", i, percent);
       print_wrapped_coords(i, TrajFileProdWr);
       print_unwrapped_coords(i, TrajFileProdUnwr);
       print_energy(i * TIMESTEP, energy_sys, EnergyValues);// save energy

       RDF_whole();
       size_distrib_whole();
   }
  }
  fraction_vanish = sigma_vanish / N_PROD;
  printf("Fraction of time when sigma becomes less than zero: %g\n", fraction_vanish);

  print_Property(NUM_VAL_G, r_values, g, rdf);
  print_Property(NUM_VAL_N, sigma_val, N_sigma, SigmaDistrib);

  return 0;
}
