#include <stdio.h>
#include <math.h>
#include "initialize.h"
#include "main.h"
#include "evolve.h"
#include "mt19937ar.h"
#include "force.h"

// Calculate the differencies between x, y and z coordinates for 1 particle
double pbc(int i, int j, int coord){

  double delta_coord;
  double coord1_wr, coord2_wr; // wraped coordinates of 2 particles

// choose what coordinate is it x, y or z
  switch (coord){
    case 1:
      coord1_wr = x[i];
      coord2_wr = x[j];
      break;
    case 2:
      coord1_wr = y[i];
      coord2_wr = y[j];
      break;
    case 3:
      coord1_wr = z[i];
      coord2_wr = z[j];
      break;
    default:
      printf("Error\n");
  }
// obtain coordinates in 1 box (0, L)
  coord1_wr = fmod(coord1_wr, box_size);
  coord2_wr = fmod(coord2_wr, box_size); // now they are in 2 boxes(0, L) and (-L, 0)
  if (coord1_wr < 0.){  // now get rid of negative coordnates
    coord1_wr += box_size;
  }
  if (coord2_wr < 0.){
    coord2_wr += box_size;
  }

  delta_coord = coord1_wr - coord2_wr;
  if (delta_coord > 0.5 * box_size){
    delta_coord -= box_size;
  }
  if (delta_coord < - 0.5 * box_size){
    delta_coord += box_size;
  }
  return delta_coord;
}

void force_pair(){
  int i, j;
  double delta_x, delta_y, delta_z;
  double sigmaij;
  double c, b; //parts of expressions for force calculation

  energy_sys = 0;
  for (i=0; i<N_PART; i++){
    fpair_x[i] = fpair_y[i] = fpair_z[i] = fpair_sigma[i] = 0;
  }

  for (i=0; i<N_PART-1; i++){

    for (j=(i+1); j<N_PART; j++){
		// Change UB : see above
      delta_x = pbc(i, j, 1); // calc coordinate differences
      delta_y = pbc(i, j, 2);
      delta_z = pbc(i, j, 3);
      // distance between particle i and j
	  // Comment UB : sqrt is a VERY EXPENSIVE calculation. Use it only when necessary.
	  // As you are using it co calculate g(r), this is ok. But generally speaking
	  // first compare squares and then if the condition holds, calculate sqrt
      rij[i][j]  = sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);
      sigmaij = 0.5 * (sigma[i] + sigma[j]);

	  // Change UB : because of j=(i+1) above, the i!=j condition is redundant
      //if (rij[i][j] <= sigmaij && i!=j){
      if (rij[i][j] <= sigmaij){
        // Calc forces for coordinates
        b = 1.0 - rij[i][j] / sigmaij;
        c =  2.5 * EPS * pow(b, 1.5) / (sigmaij * rij[i][j]);
        fpair_x[i] += c * delta_x;
        fpair_y[i] += c * delta_y;
        fpair_z[i] += c * delta_z;

        fpair_x[j] -= c * delta_x;
        fpair_y[j] -= c * delta_y;
        fpair_z[j] -= c * delta_z;

        fpair_sigma[i] -= 0.5 * c * rij[i][j] * rij[i][j] / sigmaij;
        fpair_sigma[j] -= 0.5 * c * rij[i][j] * rij[i][j] / sigmaij;

        energy_sys += EPS * pow(b, 2.5); // calc of energy of the whole system at each step
      }
    }
 }
}

void force_property(){
  int i;
  for(i=0; i<N_PART; i++){
    fsingle_sigma[i] = - K * T_0 * (sigma[i] - SIGMA_0) / (SIGMA_STD * SIGMA_STD); // * LAMBDA (in previous siulations, when I scaled psi with T_sigma)
  }

}
