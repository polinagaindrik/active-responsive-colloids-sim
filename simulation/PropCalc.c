#include <stdio.h>
#include <math.h>
#include "PropCalc.h"
#include "main.h"

void initialize_RDF_N(){
  int n, l;
  // For g(r) create bins
  delta_r = 0.05;
  for (n=0; n<NUM_VAL_G; n++){
    r_values[n] = delta_r * n + 0.5 * delta_r;
    g[n] = 0;
  }
  // For N(sigma) create bins
  delta_sigma = 0.03;
  for (l=0; l<NUM_VAL_N; l++){
    sigma_val[l] = delta_sigma * l + 0.5 * delta_sigma;
    N_sigma[l] = 0;
  }
}

void RDF_onestep(){
  int i, j, k;
  double normalization, r_next, N_count;

  for (k=0; k<NUM_VAL_G; k++){
    N_count = 0.;
    r_next = r_values[k] + delta_r;
    normalization = (1.0 * N_PART) * DENSITY * 3.14 * 4. * ((r_next * r_next * r_next - r_values[k] * r_values[k] * r_values[k])) / 3.;
    for (i=0; i<N_PART; i++){
      for (j=(i+1); j<N_PART; j++){
        if (rij[i][j] != 0 && rij[i][j] >=  r_values[k] && rij[i][j] < r_values[k] + delta_r){
          N_count += 2.0;
        }
      }
    }
    g_perstep[k] = N_count / normalization;
  }
}

void RDF_whole(){
  int k;
  RDF_onestep();
  for (k=0; k<NUM_VAL_G; k++){
    g[k] += g_perstep[k] / num_steps; // mean over all saved steps
}
}

// calc size distrirbution N(sigma) for 1 step
void size_distrib(){
  int i, k;
  double N_count;
  double normalization = 1.0 * N_PART;
  for (k=0; k<NUM_VAL_N; k++){
    N_count = 0.;
    for (i=0; i<N_PART; i++){
        if (sigma[i] >= sigma_val[k] && sigma[i] < sigma_val[k] + delta_sigma){
          N_count += 1.0;
        }
    }
    N_sigma_perstep[k] = N_count / normalization;
  }
}

void size_distrib_whole(){
  int jj;
  // N(sigma) calc
  size_distrib();
  for (jj=0; jj<NUM_VAL_N; jj++){
    N_sigma[jj] += N_sigma_perstep[jj] / num_steps; // mean over all saved steps
  }
}
