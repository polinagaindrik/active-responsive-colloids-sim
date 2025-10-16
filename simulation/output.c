#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include "main.h"
#include "initialize.h"
#include "force.h"
#include "evolve.h"
#include "output.h"
#include "mt19937ar.h"

void print_wrapped_coords(long Tstep, char Filename[600])
{
	int i;
	double xr,yr,zr;
	FILE *fp;

	fp = fopen(Filename, "a");

	fprintf(fp,"ITEM: TIMESTEP\n %ld \nITEM: NUMBER OF ATOMS\n %d \nITEM: BOX BOUNDS pp pp pp\n", Tstep, N_PART);
	fprintf(fp,"0.000 \t %.3f \n", box_size);
	fprintf(fp,"0.000 \t %.3f \n", box_size);
	fprintf(fp,"0.000 \t %.3f \n", box_size);
	fprintf(fp,"ITEM: ATOMS id type x y z mass radius\n");	// mass simply for VMD visual coloring. Actually radius

	for(i=0;i<N_PART;i++)
	{
		/* projecting onto [0:boxL] */
		xr = fmod(x[i], box_size);
		yr = fmod(y[i], box_size);
		zr = fmod(z[i], box_size);

		if(xr < 0.)
			xr += box_size;
		if(yr < 0.)
			yr += box_size;
		if(zr < 0.)
			zr += box_size;

		fprintf(fp,"%d \t %d \t %.6f \t %.6f \t %.6f \t %.3f \t %.3f \n", i, i, xr, yr, zr, sigma[i], sigma[i]);
	}

	 fflush(fp);
	 fclose(fp);

}

void print_unwrapped_coords(long Tstep, char Filename[600])
{
	int i;
	FILE *fp;

	fp = fopen(Filename, "a");

	fprintf(fp,"ITEM: TIMESTEP\n %ld \nITEM: NUMBER OF ATOMS\n %d \nITEM: BOX BOUNDS pp pp pp\n", Tstep, N_PART);
	fprintf(fp,"0.000 \t %.3f \n", box_size);
	fprintf(fp,"0.000 \t %.3f \n", box_size);
	fprintf(fp,"0.000 \t %.3f \n", box_size);
	fprintf(fp,"ITEM: ATOMS id type x y z mass radius\n");

	for(i=0;i<N_PART;i++)
		fprintf(fp,"%d \t %d \t %.6f \t %.6f \t %.6f \t %.3f \t %.3f \n", i, i, x[i], y[i], z[i], sigma[i], sigma[i]);

	 fflush(fp);
	 fclose(fp);

}

void print_Property(int num_val, double x_val[100], double property_val[100], char Filename[600]){
	int i;
	FILE *fp;
	fp = fopen(Filename, "a");
	for (i=0; i<num_val; i++){
		fprintf(fp," %g \t %g \n", x_val[i], property_val[i]);
	}
	fflush(fp);
	fclose(fp);

}

void print_energy(double timestep_val, double energy_val, char Filename[600]){
	FILE *fp;
	fp = fopen(Filename, "a");
	fprintf(fp," %g \t %g \n", timestep_val, energy_val);
	fflush(fp);
	fclose(fp);
}
