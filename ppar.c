#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpi.h"

// numerical parameters
#define N 100
#define STEP 1.0/(N+1)
#define STEP2 STEP*STEP
#define TOL 1e-6
#define MAXITER 10000


int main(int argc, char *argv[])
{
	//if (x % 2) { /* x is odd */ }

	/* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);

    if (N < P) {
		fprintf(stdout, "Too few discretization points...\n");
		exit(1);
    }



	double f[N+1], r[N+1], u[N+1], prev_u[N+1];
	double max_diff = TOL*2;
	int i;

	// define function values and give initial guess for u:
	for (i=0; i<N+1; i++){
		f[i] = -1.0*pow(i*STEP, -1.5);
		r[i] = pow(i*STEP, 3.0);
		u[i] = 1.0;
	}

	// boundary values:
	u[0]=0.0;
	u[N]=0.0;

	int iter = 0;
	while(max_diff > TOL && iter < MAXITER){

		// save old values 
		for (i=0; i<N+1; i++){
			prev_u[i] = u[i];
		}

		double biggest = 0.0;

		// update u
		for(i=1; i<N; i++){
			u[i] = (u[i-1] + u[i+1] - STEP2*f[i])/(2-STEP2*r[i]);
			double current_diff = fabs(u[i] - prev_u[i]);
			if(current_diff > biggest){
				biggest = current_diff;
			}
		}

		max_diff = biggest;
		++iter;
	}

    // write to file
 	FILE *fp;

    fp = fopen("u_values.txt", "w");
    for (i = 0; i < N+1; i++){
    	fprintf(fp, "%f", u[i]);
        fprintf(fp, "\n");
    }

    fclose(fp);

    // show the maximum difference:
	printf( "max_diff: %8f\n",max_diff );
	printf("iterations: %d\n", iter);
	
	// over and out
	MPI_Finalize();
    return 0;	
}