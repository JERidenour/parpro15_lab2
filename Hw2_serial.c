#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// numerical parameters
#define N 1000
#define h 1.0/(N+1)
#define tol 1e-6
#define MAXITER 100000

int main(int argc, char *argv[])
{

	double f[N+1], r[N+1], u[N+1], prev_u[N+1];
	double max_diff = tol*2;
	int i;

	// define function values and give initial guess for u:
	for (i=0; i<N+1; i++){
		f[i] = -1.0*pow(i*h, -1.5);
		r[i] = pow(i*h, 3.0);
		u[i] = 1.0;
	}

	// boundary values:
	u[0]=0.0;
	u[N]=0.0;

	int iter = 0;
	while(iter < MAXITER){

		// save old values 
		for (i=0; i<N+1; i++){
			prev_u[i] = u[i];
		}

		double biggest = 0.0;

		// update u
		for(i=1; i<N; i++){
			u[i] = (u[i-1] + u[i+1] - h*h*f[i])/(2-h*h*r[i]);
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

    fp = fopen("u_values_jr.txt", "w");
    for (i = 0; i < N+1; i++){
    	fprintf(fp, "%f", u[i]);
        fprintf(fp, "\n");
    }

    fclose(fp);

    // show the maximum difference:
	printf( "max_diff:\n" );
	printf("%8f", max_diff);
	printf("\n");

    return 0;
}