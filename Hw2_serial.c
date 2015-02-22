#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// numerical parameters
#define N 100
#define h 1.0/(N+1)
#define tol 1e-6

// compare is used by qsort() to compare values in an array
static int compare (const void * a, const void * b)
{
  if (*(double*)a > *(double*)b) return 1;
  else if (*(double*)a < *(double*)b) return -1;
  else return 0;  
}

int main(int argc, char *argv[])
{

	double f[N+1], r[N+1], u[N+1], prev_u[N+1], u_diff[N+1];
	double max_diff = 1.0;
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

	while(max_diff > tol){

		// save old values 
		for (i=0; i<N+1; i++){
			prev_u[i] = u[i];
		}

		// update u
		for(i=1; i<N; i++){
			u[i] = (u[i-1] + u[i+1] - h*h*f[i])/(2-h*h*r[i]);
		}

		// get the difference between iterations
		for (i=0; i<N+1; i++){
			u_diff[i] = fabs(u[i]-prev_u[i]);
		}		

		// sort u_diff to find maximum value
		qsort(u_diff, N+1, sizeof(u_diff[0]), compare);

		// set max value to check against tolerance
		max_diff = u_diff[N];
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
	printf( "max_diff:\n" );
	printf("%8f", max_diff);
	printf("\n");

    return 0;
}