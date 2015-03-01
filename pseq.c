#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// numerical parameters
#define N 100
#define STEP 1.0/(N-1)
#define STEP2 STEP*STEP
#define TOL 1e-6
#define MAXITER 10000

// the f and r coeff functions
double ffun(const double x) {
	return -1.0*pow(x, -1.5);
}
double rfun(const double x) {
	return pow(x, 3.0);
}


int main(int argc, char *argv[])
{
	//double *f, *r, *u1, *u2;
	double *f = (double *) malloc( (N)*sizeof(double) );
	double *r = (double *) malloc( (N)*sizeof(double) );
	double *u1 = (double *) malloc( (N)*sizeof(double) );
	double *u2 = (double *) malloc( (N)*sizeof(double) );

	// define function values and give initial guess for u:
	for (int i=0; i<N; i++) 	{f[i] = ffun(i*STEP);}
	for (int i=0; i<N; i++) 	{r[i] = rfun(i*STEP);}
	for (int i=1; i<N; i++) 	{u1[i] = 1.0;}

	// boundary values
	u1[0] = u2[0] = u1[N-1] = u2[N-1] = 0.0;
	
	// report
	printf( "----------------------\n");
	printf( "N: %d, x0: %2f, x(N-1): %2f, h: %8f.\n", N, 0*STEP, (N-1)*STEP, STEP );

	double max_diff = TOL*2;
	int iter = 0;
	while(/*max_diff > TOL &&*/ iter < MAXITER) 
	{
		double biggest = 0.0;
		for(int i=1; i<(N-1); i++)	// iterate inner elements
		{
			u2[i] = (u1[i-1] + u1[i+1] - STEP2*f[i])/(2.0-STEP2*r[i]);

			double current_diff = fabs(u1[i] - u2[i]);
			if(current_diff > biggest){
				biggest = current_diff;
			}
		}
		max_diff = biggest;
		++iter;

		// swap buffers
		double *utemp = u1;
		u1 = u2;		// now u1 points to the new values
		u2 = utemp;		// u2 will be overwritten next iteration
	}

    // write to file
 	FILE *fp;

    fp = fopen("u_values_mh.txt", "w");
    for (int i = 0; i < N; i++){
    	fprintf(fp, "%f", u1[i]);
        fprintf(fp, "\n");
    }

    fclose(fp);

    // show the maximum difference:
	printf( "max_diff: %8f\n",max_diff );
	printf("iterations: %d\n", iter);
	

    return 0;
}