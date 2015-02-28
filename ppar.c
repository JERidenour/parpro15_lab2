#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// numerical parameters
#define N 100
#define STEP 1.0/(N+1)
#define STEP2 STEP*STEP
#define TOL 1e-6
#define MAXITER 10000

#define MIN(a,b) ((a) < (b) ? (a) : (b))

// the f and r coeff functions
double ffun(const double x) {
	return -1.0*pow(x, -1.5);
}
double rfun(const double x) {
	return pow(x, 3.0);
}

int main(int argc, char *argv[])
{
	int P, rank;

	/* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (N < P) {
		fprintf(stdout, "Too few discretization points...\n");
		exit(1);
    }

    bool black = false; 			// are we black?
    if (rank % 2) {black = true;} 	// rank is odd if there exists a remainder

    int L = N/P;				// trunc
    int R = N%P;
    int I = (N+P-rank-1)/P;    	// number of local elements
    //n = p*L+MIN(rank,R)+i; 	// global index for given (p,i)

    // TODO: parallellize below
    
	//double *f, *r, *u1, *u2;
	double *f = (double *) malloc( (N+1)*sizeof(double) );
	double *r = (double *) malloc( (N+1)*sizeof(double) );
	double *u1 = (double *) malloc( (N+1)*sizeof(double) );
	double *u2 = (double *) malloc( (N+1)*sizeof(double) );

	// define function values and give initial guess for u:
	for (int i=0; i<N+1; i++) 	{f[i] = ffun(i*STEP);}
	for (int i=0; i<N+1; i++) 	{r[i] = rfun(i*STEP);}
	for (int i=1; i<N; i++) 	{u1[i] = 1.0;}

	// boundary values
	u1[0] = u2[0] = u1[N] = u2[N] = 0.0;

	double max_diff = TOL*2;
	int iter = 0;
	while(max_diff > TOL && iter < MAXITER) 
	{
		double biggest = 0.0;
		for(int i=1; i<N; i++)	// iterate inner elements
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
    for (int i = 0; i < N+1; i++){
    	fprintf(fp, "%f", u1[i]);
        fprintf(fp, "\n");
    }

    fclose(fp);

    // show the maximum difference:
	printf( "max_diff: %8f\n",max_diff );
	printf("iterations: %d\n", iter);
	

    return 0;
}