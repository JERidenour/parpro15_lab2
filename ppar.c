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
	return -1.0 * pow(x, -1.5);
}
double rfun(const double x) {
	return pow(x, 3.0);
}

/*The idea is a red-black (chequerboard) coloring:
• Even p: assign red
• Odd p: assign black
Communication appears in two steps: red/black and black red:
     if mycolor == red
        send(u[Ip_2],p+1);
        receive(u[Ip-1],p+1);
        send(u[1],p-1);
        receive(u[0],p-1);
else
        receive(u[0],p-1);
        send(u[1],p-1);
        receive(u[Ip-1],p+1);
        send(u[Ip-2],p+1);
end
Communication time is only doubled compared to the previous version.*/

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

	bool red = true; 				// are we red?
	if (rank % 2) {red = false;}	// rank is odd if there exists a remainder

	int L = N / P;				// trunc
	int R = N % P;
	int I = (N + P - rank - 1) / P;    	// number of local elements
	//n = p*L+MIN(rank,R)+i; 	// global index for given (p,i)

	// TODO: set correct values for array sizes

	double *f = (double *) malloc( (I) * sizeof(double) );
	double *r = (double *) malloc( (I) * sizeof(double) );
	double *u1 = (double *) malloc( (I) * sizeof(double) );
	double *u2 = (double *) malloc( (I) * sizeof(double) );

	// compute global indexes g for each local i
	int *g = (int*) malloc( (I) * sizeof(int) );	// global index table
	for (int i = 0; i < I; ++i) {
		g[i] = rank * L + MIN(rank, R) + i;	// global index for given (p,i)
	}

	// define function values using g[i] * STEP
	for (int i = 0; i < I; ++i) 	{f[i] = ffun(g[i] * STEP);}
	for (int i = 0; i < I; ++i) 	{r[i] = rfun(g[i] * STEP);}

	// initial guess for u
	for (int i = 1; i < I; ++i) 	{u1[i] = 1.0;}

	// boundary values
	if (rank == 0) {u1[0] = u2[0] = 0.0;}			// leftmost chunk
	if (rank == (P-1)) {u1[I-1] = u2[I-1] = 0.0;}	// rightmost chunk

	double max_diff = TOL * 2;
	int iter = 0;
	while (max_diff > TOL && iter < MAXITER)
	{

		// TODO: RB communication of overlap here
		if (red) {
			// send(u[Ip-2],p+1);
			MPI_Send(	u1 + I - 2,			// void* data
			            1,					// int count
			            MPI_DOUBLE,			// MPI_Datatype datatype
			            rank + 1,			// int destination
			            666,				// int tag
			            MPI_COMM_WORLD);	// MPI_Comm communicator
			// receive(u[Ip-1],p+1);
			MPI_Recv(	u1 + I - 1,			// void* data
			            1,					// int count
			            MPI_DOUBLE,			// MPI_Datatype datatype
			            rank + 1,			// int source
			            666,				// int tag
			            MPI_COMM_WORLD,		// MPI_Comm communicator
			            MPI_STATUS_IGNORE);	// MPI_Status * status
			// send(u[1],p-1);
			MPI_Send(	u1 + 1,
			            1,
			            MPI_DOUBLE,
			            rank - 1,
			            666,
			            MPI_COMM_WORLD);
			//receive(u[0],p-1);
			MPI_Recv(	u1,
			            1,
			            MPI_DOUBLE,
			            rank - 1,
			            666,
			            MPI_COMM_WORLD,
			            MPI_STATUS_IGNORE);
		} else // black
		{
			// receive(u[0], p - 1);
			MPI_Recv(	u1,
			            1,
			            MPI_DOUBLE,
			            rank - 1,
			            666,
			            MPI_COMM_WORLD,
			            MPI_STATUS_IGNORE);
			// send(u[1], p - 1);
			MPI_Send(	u1 + 1,
			            1,
			            MPI_DOUBLE,
			            rank - 1,
			            666,
			            MPI_COMM_WORLD);
			// receive(u[Ip - 1], p + 1);
			MPI_Recv(	u1 + I - 1,
			            1,
			            MPI_DOUBLE,
			            rank + 1,
			            666,
			            MPI_COMM_WORLD,
			            MPI_STATUS_IGNORE);
			// send(u[Ip - 2], p + 1);
			MPI_Send(	u1 + I - 2,
			            1,
			            MPI_DOUBLE,
			            rank + 1,
			            666,
			            MPI_COMM_WORLD);
		}

		// TODO: rewrite below as the local iteration step
		double biggest = 0.0;
		for (int i = 1; i < N; i++)	// iterate inner elements
		{
			// unlike this, Hanke uses i,i+2 and assumes shifted f and r arrays
			u2[i] = (u1[i - 1] + u1[i + 1] - STEP2 * f[i]) / (2.0 - STEP2 * r[i]);

			double current_diff = fabs(u1[i] - u2[i]);
			if (current_diff > biggest) {
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

	// TODO: rewrite

	// write to file
	FILE *fp;

	fp = fopen("u_values_mh.txt", "w");
	for (int i = 0; i < N + 1; i++) {
		fprintf(fp, "%f", u1[i]);
		fprintf(fp, "\n");
	}

	fclose(fp);

	// show the maximum difference:
	printf( "max_diff: %8f\n", max_diff );
	printf("iterations: %d\n", iter);


	return 0;
}