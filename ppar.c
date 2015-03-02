#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>

// numerical parameters
#define N 1000   			// n is [0,N-1]
//#define STEP 1.0/(N-1)
#define STEP 1.0/(N+1)		// x is [h,1.0-h] because we expanded u with 2
#define STEP2 STEP*STEP
//#define TOL 1e-6
#define MAXITER 1000000

#define MIN(a,b) ((a) < (b) ? (a) : (b))

// the f and r coeff functions
double ffun(const double x) {
	return -1.0 * pow(x, -1.5);
}
double rfun(const double x) {
	return pow(x, 3.0);
}

int main(int argc, char *argv[])
{
	int P, rank, rc;

	/* Initialize MPI */
	rc = MPI_Init(&argc, &argv);
	rc = MPI_Comm_size(MPI_COMM_WORLD, &P);
	rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (N < P) {
		fprintf(stdout, "Too few discretization points...\n");
		exit(1);
	}

	// are we red (if not, we are black)?
	bool red = (rank % 2) ? false : true;
	//bool red = true;
	//if (rank % 2) {red = false;}	// rank is odd if there exists a remainder

	// are we leftmost or rightmost chunk?
	bool leftbound = rank == 0 ? true : false;
	/*bool leftbound = false;
	if (rank == 0) {leftbound = true;}*/
	bool rightbound = rank == (P - 1) ? true : false;
	/*bool rightbound = false;
	if (rank == P - 1) {rightbound = true;}*/
	//int expand = (leftbound || rightbound ? 1 : 2);
	int expand = 2;

	int L = N / P;				// trunc
	int R = N % P;
	int I = ((N + P - rank - 1) / P);    	// number of inner local elements
	int Ie = I + expand;
	//n = p*L+MIN(rank,R)+i; 	// global index for given (p,i)

	// report to console
	printf( "----------------------\n");
	printf( "Started process # %d, %s, %s, %s.\n", rank, (red ? "red" : "black"),
	        (leftbound ? "leftbound" : "-"),
	        (rightbound ? "rightbound" : "-")  );
	printf( "Process # %d, P: %d, N: %d, h: %8f L: %d, R: %d, I: %d, Ie: %d\n",rank, P, N, STEP, L, R, I, Ie );


	double *f = (double *) malloc( (Ie) * sizeof(double) );
	double *r = (double *) malloc( (Ie) * sizeof(double) );
	double *u1 = (double *) malloc( (Ie) * sizeof(double) );
	double *u2 = (double *) malloc( (Ie) * sizeof(double) );

	// compute global indexes g for each local i
	int *g = (int*) malloc( (Ie) * sizeof(int) );	// global index lookup table
	for (int i = 0; i < Ie; ++i) {
		g[i] = rank * L + MIN(rank, R) + i;			// global index for given (p,i)
	}

	printf( "Process # %d, g[0]: %d, g[1]: %d, g[end-1]: %d, g[end]: %d\n", rank, g[0], g[1], g[Ie-2], g[Ie-1]);

	// define function values using x = g[i] * STEP
	for (int i = 0; i < Ie; ++i) 	{f[i] = ffun(g[i] * STEP);}
	for (int i = 0; i < Ie; ++i) 	{r[i] = rfun(g[i] * STEP);}

	// initial guess for u
	for (int i = 0; i < Ie; ++i) 	{u1[i] = 1.0;}

	// boundary values
	if (leftbound) {u1[0] = u2[0] = 0.0;}
	if (rightbound) {u1[Ie - 1] = u2[Ie - 1] = 0.0;}

	//double max_diff = TOL * 2;	// bogus value
	int iter = 0;				// keep track of iterations
	while (/*max_diff > TOL &&*/ iter < MAXITER)	// each process needs to iterate same amount
	{
		//printf( "-process # %d, iter %d.\n", rank, iter);
		if (red) {
			if (!rightbound) {
				// send(u[Ip-2],p+1);
				MPI_Send(	u1 + Ie - 2,			// void* data
				            1,					// int count
				            MPI_DOUBLE,			// MPI_Datatype datatype
				            rank + 1,			// int destination
				            1,				// int tag
				            MPI_COMM_WORLD);	// MPI_Comm communicator
				// receive(u[Ip-1],p+1);
				MPI_Recv(	u1 + Ie - 1,			// void* data
				            1,					// int count
				            MPI_DOUBLE,			// MPI_Datatype datatype
				            rank + 1,			// int source
				            2,				// int tag
				            MPI_COMM_WORLD,		// MPI_Comm communicator
				            MPI_STATUS_IGNORE);	// MPI_Status * status
			}
			if (!leftbound) {
				// send(u[1],p-1);
				MPI_Send(	u1 + 1,
				            1,
				            MPI_DOUBLE,
				            rank - 1,
				            3,
				            MPI_COMM_WORLD);
				//receive(u[0],p-1);
				MPI_Recv(	u1,
				            1,
				            MPI_DOUBLE,
				            rank - 1,
				            4,
				            MPI_COMM_WORLD,
				            MPI_STATUS_IGNORE);
			}
		} else // black
		{
			if (!leftbound) {
				// receive(u[0], p - 1);
				MPI_Recv(	u1,
				            1,
				            MPI_DOUBLE,
				            rank - 1,
				            1,
				            MPI_COMM_WORLD,
				            MPI_STATUS_IGNORE);
				// send(u[1], p - 1);
				MPI_Send(	u1 + 1,
				            1,
				            MPI_DOUBLE,
				            rank - 1,
				            2,
				            MPI_COMM_WORLD);
			}
			if (!rightbound) {
				// receive(u[Ip - 1], p + 1);
				MPI_Recv(	u1 + Ie - 1,
				            1,
				            MPI_DOUBLE,
				            rank + 1,
				            3,
				            MPI_COMM_WORLD,
				            MPI_STATUS_IGNORE);
				// send(u[Ip - 2], p + 1);
				MPI_Send(	u1 + Ie - 2,
				            1,
				            MPI_DOUBLE,
				            rank + 1,
				            4,
				            MPI_COMM_WORLD);
			}
		}

		//double biggest = 0.0;
		for (int i = 1; i < (Ie - 1); i++)	// iterate inner elements
		{
			// unlike this, Hanke uses i,i+2 and assumes shifted f and r arrays
			u2[i] = (u1[i - 1] + u1[i + 1] - STEP2 * f[i]) / (2.0 - STEP2 * r[i]);

			/*double current_diff = fabs(u1[i] - u2[i]);
			if (current_diff > biggest) {
				biggest = current_diff;
			}*/
		}
		//max_diff = biggest;
		++iter;

		// swap buffers
		double *utemp = u1;
		u1 = u2;		// now u1 points to the new values
		u2 = utemp;		// u2 will be overwritten next iteration
	}

	// report to console
	printf( "Process # %d done iteration.\n", rank);
	//printf( "max_diff: %8f\n", max_diff );
	//printf("iterations: %d\n", iter);
	//printf( "----------------------\n");

	MPI_Barrier(MPI_COMM_WORLD);

	if (leftbound) {
		FILE *fp;
		fp = fopen("u_values_par.txt", "w");
		int c = 0;
		for (int i = 0; i < (Ie - 1); ++i) {	// write left boundary but not right ghost
			fprintf(fp, "%f\n", u1[i]);
			++c;
		}
		printf("Leftbound wrote %d own elements.\n", c);
		// we expect at most I+2 doubles to arrive, but could be less
		MPI_Status rstat;	// indata info
		double *rbuf = (double *) malloc( (Ie + 2) * sizeof(double) ); // inbuffer

		for (int source = 1; source < P; ++source)
		{
			// receive all the middle chunks in order, blocking
			MPI_Recv(	rbuf,
			            Ie + 2,
			            MPI_DOUBLE,
			            source,
			            source,	// tag
			            MPI_COMM_WORLD,
			            &rstat);

			int rsize = 0;
			MPI_Get_count(&rstat, MPI_DOUBLE, &rsize);	// how many doubles arrived?
			// we should get more than 0, report if this happens
			if (rsize == 0) {fprintf(stdout, "Leftbound received 0 elements from process %d!\n", source);}
			printf("Leftbound received %d elements from %d.\n", rsize, source);
			// include right boundary if this is last chunk
			int end = rsize - 1;
			if (source == (P - 1)) {end = rsize;}

			// append to file
			int cc = 0;
			for (int i = 1; i < end; ++i) {
				fprintf(fp, "%f\n", rbuf[i]);
				++cc;
			}
			printf("Leftbound wrote %d elements from %d.\n", cc, source);
		}

		// cleanup
		fclose(fp);
		free(rbuf);

	} else // we are one of the others
	{
		// send what you got to leftbound
		MPI_Send(	u1,
		            Ie,
		            MPI_DOUBLE,
		            0,
		            rank,
		            MPI_COMM_WORLD);
	}

	// cleanup
	free(f);
	free(r);
	free(u1);
	free(u2);
	free(g);
	// over and out
	MPI_Finalize();
	return 0;
}